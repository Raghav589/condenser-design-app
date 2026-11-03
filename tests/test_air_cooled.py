"""Air-cooled condenser module.

BENCHMARK tests reproduce printed worked examples from their inputs: HSB Ex. 9.1, 19.1, 19.2, Ex. 2.11 (air-side pressure drop, which Ex. 9.1 does not work)
and the HEDH 3.8.9 worked problem (fan power and noise). None of these is a vertical downflow case; see docs/algorithms/air_cooled.md. The gravity
film term and the static head of vertical downflow are tested only by limits, hand values, closure and trends, which are never presented as benchmarks.

Tolerances: the book rounds to 3 or 4 s.f. and to rounded intermediates, so 0.3% unless a comment says otherwise.
"""
import copy
import math
from dataclasses import replace

import pytest

from condenser_app.common.condensation import boyko_kruzhilin_local, boyko_kruzhilin_mean, mean_coefficient_over_quality
from condenser_app.common.properties import air as air_table
from condenser_app.common.thermal import overall_u_outer
from condenser_app.modules import MODULES
from condenser_app.modules.air_cooled import (
    METHOD_BK,
    METHOD_VERTICAL,
    SIZE_TUBES,
    SPEC_BOOK,
    SPEC_FLOW,
    SPEC_OUTLET,
    SPEC_VELOCITY,
    WALL_CORBELS,
    WALL_EX91,
    WALL_NONE,
    acceleration_constant,
    air_bundle_pressure_drop,
    alpha_on_bare_area,
    area_bare,
    area_between_fins,
    area_fins,
    area_ratio,
    diagonal_pitch,
    esdu_nusselt,
    fan_dynamic_pressure,
    fan_power,
    fan_section,
    fin_efficiency,
    fin_height,
    fin_psi,
    friction_constant_high_fin,
    hedh_film_thickness,
    min_flow_area,
    natural_draft_height,
    noise_table_at_1m,
    regime_coefficient,
    regime_of,
    row_duty_fractions,
    row_factor,
    shear_parameter,
    sonic_velocity,
    sound_power_level,
    sound_pressure_level,
    total_sound_power_level,
    _duty,
    _evaluate,
)
from condenser_app.reporting import to_excel, to_pdf
from condenser_app.schemas import Category, CheckState, Status
from condenser_app.service import run_design

MODULE = MODULES["air_cooled"]
BOOK = 3e-3  # 0.3%: book rounding


def example():
    return copy.deepcopy(MODULE.example_inputs())


def consistent(spec=SPEC_OUTLET, **air):
    """The Ex. 19.2 bundle with the air balance closed and the module's default method."""
    i = example()
    fields = {"spec": spec, "T_out": None, "v_approach": None, "m_dot": None}
    fields.update(air)
    i.air = replace(i.air, **fields)
    i.tube_side.method = METHOD_VERTICAL
    return i


def run(mutate=None, inputs=None):
    i = inputs if inputs is not None else example()
    if mutate:
        mutate(i)
    return run_design("air_cooled", i)


def step(result, symbol):
    return next(s for s in result.steps if s.symbol == symbol).value


def check(result, name):
    return next(c for c in result.checks if c.name == name)


# --------------------------------------------------------------------- BENCHMARK: HSB Example 9.1, air side, p. 381-386
D_R, D_I, D_F, W, S, P1, P2 = 25.4e-3, 18.6e-3, 57.2e-3, 0.406e-3, 1.903e-3, 60.33e-3, 52.07e-3
RHO, CP, MU, K_AIR, PR_BOOK = 1.161, 1007.0, 1.85e-5, 0.0263, 0.707
N_91, L_91 = 80, 3.0


def test_ex91_areas():
    A_F, A_W = area_fins(N_91, L_91, D_F, D_R, W, S), area_between_fins(N_91, L_91, D_R, W, S)
    assert A_F == pytest.approx(436.4, rel=BOOK)
    assert A_W == pytest.approx(15.8, rel=BOOK)
    assert A_F + A_W == pytest.approx(452.2, rel=BOOK)
    assert area_bare(N_91, L_91, D_R) == pytest.approx(19.15, rel=BOOK)
    assert area_ratio(D_F, D_R, W, S) == pytest.approx(23.61, rel=BOOK)
    assert fin_height(D_R, D_F) == pytest.approx(15.9e-3, rel=BOOK)


def test_ex91_flow_area_and_velocity():
    S_min, plane = min_flow_area(20, L_91, P1, P2, D_R, 15.9e-3, W, S)
    assert plane == "transverse"
    assert S_min == pytest.approx(1.760, rel=BOOK)
    S_0 = 20.5 * P1 * L_91
    assert S_0 == pytest.approx(3.711, rel=BOOK)
    assert S_min / S_0 == pytest.approx(0.474, rel=BOOK)
    # Book's rounded sigma gives 10.55 m/s; exact sigma gives 10.54.
    v_max = 5.0 * S_0 / S_min
    assert v_max == pytest.approx(10.55, rel=BOOK)
    assert v_max * S_min * RHO == pytest.approx(21.56, rel=BOOK)
    assert RHO * v_max * D_R / MU == pytest.approx(1.68e4, rel=BOOK)


def test_ex91_coefficient_chain():
    # The book's own 1.68e4 and Pr = 0.707, so this compares the equations and not rounding upstream.
    nu = esdu_nusselt(1.68e4, S, 15.9e-3, P1, P2, PR_BOOK)
    assert nu == pytest.approx(68.26, rel=1e-3)
    alpha = nu * K_AIR / D_R
    assert alpha == pytest.approx(70.68, rel=1e-3)
    psi = fin_psi(D_R, D_F)
    assert psi == pytest.approx(20.42e-3, rel=1e-3)
    # The book's text says k_fin = 205 but it calculates with 206: 206 reproduces 0.839, while 205 gives 0.8416.
    x = math.sqrt(2.0 * 70.68 / (W * 206.0)) * psi
    assert x == pytest.approx(0.839, rel=1e-3)
    eta = fin_efficiency(70.68, W, 206.0, psi)
    assert eta == pytest.approx(0.817, rel=1e-3)
    A_F, A_W, A_T = area_fins(N_91, L_91, D_F, D_R, W, S), area_between_fins(N_91, L_91, D_R, W, S), area_bare(N_91, L_91, D_R)
    alpha_r = alpha_on_bare_area(70.68, eta, A_F, A_W, A_T)
    assert alpha_r == pytest.approx(1374.0, rel=1e-3)
    assert alpha_r / (A_F + A_W) * A_T == pytest.approx(58.2, rel=1e-3)


def test_ex91_overall_and_water_outlet():
    m_w = 5.0 / (20 * math.pi / 4.0 * D_I**2)
    v_w = m_w / 995.0
    assert v_w == pytest.approx(0.925, rel=BOOK)
    re_w = m_w * D_I / 0.769e-3
    assert re_w == pytest.approx(2.226e4, rel=BOOK)
    nu_w = 0.023 * re_w**0.8 * 5.20**0.4  # HSB Eq. 2.300 as used in Ex. 9.1
    assert nu_w == pytest.approx(133.7, rel=BOOK)
    alpha_w = nu_w * 0.620 / D_I
    assert alpha_w == pytest.approx(4457.0, rel=BOOK)
    U = overall_u_outer(1374.0, 4457.0, D_R, D_I, 60.0)
    assert U == pytest.approx(909.0, rel=BOOK)
    ntu = 19.15 * 909.0 / (21.55 * 1007.0)
    assert ntu == pytest.approx(0.802, rel=BOOK)
    R = 21.55 * 1007.0 / (5.0 * 4178.0)
    assert R == pytest.approx(1.039, rel=BOOK)
    # The book reads Fig. 9.17. HEDH 3.8.5 Eq. (10) for four or more passes reproduces the reading: P = 0.44.
    e = math.exp(-(1.0 - R) * ntu)
    P = (1.0 - e) / (1.0 - R * e)
    assert P == pytest.approx(0.44, rel=5e-3)  # book reads a chart to 2 s.f.
    Q = 21.55 * 1007.0 * 25.0 * P
    assert Q == pytest.approx(2.394e5, rel=BOOK)
    assert 40.0 - Q / (5.0 * 4178.0) == pytest.approx(28.5, abs=0.1)
    assert 15.0 + Q / (21.55 * 1007.0) == pytest.approx(26.0, abs=0.1)


def test_ex91_through_the_module_evaluation():
    """Same air side through _evaluate at L = 3 m. Pr = 0.7083 from the stated properties, not the book's 0.707, hence Nu 0.07% high."""
    i = example()
    ev = _evaluate(i, _duty(i), 20, 3.0, METHOD_BK, "kosky_staub")
    assert ev.A == pytest.approx(452.2, rel=BOOK)
    assert ev.S_min == pytest.approx(1.760, rel=BOOK)
    assert ev.M_a == pytest.approx(21.56, rel=BOOK)
    assert ev.nu == pytest.approx(68.26, rel=BOOK)
    assert ev.alpha_a == pytest.approx(70.68, rel=BOOK)
    assert ev.eta == pytest.approx(0.817, rel=BOOK)
    assert ev.alpha_ar == pytest.approx(1374.0, rel=BOOK)


# --------------------------------------------------------------------- BENCHMARK: HSB Example 2.11, air-side pressure drop, p. 90-93
def test_ex211_finned_bank():
    D_r, D_f, h, w, s, p1, p2, L = 16.4e-3, 24.6e-3, 4.1e-3, 1e-3, 2e-3, 31.3e-3, 27.1e-3, 0.5
    S_min, plane = min_flow_area(8.5, L, p1, p2, D_r, h, w, s)
    assert plane == "transverse"
    assert S_min == pytest.approx(5.172e-2, rel=BOOK)
    v_max = 0.914 / (1.217 * S_min)
    assert v_max == pytest.approx(14.52, rel=BOOK)
    assert 0.5 * 1.217 * 14.52**2 == pytest.approx(128.3, rel=BOOK)
    assert area_between_fins(32, L, D_r, w, s) == pytest.approx(0.550, rel=BOOK)
    assert area_bare(32, L, D_r) == pytest.approx(0.824, rel=BOOK)
    # SOURCE ERROR, not a pass: the book prints D_r w for D_f w in A_F and then 2.542 m2, total 3.092 m2. Eq. 2.225 and 2.224 give 3.229 and 3.778.
    assert area_fins(32, L, D_f, D_r, w, s) == pytest.approx(3.229, rel=BOOK)
    assert area_ratio(D_f, D_r, w, s) == pytest.approx(3.778 / 0.824, rel=BOOK)
    # The book prints 1.6e4 in the substitution; its 100.9 needs its own Re = 1.61e4 (1.6e4 gives 100.5).
    assert esdu_nusselt(1.61e4, s, h, p1, p2, 0.71) == pytest.approx(100.9, rel=BOOK)
    assert fin_psi(D_r, D_f) == pytest.approx(4.68e-3, rel=BOOK)
    assert fin_efficiency(127.3, w, 15.0, fin_psi(D_r, D_f)) == pytest.approx(0.89, rel=5e-3)  # book 2 s.f.
    sigma = (p1 - D_r - 2.0 * h * w / (w + s)) / p1
    assert sigma == pytest.approx(0.39, rel=5e-3)  # book 2 s.f.
    assert acceleration_constant(sigma) == pytest.approx(1.15, rel=BOOK)
    # Book's note: Eq. 2.252 gives 0.50 (2 s.f.), and that rests on its A/A_T = 3.092/0.824. With the correct ratio K_f = 0.563.
    assert friction_constant_high_fin(1.61e4, 3.092 / 0.824, p1, p2, D_r) == pytest.approx(0.50, abs=0.01)
    assert friction_constant_high_fin(1.61e4, area_ratio(D_f, D_r, w, s), p1, p2, D_r) == pytest.approx(0.563, abs=0.002)
    # Book total, with its low-fin K_f = 0.52 and the rounded head.
    assert (1.15 + 4 * 0.52) * 128.3 == pytest.approx(414.0, rel=BOOK)
    assert air_bundle_pressure_drop(1.15, 4, 0.52, 1.217, 14.52) == pytest.approx(414.0, rel=BOOK)


# --------------------------------------------------------------------- BENCHMARK: HSB Examples 19.1 and 19.2, p. 681-683
def test_ex191_estimate_is_the_first_pass():
    r = run()
    assert step(r, "Q") == pytest.approx(269.8e3, rel=BOOK)
    assert step(r, "dT_LM") == pytest.approx(32.93, rel=1e-3)
    assert step(r, "A_est") == pytest.approx(16.38, rel=BOOK)
    assert step(r, "a_T") == pytest.approx(6.38, rel=BOOK)
    assert step(r, "L_est") == pytest.approx(2.57, rel=BOOK)
    assert r.convergence[0].rows[0][1] == pytest.approx(2.57, rel=BOOK)


def test_ex192_tube_side_chain():
    r = run()
    assert r.get("mass_flux") == pytest.approx(46.0, rel=BOOK)
    assert step(r, "Re_LO") == pytest.approx(5381.0, rel=BOOK)
    assert step(r, "Pr_L") == pytest.approx(4.09, rel=BOOK)
    assert step(r, "alpha_LO") == pytest.approx(195.7, rel=BOOK)
    assert step(r, "alpha_in") == pytest.approx(1729.0, rel=BOOK)
    assert step(r, "alpha_out") == pytest.approx(460.0, rel=BOOK)
    assert r.get("alpha_i") == pytest.approx(1095.0, rel=BOOK)


def test_ex192_result():
    r = run()
    assert r.status is Status.OK_WITH_WARNINGS
    assert r.get("alpha_a_r") == pytest.approx(1374.0, rel=BOOK)
    assert r.get("U") == pytest.approx(490.0, rel=BOOK)
    assert r.get("A_required") == pytest.approx(16.74, rel=BOOK)  # computed 16.72: the book's last step uses 32.90 where its Ex. 19.1 has 32.93
    assert r.get("L") == pytest.approx(2.62, rel=BOOK)
    assert r.get("n_t") == 20 and r.get("N") == 80
    assert r.get("width") == pytest.approx(1.237, rel=BOOK)
    assert r.convergence[0].converged and r.convergence[0].iterations == 2


def test_ex192_air_balance_disagreement_is_reported_not_hidden():
    """Not a benchmark: the book fixes both velocity and outlet. The module keeps its L and reports what closing the balance does."""
    r = run()
    assert check(r, "Air heat balance closes").state is CheckState.FAIL
    assert any("does not close" in w for w in r.warnings)
    assert r.get("m_air") > 1.5 * r.get("m_air_balance")
    assert r.get("T_air_out_balance") < 30.0  # 5 m/s cannot heat the air to the stated 40 degC
    # Both ways of closing it agree with an independent run of those specifications.
    assert r.get("L_outlet_kept") == pytest.approx(run(inputs=_bk(SPEC_OUTLET, T_out=40.0)).get("L"), rel=1e-9)
    assert r.get("L_velocity_kept") == pytest.approx(run(inputs=_bk(SPEC_VELOCITY, v_approach=5.0)).get("L"), rel=1e-9)
    assert r.get("L_velocity_kept") < r.get("L") < r.get("L_outlet_kept")


def _bk(spec, **air):
    i = consistent(spec, **air)
    i.tube_side.method = METHOD_BK
    return i


# --------------------------------------------------------------------- BENCHMARK: HEDH 3.8.9 worked problem, p. 3.8.9-3, -4
def test_hedh_fan_and_noise():
    assert fan_power(68.0, 110.0, 30.0, 0.65) == pytest.approx(14.7e3, abs=60.0)  # 14.65 kW exact; the book prints 14.7
    # Its 30 Pa dynamic pressure: 68 m3/s over the annulus at the air outlet, 25 + 16.8 degC. The full disc would give 23 Pa, so the annulus is the right reading.
    rho = air_table(41.8).rho
    assert fan_dynamic_pressure(68.0, fan_section(3.65, 1.2), rho) == pytest.approx(30.0, abs=1.0)
    assert fan_dynamic_pressure(68.0, fan_section(3.65), rho) < 25.0
    # Tip speed: pi D rpm/60 = 55 m/s as the book uses.
    assert math.pi * 3.65 * 288.0 / 60.0 == pytest.approx(55.0, rel=BOOK)
    l_w = sound_power_level(44.0, 55.0, 110.0, 68.0, 3.65)
    assert l_w == pytest.approx(102.0, abs=0.5)  # book to whole dB
    l_4 = total_sound_power_level(102.0, 4)
    assert l_4 == pytest.approx(108.0, abs=0.5)
    assert sound_pressure_level(108.0, 300.0) == pytest.approx(50.0, abs=0.6)  # book rounds 10 log(2 pi 300^2) = 57.5


def test_noise_falls_with_distance():
    """HSB Eq. 19.3 prints a plus; HEDH Eq. (2) and its worked problem subtract, as the level must."""
    assert sound_pressure_level(100.0, 100.0) < sound_pressure_level(100.0, 10.0) < 100.0


# ------------------------------------------------------------- not benchmarks: identities, closures, limits, trends
def test_boyko_kruzhilin_end_mean_is_the_area_mean():
    """Analytic: over x linear in area, the harmonic mean of alpha_LO (1 + k x)^(1/2) is alpha_LO (a + b)/2. So HSB's arithmetic mean is exact."""
    i = example()
    fluid = i.fluid
    harmonic = mean_coefficient_over_quality(lambda x: boyko_kruzhilin_local(fluid, 46.0, x, D_I), 0.85, 0.05, 4000, "area")
    assert harmonic == pytest.approx(boyko_kruzhilin_mean(fluid, 46.0, D_I, 0.85, 0.05), rel=1e-5)


def test_outlet_spec_closes_the_balance_and_ntu():
    r = run(inputs=consistent(T_out=40.0))
    assert r.convergence[0].iterations > 2  # the velocity moves with L, so it really iterates
    assert check(r, "Air heat balance closes").state is CheckState.PASS
    assert check(r, "Effectiveness closure").state is CheckState.PASS
    assert r.get("m_air") == pytest.approx(269840.0 / (1007.0 * 25.0))
    P = 25.0 / 47.0
    assert r.get("NTU_air") == pytest.approx(-math.log(1.0 - P), rel=2e-3)
    assert r.get("A_bare") == pytest.approx(r.get("A_required"), rel=2e-3)
    assert r.get("v_0") == pytest.approx(r.get("m_air") / (RHO * r.get("S_0")))


def test_mass_flow_spec_matches_outlet_spec():
    a = run(inputs=consistent(T_out=40.0))
    b = run(inputs=consistent(SPEC_FLOW, m_dot=a.get("m_air")))
    assert b.get("T_air_out") == pytest.approx(40.0)
    assert b.get("L") == pytest.approx(a.get("L"), rel=1e-3)


def test_velocity_spec_keeps_the_velocity_and_closes():
    r = run(inputs=consistent(SPEC_VELOCITY, v_approach=3.0))
    assert r.get("v_0") == pytest.approx(3.0)
    Q = r.get("m_air") * CP * (r.get("T_air_out") - 15.0)
    assert Q == pytest.approx(269840.0)
    assert r.get("A_bare") == pytest.approx(r.get("A_required"), rel=2e-3)


def test_air_side_controls_length_trend():
    slow = run(inputs=consistent(SPEC_VELOCITY, v_approach=2.5))
    fast = run(inputs=consistent(SPEC_VELOCITY, v_approach=4.0))
    assert fast.get("U") > slow.get("U")
    assert fast.get("dp_air_bundle") > slow.get("dp_air_bundle")
    assert fast.get("L") < slow.get("L")


def test_tube_count_sizing_finds_smallest_count():
    i = consistent(T_out=40.0)
    i.solver.sizing = SIZE_TUBES
    i.geometry.L = 3.0
    i.geometry.n_t = None
    r = run(inputs=i)
    n = r.get("n_t")
    assert r.get("A_bare") >= r.get("A_required")
    assert check(r, "Area provided meets the area required").state is CheckState.PASS
    below = _evaluate(i, _duty(i), n - 1, 3.0, METHOD_VERTICAL, "kosky_staub")
    assert below.A_T < below.A_req


def test_tube_count_agrees_with_length_mode():
    i = consistent(T_out=40.0)
    i.solver.sizing = SIZE_TUBES
    i.geometry.L = 3.0
    i.geometry.n_t = None
    n = run(inputs=i).get("n_t")
    j = consistent(T_out=40.0)
    j.geometry.n_t = n
    assert run(inputs=j).get("L") <= 3.0
    j.geometry.n_t = n - 1
    assert run(inputs=j).get("L") > 3.0


def test_gravity_term_raises_the_coefficient_at_low_mass_flux():
    """What vertical downflow adds. Both are computed, each cross-checks the other, and Ex. 19.2's Boyko-Kruzhilin runs well below Eq. 16.56 at m = 46."""
    book = run()
    vert = run(inputs=consistent(T_out=40.0))
    assert book.get("alpha_i_cross_check") > 1.5 * book.get("alpha_i")
    assert vert.get("alpha_i_cross_check") == pytest.approx(book.get("alpha_i"))
    assert vert.get("L") < vert.get("L_cross_check")


def test_inclination_lowers_gravity_term_and_static_head():
    def at(angle, method=METHOD_VERTICAL):
        i = consistent(T_out=40.0)
        i.tube_side.tube_angle = angle
        i.tube_side.method = method
        return run(inputs=i)

    v, tilted = at(90.0), at(30.0)
    assert v.get("alpha_i") > tilted.get("alpha_i")
    assert v.get("sin_alpha") == pytest.approx(-1.0)
    assert tilted.get("sin_alpha") == pytest.approx(-0.5)
    assert abs(v.get("dp_tube_gravity")) > abs(tilted.get("dp_tube_gravity"))
    # Boyko-Kruzhilin has no gravity term, so its length does not depend on the angle.
    assert at(90.0, METHOD_BK).get("L") == pytest.approx(at(10.0, METHOD_BK).get("L"))
    flat = at(0.0, METHOD_BK)
    assert flat.get("dp_tube_gravity") == pytest.approx(0.0, abs=1e-9)


def test_static_head_is_a_gain_in_downflow_and_sum_closes():
    r = run()
    assert r.get("dp_tube_gravity") < 0.0 and r.get("dp_tube_acceleration") < 0.0 and r.get("dp_tube_friction") > 0.0
    assert r.get("dp_tube") == pytest.approx(r.get("dp_tube_friction") + r.get("dp_tube_acceleration") + r.get("dp_tube_gravity"))
    # Hand bound: the head cannot exceed a full liquid column or fall below a full vapour column.
    L = r.get("L")
    assert 6.36 * 9.81 * L < -r.get("dp_tube_gravity") < 582.9 * 9.81 * L


def test_air_bundle_drop_by_hand():
    r = run()
    i = example()
    sigma = r.get("S_min") / r.get("S_0")
    k_f = friction_constant_high_fin(r.get("Re_air"), r.get("area_ratio"), P1, P2, D_R)
    hand = (1.0 + sigma**2 + 4.0 * k_f) * 0.5 * RHO * r.get("v_max") ** 2
    assert r.get("dp_air_bundle") == pytest.approx(hand)
    assert check(r, "Static pressure against conventional fans").state is CheckState.FAIL  # 5 m/s, above the book's own 2-4 m/s
    assert i.fan.dp_other == 0.0


def test_diagonal_plane_governs_only_when_the_rows_are_close():
    """The printed test of Eq. 2.229 is reversed: at Ex. 9.1's pitches it sends the flow to the diagonal, and its own working does not."""
    blockage = D_R + 2.0 * W * 15.9e-3 / (W + S)
    p3 = diagonal_pitch(P1, P2)
    assert P1 < 2.0 * p3 - blockage  # the printed condition holds at Ex. 9.1
    assert min_flow_area(20, 3.0, P1, P2, D_R, 15.9e-3, W, S)[1] == "transverse"  # and yet Ex. 9.1 uses this plane, the smaller one
    wide = 4.0 * D_F
    near = D_F * 0.87
    area, plane = min_flow_area(20, 3.0, wide, near, D_R, 15.9e-3, W, S)
    assert plane == "diagonal"
    assert area == pytest.approx(2.0 * 20 * 3.0 * (diagonal_pitch(wide, near) - blockage))


@pytest.mark.parametrize("wall,flow,face", [(WALL_EX91, 20.0, 20.5), (WALL_CORBELS, 20.5, 20.5), (WALL_NONE, 21.0, 21.0)])
def test_wall_treatment_counts(wall, flow, face):
    i = example()
    i.geometry.wall_treatment = wall
    ev = _evaluate(i, _duty(i), 20, 3.0, METHOD_BK, "kosky_staub")
    assert (ev.n_flow, ev.n_face) == (flow, face)


def test_row_factor_and_duty_split():
    assert [row_factor(n) for n in (1, 2, 3, 4, 8)] == [0.76, 0.84, 0.92, 1.0, 1.0]
    fractions = row_duty_fractions(0.8, 4)
    assert sum(fractions) == pytest.approx(1.0)
    assert fractions == sorted(fractions, reverse=True)
    assert fractions[0] / fractions[-1] == pytest.approx(math.exp(0.8 * 3 / 4))
    assert row_duty_fractions(1e-6, 4)[0] == pytest.approx(0.25, rel=1e-5)


def test_rows_below_four_lower_the_coefficient():
    i = consistent(T_out=40.0)
    i.geometry.n_rows = 3
    i.geometry.n_t = 27
    r = run(inputs=i)
    assert step(r, "Nu") == pytest.approx(esdu_nusselt(r.get("Re_air"), S, 15.9e-3, P1, P2, CP * MU / K_AIR, 0.92))
    assert check(r, "Number of rows").state is CheckState.FAIL


def test_fouling_and_contact_raise_the_area():
    clean = run(inputs=consistent(T_out=40.0))
    dirty = consistent(T_out=40.0)
    dirty.geometry.R_f_tube = 2e-4
    dirty.geometry.R_contact = 1e-4
    r = run(inputs=dirty)
    assert r.get("U") < clean.get("U")
    assert r.get("L") > clean.get("L")


def test_suction_correction_raises_the_coefficient():
    base = run(inputs=consistent(T_out=40.0)).get("alpha_i")
    i = consistent(T_out=40.0)
    i.tube_side.suction_correction = "yes"
    assert run(inputs=i).get("alpha_i") > base


@pytest.mark.parametrize("shear", ["boyko_kruzhilin", "butterworth", "chen"])
def test_other_shear_terms_run(shear):
    i = consistent(T_out=40.0)
    i.tube_side.shear_term = shear
    r = run(inputs=i)
    assert r.status in (Status.OK, Status.OK_WITH_WARNINGS)
    assert r.get("L") > 0.0


# --------------------------------------------------------------------------------- fan, natural draft, checks
def test_induced_draft_moves_more_volume():
    forced = run(inputs=consistent(T_out=40.0))
    i = consistent(T_out=40.0)
    i.fan.draft = "induced"
    induced = run(inputs=i)
    assert induced.get("V_air") > forced.get("V_air")
    assert induced.get("fan_power") > forced.get("fan_power")  # HSB Sec. 9.3.1.1, advantage (a) of forced draft
    assert forced.get("V_air") == pytest.approx(forced.get("m_air") / air_table(15.0).rho)
    assert induced.get("V_air") == pytest.approx(induced.get("m_air") / air_table(40.0).rho)


def test_fan_count_from_the_area_rule_and_noise_through_design():
    i = consistent(T_out=40.0)
    i.fan.fan_diameter = 1.0
    i.fan.hub_diameter = 0.3
    i.fan.tip_speed = 40.0
    i.fan.noise_distance = 50.0
    r = run(inputs=i)
    S_v = math.pi / 4.0
    n = r.get("n_fans")
    assert r.get("S_0") / (n * S_v) < 2.6
    assert r.get("S_0") / ((n - 1) * S_v) > 2.6 if n > 1 else True
    assert r.get("dp_dynamic") == pytest.approx(fan_dynamic_pressure(r.get("V_air") / n, fan_section(1.0, 0.3), air_table(15.0).rho))
    assert r.get("l_n") == pytest.approx(total_sound_power_level(r.get("l_w"), n))
    assert r.get("l_p") == pytest.approx(sound_pressure_level(r.get("l_n"), 50.0))
    assert check(r, "Noise").state is CheckState.NOT_ASSESSED


def test_no_fan_diameter_uses_the_source_range():
    r = run()
    assert r.get("dp_dynamic") == pytest.approx(50.0)
    assert r.get("n_fans") == 1
    assert check(r, "Face area to fan section").state is CheckState.PASS
    assert "by construction" in check(r, "Face area to fan section").detail
    assert r.get("fan_power") == pytest.approx(r.get("V_air") * (r.get("dp_air_static") + 50.0) / 0.65)


def test_noise_reports_a_number_from_the_derived_diameter_and_mid_tip_speed():
    # Hand check: D from S/S_v = 2.2, the middle of 1.8 to 2.6; w = 48.5, the middle of 46 to 51 m/s (HEDH 4.4.1-6).
    r = run()
    D = math.sqrt(4.0 * r.get("S_0") / (2.2 * math.pi))
    assert r.get("fan_diameter") == pytest.approx(D)
    assert r.get("tip_speed") == pytest.approx(48.5)
    v_f = r.get("V_air")
    assert r.get("l_w") == pytest.approx(sound_power_level(44.0, 48.5, r.get("dp_air_static"), v_f, D))
    assert r.get("l_1m") == pytest.approx(r.get("l_w") - 10.0 * math.log10(2.0 * math.pi))
    assert "HEDH 4.4.1-6" in check(r, "Noise").detail


def test_noise_table_hand_check_and_hedh_worked_problem():
    assert noise_table_at_1m(60.0) == 85.0
    assert noise_table_at_1m(30.0) == 75.0
    assert noise_table_at_1m(70.0) == 85.0
    assert noise_table_at_1m(20.0) == 75.0
    # HEDH 3.8.9's worked problem: w = 55 m/s, L_w = 102 dB(A). Table 8 between 45 m/s (81) and 60 m/s (85): 83.7. Eq. (2) at 1 m: 94.
    assert noise_table_at_1m(55.0) == pytest.approx(81.0 + 4.0 * 10.0 / 15.0)
    assert sound_pressure_level(sound_power_level(44.0, 55.0, 110.0, 68.0, 3.65), 1.0) == pytest.approx(94.1, abs=0.1)


def test_fan_diameter_and_tip_speed_checks():
    i = consistent(T_out=40.0)
    i.fan.fan_diameter, i.fan.tip_speed = 1.0, 70.0
    r = run(inputs=i)
    assert check(r, "Fan diameter").state is CheckState.FAIL
    assert check(r, "Fan tip speed").state is CheckState.FAIL
    i.fan.fan_diameter, i.fan.tip_speed = 3.0, 50.0
    r = run(inputs=i)
    assert check(r, "Fan diameter").state is CheckState.PASS
    assert check(r, "Fan tip speed").state is CheckState.PASS
    assert r.get("fan_diameter") == 3.0 and r.get("tip_speed") == 50.0


def test_large_face_area_splits_the_fans_to_stay_within_the_diameter_range():
    i = consistent(T_out=40.0)
    i.stream.m_dot = 40.0
    i.geometry.n_t = 400
    r = run(inputs=i)
    n = r.get("n_fans")
    assert n > 1
    assert r.get("fan_diameter") <= 4.1 + 1e-9
    assert math.sqrt(4.0 * r.get("S_0") / ((n - 1) * 2.2 * math.pi)) > 4.1


def test_motor_margin_check():
    i = consistent(T_out=40.0)
    assert check(run(inputs=i), "Motor margin").state is CheckState.PASS
    power = run(inputs=i).get("fan_power")
    i.fan.motor_power = 1.20 * power
    assert check(run(inputs=i), "Motor margin").state is CheckState.FAIL  # fixed blades: 25-30%
    i.fan.blade_control = "seasonal"
    assert check(run(inputs=i), "Motor margin").state is CheckState.PASS  # reset blades: 15%
    i.fan.blade_control = "autovariable"
    assert check(run(inputs=i), "Motor margin").state is CheckState.NOT_ASSESSED


# --------------------------------------------------------------------------------- regime criteria, HEDH 3.4.6 B
def test_shear_parameter_hand_check():
    # Re_l = 1000: s+ = (1000/2)^(1/2) = 22.36. tau = 1, rho_l = 582.9, eta_l = 159e-6: s = 22.36 (159e-6)/(582.9)^(1/2) = 1.4726e-4 m.
    s = hedh_film_thickness(1000.0, 1.0, 582.9, 159e-6)
    assert s == pytest.approx(math.sqrt(500.0) * 159e-6 / math.sqrt(582.9))
    # Above Re_l = 1250: s+ = 0.0504 Re^(7/8).
    assert hedh_film_thickness(2000.0, 1.0, 582.9, 159e-6) == pytest.approx(0.0504 * 2000.0**0.875 * 159e-6 / math.sqrt(582.9))
    fluid = example().fluid
    assert shear_parameter(fluid, 1000.0, 1.0) == pytest.approx(1.0 / ((fluid.rho_l - fluid.rho_v) * 9.81 * hedh_film_thickness(1000.0, 1.0, fluid.rho_l, fluid.mu_l)), rel=1e-3)


@pytest.mark.parametrize("value,limits,expected", [
    (5.0, (5.0, 29.0), "no-shear"), (29.0, (5.0, 29.0), "shear"), (10.0, (5.0, 29.0), "transition"),
    (0.5, (0.5, 1.5), "no-shear"), (1.5, (0.5, 1.5), "shear"), (1.0, (0.5, 1.5), "transition"),
])
def test_regime_limits(value, limits, expected):
    assert regime_of(value, limits) == expected


def test_regime_rule_takes_each_limit_and_interpolates():
    assert regime_coefficient(2.0, 800.0, 1500.0) == 800.0
    assert regime_coefficient(0.3, 800.0, 1500.0) == 1500.0
    assert regime_coefficient(1.0, 800.0, 1500.0) == pytest.approx(1150.0)


def regime_rows(result):
    return next(c for c in result.convergence if c.name.startswith("Tube-side regime")).rows


def test_regime_table_along_the_tube():
    r = run()  # Ex. 19.2: m = 46 kg/m2.s
    fluid = example().fluid
    c = 46.0 / math.sqrt(9.81 * 18.6e-3 * fluid.rho_v * (fluid.rho_l - fluid.rho_v))
    assert r.get("x_shear_limit") == pytest.approx(1.5 / c, rel=1e-3)
    assert r.get("x_no_shear_limit") == pytest.approx(0.5 / c, rel=1e-3)
    rows = regime_rows(r)
    assert len(rows) == 10
    for row in rows:
        assert row[1] == pytest.approx(c * row[0], rel=1e-3)
        assert row[3] == regime_of(row[1], (0.5, 1.5))
        assert row[4] == regime_of(row[2], (5.0, 29.0))
        assert row[3] in row[5]
    assert rows[-1][3] == "no-shear" and "does not apply" in rows[-1][5]
    assert r.get("alpha_i_regime_rule") > 0


def test_method_choice_is_reported_against_the_regime():
    # Ex. 19.2's Boyko-Kruzhilin runs into no-shear stations at the low-quality end; Eq. 16.56 holds through all of them.
    assert check(run(), "Tube-side method fits the regime along the tube").state is CheckState.FAIL
    r = run(inputs=consistent(T_out=40.0))
    assert check(r, "Tube-side method fits the regime along the tube").state is CheckState.PASS
    assert any("HEDH 3.4.6 B" in a and "Eq. 16.56" in a for a in r.design_basis)
    assert r.get("L_cross_check") > r.get("L")  # the method spread stays a reported fact


def test_regime_table_without_surface_tension_keeps_j_g():
    i = example()
    i.fluid = replace(i.fluid, sigma=None)
    r = run(inputs=i)
    rows = regime_rows(r)
    assert rows[0][1] > rows[-1][1] and rows[0][2] is None


# --------------------------------------------------------------------------------- vacuum, HEDH 3.8.9 B
def test_sonic_limit_hand_check():
    # c = (gamma P/rho)^(1/2). 3000 Pa, 6.36 kg/m3: 21.72 m/s. v_in = 46 (0.85)/6.36 = 6.15 m/s, below 50% of c.
    assert sonic_velocity(3000.0, 6.36) == pytest.approx(21.72, abs=0.01)
    assert sonic_velocity(3000.0, 6.36, 1.4) > sonic_velocity(3000.0, 6.36)
    i = example()
    i.stream.pressure = 3000.0
    r = run(inputs=i)
    assert r.get("v_vapour_inlet") == pytest.approx(step(r, "m") * 0.85 / i.fluid.rho_v)
    assert r.get("sonic_velocity") == pytest.approx(sonic_velocity(3000.0, i.fluid.rho_v))
    assert check(r, "Sonic limit at the tube inlet").state is CheckState.PASS
    assert check(r, "Vacuum limit").state is CheckState.FAIL
    i.stream.pressure = 300.0
    r = run(inputs=i)
    assert check(r, "Sonic limit at the tube inlet").state is CheckState.FAIL
    assert any("short tubes" in w for w in r.warnings)
    i.stream.pressure = 1e5
    assert check(run(inputs=i), "Vacuum limit").state is CheckState.PASS


def test_vacuum_checks_need_the_pressure():
    r = run()
    assert check(r, "Vacuum limit").state is CheckState.NOT_ASSESSED
    assert not any(c.name == "Sonic limit at the tube inlet" for c in r.checks)


def test_api_661_values_are_guidance_citing_the_edition():
    detail = check(run(), "API 661 and other air-cooled exchanger standards").detail
    for text in ("2nd ed.", "0.5%", "19 mm", "9.5 mm", "1.9 mm", "1.5 mm", "2.7 mm", "6.4 mm", "no compliance is claimed"):
        assert text in detail


def test_result_wording_is_decisive():
    banned = ("supplied", "project drive", "replace or confirm", "when it arrives", "not available", "unconfirmed", "should confirm")
    for inputs in (example(), consistent(T_out=40.0)):
        r = run(inputs=inputs)
        text = " ".join(r.design_basis + r.scope_limits + r.warnings + [c.detail for c in r.checks]).lower()
        assert not [b for b in banned if b in text]


def test_natural_draft_height():
    i = consistent(T_out=40.0)
    i.fan.draft = "natural"
    r = run(inputs=i)
    delta = air_table(15.0).rho - air_table(40.0).rho
    assert r.get("H_eff") == pytest.approx(natural_draft_height(r.get("dp_air_static"), r.get("dp_air_static"), delta))
    assert check(r, "Natural draft economics").state is CheckState.FAIL  # far below 100 MW
    with pytest.raises(KeyError):
        r.get("fan_power")


def test_winter_check():
    i = consistent(T_out=40.0)
    i.winter.T_air_min, i.winter.T_freeze = -20.0, 0.0
    assert check(run(inputs=i), "Freezing in winter").state is CheckState.FAIL
    i.winter.T_air_min = 5.0
    assert check(run(inputs=i), "Freezing in winter").state is CheckState.PASS
    assert check(run(), "Freezing in winter").state is CheckState.NOT_ASSESSED
    assert "Lecture 23" in check(run(), "Freezing in winter").source


def test_warnings_and_contact_check():
    i = consistent(T_out=40.0)
    i.stream.noncondensables = "yes"
    i.stream.pressure = 3000.0
    r = run(inputs=i)
    assert any("Non-condensables" in w for w in r.warnings)
    assert any("30 torr" in w for w in r.warnings)
    hot = consistent(T_out=40.0)
    hot.fluid = replace(hot.fluid, T_sat=160.0)
    assert check(run(inputs=hot), "Fin contact at temperature").state is CheckState.FAIL


def test_allowable_pressure_drops():
    i = consistent(T_out=40.0)
    i.limits.dp_allow_coolant = 50.0
    i.limits.dp_allow_condensing = 1e5
    r = run(inputs=i)
    assert check(r, "Air-side allowable pressure drop").state is CheckState.FAIL
    assert check(r, "Tube-side allowable pressure drop").state is CheckState.PASS


def test_codes_are_not_assessed_and_no_compliance_claimed():
    r = run()
    for name in ("API 661 and other air-cooled exchanger standards", "Header box, tube and pressure design", "Vent and non-condensable placement", "Condensate backup"):
        assert check(r, name).state is CheckState.NOT_ASSESSED
    assert "no compliance is claimed" in check(r, "API 661 and other air-cooled exchanger standards").detail


def test_traceability():
    r = run(inputs=consistent(T_out=40.0))
    assert all(s.formula and s.substitution for s in r.steps)
    assert any("Why the length iterates" in a for a in r.design_basis)
    assert any("F = 1" in a for a in r.design_basis)
    assert any("BENCHMARKS EXTEND" in g for g in r.scope_limits)
    assert {Category.THERMAL, Category.HYDRAULIC, Category.GEOMETRIC, Category.MECHANICAL} <= {s.category for s in r.steps}
    assert r.convergence[0].columns and r.convergence[0].rows


def test_report_exports():
    r = run()
    assert to_excel(r, example())
    assert to_pdf(r, example())


# --------------------------------------------------------------------------------- validation
def issues(mutate):
    result = run(mutate)
    assert result.status is Status.INVALID_INPUT
    return " ".join(i.message for i in result.input_issues)


def test_empty_form_reports_issues():
    r = run_design("air_cooled", MODULE.new_inputs())
    assert r.status is Status.INVALID_INPUT and r.input_issues


def test_missing_required_value():
    def missing(i):
        i.stream.m_dot = None

    assert "Total process mass flow" in issues(missing)


@pytest.mark.parametrize("mutate,text", [
    (lambda i: setattr(i.stream, "x_out", 0.9), "below the inlet quality"),
    (lambda i: setattr(i.air, "T_in", 70.0), "below the saturation"),
    (lambda i: setattr(i.air, "T_out", 70.0), "temperature cross"),
    (lambda i: setattr(i.air, "T_out", None), "outlet temperature is required"),
    (lambda i: setattr(i.air, "v_approach", None), "approach velocity is required"),
    (lambda i: setattr(i.geometry, "D_i", 30e-3), "inside diameter"),
    (lambda i: setattr(i.geometry, "D_f", 20e-3), "fin diameter"),
    (lambda i: setattr(i.geometry, "p_1", 50e-3), "overlap"),
    (lambda i: setattr(i.geometry, "p_2", 20e-3), "overlap"),
    (lambda i: setattr(i.geometry, "n_t", None), "Tubes per row"),
    (lambda i: setattr(i.geometry, "n_t", 20.5), "whole number"),
    (lambda i: setattr(i.fan, "hub_diameter", 0.5), "hub diameter"),
])
def test_invalid_values(mutate, text):
    assert text in issues(mutate)


def test_air_flow_too_small_is_a_cross():
    i = consistent(SPEC_FLOW, m_dot=2.0)
    r = run(inputs=i)
    assert r.status is Status.INVALID_INPUT and "at or above T_sat" in r.input_issues[0].message


def test_vertical_method_needs_sigma_and_a_slope():
    i = consistent(T_out=40.0)
    i.fluid = replace(i.fluid, sigma=None)
    assert "Surface tension" in run(inputs=i).input_issues[0].message
    j = consistent(T_out=40.0)
    j.tube_side.tube_angle = 0.0
    assert "horizontal tube" in run(inputs=j).input_issues[0].message


def test_boyko_kruzhilin_without_sigma_falls_back_as_the_source_says():
    i = example()
    i.fluid = replace(i.fluid, sigma=None)
    r = run(inputs=i)
    assert r.status is Status.OK_WITH_WARNINGS
    assert r.get("L") == pytest.approx(2.62, rel=BOOK)
    assert r.get("dp_tube_friction") > 0.0  # Chisholm, HSB Sec. 10.5
    assert any("cross-check was not computed" in g for g in r.scope_limits)


def test_non_convergence_reports_with_its_history():
    def one_iteration(i):
        i.solver.max_iterations = 1

    result = run(one_iteration)
    assert result.status is Status.NOT_CONVERGED
    assert result.convergence and not result.convergence[0].converged
    assert result.convergence[0].rows
    assert "no convergence" in result.message


def test_relevant_fields_hide_what_does_not_apply():
    new = MODULE.new_inputs()
    shown = MODULE.relevant_fields(new)
    assert "air.T_out" in shown and "air.m_dot" not in shown and "air.v_approach" not in shown
    assert "geometry.n_t" in shown and "geometry.L" not in shown
    new.air.spec = SPEC_BOOK
    new.solver.sizing = SIZE_TUBES
    new.tube_side.method = METHOD_BK
    new.fan.draft = "natural"
    shown = MODULE.relevant_fields(new)
    assert {"air.T_out", "air.v_approach", "geometry.L", "fan.dp_tower"} <= shown
    assert not {"geometry.n_t", "tube_side.shear_term", "fan.efficiency", "fan.tip_speed"} & shown


def test_module_registered():
    assert MODULE.key == "air_cooled"
    assert MODULE.example_inputs().air.spec == SPEC_BOOK
