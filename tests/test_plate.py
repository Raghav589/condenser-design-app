"""Plate condenser module.

No source works a condensing plate example. The BENCHMARK tests are single-phase PHE examples of the coolant-side and
sizing path (HSB Ex. 8.1-8.4, HEDH 3.7.10), run through the module's own functions; the two-phase relations are
benchmarked only where HSB Ch. 10 works numbers. Everything else is a hand check, limit, closure or trend, not a benchmark.
"""
import copy
import math

import pytest

from condenser_app.common.thermal import lmtd, lmtd_counterflow
from condenser_app.modules import MODULES
from condenser_app.modules.plate import (
    MAX_PLATES,
    akers_mass_flux,
    channel_dp,
    channel_mass_flow,
    chisholm_c,
    effectiveness_counterflow,
    end_effect_factor,
    equivalent_diameter,
    friction_factor,
    ntu_counterflow,
    nusselt,
    overall_u,
    phi_l_squared,
    plate_area,
    port_dp,
    property_correction_factors,
    smallest_plate_count,
    void_fraction_lm,
    _duty,
    _evaluate,
    _hydraulics,
    _snap_plates,
)
from condenser_app.reporting import to_excel, to_pdf
from condenser_app.schemas import Category, CheckState, Status
from condenser_app.service import run_design

MODULE = MODULES["plate"]
C_H, Y, N_PR, K_P, Z = 0.4, 0.64, 0.4, 2.78, 0.18  # HSB Ex. 8.1 and 8.2: fits to Fig. 8.8
L, W, B = 1.0, 0.25, 5e-3  # HSB Ex. 8.1, 8.3 plate


def example():
    return copy.deepcopy(MODULE.example_inputs())


def run(mutate=None):
    inputs = example()
    if mutate:
        mutate(inputs)
    return run_design("plate", inputs)


def check(result, name):
    return next(c for c in result.checks if c.name == name)


def h_channel(M, N, rho, mu, k, cp, passes=1, enh=1.0):
    """Lecture 17 rating Steps 1-5 on one stream, with the module's functions."""
    de = equivalent_diameter(B)
    G = channel_mass_flow(M, N, passes) / (B * W)
    re = G * de / mu
    nu = nusselt(re, cp * mu / k, C_H, Y, N_PR)
    return enh * nu * k / de, re, nu, G


# --------------------------------------------------------------------- BENCHMARK: HSB Example 8.1
WATER_55 = dict(rho=985.2, mu=5.05e-4, k=0.648, cp=4183.0)  # the book's Pr = 3.26 agrees with these


def ex81_ntu(N, passes=1):
    h_h, *_ = h_channel(10.0, N, **WATER_55, passes=passes)
    h_c, *_ = h_channel(5.0, N, **WATER_55, passes=passes)
    return overall_u(h_h, h_c) * N * plate_area(L, W) / (5.0 * WATER_55["cp"])


def test_benchmark_hsb_example_8_1_intermediates():
    # G and Re: 0.5%, the book prints 3 s.f. Nu onward: 1%. The book prints Nu = 111, yet 0.4 (3.26)^0.4 (3160)^0.64 with
    # its own rounded numbers is 111.5 (unrounded inputs give 111.65), so its 111 is truncated and carries 0.6% downstream.
    h_h, re, nu, G = h_channel(10.0, 99, **WATER_55)
    assert G / WATER_55["rho"] == pytest.approx(0.162, rel=5e-3)
    assert re == pytest.approx(3.16e3, rel=5e-3)
    assert nu == pytest.approx(111, rel=1e-2)
    assert h_h == pytest.approx(7.19e3, rel=1e-2)
    h_c, *_ = h_channel(5.0, 99, **WATER_55)
    assert h_c == pytest.approx(4.61e3, rel=1e-2)
    assert overall_u(h_h, h_c) == pytest.approx(2.81e3, rel=1e-2)
    assert ex81_ntu(99) == pytest.approx(3.33, rel=1e-2)


@pytest.mark.parametrize(
    "N,ntu,E,T_out",
    # HSB Ex. 8.1 table. The N = 49 and 25 rows are garbled in the text extraction; values read from the page and
    # recalculated from the book's own scaling rule: 2.57, 0.839, 82.1 and 1.99, 0.773, 76.8.
    [(99, 3.33, 0.896, 86.7), (75, 3.01, 0.875, 85.0), (49, 2.57, 0.839, 82.1), (25, 1.99, 0.773, 76.8)],
)
def test_benchmark_hsb_example_8_1_table(N, ntu, E, T_out):
    # 1% on NTU (carries the truncated Nu above); 0.5% on E; 0.25 K on the outlet, which is that same carry (0.3% of 80 K).
    computed = ex81_ntu(N)
    assert computed == pytest.approx(ntu, rel=1e-2)
    e = effectiveness_counterflow(computed, 0.5)
    assert e == pytest.approx(E, rel=5e-3)
    assert 15 + e * 80 == pytest.approx(T_out, abs=0.25)


# --------------------------------------------------------------------- BENCHMARK: HSB Example 8.2
def test_benchmark_hsb_example_8_2_two_pass():
    ntu = ex81_ntu(99, passes=2)
    assert ntu == pytest.approx(5.19, rel=5e-3)
    e = effectiveness_counterflow(ntu, 0.5)
    assert e == pytest.approx(0.961, rel=2e-3)
    assert 15 + e * 80 == pytest.approx(91.9, abs=0.15)


def test_benchmark_hsb_example_8_2_pressure_drops():
    # The book prints 3.36e3 and 23.73e3 'kN/m2'. Its own numbers give N/m2: a unit error in the book, see the doc.
    # 1%: the book rounds u to 0.162 m/s before squaring (the unrounded 0.1624 gives 3386 N/m2, +0.8%).
    rho = WATER_55["rho"]
    drops = []
    for passes in (1, 2):
        _, re, _, G = h_channel(10.0, 99, **WATER_55, passes=passes)
        drops.append(channel_dp(friction_factor(re, K_P, Z), L * passes, equivalent_diameter(B), G, rho))
    assert drops[0] == pytest.approx(3.36e3, rel=1e-2)
    assert drops[1] == pytest.approx(23.73e3, rel=1e-2)  # book scales by 7.06 from the rounded single-pass value
    assert drops[1] / drops[0] == pytest.approx(7.06, rel=2e-3)


# --------------------------------------------------------------------- BENCHMARK: HSB Example 8.3, 8.4
DOWTHERM = dict(rho=1036.0, mu=2.15e-3, k=0.137, cp=1660.0)
WATER_77 = dict(rho=973.7, mu=3.65e-4, k=0.668, cp=4195.0)
M_DOWTHERM = 5.0


def ex83_heat_balance():
    M_h = M_DOWTHERM * DOWTHERM["cp"] * 70 / (WATER_77["cp"] * 35)
    R = (M_DOWTHERM * DOWTHERM["cp"]) / (M_h * WATER_77["cp"])
    return M_h, R, 70 / 80


def ex83_size(passes=1, enh_h=1.0):
    M_h, R, E = ex83_heat_balance()
    ntu = ntu_counterflow(E, R)
    C_c = M_DOWTHERM * DOWTHERM["cp"]

    def required(n):
        h_h, *_ = h_channel(M_h, n, **WATER_77, passes=passes, enh=enh_h)
        h_c, *_ = h_channel(M_DOWTHERM, n, **DOWTHERM, passes=passes, enh=enh_h)
        return ntu * C_c / overall_u(h_h, h_c) / plate_area(L, W)

    N, info = smallest_plate_count(required, start=1, snap=math.ceil, name="Ex 8.3")
    return int(N), info


def ex83_dp(N, passes=1, enh_f=1.0):
    M_h, *_ = ex83_heat_balance()
    _, re, _, G = h_channel(M_h, N, **WATER_77, passes=passes)
    f = friction_factor(re, K_P * enh_f, Z)
    return re, f, channel_dp(f, L * passes, equivalent_diameter(B), G, WATER_77["rho"])


def test_benchmark_hsb_example_8_3_required_ntu():
    M_h, R, E = ex83_heat_balance()
    assert M_h == pytest.approx(3.957, rel=1e-3)
    assert R == pytest.approx(0.5, rel=1e-9)
    assert E == pytest.approx(0.875)
    assert ntu_counterflow(E, R) == pytest.approx(3.008, rel=1e-3)


def test_benchmark_hsb_example_8_3_one_pass_size():
    # Book: 230 plates. Computed: 231. Within 1 plate: with the book's rounded constants n = 231 passages sits on the line
    # (7.0634 against 7.063). With unrounded values it falls just short, and the module requires the NTU to be met.
    N, info = ex83_size()
    assert abs(N - 230) <= 1
    assert info.converged
    re, f, dp = ex83_dp(230)
    assert re == pytest.approx(751, rel=5e-3)
    assert f == pytest.approx(0.844, rel=5e-3)
    assert dp == pytest.approx(130, rel=1e-2)


def test_benchmark_hsb_example_8_3_two_pass_size():
    # Book: 69 plates, 8.05e3 N/m2 scaled from the rounded 130 N/m2; computed directly, hence 1%.
    N, _ = ex83_size(passes=2)
    assert abs(N - 69) <= 1
    assert ex83_dp(69, passes=2)[2] == pytest.approx(8.05e3, rel=1e-2)


@pytest.mark.parametrize("enh_h,enh_f,n_passages,re,f,dp", [(2.0, 10.0, 36, 4.819e3, 6.04, 3.85e4), (1.4, 3.0, 92, 1.886e3, 2.146, 2.08e3)])
def test_benchmark_hsb_example_8_4_chevron_angle(enh_h, enh_f, n_passages, re, f, dp):
    # Enhancement factors are HSB's readings of Figs. 8.9 and 8.10. Within 1 passage, for the reason in the one-pass test:
    # at n = 36 the book's rounded expression gives 3.5321 against 3.532. 1% on values.
    N, _ = ex83_size(enh_h=enh_h)
    assert abs((N + 1) - n_passages) <= 1
    computed_re, computed_f, computed_dp = ex83_dp(n_passages - 1, enh_f=enh_f)
    assert computed_re == pytest.approx(re, rel=1e-2)
    assert computed_f == pytest.approx(f, rel=1e-2)
    assert computed_dp == pytest.approx(dp, rel=1e-2)


# ------------------------------------------------------------------- BENCHMARK: HEDH Sec. 3.7.10
def test_benchmark_hedh_example_1_property_corrections():
    # Book rounds to 0.17 and 1.34. The chart readings (a, alpha, q) are not reproduced; the curves are not digitised.
    psi, phi_dp = property_correction_factors(rho=880.0, cp=1926.0, k=0.125, mu=3e-3)
    assert psi == pytest.approx(0.17, abs=0.005)
    assert phi_dp == pytest.approx(1.34, abs=0.005)


def test_benchmark_hedh_example_2_reciprocal_and_lmtd():
    # alpha_1, alpha_2 are the book's chart readings, entered as inputs.
    assert overall_u(9000.0, 10200.0) == pytest.approx(4780, rel=1e-3)
    T_hot_out = 368.0 - 22.7 * 78.0 / 34.0
    assert lmtd_counterflow(368.0, T_hot_out, 283.0, 361.0) == pytest.approx(16.77, rel=5e-3)  # book 16.77
    assert 4780 * 185 * 0.52 * 16.77 == pytest.approx(7710e3, rel=1e-3)


# ---------------------------------------------- CHART CHECK: Lecture 15 Fig. 2.25, F at C = 0
@pytest.mark.parametrize(
    "n_channels,E,F_chart",
    [
        (4, 0.5, 0.96), (4, 0.9, 0.867), (4, 0.99, 0.78), (4, 0.9999, 0.72),
        # Right-hand edge, n = 18: the scan is tilted, with F = 1 about 0.008 higher there; readings corrected for it.
        (18, 0.9, 0.955), (18, 0.95, 0.935), (18, 0.99, 0.87), (18, 0.999, 0.768), (18, 0.9999, 0.695),
    ],
)
def test_end_effect_factor_reproduces_lecture_15_fig_2_25(n_channels, E, F_chart):
    # Read from the chart, so 0.012: this is a chart check of a derived formula, not a printed benchmark.
    assert end_effect_factor(E, n_channels - 1) == pytest.approx(F_chart, abs=0.012)


def test_end_effect_factor_limits():
    assert end_effect_factor(0.9, 1) == pytest.approx(1.0)  # one plate, one coolant channel
    assert end_effect_factor(0.7, 999) == pytest.approx(1.0, abs=1e-3)  # HSB Eq. 8.26: F = 1 for large N
    assert end_effect_factor(0.7, 999, passes=2) == pytest.approx(1.0, abs=2e-3)  # isothermal: any arrangement


# ---------------------------------------------- BENCHMARK: HSB Ch. 10 worked numbers
def test_benchmark_hsb_example_10_1_multiplier():
    # Book: X = 1.796, C = 20, phi_L^2 = 1 + 11.13 + 0.310 = 12.44 (it rounds; 0.1%).
    assert phi_l_squared(1.796, 20.0) == pytest.approx(12.44, rel=1e-3)


def test_benchmark_hsb_example_10_3_void_fraction():
    # Book: phi_L^2 = 12.44 gives eps_G = 0.716 by Eq. 10.58.
    assert void_fraction_lm(math.sqrt(12.44)) == pytest.approx(0.716, abs=1e-3)


def test_chisholm_table():
    assert [chisholm_c(1e4, 1e4, 150), chisholm_c(10, 1e4, 150), chisholm_c(1e4, 10, 150), chisholm_c(10, 10, 150)] == [20, 12, 10, 5]


# ---------------------------------------------- HAND CHECKS on the plate relations (not benchmarks)
def test_diameter_and_area_bases():
    assert equivalent_diameter(4e-3) == pytest.approx(8e-3)
    assert equivalent_diameter(4e-3, 1.25, "developed") == pytest.approx(6.4e-3)
    assert plate_area(1.0, 0.5) == pytest.approx(0.5)
    assert plate_area(1.0, 0.5, 1.2, "developed") == pytest.approx(0.6)


def test_port_loss_by_hand():
    # G_p = 1 / (pi/4 0.1^2) = 127.32; 1.4 (127.32^2 / 2000) (2) = 22.70
    assert port_dp(1.0, 0.1, 1000.0, 2) == pytest.approx(22.70, abs=0.01)


def test_akers_flux_limits_and_correction_of_eq_20_7():
    p = dict(rho_l=958.3, rho_g=0.597, mu_l=277.53e-6, mu_g=12.55e-6, z=0.18)
    G = 16.0
    # x -> 0: the actual liquid flux. x -> 1: all effective flux from the vapour.
    assert akers_mass_flux(G, 0.0, **p) == pytest.approx(G)
    expected = G * math.sqrt(958.3 / 0.597) * math.sqrt((12.55e-6 / 277.53e-6) ** 0.18)
    assert akers_mass_flux(G, 1.0, **p) == pytest.approx(expected)
    assert akers_mass_flux(G, 0.8, **p) > akers_mass_flux(G, 0.2, **p)
    # With z = 0 it is the classical (1 - x) + x (rho_L/rho_G)^1/2.
    assert akers_mass_flux(G, 0.5, 958.3, 0.597, 1.0, 1.0, 0.0) == pytest.approx(G * (0.5 + 0.5 * math.sqrt(958.3 / 0.597)))
    # Eq. 20.7 as printed, with actual-flow gradients and constant f, reduces to 2 x G, so it vanishes as x -> 0.
    x = 1e-3
    ratio = (x / (1 - x)) ** 2 * (958.3 / 0.597)
    printed = G * (x + (1 - x) * math.sqrt(0.597 / 958.3) * math.sqrt(ratio))
    assert printed == pytest.approx(2 * x * G, rel=1e-9)


def test_allowed_plate_counts():
    assert [_snap_plates(n, 1) for n in (1, 2, 3, 4)] == [1, 3, 3, 5]
    assert [_snap_plates(n, 2) for n in (1, 3, 4, 8)] == [3, 3, 7, 11]
    assert _snap_plates(99, 2) == 99  # HSB Ex. 8.2: 25 channels per pass


# --------------------------------------------------------- THE CONDENSING PATH (no benchmark exists)
def test_example_is_cited_as_illustrative_and_runs():
    result = run()
    assert result.status is Status.OK_WITH_WARNINGS
    assert result.get("n_plates") % 2 == 1
    assert any("No source works a condensing plate example" in g for g in result.source_gaps)
    assert any("proprietary" in a for a in result.assumptions)


def test_every_output_requested_is_present():
    result = run()
    for key in ("Q", "m_dot_coolant", "u_coolant", "Re_coolant", "G_condensing", "h_coolant", "h_condensing", "U", "lmtd", "F", "area", "area_required", "n_plates",
                "dp_coolant_channel", "dp_coolant_port", "dp_coolant_static", "dp_coolant",
                "dp_condensing_friction", "dp_condensing_acceleration", "dp_condensing_gravity", "dp_condensing_port", "dp_condensing"):
        assert isinstance(result.get(key), (int, float)), key


def test_energy_balance_closure():
    i = example()
    result = run_design("plate", i)
    Q = i.vapour.m_dot * i.condensing.h_fg
    assert result.get("Q") == pytest.approx(Q)
    assert result.get("m_dot_coolant") * i.coolant_props.cp * (result.get("T_c_out") - i.coolant.T_in) == pytest.approx(Q)


def test_zones_sum_to_the_duty_area_and_plate_length():
    result = run()
    table = next(c for c in result.convergence if c.name.startswith("Zone table"))
    col = {name: n for n, name in enumerate(table.columns)}
    assert sum(row[col["dA (m2)"]] for row in table.rows) == pytest.approx(result.get("area_required"))
    assert sum(row[col["l (m)"]] for row in table.rows) == pytest.approx(example().geometry.L)
    assert sum(row[col["dP_F (Pa)"]] for row in table.rows) == pytest.approx(result.get("dp_condensing_friction"))
    assert table.rows[0][col["x_start"]] == 1.0 and table.rows[-1][col["x_end"]] == 0.0


def test_constant_u_reduces_to_q_over_u_f_lmtd():
    # Equal phase properties make G_Le = G at every x, so U is constant and the zones must sum to Q/(U F LMTD).
    def equal_phases(i):
        c = i.condensing
        c.rho_v, c.mu_v = c.rho_l, c.mu_l

    result = run(equal_phases)
    i = example()
    equal_phases(i)
    assert result.get("h_condensing_inlet") == pytest.approx(result.get("h_condensing_outlet"))
    u = overall_u(result.get("h_condensing"), result.get("h_coolant"), i.geometry.t, i.geometry.k_plate)
    assert result.get("U") == pytest.approx(u, rel=1e-9)
    assert result.get("area_required") == pytest.approx(result.get("Q") / (u * result.get("F") * result.get("lmtd")), rel=1e-9)


def test_h_by_hand_at_the_first_zone():
    # x = 0.975, steam at 1 atm: see the formula chain in docs/algorithms/plate.md.
    i = example()
    result = run_design("plate", i)
    n = result.get("n_channels_per_side")
    G = i.vapour.m_dot / n / (i.geometry.b * i.geometry.W)
    assert result.get("G_condensing") == pytest.approx(G)
    G_le = G * (0.025 + 0.975 * math.sqrt(958.3 / 0.597) * math.sqrt((12.55 / 277.53) ** 0.18))
    re = G_le * 0.01 / 277.53e-6
    pr = 4220 * 277.53e-6 / 0.679
    h = 0.4 * re**0.64 * pr**0.4 * 0.679 / 0.01
    assert result.get("h_condensing_inlet") == pytest.approx(h, rel=1e-9)


def test_acceleration_term_telescopes_to_the_end_states():
    # Total condensation: dP_a = G^2 (1/rho_L - 1/rho_G) exactly, whatever the void fraction in between. Negative: recovery.
    i = example()
    result = run_design("plate", i)
    G = result.get("G_condensing")
    assert result.get("dp_condensing_acceleration") == pytest.approx(G**2 * (1 / 958.3 - 1 / 0.597), rel=1e-9)
    assert result.get("dp_condensing_acceleration") < 0


def test_gravity_head_is_bounded():
    result = run()
    assert 0 < result.get("dp_condensing_gravity") < 958.3 * 9.81 * 1.0


def test_condensing_total_adds_up():
    r = run()
    total = r.get("dp_condensing_friction") + r.get("dp_condensing_acceleration") - r.get("dp_condensing_gravity") + r.get("dp_condensing_port")
    assert r.get("dp_condensing") == pytest.approx(total)
    total_c = r.get("dp_coolant_channel") + r.get("dp_coolant_port") + r.get("dp_coolant_static")
    assert r.get("dp_coolant") == pytest.approx(total_c)


def test_static_head_follows_the_coolant_path():
    up = run()
    assert up.get("dp_coolant_static") == pytest.approx(985.2 * 9.81 * 1.0)

    def two_pass(i):
        i.coolant.passes = 2

    assert run(two_pass).get("dp_coolant_static") == 0.0

    def cocurrent(i):
        i.coolant.direction = "cocurrent"

    assert run(cocurrent).get("dp_coolant_static") == pytest.approx(-985.2 * 9.81 * 1.0)


def test_zone_refinement_converges():
    def zones(n):
        def mutate(i):
            i.solver.n_zones = n
            i.limits.dp_allow_condensing = None
            i.limits.dp_allow_coolant = None

        return run(mutate)

    assert zones(20).get("area_required") == pytest.approx(zones(160).get("area_required"), rel=0.01)
    assert zones(20).get("dp_condensing_friction") == pytest.approx(zones(160).get("dp_condensing_friction"), rel=0.03)


# ------------------------------------------------------------------------- iteration and its record
def no_limits(i):
    i.limits.dp_allow_condensing = None
    i.limits.dp_allow_coolant = None


def test_iteration_is_recorded_and_stops_at_the_smallest_count():
    result = run(no_limits)
    info = result.convergence[0]
    assert info.converged and info.rows and len(info.rows[0]) == len(info.columns)
    assert any("First guess" in a for a in result.assumptions)
    N = result.get("n_plates")
    assert result.get("sized_by") == "duty"
    assert result.get("area") >= result.get("area_required")
    # the next smaller allowed count must fail on area
    below = [row for row in info.rows if row[0] == N - 2]
    assert below and below[0][9] < below[0][8]


def test_multipass_count_is_a_multiple_of_passes():
    def mutate(i):
        no_limits(i)
        i.coolant.passes = 2

    result = run(mutate)
    assert result.get("n_channels_per_side") % 2 == 0
    assert check(result, "Single-pass arrangement avoids thermally induced maldistribution").state is CheckState.FAIL
    assert any("should be avoided" in w for w in result.warnings)


def test_pressure_limited_design_takes_more_plates_and_excess_area():
    free = run(no_limits)
    limited = run()
    assert limited.get("n_plates") > free.get("n_plates")
    assert limited.get("sized_by") == "allowable pressure drop"
    assert limited.get("excess_area_percent") > free.get("excess_area_percent")
    assert limited.get("dp_condensing") <= example().limits.dp_allow_condensing
    assert check(limited, "Condensing pressure drop within the allowable").state is CheckState.PASS
    assert limited.get("T_c_out_rated") > limited.get("T_c_out")
    # the search stops at the smallest count that meets the limit
    search = next(c for c in limited.convergence if "allowable" in c.name)
    assert search.converged and search.rows


def test_port_loss_cannot_be_cured_by_plates():
    def small_port(i):
        i.geometry.port_d_vapour = 0.005
        i.geometry.port_d = 0.1

    result = run(small_port)
    assert any("port loss does not depend on N" in w for w in result.warnings)
    assert check(result, "Condensing pressure drop within the allowable").state is CheckState.FAIL


def test_stops_at_the_plate_limit_when_the_allowable_is_out_of_reach():
    def tight(i):
        i.limits.dp_allow_condensing = None
        # Port 13.41 Pa plus static 9664.81 Pa leave 0.28 Pa for the channels, which 700 plates cannot reach.
        i.limits.dp_allow_coolant = 9678.5

    result = run(tight)
    assert result.status is Status.OK_WITH_WARNINGS
    assert any(f"{MAX_PLATES}-plate limit" in w for w in result.warnings)
    assert result.get("n_plates") == 699
    assert check(result, "Coolant pressure drop within the allowable").state is CheckState.FAIL


# ---------------------------------------------------------------------------------------- trends
def test_trend_fouling_needs_more_area():
    def fouled(i):
        no_limits(i)
        i.fouling.R_f_coolant = 5e-5

    clean, dirty = run(no_limits), run(fouled)
    assert dirty.get("area_required") > clean.get("area_required")
    assert dirty.get("cleanliness_factor") < 1.0
    assert clean.get("cleanliness_factor") == pytest.approx(1.0)


def test_trend_cocurrent_gives_lower_vapour_drop_than_countercurrent():
    # HEDH Sec. 3.7.12 states this trend at a given steam flow: in concurrent flow most of the condensation occurs in the
    # top half of the plate. Same plate count on both, so the same mass flux.
    counter, co = example(), example()
    co.coolant.direction = "cocurrent"
    drops = [_hydraulics(i, _duty(i), _evaluate(i, _duty(i), 15)).dp_f for i in (counter, co)]
    assert drops[1] < drops[0]


def test_coolant_passes_at_the_same_plate_count_follow_hsb_example_8_2():
    # Same N: velocity doubles, h rises by 2^0.64 and the channel drop by 2 x 2^1.82 = 7.06 (HSB Ex. 8.2).
    one, two = example(), example()
    two.coolant.passes = 2
    ev1, ev2 = _evaluate(one, _duty(one), 15), _evaluate(two, _duty(two), 15)
    assert ev2.u_c == pytest.approx(2 * ev1.u_c)
    assert ev2.h_c == pytest.approx(2**0.64 * ev1.h_c)
    hy1, hy2 = _hydraulics(one, _duty(one), ev1), _hydraulics(two, _duty(two), ev2)
    assert hy2.dp_c_channel / hy1.dp_c_channel == pytest.approx(2 * 2**1.82)
    assert hy2.dp_c_port == pytest.approx(2 * hy1.dp_c_port)  # Lecture 14: per pass


def test_double_vapour_entry_quarters_the_port_loss():
    def double(i):
        i.vapour.vapour_entry = "double"

    assert run(double).get("dp_condensing_port") == pytest.approx(run().get("dp_condensing_port") / 4.0)


def test_enhancement_factors_act_on_h_and_dp():
    def enhanced(i):
        no_limits(i)
        i.correlation.enhancement_h = 2.0
        i.correlation.enhancement_f = 10.0

    base, high = run(no_limits), run(enhanced)
    assert high.get("n_plates") < base.get("n_plates")
    assert high.get("h_coolant") > base.get("h_coolant")


def test_chevron_angle_warns_when_factors_are_left_at_one():
    def angle(i):
        i.geometry.chevron_angle = 60.0

    assert any("HSB Ex. 8.4 uses 2 on h and 10 on f" in w for w in run(angle).warnings)


def test_developed_basis_uses_its_own_diameter_and_area():
    def developed(i):
        no_limits(i)
        i.geometry.area_basis = "developed"

    result = run(developed)
    assert result.get("area") == pytest.approx(result.get("n_plates") * 0.25 * 1.17)
    assert any(s.symbol == "De" and s.value == pytest.approx(2 * 5e-3 / 1.17) for s in result.steps)


def test_wall_viscosity_shifts_nu_and_f_oppositely():
    base, wall = example(), example()
    wall.coolant.mu_wall = 2.5e-4
    wall.correlation.m_visc = 0.1
    ratio = (5.05e-4 / 2.5e-4) ** 0.1
    ev0, ev1 = _evaluate(base, _duty(base), 15), _evaluate(wall, _duty(wall), 15)
    assert ev1.h_c == pytest.approx(ev0.h_c * ratio)
    assert _hydraulics(wall, _duty(wall), ev1).dp_c_channel == pytest.approx(_hydraulics(base, _duty(base), ev0).dp_c_channel / ratio)


# -------------------------------------------------------------------------- coolant specification
def test_coolant_flow_given_gives_outlet():
    def by_flow(i):
        i.coolant.T_out = None
        i.coolant.m_dot = 2.0

    result = run(by_flow)
    assert result.get("T_c_out") == pytest.approx(30 + 225670 / (2.0 * 4183.0))


def test_assumed_effectiveness_of_lecture_17():
    def by_p(i):
        i.coolant.T_out = None
        i.coolant.assumed_effectiveness = 0.85

    result = run(by_p)
    assert result.get("T_c_out") == pytest.approx(30 + 0.85 * 70)
    assert check(result, "Assumed effectiveness in the range Lecture 17 gives").state is CheckState.PASS


# ------------------------------------------------------------------------------- limits and checks
def test_gasket_and_operating_limits():
    def hot(i):
        i.condensing.T_sat = 200.0
        i.vapour.pressure = 40e5
        i.vapour.flammable_or_toxic = "yes"

    result = run(hot)
    assert check(result, "Operating temperature limit").state is CheckState.FAIL
    assert check(result, "Operating pressure limit").state is CheckState.FAIL
    assert check(result, "Gasket material temperature limit").state is CheckState.FAIL
    assert check(result, "Flammable or toxic fluid with gaskets").state is CheckState.FAIL


def test_limits_at_180_degc_and_30_bar_pass():
    def at_limit(i):
        i.condensing.T_sat = 180.0
        i.vapour.pressure = 30e5
        i.geometry.gasket_material = None

    result = run(at_limit)
    assert check(result, "Operating temperature limit").state is CheckState.PASS
    assert check(result, "Operating pressure limit").state is CheckState.PASS
    assert check(result, "Gasket material temperature limit").state is CheckState.NOT_ASSESSED


def test_welded_unit_uses_its_own_range_and_drops_the_gasket_checks():
    def welded(i):
        i.geometry.construction = "brazed or welded"
        i.condensing.T_sat = 200.0

    result = run(welded)
    assert check(result, "Operating temperature limit").state is CheckState.PASS
    assert not [c for c in result.checks if "Gasket" in c.name or "gasket" in c.name]


def test_vacuum_is_flagged():
    def vacuum(i):
        i.vapour.pressure = 0.2e5

    assert check(run(vacuum), "Low pressure and vacuum").state is CheckState.FAIL


def test_correlation_range_is_checked_where_used():
    def slow(i):
        no_limits(i)
        i.correlation.re_min = 1e5

    result = run(slow)
    assert check(result, "Coolant Reynolds number in the correlation range").state is CheckState.FAIL
    assert check(result, "Effective condensing Reynolds number in the correlation range").state is CheckState.FAIL


def test_code_and_maldistribution_items_are_not_assessed():
    result = run()
    code = check(result, "Plate, frame and tie-bolt strength, and pressure-boundary code")
    assert code.state is CheckState.NOT_ASSESSED and code.category is Category.MECHANICAL and "ASME" in code.detail
    assert check(result, "Vapour port velocity and maldistribution along the port").state is CheckState.NOT_ASSESSED
    assert check(result, "Pressure-drop feedback on the condensing temperature").state is CheckState.NOT_ASSESSED


def test_coolant_static_head_rule_applies_only_in_downflow():
    assert check(run(), "Coolant channels run full in downflow: channel drop at least the static head").state is CheckState.NOT_ASSESSED

    def cocurrent(i):
        i.coolant.direction = "cocurrent"

    assert check(run(cocurrent), "Coolant channels run full in downflow: channel drop at least the static head").state is not CheckState.NOT_ASSESSED


# ------------------------------------------------------------------------------------ traceability
def test_every_step_is_a_full_chain_with_a_source_and_category():
    for step in run().steps:
        assert step.formula, step.label
        assert step.source, step.label
        assert step.category in (Category.THERMAL, Category.HYDRAULIC, Category.GEOMETRIC), step.label
        if step.symbol != "N":
            assert step.substitution, step.label


def test_every_check_states_its_criterion():
    for c in run().checks:
        assert c.detail, c.name


def test_reports_build_with_the_result():
    inputs = example()
    result = run_design("plate", inputs)
    assert to_excel(result, inputs)[:2] == b"PK"
    assert to_pdf(result, inputs).startswith(b"%PDF")


# ------------------------------------------------------------------------ bad input, via the service
def test_missing_inputs_are_named():
    result = run_design("plate", MODULE.new_inputs())
    assert result.status is Status.INVALID_INPUT
    fields = {issue.field for issue in result.input_issues}
    assert {"geometry.L", "geometry.W", "geometry.b", "correlation.C_h", "correlation.K_p", "condensing.T_sat", "vapour.m_dot", "coolant.T_in", "coolant_props.rho"} <= fields


def test_correlation_coefficients_have_no_defaults():
    inputs = MODULE.new_inputs()
    assert inputs.correlation.C_h is None and inputs.correlation.y is None and inputs.correlation.K_p is None and inputs.correlation.z is None


@pytest.mark.parametrize(
    "mutate,field",
    [
        (lambda i: setattr(i.coolant, "T_out", 120.0), "coolant.T_out"),
        (lambda i: setattr(i.coolant, "T_in", 100.0), "coolant.T_in"),
        (lambda i: setattr(i.coolant, "T_out", 20.0), "coolant.T_out"),
        (lambda i: setattr(i.coolant, "m_dot", 1.0), "coolant.T_out"),
        (lambda i: setattr(i.vapour, "x_out", 1.0), "vapour.x_out"),
        (lambda i: setattr(i.coolant, "passes", 1.5), "coolant.passes"),
        (lambda i: setattr(i.solver, "n_zones", 2.5), "solver.n_zones"),
        (lambda i: (setattr(i.geometry, "area_basis", "developed"), setattr(i.geometry, "phi", None)), "geometry.phi"),
        (lambda i: setattr(i.geometry, "b", 1.0), "geometry.b"),
        (lambda i: setattr(i.coolant, "mu_wall", 1e-4), "correlation.m_visc"),
        (lambda i: setattr(i.geometry, "L", -1.0), "geometry.L"),
        (lambda i: setattr(i.correlation, "y", 1.5), "correlation.y"),
        (lambda i: (setattr(i.coolant, "passes", 2), setattr(i.coolant, "direction", "cocurrent")), "coolant.direction"),
    ],
)
def test_invalid_values(mutate, field):
    result = run(mutate)
    assert result.status is Status.INVALID_INPUT
    assert field in {issue.field for issue in result.input_issues}


def test_coolant_flow_too_low_is_invalid():
    def starved(i):
        i.coolant.T_out = None
        i.coolant.m_dot = 0.5

    result = run(starved)
    assert result.status is Status.INVALID_INPUT
    assert "Raise the flow" in result.input_issues[0].message


def test_non_convergence_reports_with_its_history():
    def one_iteration(i):
        i.solver.max_iterations = 1

    result = run(one_iteration)
    assert result.status is Status.NOT_CONVERGED
    assert result.convergence and not result.convergence[0].converged
    assert result.convergence[0].rows
    assert "no convergence" in result.message


def test_module_registered():
    assert MODULE.key == "plate"
    assert lmtd(70.0, 20.0) == pytest.approx(run().get("lmtd"))
