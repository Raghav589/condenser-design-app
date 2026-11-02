"""Vertical shell-and-tube condenser.

NO WORKED VERTICAL CONDENSER EXAMPLE exists in HSB Ch. 16 to 18, HEDH 3.3 / 3.4 or Lectures 9 to 12, 23, 24, so the whole path has no benchmark.
What is benchmarked are SUB-PATHS, labelled as such: shell-side Bell-Delaware (HEDH 3.3.9, HSB Ex. 6.2), Kern (HSB Ex. 6.1), tube-side coolant (HSB Ex. 6.4), the
tube count (HEDH 3.3.9), the shear term of the condensing coefficient (HSB Ex. 16.2) and flooding (HSB Ex. 10.5). The rest is tested by hand-checked values,
limiting cases, closure and trends, which are not benchmarks.
"""
import math

import pytest

from condenser_app.common.condensation import gravity_local, gravity_mean, nusselt_mean, vertical_downflow_local
from condenser_app.common.properties import pentane_hsb_example_16_2, saturated, saturated_vapour, saturation_pressure, water_liquid
from condenser_app.modules import MODULES
from condenser_app.modules import vertical_shell_tube as vst
from condenser_app.modules.vertical_shell_tube import (
    DOWNFLOW, LAYOUT_30, LAYOUT_90, MATERIAL_A, MATERIAL_B, REFLUX, SHELL_SIDE,
    BundleGeometry, CoolantStream, SolverOptions, VapourStream, VelocityLimits, VerticalShellTubeInputs,
)
from condenser_app.reporting import to_excel, to_pdf
from condenser_app.schemas import CheckState, FoulingFactors, PressureDropLimits, Status
from condenser_app.service import run_design

KEY = "vertical_shell_tube"
T_SAT = 61.85


def inputs(variant=DOWNFLOW, geometry=None, vapour=None, coolant=None, fouling=None, limits=None, velocity=None, solver=None):
    fluid = saturated("n_pentane", T_SAT)
    gas = saturated_vapour("n_pentane", T_SAT)
    return VerticalShellTubeInputs(
        variant=variant,
        geometry=BundleGeometry(**{**dict(D_o=25.4e-3, t_wall=2.11e-3, k_wall=18.0, L=4.0, pitch=31.75e-3, layout=LAYOUT_90, material_group=MATERIAL_A, shell_fouling="yes", L_bb=15e-3, L_tb=0.8e-3, N_ss=0, baffle_cut=25.0, baffle_spacing_ratio=0.4), **(geometry or {})}),
        condensing=fluid,
        vapour=VapourStream(**{**dict(m_dot=5.0, x_in=1.0, x_out=0.0, pressure=saturation_pressure("n_pentane", T_SAT), cp_v=gas.cp, k_v=gas.k, condensation_curve=""), **(vapour or {})}),
        coolant_props=water_liquid(32.5),
        coolant=CoolantStream(**{**dict(T_in=25.0, T_out=40.0), **(coolant or {})}),
        fouling=FoulingFactors(**{**dict(R_f_condensing=0.0, R_f_coolant=0.0), **(fouling or {})}),
        limits=PressureDropLimits(**(limits or {})),
        velocity=VelocityLimits(**(velocity or {})),
        solver=SolverOptions(**(solver or {})),
    )


def run(**kw):
    return run_design(KEY, inputs(**kw))


def check_state(result, name_part):
    return next(c.state for c in result.checks if name_part in c.name)


# ------------------------------------------------------------------------ SUB-PATH benchmark: Bell-Delaware, HEDH 3.3.9
def hedh_case(**kw):
    # HEDH 3.3.9B data (Table 1). Wall viscosity 0.405 cP at T_w = 169 degC as the book estimates it; mean viscosity 0.494 cP; two tube passes, so L_pl = 6.35 mm in the bypass area.
    args = dict(M=22.4, rho=740.0, cp=2407.0, mu=0.494e-3, k=0.105, D_s=0.591, D_t=0.016, pitch=20.8e-3, layout=LAYOUT_30, baffle_cut=24.0, L_bc=0.250, L_ti=2.276, N_tt=652, L_bb=0.010, L_tb=0.5e-3, L_sb=5.5e-3, N_ss=0, mu_wall=0.405e-3, L_pl=6.35e-3)
    return vst.bell_delaware(**{**args, **kw})


@pytest.mark.parametrize("attr, expected, tol", [
    # The book prints 3 or 4 significant figures and reads j, f, J and R from charts, which it says match the fits used here for computers. 0.5% unless noted.
    ("S_m", 35096e-6, 0.005), ("m_s", 638.0, 0.005),
    ("Re", 20695.0, 0.005),  # the book prints both 20 695 and 20 659; 16 (638)/0.494 = 20 664
    ("Pr", 11.29, 0.005), ("theta_ds", 117.3, 0.005), ("theta_ctl", 114.0, 0.005),
    ("F_w", 0.172, 0.005), ("F_c", 0.657, 0.005), ("S_w", 28117e-6, 0.005), ("N_tcc", 17.06, 0.005),
    ("N_tcw", 5.7, 0.005), ("N_b", 8, 0.0), ("S_b", 4088e-6, 0.005), ("F_sbp", 0.116, 0.005),
    ("S_sb", 3442e-6, 0.005), ("S_tb", 6893e-6, 0.005), ("r_lm", 0.294, 0.005), ("r_s", 0.333, 0.005),
    ("J_c", 1.02, 0.005), ("J_l", 0.665, 0.005), ("R_l", 0.427, 0.005), ("J_b", 0.865, 0.005), ("R_b", 0.65, 0.005),
    ("j_i", 0.0068, 0.005), ("f_i", 0.11, 0.005), ("m_w", 713.0, 0.005), ("alpha_i", 2134.0, 0.005),
    ("alpha_s", 1252.0, 0.005), ("J_tot", 0.59, 0.01),  # book 2 s.f.
    ("dp_bi", 2.0e3, 0.015), ("dp_c", 3.9e3, 0.015), ("dp_w", 6.4e3, 0.015),  # the book's kPa values have 2 s.f.
    ("dp_e", 3.4e3, 0.04),  # the book carries its rounded 2.0 kPa; with 2.017 this is 3.50
    ("dp_s", 13.7e3, 0.01),
])
def test_hedh_3_3_9_end_to_end(attr, expected, tol):
    assert getattr(hedh_case(), attr) == pytest.approx(expected, rel=tol if tol else 1e-12)


def test_hedh_3_3_9_overall_coefficient_and_duty():
    # Book U_o = 971 W/m2.K, Q = 3 947 794 W, with alpha_t = 6000 and its plain-wall R_w = t/k; the cylinder wall used here adds 0.3%. HEDH measured alpha_s = 1335, dP_s = 12.5 kPa.
    s = hedh_case()
    u = vst.overall_u_outer(s.alpha_s, 6000.0, 0.016, 0.01352, 42.0)
    assert u == pytest.approx(971.0, rel=0.005)
    area = math.pi * 0.016 * 2.276 * 652
    assert area == pytest.approx(74.6, rel=0.001)
    assert area * u * vst.lmtd(88.2, 30.7) == pytest.approx(3947794.0, rel=0.01)


def test_no_pass_lane_means_less_bypass():
    one_pass, two_pass = hedh_case(L_pl=0.0), hedh_case()
    assert one_pass.F_sbp < two_pass.F_sbp
    assert one_pass.alpha_s > two_pass.alpha_s
    assert one_pass.dp_s > two_pass.dp_s


def test_hedh_tube_count_estimate():
    # HEDH 3.3.9B(c) item 15: 0.78 (565)^2 / ((0.866)(20.8)^2) (1 - 0.04) = 638. Book 3 s.f.
    assert vst.tube_count(0.565, 20.8e-3, LAYOUT_30, 0.04) == pytest.approx(638.0, rel=0.002)
    d_ctl, d_otl, d_s = vst.bundle_diameters(vst.tube_count(0.565, 20.8e-3, LAYOUT_30), 0.016, 20.8e-3, LAYOUT_30, 0.010)
    assert (d_ctl, d_otl, d_s) == pytest.approx((0.565, 0.581, 0.591), rel=1e-9)


# ------------------------------------------------------------------- SUB-PATH benchmark: Kern, HSB Ex. 6.1
def test_hsb_example_6_1_kern():
    k = vst.kern_shell_side(5.5188, 730.0, 2470.0, 401e-6, 0.1324, 0.5398, 0.0254, 0.03175, LAYOUT_90, 0.1270)
    assert k.S_s == pytest.approx(0.01371, rel=0.001)
    assert k.D_e == pytest.approx(0.02513, rel=0.001)
    assert k.m_s == pytest.approx(402.5, rel=0.001)
    assert k.Re == pytest.approx(25224.0, rel=0.001)
    assert k.Pr == pytest.approx(7.481, rel=0.001)
    assert k.alpha == pytest.approx(977.94, rel=0.002)  # book 5 s.f. but rounded inputs
    # f = 0.063 read by the book from Fig. 6.6, N = 36. Book 22.2 kPa; its own printed line gives 22,224 Pa.
    assert vst.kern_pressure_drop(0.063, k.m_s, 0.5398, 36, 730.0, k.D_e) == pytest.approx(22224.0, rel=0.002)


# ------------------------------------------------------------------- SUB-PATH benchmark: Bell-Delaware, HSB Ex. 6.2
@pytest.fixture
def hsb_62():
    # Clearances are diametral: 0.8, 5, 35 mm. N_ss/N_c = 1/5, so N_ss = 8.5/5. N = 36 is the book's count, with the baffle thickness.
    return vst.bell_delaware(M=5.5188, rho=730.0, cp=2470.0, mu=401e-6, k=0.1324, D_s=0.5398, D_t=0.0254, pitch=0.03175, layout=LAYOUT_90, baffle_cut=25.0, L_bc=0.1270, L_ti=4.8768, N_tt=158, L_bb=0.035, L_tb=0.8e-3, L_sb=5e-3, N_ss=8.5 / 5.0, N_b=36)


def test_hsb_example_6_2_geometry_and_factors(hsb_62):
    s = hsb_62
    # S_m: book 0.01667 m2, but its own Eq. 6.9 and numbers give 0.01662; 0.5% covers this.
    assert s.S_m == pytest.approx(0.01667, rel=0.005)
    assert s.u_s == pytest.approx(0.4535, rel=0.005)
    assert s.Re == pytest.approx(2.097e4, rel=0.005)
    assert s.F_sbp == pytest.approx(0.2666, rel=0.005)
    assert s.N_tcc == pytest.approx(8.5, rel=0.002)
    assert s.S_sb == pytest.approx(0.002827, rel=0.005)
    # Book chart readings against the HEDH fits: R_l 0.33, R_b 0.76, J_b 0.92; J_l 0.58 read from Fig. 6.8 (the fit gives 4% lower).
    assert s.R_l == pytest.approx(0.33, abs=0.01)
    assert s.R_b == pytest.approx(0.76, abs=0.015)
    assert s.J_b == pytest.approx(0.92, abs=0.01)
    assert s.J_l == pytest.approx(0.58, rel=0.05)
    # HSB Eq. 6.13 and N_cw = 0.8 L_c/P use D_OTL and the whole cut; HEDH Eqs. 3.3.6(2), (13) use D_ctl and subtract the bypass gap. That moves F_c from 0.647 to 0.677, J_c from 1.016 to 1.037, N_tcw from 3.4 to 2.6.
    assert s.F_c == pytest.approx(0.647, rel=0.05)
    assert s.J_c == pytest.approx(1.016, rel=0.025)
    assert s.N_tcw == pytest.approx(3.4, rel=0.25)


def test_hsb_example_6_2_results_and_cause_of_difference(hsb_62):
    s = hsb_62
    # HSB's ideal coefficient is 1418 W/m2.K (Zukauskas, Table 2.5). HEDH's j_i fit gives about 1570, 11% higher: two different ideal tube-bank data sets. Not tuned.
    assert s.alpha_i == pytest.approx(1418.0, rel=0.13)
    # Book result 769 W/m2.K. The difference is that 11%, J_c +2%, J_l -5%.
    assert s.alpha_s == pytest.approx(769.0, rel=0.10)
    # Book Delta p_bi = 233 Pa with Zukauskas K_f; the f_i fit gives about 242 Pa.
    assert s.dp_bi == pytest.approx(233.0, rel=0.06)
    # Book total 4501 Pa, with 165 Pa per window, and a smaller N_tcw: window terms and R_s = 2 differ.
    assert s.dp_s == pytest.approx(4501.0, rel=0.12)


# ------------------------------------------------------------------- SUB-PATH benchmark: tube side, HSB Ex. 6.4
def test_hsb_example_6_4_tube_side():
    # 158 tubes in 4 passes: 39.5 per pass; four passes of 4.8768 m. Book friction only, with Blasius, 38.76 kPa; Dittus-Boelter 985 W/m2.K.
    t = vst.tube_side_coolant(18.77, 158 / 4, 0.020574, 4.8768, 830.0, 1971.0, 3596e-6, 0.1333, passes=4)
    assert t.Re == pytest.approx(8175.0, rel=0.002)
    assert t.Pr == pytest.approx(53.2, rel=0.002)
    assert t.dp_friction == pytest.approx(38760.0, rel=0.003)
    assert vst.dittus_boelter(t.Re, t.Pr, True) * 0.1333 / 0.020574 == pytest.approx(985.0, rel=0.002)
    # Below Re = 1e4 this code selects Gnielinski; the book's Dittus-Boelter value is used above only as the correlation check.
    assert t.correlation == "Gnielinski"


# ------------------------------------------------------------------- SUB-PATH benchmark: condensing shear term and flooding
def test_hsb_example_16_2_is_the_shear_term_used():
    # n-pentane, m = 500, x = 0.85, D = 0.03. Book alpha_f = 11,104 W/m2.K, tau rounded (library test: -0.04%). The module calls this function.
    local = vertical_downflow_local(pentane_hsb_example_16_2(), 500.0, 0.85, 0.03)
    assert local.alpha_shear == pytest.approx(11104.0, rel=0.002)


def test_hsb_example_10_5_flooding_in_the_reflux_variant():
    # Ex. 10.5: pentane reflux, D = 3 cm, rho_L = 529, rho_G = 19, C = 0.75. Book m_crit = 14.6 kg/m2.s, 3 s.f. (library test: -0.13%). Through the module's output here.
    i = inputs(variant=REFLUX, geometry={"L": 2.5, "D_o": 34e-3, "t_wall": 2e-3, "pitch": 42.5e-3, "layout": LAYOUT_30, "shell_fouling": "no"}, vapour={"m_dot": 1.0})
    i.condensing.rho_l, i.condensing.rho_v = 529.0, 19.0
    assert run_design(KEY, i).get("G_flood_hsb") == pytest.approx(14.6, abs=0.05)


# ---------------------------------------------------------------------------------- whole module
def test_registered_and_example_runs():
    assert KEY in MODULES
    module = MODULES[KEY]
    r = run_design(KEY, module.example_inputs())
    assert r.status in (Status.OK, Status.OK_WITH_WARNINGS)
    assert r.get("L_required") <= r.get("L_available")
    assert r.get("dp_coolant") <= module.example_inputs().limits.dp_allow_coolant
    assert r.get("sized_by") == "required length"
    assert r.get("n_tubes") >= 8


@pytest.mark.parametrize("make", [lambda: MODULES[KEY].example_inputs(), lambda: shell_inputs(), lambda: reflux_inputs(), lambda: dry_inputs()])
def test_reports_build(make):
    i = make()
    r = run_design(KEY, i)
    assert to_pdf(r, i).startswith(b"%PDF")
    assert to_excel(r, i)[:2] == b"PK"


def test_every_step_is_a_chain_and_gaps_are_stated():
    r = run_design(KEY, MODULES[KEY].example_inputs())
    for s in r.steps:
        assert s.label and s.symbol and s.formula and s.substitution and s.source, s
        assert "=" in s.formula or s.symbol in ("N",), s
    assert any("NO WORKED VERTICAL CONDENSER EXAMPLE" in g for g in r.source_gaps)
    assert any("Why it iterates" in a and "First guess" in a and "Update" in a and "Stop" in a for a in r.assumptions)
    assert any("Bell-Delaware" in s.label and "Kern" in s.label for s in r.steps)
    assert any("Kosky-Staub" in a and "Boyko-Kruzhilin" in a for a in r.assumptions)
    assert any("Suction" in a for a in r.assumptions)
    for name in ("Tube vibration", "Pressure boundary", "Coolant nozzle size", "Venting"):
        assert check_state(r, name) is CheckState.NOT_ASSESSED
    assert any("ASME" in c.detail and "TEMA" in c.detail for c in r.checks)
    assert {s.category.value for s in r.steps} >= {"thermal", "hydraulic", "geometric"}


def test_energy_closure_and_coolant_profile():
    r = run(vapour={"T_in": 110.0, "T_out": 45.0})
    module_inputs = inputs(vapour={"T_in": 110.0, "T_out": 45.0})
    d = vst._duty(module_inputs)
    ev = vst._evaluate(module_inputs, d, r.get("n_tubes"))
    q_zones = (ev.dry.Q if ev.dry else 0.0) + sum(s.dQ for s in ev.steps) + ev.sub.Q
    assert q_zones == pytest.approx(d.Q, rel=1e-9)
    assert ev.steps[0].t1 == pytest.approx(ev.dry.t_c_star)
    assert ev.steps[-1].t2 == pytest.approx(ev.sub.t_c_prime)
    assert d.C_c * (d.T_c_out - 25.0) == pytest.approx(d.Q, rel=1e-9)


def test_equation_18_13_solves_18_10_and_18_12():
    t_sat, t_c_out, t_h_in, r_u, R = 60.0, 45.0, 120.0, 0.4, 0.05
    t_h = vst.wet_wall_vapour_temperature(t_sat, t_c_out, t_h_in, r_u, R)
    t_c = t_c_out - (t_h_in - t_h) * R  # Eq. 18.12
    assert t_h - r_u * (t_h - t_c) == pytest.approx(t_sat, rel=1e-12)  # Eq. 18.10


def dry_inputs():
    # A hot vapour, warm coolant and a coolant-side resistance, so the wall stays above T_sat near the inlet.
    i = inputs(vapour={"T_in": 300.0, "cp_v": 2500.0, "k_v": 0.03}, coolant={"T_in": 50.0, "T_out": 60.0}, fouling={"R_f_coolant": 0.001})
    i.condensing = saturated("n_pentane", 76.85)
    return i


def test_no_dry_zone_when_the_wall_is_below_saturation():
    wet = run(vapour={"T_in": 66.0})
    assert wet.get("A_dsh") == 0.0
    assert wet.get("T_h_star") == 66.0
    assert any("no dry zone" in s.label for s in wet.steps)


def test_dry_zone_area_against_hand_values():
    i = dry_inputs()
    d = vst._duty(i)
    ev = vst._evaluate(i, d, 100)
    z = ev.dry
    assert z.has_zone and i.condensing.T_sat < z.t_h_star < 300.0
    hand_q = i.vapour.m_dot * 2500.0 * (300.0 - z.t_h_star)
    assert z.Q == pytest.approx(hand_q, rel=1e-9)
    hand_mtd = ((300.0 - d.T_c_out) - (z.t_h_star - z.t_c_star)) / math.log((300.0 - d.T_c_out) / (z.t_h_star - z.t_c_star))
    assert z.A == pytest.approx(hand_q / (z.U * hand_mtd), rel=1e-9)
    assert ev.Q_con == pytest.approx(i.vapour.m_dot * (2500.0 * (z.t_h_star - i.condensing.T_sat) + i.condensing.h_fg), rel=1e-9)  # HSB Eq. 18.17


def test_dry_zone_in_a_full_run_and_hsbs_conservative_comparison():
    r = run_design(KEY, dry_inputs())
    assert r.get("A_dsh") > 0.0
    conservative = next(s.value for s in r.steps if s.symbol == "A_dsh,cons")
    assert conservative > r.get("A_dsh")  # HSB Sec. 18.3.2: no condensation until T_sat is the conservative extreme
    assert any("Dittus-Boelter range" in c.name and "Re_v" in c.detail for c in r.checks)  # recorded where it is used


def test_subcooling_coefficient_hand_values_and_lecture_24_form():
    p = pentane_hsb_example_16_2()
    gamma, length = 0.02, 0.5
    alpha, regime, re_tr = vst.subcooling_film_coefficient(p, gamma, length)
    pr = p.cp_l * p.mu_l / p.k_l
    assert re_tr == pytest.approx(2460.0 * pr**-0.65)
    assert regime.startswith("laminar")
    x = p.mu_l**2 / (9.81 * p.rho_l**2)
    # Lecture 24: h_sc = 0.78 k (cp mu/(k L))^(1/3) (mu^2/(g rho^2))^(-2/9) (4 Gamma/mu)^(1/9)
    lecture = 0.78 * p.k_l * (pr / length) ** (1 / 3) * x ** (-2 / 9) * (4 * gamma / p.mu_l) ** (1 / 9)
    assert alpha == pytest.approx(lecture, rel=1e-6)
    alpha2, regime2, _ = vst.subcooling_film_coefficient(p, 0.15, length)
    re = 4 * 0.15 / p.mu_l
    assert regime2.startswith("transition")
    assert alpha2 == pytest.approx(0.032 * p.k_l * re**0.2 * pr**0.34 * x ** (-1 / 3), rel=1e-6)
    alpha3, regime3, _ = vst.subcooling_film_coefficient(p, 1.0, length)
    re3 = 4.0 / p.mu_l
    assert regime3.startswith("turbulent")
    assert alpha3 == pytest.approx(0.0057 * p.k_l * re3**0.4 * pr**0.34 * x ** (-1 / 3), rel=1e-6)


def test_constant_u_limit_recovers_lmtd():
    # A large fouling resistance dominates, so U is nearly constant and HSB Eq. 18.9 must give the same area.
    i = inputs(fouling={"R_f_coolant": 0.05}, geometry={"n_tubes": 300})
    r = run_design(KEY, i)
    d = vst._duty(i)
    dt = vst.lmtd(T_SAT - 25.0, T_SAT - 40.0)
    assert r.get("area_required") == pytest.approx(d.Q / (r.get("U_mean") * dt), rel=0.001)


def test_zero_length_like_limits_and_shear_goes_to_gravity_at_low_flow():
    # Many tubes: low shear, so alpha tends to the gravity term alone.
    i = inputs(geometry={"n_tubes": 5000})
    d = vst._duty(i)
    ev = vst._evaluate(i, d, 5000)
    s = ev.steps[-1]
    assert s.alpha == pytest.approx(s.alpha_gravity, rel=0.02)
    # Few tubes: high shear, so the shear term dominates near the inlet.
    i = inputs(geometry={"n_tubes": 20})
    d = vst._duty(i)
    ev = vst._evaluate(i, d, 20)
    assert ev.steps[0].alpha == pytest.approx(ev.steps[0].alpha_shear, rel=0.05)


def test_iteration_finds_smallest_count():
    i = inputs(limits={"dp_allow_condensing": 20e3, "dp_allow_coolant": 35e3})
    r = run_design(KEY, i)
    n = r.get("n_tubes")
    assert r.convergence[0].name == "Number of tubes"
    assert r.convergence[0].columns[0] == "N"
    assert r.get("L_required") <= 4.0
    d = vst._duty(i)
    assert max(vst._ratios(i, vst._evaluate(i, d, n - 1)).values()) > 1.0
    assert max(vst._ratios(i, vst._evaluate(i, d, n)).values()) <= 1.0


def test_first_guess_does_not_change_the_answer():
    a = run().get("n_tubes")
    b = run(solver={"U_estimate": 800.0})
    assert b.get("n_tubes") == a
    assert b.convergence[0].rows[0][0] > 8  # it did start from A = Q/(U MTD)


def test_limits_drive_the_count():
    base = run().get("n_tubes")
    r = run(velocity={"u_coolant_max": 0.5})
    assert r.get("u_coolant") <= 0.5 * (1 + 1e-9)
    assert r.get("n_tubes") > base
    assert r.get("sized_by") == "coolant velocity"
    r = run(limits={"dp_allow_condensing": 200.0})
    assert r.get("dp_condensing") <= 200.0 * (1 + 1e-9)
    assert r.get("n_tubes") > base
    assert r.get("sized_by") == "condensing pressure drop"
    r = run(velocity={"u_vapour_max": 3.0})
    assert r.get("u_vapour_inlet") <= 3.0 * (1 + 1e-9)


def test_infeasible_limit_warns_at_the_largest_count():
    r = run(limits={"dp_allow_coolant": 1.0}, solver={"n_tubes_max": 500})
    assert r.get("n_tubes") == 500
    assert any("largest tube count" in w for w in r.warnings)
    assert check_state(r, "Coolant pressure drop") is CheckState.FAIL


def test_non_convergence_reported():
    r = run(solver={"max_iterations": 1})
    assert r.status is Status.NOT_CONVERGED
    assert r.convergence and r.convergence[0].rows


def test_rating_mode_reports_length():
    r = run(geometry={"n_tubes": 100})
    assert r.get("n_tubes") == 100
    assert r.get("sized_by") == "fixed by the user"
    assert r.get("L_required") > 4.0
    assert check_state(r, "Tube length covers") is CheckState.FAIL


def test_trends():
    base = run(geometry={"n_tubes": 250})
    fouled = run(geometry={"n_tubes": 250}, fouling={"R_f_coolant": 0.0005})
    assert fouled.get("area_required") > base.get("area_required")
    longer = run(geometry={"L": 6.0})
    assert longer.get("n_tubes") < run().get("n_tubes")
    more_duty = run(geometry={"n_tubes": 250}, vapour={"m_dot": 6.0})
    assert more_duty.get("area_required") > base.get("area_required")
    subcool = run(geometry={"n_tubes": 250}, vapour={"T_out": 45.0})
    assert subcool.get("A_sub") > 0.0
    assert subcool.get("Q") > base.get("Q")
    warmer = run(geometry={"n_tubes": 250}, coolant={"T_in": 30.0, "T_out": 45.0})
    assert warmer.get("area_required") > base.get("area_required")


def test_suction_is_a_choice_and_raises_the_coefficient():
    off = run(geometry={"n_tubes": 60})
    on = run(geometry={"n_tubes": 60}, solver={"suction": "on"})
    assert on.get("alpha_condensing_first") > off.get("alpha_condensing_first")
    assert on.get("area_required") < off.get("area_required")
    assert off.get("area_other_suction") == pytest.approx(on.get("area_required"), rel=1e-9)
    assert on.get("area_other_suction") == pytest.approx(off.get("area_required"), rel=1e-9)
    assert any("OFF, the default" in a for a in off.assumptions)
    assert any("ON by the user" in a for a in on.assumptions)


def test_boyko_kruzhilin_and_chen_are_recorded():
    r = run(geometry={"n_tubes": 60})
    assert r.get("area_boyko_kruzhilin") == pytest.approx(r.get("area_required"), rel=0.15)
    assert r.get("area_chen") > r.get("area_required")  # Chen runs lowest in alpha, library table
    assert check_state(r, "Boyko-Kruzhilin") is CheckState.PASS


def test_pressure_drop_uses_the_library_and_recovery_is_not_credited():
    r = run(geometry={"n_tubes": 100})
    assert r.get("dp_condensing_acceleration") < 0.0
    assert r.get("dp_condensing_gravity") < 0.0
    assert r.get("dp_condensing") == pytest.approx(r.get("dp_condensing_friction"), rel=1e-9)
    net = r.get("dp_condensing_friction") + r.get("dp_condensing_acceleration") + r.get("dp_condensing_gravity")
    assert r.get("dp_condensing_net") == pytest.approx(net, rel=1e-9)


# ---------------------------------------------------------------------------------- geometry rules
def test_tema_span_equations_against_hsb_table_6_1():
    assert vst.max_unsupported_span(19e-3, MATERIAL_A) == pytest.approx(1.520, abs=1e-9)  # Table 6.1: 1520
    assert vst.max_unsupported_span(19.05e-3, MATERIAL_A) / 2 == pytest.approx(0.762, abs=0.002)  # Lecture 9: 762 mm segmental
    assert vst.max_unsupported_span(19e-3, MATERIAL_B) == pytest.approx(1.310, abs=1e-9)  # Table 6.1 prints 1321; HEDH's fit is within 1%
    assert vst.max_unsupported_span(25e-3, MATERIAL_A) == pytest.approx(1.880, rel=0.03)
    assert vst.max_unsupported_span(100e-3, MATERIAL_A) is None


def test_geometry_checks_fire():
    r = run(geometry={"pitch": 30.0e-3, "n_tubes": 250})
    assert check_state(r, "Pitch at least") is CheckState.FAIL
    assert check_state(r, "Lane between tubes") is CheckState.FAIL
    r = run(geometry={"baffle_spacing": 0.03, "n_tubes": 250})
    assert check_state(r, "Baffle spacing at least") is CheckState.FAIL
    r = run(geometry={"baffle_spacing": 1.0, "n_tubes": 250})
    assert check_state(r, "Baffle spacing at most") is CheckState.FAIL
    assert check_state(r, "Number of baffles") is CheckState.FAIL
    r = run(geometry={"baffle_cut": 10.0, "n_tubes": 250})
    assert check_state(r, "Baffle cut in the range") is CheckState.FAIL


def test_baffle_spacing_follows_the_rules():
    i = inputs()
    assert vst._baffle_spacing(i, 0.05) == pytest.approx(50.8e-3)  # 50.8 mm floor
    assert vst._baffle_spacing(i, 1.0) == pytest.approx(0.4)
    assert vst._baffle_spacing(i, 10.0) == pytest.approx(vst.max_unsupported_span(25.4e-3, MATERIAL_A) / 2)


def test_window_angle_limited_when_the_cut_misses_the_tube_field():
    s = vst.bell_delaware(5.0, 1000.0, 4180.0, 1e-3, 0.6, 0.20, 0.025, 0.03125, LAYOUT_30, 25.0, 0.1, 3.0, 8, 0.14, 0.8e-3)
    assert s.theta_ctl == 0.0
    assert s.F_w == 0.0
    assert s.N_tcw == 0.0


def test_laminar_branches():
    s = vst.bell_delaware(5.0, 900.0, 2000.0, 0.5, 0.13, 0.5, 0.0254, 0.03175, LAYOUT_90, 25.0, 0.2, 4.0, 150, 0.015, 0.8e-3)
    assert s.Re < 100.0
    assert s.J_r < 1.0
    assert s.J_b == pytest.approx(math.exp(-1.35 * s.F_sbp))
    assert s.R_b == pytest.approx(math.exp(-4.5 * s.F_sbp))
    hand = s.N_b * (26.0 * s.m_w * 0.5 / 900.0 * (s.N_tcw / (0.03175 - 0.0254) + 0.2 / s.D_w**2) + s.m_w**2 / 900.0) * s.R_l
    assert s.dp_w == pytest.approx(hand)


def test_sealing_strips_limit_at_one_half():
    a = vst.bell_delaware(5.0, 1000.0, 4180.0, 1e-3, 0.6, 0.5, 0.0254, 0.03175, LAYOUT_90, 25.0, 0.2, 4.0, 150, 0.05, 0.8e-3, N_ss=0)
    b = vst.bell_delaware(5.0, 1000.0, 4180.0, 1e-3, 0.6, 0.5, 0.0254, 0.03175, LAYOUT_90, 25.0, 0.2, 4.0, 150, 0.05, 0.8e-3, N_ss=100)
    assert b.J_b == 1.0 and b.R_b == 1.0
    assert a.J_b < 1.0 and a.R_b < 1.0


# ---------------------------------------------------------------------------------- shell-side variant
def shell_inputs(**kw):
    defaults = dict(variant=SHELL_SIDE, geometry={"L_tb": None, "baffle_spacing": 0.8, "shell_fouling": "no"}, velocity={})
    return inputs(**{**defaults, **kw})


def test_shell_side_runs_without_baffle_data():
    r = run_design(KEY, shell_inputs())
    assert r.status in (Status.OK, Status.OK_WITH_WARNINGS), r.message
    assert r.get("L_required") <= 4.0
    assert any("unbaffled" in a for a in r.assumptions)
    assert not any(s.symbol == "S_m" for s in r.steps)


def test_shell_side_matches_the_mean_film_coefficient_at_uniform_conditions():
    # Dominant film, nearly uniform dT: A from the steps must equal Q/(alpha_mean dT) with HSB's mean for the exit film, gravity_mean (HEDH Eq. 18 links them).
    i = shell_inputs(geometry={"L_tb": None, "n_tubes": 400, "k_wall": 1e6, "t_wall": 0.2e-3}, coolant={"T_in": 25.0, "T_out": 25.5})
    r = run_design(KEY, i)
    p = i.condensing
    re_l = 4.0 * i.vapour.m_dot / (400 * math.pi * 25.4e-3 * p.mu_l)
    alpha = gravity_mean(p, re_l)
    dT = T_SAT - 25.25
    h_c = r.get("h_coolant")
    # the wall and the coolant remain, so include them; thin wall, and D_o/D_i near 1
    u = 1.0 / (1.0 / alpha + (25.4 / 25.0) / h_c)
    assert r.get("area_required") == pytest.approx(vst._duty(i).Q / (u * dT), rel=0.03)


def test_gravity_steps_reduce_to_nusselt_when_wave_free():
    # Exit film Re below 30: the mean must be HSB Eq. 16.14 with the dT that gives this film flow by Eq. 16.15. Printed 1.47 and 0.943 differ in rounding only.
    p = pentane_hsb_example_16_2()
    gamma_exit = 0.001
    mean = gravity_mean(p, 4 * gamma_exit / p.mu_l)
    dT = gamma_exit * p.h_fg / (mean * 1.0)
    assert nusselt_mean(p, dT, 1.0) == pytest.approx(mean, rel=0.01)
    assert gravity_local(p, 10.0) > gravity_local(p, 25.0)


def test_shell_side_coolant_in_tubes():
    i = shell_inputs()
    r = run_design(KEY, i)
    t = vst.tube_side_coolant(vst._duty(i).m_c, r.get("n_tubes"), 25.4e-3 - 2 * 2.11e-3, 4.0, i.coolant_props.rho, i.coolant_props.cp, i.coolant_props.mu, i.coolant_props.k)
    assert r.get("h_coolant") == pytest.approx(t.h)
    assert r.get("dp_coolant") == pytest.approx(t.dp)
    assert t.dp == pytest.approx(t.f_darcy * (4.0 / (25.4e-3 - 2 * 2.11e-3)) * i.coolant_props.rho * t.u**2 / 2 + 2.5 * i.coolant_props.rho * t.u**2 / 2)


def test_shell_side_subcooling_allowed_superheat_refused():
    r = run_design(KEY, shell_inputs(vapour={"T_out": 50.0}))
    assert r.get("A_sub") > 0.0
    r = run_design(KEY, shell_inputs(vapour={"T_in": 90.0}))
    assert r.status is Status.INVALID_INPUT


# ---------------------------------------------------------------------------------- reflux variant
def reflux_inputs(**kw):
    defaults = dict(variant=REFLUX, geometry={"L": 2.5, "D_o": 38e-3, "t_wall": 2.5e-3, "pitch": 47.5e-3, "layout": LAYOUT_30, "shell_fouling": "no"}, vapour={"m_dot": 2.0})
    return inputs(**{**defaults, **kw})


def test_reflux_reports_both_limits_and_stays_below_75_percent():
    r = run_design(KEY, reflux_inputs())
    assert r.status in (Status.OK, Status.OK_WITH_WARNINGS), r.message
    hsb, hedh = r.get("G_flood_hsb"), r.get("G_flood_hedh")
    assert hsb != hedh
    # pentane: HEDH's surface tension factor makes its limit much lower (condensation_library.md: 63% below at 11.3 mN/m)
    assert hedh == pytest.approx(hsb * (0.53 * (saturated("n_pentane", T_SAT).sigma / 0.072) ** 0.5) / 0.75**2, rel=1e-9)
    assert hedh < hsb
    assert r.get("G_condensing") <= 0.75 * min(hsb, hedh) * (1 + 1e-9)
    assert check_state(r, "Below 75%") is CheckState.PASS
    assert r.get("G_flood_lecture") == pytest.approx(hsb * (0.8 / 0.75) ** 2, rel=1e-9)
    assert any("English" in g and "unit" in g for g in r.source_gaps)


def test_reflux_flooding_sets_the_count_when_the_tube_is_long_enough():
    i = reflux_inputs(geometry={"L": 2.5, "D_o": 25.4e-3, "t_wall": 2.11e-3, "pitch": 31.75e-3, "layout": LAYOUT_30, "shell_fouling": "no"}, vapour={"m_dot": 5.0})
    r = run_design(KEY, i)
    assert r.get("sized_by").startswith("flooding")
    d_i = 25.4e-3 - 2 * 2.11e-3
    assert r.get("n_tubes") == math.ceil(5.0 / (0.75 * r.get("G_flood_hedh") * math.pi * d_i**2 / 4))
    assert r.get("L_required") < 2.5


def test_reflux_length_sets_the_count_when_the_tube_is_short():
    i = reflux_inputs(geometry={"L": 0.2, "D_o": 25.4e-3, "t_wall": 2.11e-3, "pitch": 31.75e-3, "layout": LAYOUT_30, "shell_fouling": "no"}, vapour={"m_dot": 5.0})
    r = run_design(KEY, i)
    assert r.get("sized_by") == "required length"
    assert r.get("L_required") <= 0.2 * (1 + 1e-9)
    assert check_state(r, "Tube length in the usual range") is CheckState.FAIL


def test_reflux_film_is_thin_at_the_top_and_thick_at_the_bottom():
    i = reflux_inputs(geometry={"L": 2.5, "D_o": 38e-3, "t_wall": 2.5e-3, "pitch": 47.5e-3, "layout": LAYOUT_30, "shell_fouling": "no", "n_tubes": 120})
    d = vst._duty(i)
    ev = vst._evaluate(i, d, 120)
    assert ev.steps[0].x > ev.steps[-1].x
    assert ev.steps[0].gamma > ev.steps[-1].gamma  # bottom, where the vapour enters, has the most condensate
    assert ev.steps[0].t1 == pytest.approx(25.0)  # HSB Fig. 17.2: coolant enters at the vapour inlet end
    assert ev.steps[-1].t2 == pytest.approx(40.0)
    e9 = vst._evaluate(i, d, 120, alpha_factor=0.9)
    assert e9.A_req > ev.A_req


@pytest.mark.parametrize("change, field", [
    (dict(vapour={"m_dot": 2.0, "x_out": 0.1}), "vapour.x_out"),
    (dict(vapour={"m_dot": 2.0, "T_out": 50.0}), "vapour.T_out"),
    (dict(vapour={"m_dot": 2.0, "T_in": 90.0}), "vapour.T_in"),
])
def test_reflux_refusals(change, field):
    r = run_design(KEY, reflux_inputs(**change))
    assert r.status is Status.INVALID_INPUT
    assert field in {issue.field for issue in r.input_issues}


# ---------------------------------------------------------------------------------- mixtures
def curve_text(t_in, t_out, q_total, m_v, rows=5, m_out_fraction=0.02):
    parts = []
    for k in range(rows):
        s = k / (rows - 1)
        parts.append(f"{q_total * s}, {t_in + (t_out - t_in) * s}, {m_v * (1 - (1 - m_out_fraction) * s)}")
    return "; ".join(parts)


def test_sbg_constant_temperature_limit_matches_the_pure_path():
    # Z = 0 when T_E does not change, so SBG reduces to integral dQ/(U' dT). Same quality range as the pure run. 2%: U' is one mean, and five increments only.
    q_total = 5.0 * saturated("n_pentane", T_SAT).h_fg * 0.98
    pure = run(geometry={"n_tubes": 250}, vapour={"x_out": 0.02})
    sbg = run(geometry={"n_tubes": 250}, vapour={"condensation_curve": curve_text(T_SAT, T_SAT, q_total, 5.0)})
    assert sbg.status in (Status.OK, Status.OK_WITH_WARNINGS), sbg.message
    assert sbg.get("Q") == pytest.approx(q_total)
    assert sbg.get("A_con") == pytest.approx(pure.get("A_con"), rel=0.02)
    assert any("Silver-Bell-Ghaly" in a for a in sbg.assumptions)


def test_sbg_gliding_mixture_needs_more_area():
    q_total = 5.0 * saturated("n_pentane", T_SAT).h_fg * 0.98
    flat = run(geometry={"n_tubes": 250}, vapour={"condensation_curve": curve_text(T_SAT, T_SAT, q_total, 5.0)})
    glide = run(geometry={"n_tubes": 250}, vapour={"condensation_curve": curve_text(T_SAT, T_SAT - 12.0, q_total, 5.0)})
    assert glide.get("A_con") > flat.get("A_con")


@pytest.mark.parametrize("text", ["1, 2", "0, 60, 5; 100, 50, 2", "10, 60, 5; 20, 50, 4; 30, 45, 3", "0, 60, 5; 10, 65, 4; 20, 60, 3", "0, 60, 5; 10, 55, 4; 20, 50, 0"])
def test_bad_curves_are_refused(text):
    r = run(vapour={"condensation_curve": text})
    assert r.status is Status.INVALID_INPUT
    assert "vapour.condensation_curve" in {issue.field for issue in r.input_issues}


def test_curve_refused_in_other_variants_and_with_zones():
    good = curve_text(T_SAT, T_SAT - 5.0, 1e6, 5.0)
    assert run_design(KEY, shell_inputs(vapour={"condensation_curve": good})).status is Status.INVALID_INPUT
    assert run(vapour={"condensation_curve": good, "T_out": 40.0}).status is Status.INVALID_INPUT


# ---------------------------------------------------------------------------------- inputs and form
@pytest.mark.parametrize("path", ["geometry.D_o", "geometry.L", "geometry.pitch", "geometry.L_bb", "vapour.m_dot", "vapour.pressure", "coolant.T_in", "condensing.T_sat", "condensing.h_fg", "coolant_props.rho"])
def test_missing_required_inputs(path):
    i = inputs()
    owner, name = path.split(".")
    setattr(getattr(i, owner), name, None)
    r = run_design(KEY, i)
    assert r.status is Status.INVALID_INPUT
    assert path in {issue.field for issue in r.input_issues}


@pytest.mark.parametrize("change, field", [
    (dict(geometry={"pitch": 20e-3}), "geometry.pitch"),
    (dict(geometry={"t_wall": 13e-3}), "geometry.t_wall"),
    (dict(geometry={"L_tb": None}), "geometry.L_tb"),
    (dict(geometry={"n_tubes": 10.5}), "geometry.n_tubes"),
    (dict(vapour={"x_out": 1.0}), "vapour.x_out"),
    (dict(vapour={"T_in": 40.0}), "vapour.T_in"),
    (dict(vapour={"T_out": 80.0}), "vapour.T_out"),
    (dict(vapour={"T_out": 20.0}), "vapour.T_out"),
    (dict(vapour={"T_out": 50.0, "x_out": 0.1}), "vapour.x_out"),
    (dict(vapour={"T_in": 90.0, "cp_v": None}), "vapour.cp_v"),
    (dict(vapour={"T_in": 90.0, "x_in": 0.9}), "vapour.x_in"),
    (dict(coolant={"T_in": 70.0}), "coolant.T_in"),
    (dict(coolant={"T_out": 20.0}), "coolant.T_out"),
    (dict(coolant={"T_out": 70.0}), "coolant.T_out"),
    (dict(coolant={"T_out": 40.0, "m_dot": 10.0}), "coolant.T_out"),
    (dict(coolant={"T_out": None, "m_dot": 1.0}), "coolant.m_dot"),
])
def test_invalid_values(change, field):
    r = run(**change)
    assert r.status is Status.INVALID_INPUT, r.message
    assert field in {issue.field for issue in r.input_issues}


def test_surface_tension_needed():
    i = inputs()
    i.condensing.sigma = None
    r = run_design(KEY, i)
    assert r.status is Status.INVALID_INPUT
    assert "condensing.sigma" in {issue.field for issue in r.input_issues}


def test_coolant_flow_gives_outlet_temperature():
    r = run(coolant={"T_out": None, "m_dot": 30.0}, geometry={"n_tubes": 250})
    assert r.get("T_c_out") == pytest.approx(25.0 + r.get("Q") / (30.0 * water_liquid(32.5).cp))


def test_relevant_fields_hide_what_a_variant_does_not_use():
    module = MODULES[KEY]
    downflow = module.relevant_fields(inputs())
    assert {"vapour.T_in", "solver.suction", "geometry.L_tb", "vapour.condensation_curve"} <= downflow
    assert "solver.flooding_f2" not in downflow
    shell = module.relevant_fields(shell_inputs())
    assert not {"vapour.T_in", "solver.suction", "geometry.L_tb", "geometry.baffle_cut"} & shell
    assert "geometry.baffle_spacing" in shell and "vapour.T_out" in shell
    reflux = module.relevant_fields(reflux_inputs())
    assert "solver.flooding_f2" in reflux
    assert not {"vapour.T_in", "vapour.T_out", "solver.suction"} & reflux
    # An empty form must not crash the hook.
    assert module.relevant_fields(module.new_inputs())
