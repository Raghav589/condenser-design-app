"""Direct-contact condenser.

Worked examples exist for four types: HSB Examples 22.1 (drops), 22.2 (sheet), 22.3 (packed bed) and 22.4 (bubble).
There is no worked example for jets, fans, baffle columns, pools, or for the balance with superheat or a vent:
those are tested by hand-checked values, limiting cases, closure and trends, and are not benchmarks.
"""
import math

import pytest

from condenser_app.modules import MODULES, direct_contact as dc
from condenser_app.modules.direct_contact import (
    BubbleData, CoolantStream, DirectContactInputs, NonCondensable, PackedData, SheetJetData, SprayData, VapourStream,
)
from condenser_app.schemas import CheckState, Status
from condenser_app.service import run_design


# Properties of HSB Example 22.1: steam at atmospheric pressure, water at 92 degC.
def steam(**kw):
    return VapourStream(**{**dict(T_sat=100.0, h_fg=2.257e6, rho_v=0.596, mu_v=1.20e-5, cp_v=1930.0), **kw})


def water(**kw):
    return CoolantStream(**{**dict(T_in=85.0, rho=963.0, cp=4209.0, k=0.677), **kw})


def spray_inputs(**kw):
    defaults = dict(
        condenser_type="spray", coolant_relation="same_fluid", vapour=steam(), coolant=water(m_coolant=1.0, utilization=0.95),
        spray=SprayData(drop_diameter=1e-3, injection_velocity=5.0),
    )
    return DirectContactInputs(**{**defaults, **kw})


def sheet_inputs(kind="sheet", **kw):
    defaults = dict(
        condenser_type=kind, coolant_relation="same_fluid", vapour=steam(), coolant=water(m_coolant=1.0, utilization=0.95),
        film=SheetJetData(liquid_dimension=1e-3, liquid_velocity=1.0),
    )
    return DirectContactInputs(**{**defaults, **kw})


def packed_inputs(**kw):
    defaults = dict(
        condenser_type="packed_bed", coolant_relation="same_fluid", column_area=1.0,
        vapour=steam(m_vapour=0.1), coolant=water(m_coolant=5.0),
        packed=PackedData(packing_size=3.2e-2, voidage=0.4),
    )
    return DirectContactInputs(**{**defaults, **kw})


def bubble_inputs(**kw):
    defaults = dict(
        condenser_type="bubble_column", coolant_relation="same_fluid",
        vapour=steam(m_vapour=0.1), coolant=water(m_coolant=5.0),
        bubble=BubbleData(bubble_diameter=1e-2, T_pool=92.0),
    )
    return DirectContactInputs(**{**defaults, **kw})


def pool_inputs(**kw):
    defaults = dict(condenser_type="pool", coolant_relation="same_fluid", vapour=steam(m_vapour=1.0), coolant=CoolantStream(T_in=85.0, T_out=95.0, rho=963.0, cp=4209.0))
    return DirectContactInputs(**{**defaults, **kw})


def check_state(result, name):
    return next(c.state for c in result.checks if c.name == name)


# --------------------------------------------------------------------------------- benchmark: HSB Example 22.1
def test_registered_and_example_runs():
    assert "direct_contact" in MODULES
    result = run_design("direct_contact", MODULES["direct_contact"].example_inputs())
    assert result.status is Status.OK_WITH_WARNINGS
    assert result.get("column_height") > 0


def test_example_22_1_balance_and_groups():
    result = run_design("direct_contact", spray_inputs())
    # Book prints 3 s.f.
    assert result.get("jakob") == pytest.approx(35.7, rel=0.005)
    # Book: M_G = 0.95/35.7 = 0.0266, printed 0.027 (2 s.f.), hence 3%.
    assert result.get("m_condensed") == pytest.approx(0.027, rel=0.03)
    assert result.get("utilization") == pytest.approx(0.95)


def test_example_22_1_book_intermediates_from_its_own_time():
    # The book takes t = 0.374 s from a chart reading and rounds to t = 0.37 in the velocity. Gamma is printed to 2 s.f., hence 1%.
    gamma = dc.gamma_drag(0.596, 1.20e-5, 963.0)
    assert gamma == pytest.approx(7.0e-8, rel=0.01)
    v_end = dc.drop_velocity(5.0, 7.0e-8, 0.37, 1e-3)
    assert v_end == pytest.approx(4.0, rel=0.01)
    assert dc.drop_distance(5.0, 4.0, 7.0e-8, 1e-3) == pytest.approx(1.71, rel=0.005)


def test_example_22_1_fourier_and_height_differ_from_the_chart_reading_for_a_known_reason():
    # The book reads log10 Fo = -0.6 (Fo = 0.25, t = 0.374 s, H = 1.71 m) from its figure. Boxed Eq. 22.9 itself gives Fo = 0.236.
    # The exact series of HEDH Eq. (12) gives 0.253, which is what the book's reading agrees with.
    result = run_design("direct_contact", spray_inputs())
    assert result.get("fourier") == pytest.approx(0.2359, rel=0.001)
    assert result.get("residence_time") == pytest.approx(0.374, rel=0.06)
    assert result.get("residence_time_series") == pytest.approx(0.374, rel=0.02)
    assert result.get("height_no_gravity") == pytest.approx(1.71, rel=0.06)


def test_gravity_is_not_negligible_in_example_22_1():
    # At 5 m/s the drag deceleration is about 3 m/s2 against 9.8 m/s2 for gravity, so the book's gravity-free height is short.
    result = run_design("direct_contact", spray_inputs())
    assert check_state(result, "Gravity term below drag term") is CheckState.FAIL
    assert result.get("height_with_gravity") > result.get("height_no_gravity")
    assert result.get("column_height") == result.get("height_with_gravity")
    assert any("Gravity is not negligible" in w for w in result.warnings)


def test_gravity_integration_without_gravity_recovers_the_closed_form(monkeypatch):
    # Hand check: with g = 0 Eq. 22.1 is Eq. 22.2. The book rounds 1/(0.84 x 20.25) = 0.0588 to 0.06, so allow 3%.
    monkeypatch.setattr(dc, "GRAVITY", 0.0)
    t, d = 0.37, 1e-3
    x, v = dc.drop_path_with_gravity(5.0, t, d, 0.596, 1.20e-5, 963.0, 2000)
    gamma = dc.gamma_drag(0.596, 1.20e-5, 963.0)
    assert v == pytest.approx(dc.drop_velocity(5.0, gamma, t, d), rel=0.002)
    assert x == pytest.approx(dc.drop_distance(5.0, v, gamma, d), rel=0.03)


# --------------------------------------------------------------------------------- benchmark: HSB Example 22.2
def test_example_22_2_sheet():
    result = run_design("direct_contact", sheet_inputs())
    # Book prints 3 s.f.
    assert result.get("graetz") == pytest.approx(14.2, rel=0.005)
    assert result.get("column_height") == pytest.approx(1.69, rel=0.005)


# --------------------------------------------------------------------------------- benchmark: HSB Example 22.3
def test_example_22_3_packed_bed():
    result = run_design("direct_contact", packed_inputs())
    # Book prints 3 s.f.
    assert result.get("packing_group") == pytest.approx(3.6)
    assert result.get("capacity_ratio") == pytest.approx(109, rel=0.005)
    assert result.get("height_ratio") == pytest.approx(6.43, rel=0.005)
    assert result.get("column_height") == pytest.approx(0.206, rel=0.005)
    assert result.get("duty") == pytest.approx(2.26e5, rel=0.005)
    assert result.get("T_coolant_out") == pytest.approx(95.7, abs=0.05)
    # The book rounds T_avg to 90.4 (exact 90.36) before dividing by a small difference, so these move by about 0.5%: allow 1.5%.
    assert result.get("stanton") == pytest.approx(5.24, rel=0.015)
    assert result.get("alpha") == pytest.approx(1011, rel=0.015)
    assert result.get("packing_surface") == pytest.approx(23.2, rel=0.005)
    assert check_state(result, "Correlation closes on the latent duty") is CheckState.PASS


# --------------------------------------------------------------------------------- benchmark: HSB Example 22.4
def test_example_22_4_with_the_printed_jakob_number():
    # Book's own intermediates. Printed to 3 s.f.
    assert dc.bubble_collapse_fourier(27.0, 2.69e4) == pytest.approx(2.67e-4, rel=0.005)
    assert 2.67e-4 * (1e-2) ** 2 / (4 * 1.67e-7) == pytest.approx(0.040, rel=0.005)
    assert dc.stationary_bubble_fourier(27.0) == pytest.approx(6.10e-3, rel=0.005)
    assert 6.10e-3 * (1e-2) ** 2 / (4 * 1.67e-7) == pytest.approx(0.91, rel=0.005)


def test_example_22_4_source_error_in_specific_heat():
    # The book's Ja' = 27.0 uses c_p = 4709, a misprint of the 4209 it says to use. With 4209, Ja' = 24.1 and t = 0.045 s:
    # that is the value the book itself quotes in Sec. 22.4.5. So the benchmark here is 0.045 s, not 0.040 s.
    result = run_design("direct_contact", bubble_inputs())
    assert result.get("jakob_bubble") == pytest.approx(24.1, rel=0.005)
    assert result.get("rise_velocity") == pytest.approx(0.45, rel=0.015)  # book prints 2 s.f.; exact 0.455
    assert result.get("peclet") == pytest.approx(2.69e4, rel=0.015)  # follows the velocity
    assert result.get("collapse_time") == pytest.approx(0.045, rel=0.02)
    assert result.get("collapse_height_single") == pytest.approx(0.02, rel=0.05)  # book: "about 2 cm"
    # Stationary bubble with the corrected Ja': the book's 0.91 s scaled by (27.0/24.1)^2.
    assert result.get("collapse_time_stationary") == pytest.approx(0.91 * (27.0 / 24.1) ** 2, rel=0.01)


def test_stepped_bubble_matches_the_closed_form_and_is_slower_than_constant_velocity():
    # Hand derivation, not a benchmark: v = c sqrt(D) in Eq. 22.41 gives dD/dt = -A D^(-1/4), so t = D0^1.25/(1.25 A), z = (c/A) D0^1.75/1.75.
    result = run_design("direct_contact", bubble_inputs())
    d0, ja = 1e-2, result.get("jakob_bubble")
    kappa = 0.677 / (963.0 * 4209.0)
    c_v = math.sqrt(4.0 / 3.0 * 9.81 * (963.0 - 0.596) / (0.5 * 963.0))
    big_a = 4.0 / math.sqrt(math.pi) * math.sqrt(kappa * c_v) * ja
    assert result.get("collapse_time_stepped") == pytest.approx(d0**1.25 / (1.25 * big_a), rel=0.002)
    assert result.get("collapse_height_stepped") == pytest.approx(c_v / big_a * d0**1.75 / 1.75, rel=0.002)
    # The book says its constant-velocity estimate may be too small.
    assert result.get("collapse_time_stepped") > result.get("collapse_time")
    assert any(c.name == "Bubble collapse in steps" and c.rows for c in result.convergence)


# --------------------------------------------------------------------------------- balances: hand checks, not benchmarks
def test_coolant_flow_for_saturated_vapour_hand_check():
    # 1 kg/s steam, water 85 to 95 degC: M_L = 2.257e6/(4209 x 10) = 53.62 kg/s.
    result = run_design("direct_contact", pool_inputs())
    assert result.get("m_coolant") == pytest.approx(2.257e6 / (4209.0 * 10.0))
    assert result.get("coolant_ratio") == pytest.approx(2.257e6 / (4209.0 * 10.0))
    # HEDH Eq. (3)
    assert result.get("volume_ratio") == pytest.approx(0.596 * 2.257e6 / (963.0 * 4209.0 * 10.0))


def test_superheat_adds_to_the_duty_hand_check():
    # 1 x 1930 x 20 = 38600 W on top of 2.257e6 W.
    inputs = pool_inputs(vapour=steam(m_vapour=1.0, T_in=120.0))
    result = run_design("direct_contact", inputs)
    assert result.get("duty") == pytest.approx(2.257e6 + 38600.0)


def test_vent_balance_hand_check():
    # 0.01 kg/s air (M 29) with steam (M 18). P = 10 kPa, p_v at vent = 5 kPa: M_vv = 0.01(18/29)(5/5) = 0.0062069 kg/s.
    gas = NonCondensable(m_gas=0.01, M_gas=29.0, M_vapour=18.0, p_v_vent=5000.0, T_vent=33.0, cp_gas=1005.0)
    inputs = DirectContactInputs(
        condenser_type="pool", coolant_relation="same_fluid", pressure=1e4, gas=gas,
        vapour=steam(T_sat=45.0, h_fg=2.39e6, rho_v=0.065, m_vapour=1.0),
        coolant=CoolantStream(T_in=30.0, T_out=40.0, rho=995.0, cp=4180.0),
    )
    result = run_design("direct_contact", inputs)
    m_vv = 0.01 * 18.0 / 29.0
    assert result.get("m_vent_vapour") == pytest.approx(m_vv)
    assert result.get("m_condensed") == pytest.approx(1.0 - m_vv)
    expected = (1.0 - m_vv) * 2.39e6 + m_vv * 1930.0 * (45.0 - 33.0) + 0.01 * 1005.0 * (45.0 - 33.0)
    assert result.get("duty") == pytest.approx(expected)
    y = (1.0 / 18.0) / (1.0 / 18.0 + 0.01 / 29.0)
    assert result.get("p_vapour_inlet") == pytest.approx(y * 1e4)
    assert result.get("m_vent_total") == pytest.approx(0.01 + m_vv)


@pytest.mark.parametrize("given", ["coolant", "vapour", "outlet"])
def test_three_ways_of_posing_the_balance_agree(given):
    gas = NonCondensable(m_gas=0.01, M_gas=29.0, M_vapour=18.0, p_v_vent=5000.0, T_vent=33.0, cp_gas=1005.0)
    base = dict(T_sat=45.0, h_fg=2.39e6, rho_v=0.065, T_in=60.0)
    ref = dc._solve_balance(DirectContactInputs(condenser_type="pool", coolant_relation="same_fluid", pressure=1e4, gas=gas, vapour=steam(m_vapour=1.0, **base), coolant=CoolantStream(T_in=30.0, T_out=40.0, rho=995.0, cp=4180.0)))
    if given == "coolant":
        vap, cool = steam(m_vapour=1.0, **base), CoolantStream(T_in=30.0, T_out=40.0, rho=995.0, cp=4180.0)
    elif given == "vapour":
        vap, cool = steam(**base), CoolantStream(T_in=30.0, T_out=40.0, m_coolant=ref.m_coolant, rho=995.0, cp=4180.0)
    else:
        vap, cool = steam(m_vapour=1.0, **base), CoolantStream(T_in=30.0, m_coolant=ref.m_coolant, rho=995.0, cp=4180.0)
    result = run_design("direct_contact", DirectContactInputs(condenser_type="pool", coolant_relation="same_fluid", pressure=1e4, gas=gas, vapour=vap, coolant=cool))
    assert result.get("m_vapour") == pytest.approx(1.0)
    assert result.get("m_coolant") == pytest.approx(ref.m_coolant)
    assert result.get("T_coolant_out") == pytest.approx(40.0)
    assert check_state(result, "Energy balance closure") is CheckState.PASS


def test_mixed_outlet_closes_on_the_first_law():
    # Example 22.3 numbers: (5 x 85 + 0.1 x 100)/5.1 + 0.1 x 2.257e6/(5.1 x 4209) = 95.81 degC, above the coolant-only 95.72 degC.
    result = run_design("direct_contact", packed_inputs())
    assert result.get("T_mixed") == pytest.approx((5 * 85 + 0.1 * 100) / 5.1 + 0.1 * 2.257e6 / (5.1 * 4209.0))
    assert result.get("T_coolant_out") < result.get("T_mixed") < 100.0
    assert result.get("m_outlet") == pytest.approx(5.1)
    assert result.get("condensate_fraction") == pytest.approx(0.1 / 5.1)


def test_limit_of_a_large_coolant_flow():
    result = run_design("direct_contact", pool_inputs(coolant=CoolantStream(T_in=85.0, m_coolant=1e7, rho=963.0, cp=4209.0)))
    assert result.get("T_coolant_out") == pytest.approx(85.0, abs=1e-2)


# --------------------------------------------------------------------------------- spray: trends and side checks
def test_height_rises_with_utilization_and_drop_size():
    def height(theta, d):
        return run_design("direct_contact", spray_inputs(
            coolant=water(m_coolant=1.0, utilization=theta), spray=SprayData(drop_diameter=d, injection_velocity=5.0),
        )).get("column_height")
    assert height(0.8, 1e-3) < height(0.9, 1e-3) < height(0.97, 1e-3)
    assert height(0.9, 0.5e-3) < height(0.9, 1e-3)


def test_series_and_eq_22_9_are_close_at_moderate_utilization():
    fo = dc.drop_fourier_required(0.5)
    assert dc.drop_theta_series(fo) == pytest.approx(0.5, abs=0.02)


def test_terminal_velocity_iteration_matches_the_power_law_solution():
    # Hand check: with C_D = 27 Re^-0.84, V^1.16 = [4 g D (rho_L - rho_G)/(81 rho_G)] (D/nu)^0.84.
    inputs = spray_inputs(spray=SprayData(drop_diameter=1e-3, injection_velocity=5.0, smallest_drop=5e-4, velocity_fraction=0.5))
    result = run_design("direct_contact", inputs)
    d, nu = 5e-4, 1.20e-5 / 0.596
    expected = (4 * 9.81 * d * (963.0 - 0.596) / (81.0 * 0.596) * (d / nu) ** 0.84) ** (1 / 1.16)
    assert result.get("terminal_velocity") == pytest.approx(expected, rel=1e-3)
    # Area from the 50% limit and the 0.0266 kg/s of steam.
    assert result.get("column_area") == pytest.approx(result.get("m_vapour") / (0.596 * 0.5 * expected), rel=1e-3)
    info = next(c for c in result.convergence if c.name == "Terminal velocity of the smallest drop")
    assert info.converged and info.rows


def test_area_not_sized_without_a_fraction():
    result = run_design("direct_contact", spray_inputs())
    assert any("Vessel cross-section" in g for g in result.source_gaps)
    assert check_state(result, "Drop velocity much greater than vapour velocity") is CheckState.NOT_ASSESSED


def test_low_jakob_number_warns():
    # R-type data: low latent heat and a large subcooling.
    inputs = spray_inputs(vapour=steam(h_fg=2.0e5), coolant=water(m_coolant=1.0, T_in=40.0, T_out=95.0))
    result = run_design("direct_contact", inputs)
    assert result.get("jakob") < 20
    assert check_state(result, "Large Jakob number range") is CheckState.FAIL
    assert any("nonlinear" in w for w in result.warnings)
    assert any("Eqs (9)-(11)" in g for g in result.source_gaps)


# --------------------------------------------------------------------------------- jets, sheets, fans
def test_gz_order_matches_hedh_fig_12():
    # Fig. 12: to reach the same temperature the fan needs the longest length (lowest Gz), then the jet, then the sheet.
    assert dc.graetz_required("fan", 0.95)[0] < dc.graetz_required("jet", 0.95)[0] < dc.graetz_required("sheet", 0.95)[0]


def test_jet_hand_check():
    # 1 - theta = 0.6915 exp(-23.136/Gz) at theta = 0.95: Gz = 23.136/ln(0.6915/0.05) = 8.81. x = (2e-3)^2 (2)/(1.67e-7 x 8.81) = 5.44 m.
    inputs = sheet_inputs("jet", film=SheetJetData(liquid_dimension=2e-3, liquid_velocity=2.0))
    result = run_design("direct_contact", inputs)
    assert result.get("graetz") == pytest.approx(8.81, rel=0.003)
    assert result.get("column_height") == pytest.approx(5.44, rel=0.005)
    assert result.get("jet_count") == math.ceil(1.0 / (963.0 * math.pi / 4 * (2e-3) ** 2 * 2.0))


def test_fan_equation_hand_check():
    # Gz = (2 delta)^2 u/(K x) with delta = V/(phi x u) must return the Gz the length was built from.
    inputs = sheet_inputs("fan", film=SheetJetData(liquid_velocity=5.0, fan_angle=math.pi, fan_nozzles=4))
    result = run_design("direct_contact", inputs)
    x, delta = result.get("column_height"), result.get("fan_thickness_end")
    kappa = 0.677 / (963.0 * 4209.0)
    assert (2 * delta) ** 2 * 5.0 / (kappa * x) == pytest.approx(result.get("graetz"), rel=1e-9)
    assert 1.0 / (963.0 * 4) == pytest.approx(math.pi * x * delta * 5.0, rel=1e-9)


def test_sheet_height_scales_with_thickness_squared():
    # HSB Example 22.2: "the height required increases in proportion to the square of the sheet thickness".
    h1 = run_design("direct_contact", sheet_inputs(film=SheetJetData(liquid_dimension=1e-3, liquid_velocity=1.0))).get("column_height")
    h2 = run_design("direct_contact", sheet_inputs(film=SheetJetData(liquid_dimension=2e-3, liquid_velocity=1.0))).get("column_height")
    assert h2 / h1 == pytest.approx(4.0)


@pytest.mark.parametrize("kind", ["sheet", "jet", "fan"])
def test_asymptotes_meet_near_the_switching_utilization(kind):
    theta = dc.THETA_ASYMPTOTE_SWITCH
    above = dc.graetz_required(kind, theta + 1e-9)[0]
    below = dc.graetz_required(kind, theta - 1e-9)[0]
    assert above == pytest.approx(below, rel=0.05)


# --------------------------------------------------------------------------------- baffle column
def baffle_inputs(**film):
    params = {**dict(liquid_dimension=2e-3, liquid_velocity=1.0, baffle_spacing=0.3), **film}
    return DirectContactInputs(
        condenser_type="baffle_column", coolant_relation="same_fluid", vapour=steam(m_vapour=1.0),
        coolant=water(T_in=70.0, T_out=95.0), film=SheetJetData(**params),
    )


def test_baffle_count_matches_the_stage_formula():
    result = run_design("direct_contact", baffle_inputs())
    theta, theta_s = result.get("utilization"), result.get("stage_utilization")
    n = result.get("n_baffles_estimate")
    assert n == math.ceil(math.log(1 - theta) / math.log(1 - theta_s))
    assert result.get("column_height") == pytest.approx(n * 0.3)
    info = next(c for c in result.convergence if c.name == "Baffle stage march")
    assert info.rows[-1][1] >= 95.0 > info.rows[-2][1]
    assert any("No source gives a sizing procedure for a baffle column" in g for g in result.source_gaps)
    assert any("preliminary estimate" in w for w in result.warnings)


def test_baffle_count_falls_with_more_fall_height():
    few = run_design("direct_contact", baffle_inputs(baffle_spacing=1.0)).get("n_baffles_estimate")
    many = run_design("direct_contact", baffle_inputs(baffle_spacing=0.1)).get("n_baffles_estimate")
    assert few < many


def test_baffle_non_convergence_is_reported():
    result = run_design("direct_contact", baffle_inputs(baffle_spacing=1e-9, liquid_dimension=1e-2))
    assert result.status is Status.NOT_CONVERGED
    assert result.convergence and not result.convergence[0].converged
    assert "not reached" in result.message


# --------------------------------------------------------------------------------- packed bed side paths
def test_packed_area_from_velocity_limit():
    inputs = packed_inputs(column_area=None)
    result = run_design("direct_contact", inputs)
    assert result.get("column_area") == pytest.approx(0.1 / (0.596 * 4.0))
    assert result.get("vapour_velocity") == pytest.approx(4.0)


def test_packed_velocity_too_high_fails_the_check():
    result = run_design("direct_contact", packed_inputs(column_area=0.01))
    assert check_state(result, "Vapour does not strip the coolant") is CheckState.FAIL


def test_packed_hedh_eq_65_hand_check():
    inputs = packed_inputs(packed=PackedData(packing_size=3.2e-2, voidage=0.4, U_v=1e5))
    result = run_design("direct_contact", inputs)
    t_out = result.get("T_coolant_out")
    assert result.get("height_volumetric") == pytest.approx(5.0 * 4209.0 / (1e5 * 1.0) * math.log(15.0 / (100.0 - t_out)))


def test_packed_specific_area_route():
    # Da = D a must give the same bed as voidage when a = 6(1 - voidage)/D.
    a = 6 * (1 - 0.4) / 3.2e-2
    h_area = run_design("direct_contact", packed_inputs(packed=PackedData(packing_size=3.2e-2, specific_area=a))).get("column_height")
    h_void = run_design("direct_contact", packed_inputs()).get("column_height")
    assert h_area == pytest.approx(h_void)


# --------------------------------------------------------------------------------- bubble column side paths
def test_bubble_train_eq_87_hand_check():
    # R_0 = 5 mm, f = 26/s, L = 12.5 mm, so L/2R_0 = 1.25. Velocity given: 0.3 m/s. T_sat - T_in = 15 degC.
    bub = BubbleData(bubble_diameter=1e-2, T_pool=92.0, rise_velocity=0.3, release_frequency=26.0, nozzle_spacing=12.5e-3)
    result = run_design("direct_contact", bubble_inputs(bubble=bub))
    kappa = 0.677 / (963.0 * 4209.0)
    pe = 0.3 * 1e-2 / kappa
    ja_star = 4209.0 * 15.0 / 2.257e6
    ratio = 1.0415 * math.sqrt(pe) * (0.596 / 963.0) * (26.0 * 5e-3 / 0.3) ** 0.054 / (1.25**0.028 * math.sinh(ja_star))
    assert result.get("column_height") == pytest.approx(ratio * 5e-3)
    n = math.ceil(0.1 / (0.596 * 26.0 * math.pi / 6 * 1e-6))
    assert result.get("nozzle_count") == n
    assert result.get("column_area") == pytest.approx(n * (12.5e-3) ** 2)


def test_bubble_train_with_gas_needs_final_ratio():
    gas = NonCondensable(m_gas=0.001, M_gas=28.0, M_vapour=18.0, p_v_vent=9e4, T_vent=96.0, cp_gas=1040.0)
    base = BubbleData(bubble_diameter=1e-2, T_pool=92.0, rise_velocity=0.3, release_frequency=26.0, nozzle_spacing=12.5e-3)
    no_beta = run_design("direct_contact", bubble_inputs(gas=gas, pressure=1.01e5, bubble=base))
    assert any("No column height" in w for w in no_beta.warnings)
    with pytest.raises(KeyError):
        no_beta.get("column_height")
    with_beta = BubbleData(bubble_diameter=1e-2, T_pool=92.0, rise_velocity=0.3, release_frequency=26.0, nozzle_spacing=12.5e-3, beta_final=0.2)
    result = run_design("direct_contact", bubble_inputs(gas=gas, pressure=1.01e5, bubble=with_beta))
    assert result.get("column_height") > 0


def test_bubble_without_injector_data_gives_no_column_height():
    result = run_design("direct_contact", bubble_inputs())
    assert any("Column height" in g for g in result.source_gaps)
    with pytest.raises(KeyError):
        result.get("column_height")


def test_small_subcooling_outside_the_correlation():
    bub = BubbleData(bubble_diameter=1e-2, rise_velocity=0.3, release_frequency=26.0, nozzle_spacing=12.5e-3)
    inputs = bubble_inputs(coolant=water(T_in=99.0, T_out=99.5), bubble=bub, vapour=steam(m_vapour=0.01))
    result = run_design("direct_contact", inputs)
    assert check_state(result, "Subcooling range of the correlation") is CheckState.FAIL


def test_small_bubble_flagged():
    result = run_design("direct_contact", bubble_inputs(bubble=BubbleData(bubble_diameter=2e-3, T_pool=92.0)))
    assert check_state(result, "Constant-velocity bubble model") is CheckState.FAIL
    assert any("1.5 mm" in g for g in result.source_gaps)


def test_immiscible_bubble_needs_viscosity_and_lowers_the_rate():
    vap = steam(m_vapour=0.1, rho_cond=626.0, cp_cond=2300.0)
    missing = run_design("direct_contact", bubble_inputs(coolant_relation="immiscible", vapour=vap))
    assert missing.status is Status.INVALID_INPUT
    assert "coolant.mu" in {i.field for i in missing.input_issues}
    result = run_design("direct_contact", bubble_inputs(coolant_relation="immiscible", vapour=vap, coolant=water(m_coolant=5.0, mu=3e-4)))
    same = run_design("direct_contact", bubble_inputs())
    assert result.get("collapse_time") > same.get("collapse_time")


# --------------------------------------------------------------------------------- pool, barometric, warnings
def test_pool_is_balances_only_and_says_so():
    result = run_design("direct_contact", pool_inputs())
    assert result.get("balances_only") == "balances only"
    assert any("balances only" in g for g in result.source_gaps)
    with pytest.raises(KeyError):
        result.get("column_height")


def test_barometric_tail_pipe_about_ten_metres():
    # HSB: about 10 m. HEDH: about 34 ft. Hand check: 1.013e5 Pa over water is about 10.5 m at 963 kg/m3.
    inputs = pool_inputs(barometric="yes", p_atm=1.013e5, pressure=5e3, vapour=steam(T_sat=33.0, h_fg=2.42e6, rho_v=0.035, m_vapour=1.0), coolant=CoolantStream(T_in=20.0, T_out=28.0, rho=996.0, cp=4180.0))
    result = run_design("direct_contact", inputs)
    assert result.get("tail_pipe_height") == pytest.approx((1.013e5 - 5e3) / (996.0 * 9.81), rel=0.005)
    assert 9.0 < result.get("tail_pipe_height") < 10.5


def test_contamination_and_separation_warnings():
    same = run_design("direct_contact", pool_inputs())
    assert any("chemically clean" in w for w in same.warnings)
    assert any("separated only if immiscible" in w for w in same.warnings)
    miscible = run_design("direct_contact", pool_inputs(coolant_relation="miscible", vapour=steam(m_vapour=1.0, rho_cond=700.0, cp_cond=2500.0)))
    assert any("cannot separate" in w for w in miscible.warnings)
    immiscible = run_design("direct_contact", pool_inputs(coolant_relation="immiscible", vapour=steam(m_vapour=1.0, rho_cond=626.0, cp_cond=2300.0)))
    assert any("emulsions" in w for w in immiscible.warnings)
    assert immiscible.get("density_difference") == pytest.approx(626.0 - 963.0)
    assert any("Separator" in g for g in immiscible.source_gaps)


def test_solids_and_multicomponent_warnings():
    result = run_design("direct_contact", spray_inputs(vapour_solids="yes", multicomponent_total="yes"))
    assert any("block up" in w for w in result.warnings)
    assert any("Multicomponent" in w for w in result.warnings)


def test_spray_with_different_fluid_warns_the_model_is_inadequate():
    inputs = spray_inputs(coolant_relation="immiscible", vapour=steam(rho_cond=626.0, cp_cond=2300.0))
    result = run_design("direct_contact", inputs)
    assert any("not adequate" in w for w in result.warnings)


def test_unassessed_items_are_named():
    result = run_design("direct_contact", spray_inputs())
    code = next(c for c in result.checks if c.name == "Vessel and code requirements")
    assert code.state is CheckState.NOT_ASSESSED and "ASME" in code.detail and "Heat Exchange Institute" in code.detail


def test_every_step_has_formula_and_source():
    for inputs in (spray_inputs(), sheet_inputs(), packed_inputs(), bubble_inputs(), pool_inputs(), baffle_inputs()):
        result = run_design("direct_contact", inputs)
        for step in result.steps:
            assert step.formula and step.substitution and step.source, step.label


# --------------------------------------------------------------------------------- input validation
def test_missing_required_inputs():
    result = run_design("direct_contact", DirectContactInputs(condenser_type="spray"))
    assert result.status is Status.INVALID_INPUT
    fields = {i.field for i in result.input_issues}
    assert {"coolant_relation", "vapour.T_sat", "vapour.h_fg", "coolant.T_in", "coolant.cp"} <= fields


@pytest.mark.parametrize("kind,field", [
    ("spray", "spray.drop_diameter"), ("spray", "vapour.mu_v"), ("jet", "film.liquid_dimension"), ("fan", "film.fan_angle"),
    ("baffle_column", "film.baffle_spacing"), ("packed_bed", "packed.packing_size"), ("bubble_column", "bubble.bubble_diameter"),
])
def test_type_specific_inputs_are_required(kind, field):
    inputs = DirectContactInputs(condenser_type=kind, coolant_relation="same_fluid", vapour=VapourStream(T_sat=100.0, h_fg=2.257e6, rho_v=0.596, m_vapour=1.0), coolant=CoolantStream(T_in=85.0, T_out=95.0, rho=963.0, cp=4209.0, k=0.677))
    result = run_design("direct_contact", inputs)
    assert result.status is Status.INVALID_INPUT
    assert field in {i.field for i in result.input_issues}


def test_bad_type_rejected():
    result = run_design("direct_contact", spray_inputs(condenser_type="magic"))
    assert result.status is Status.INVALID_INPUT


def test_coolant_outlet_at_or_above_saturation_rejected():
    result = run_design("direct_contact", spray_inputs(coolant=water(m_coolant=1.0, T_out=100.0)))
    assert result.status is Status.INVALID_INPUT
    assert "coolant.T_out" in {i.field for i in result.input_issues}


def test_coolant_inlet_above_saturation_rejected():
    result = run_design("direct_contact", spray_inputs(coolant=water(T_in=110.0, m_coolant=1.0, T_out=120.0)))
    assert result.status is Status.INVALID_INPUT


def test_too_little_coolant_rejected():
    inputs = pool_inputs(coolant=CoolantStream(T_in=85.0, m_coolant=1.0, rho=963.0, cp=4209.0))
    result = run_design("direct_contact", inputs)
    assert result.status is Status.INVALID_INPUT
    assert any("above T_sat" in i.message for i in result.input_issues)


@pytest.mark.parametrize("coolant", [
    dict(m_coolant=1.0),  # only one of the three
    dict(m_coolant=1.0, T_out=95.0, utilization=0.5),  # outlet given twice
])
def test_balance_specification(coolant):
    result = run_design("direct_contact", pool_inputs(vapour=steam(), coolant=CoolantStream(T_in=85.0, rho=963.0, cp=4209.0, **coolant)))
    assert result.status is Status.INVALID_INPUT


def test_all_three_balance_quantities_rejected():
    result = run_design("direct_contact", pool_inputs(coolant=CoolantStream(T_in=85.0, m_coolant=60.0, T_out=95.0, rho=963.0, cp=4209.0)))
    assert result.status is Status.INVALID_INPUT


def test_superheat_needs_vapour_specific_heat():
    inputs = pool_inputs(vapour=VapourStream(T_sat=100.0, h_fg=2.257e6, rho_v=0.596, m_vapour=1.0, T_in=120.0))
    result = run_design("direct_contact", inputs)
    assert "vapour.cp_v" in {i.field for i in result.input_issues}


def test_vapour_below_saturation_rejected():
    result = run_design("direct_contact", pool_inputs(vapour=steam(m_vapour=1.0, T_in=90.0)))
    assert "vapour.T_in" in {i.field for i in result.input_issues}


def test_gas_needs_its_data():
    result = run_design("direct_contact", pool_inputs(gas=NonCondensable(m_gas=0.01)))
    fields = {i.field for i in result.input_issues}
    assert {"pressure", "gas.M_gas", "gas.M_vapour", "gas.p_v_vent", "gas.T_vent", "gas.cp_gas"} <= fields


def test_gas_vent_pressure_and_temperature_range():
    gas = NonCondensable(m_gas=0.01, M_gas=29.0, M_vapour=18.0, p_v_vent=2e5, T_vent=70.0, cp_gas=1005.0)
    result = run_design("direct_contact", pool_inputs(gas=gas, pressure=1e5))
    fields = {i.field for i in result.input_issues}
    assert {"gas.p_v_vent", "gas.T_vent"} <= fields


def test_all_vapour_vented_rejected():
    gas = NonCondensable(m_gas=10.0, M_gas=29.0, M_vapour=18.0, p_v_vent=9e4, T_vent=95.0, cp_gas=1005.0)
    result = run_design("direct_contact", pool_inputs(gas=gas, pressure=1e5))
    assert result.status is Status.INVALID_INPUT
    assert "nothing condenses" in result.input_issues[0].message


def test_barometric_needs_pressures():
    result = run_design("direct_contact", pool_inputs(barometric="yes"))
    assert {"pressure", "p_atm"} <= {i.field for i in result.input_issues}


def test_packed_needs_voidage_or_area():
    result = run_design("direct_contact", packed_inputs(packed=PackedData(packing_size=0.03)))
    assert "packed.voidage" in {i.field for i in result.input_issues}


def test_pool_liquid_above_saturation_rejected():
    result = run_design("direct_contact", bubble_inputs(bubble=BubbleData(bubble_diameter=1e-2, T_pool=101.0)))
    assert "bubble.T_pool" in {i.field for i in result.input_issues}


def test_negative_diameter_rejected_by_generic_validation():
    result = run_design("direct_contact", spray_inputs(spray=SprayData(drop_diameter=-1e-3, injection_velocity=5.0)))
    assert result.status is Status.INVALID_INPUT
