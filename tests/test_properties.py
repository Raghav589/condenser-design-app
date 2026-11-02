"""Property tables: printed values returned at nodes, no extrapolation, and checks independent of the transcription.

Independent checks used here, because the OCR of these tables is poor:
  * the book's own Pr column against c_p mu / k from the transcribed columns;
  * held-out rows of the HEDH steam table (1 degC spacing) against interpolation between the 5 K or 10 K nodes kept in the code;
  * HSB Table II.1 water against HEDH 5.5.3, and app_props.tex (Incropera) against HEDH 5.5.3.
"""
import math

import pytest

from condenser_app.common import properties as p
from condenser_app.schemas import CondensingFluidProperties, SinglePhaseFluidProperties


class TestWaterSaturated:
    def test_node_values_from_table_are_exact(self):
        # HEDH 5.5.3 Table 1 and Table 8 at 50 degC, printed values
        w = p.water_saturated(50.0)
        assert isinstance(w, CondensingFluidProperties)
        assert w.h_fg == pytest.approx(2381.9e3)
        assert w.rho_l == pytest.approx(987.99)
        assert w.rho_v == pytest.approx(0.08308)
        assert w.mu_l == pytest.approx(547.1e-6)
        assert w.mu_v == pytest.approx(10.62e-6)
        assert w.k_l == pytest.approx(0.6435)
        assert w.sigma == pytest.approx(67.95e-3)
        assert w.T_sat == 50.0

    @pytest.mark.parametrize("t_c, rho_l, rho_v, r, p_bar", [
        (37.0, 993.28, 0.043973, 2413.0, 0.062795),  # between the 35 and 40 degC nodes
        (72.0, 976.60, 0.21485, 2328.1, 0.33972),  # between 70 and 75
        (123.0, 940.72, 1.2242, 2194.0, 2.1809),  # between 120 and 130
        (237.0, 817.71, 15.865, 1780.2, 31.718),  # between 230 and 240: widest gap
    ])
    def test_held_out_rows_of_steam_table(self, t_c, rho_l, rho_v, r, p_bar):
        # These rows are NOT in the code's tables. Tolerance 0.3% for the steep ln(p), ln(rho_v), from interpolating over up to 10 K.
        w = p.water_saturated(t_c)
        assert w.rho_l == pytest.approx(rho_l, rel=1e-3)
        assert w.rho_v == pytest.approx(rho_v, rel=3e-3)
        assert w.h_fg / 1e3 == pytest.approx(r, rel=1e-3)
        assert p.saturation_pressure("water", t_c) / 1e5 == pytest.approx(p_bar, rel=3e-3)

    def test_log_rule_beats_straight_line_for_vapour_density(self):
        # HSB's preamble says rho_v is about linear in T; the held-out row at 37 degC justifies ln(rho_v) vs 1/T instead.
        true = 0.043973
        straight = 0.03965 + (0.05121 - 0.03965) * 0.4
        assert abs(p.water_saturated(37.0).rho_v / true - 1) < 0.0005
        assert abs(straight / true - 1) > 0.005

    @pytest.mark.parametrize("t_k, rho_l, rho_v, h_fg, mu_l, mu_v, k_l, sigma", [
        # HSB Table II.1 water, p. 1017, against HEDH 5.5.3 (NBS/NRC). State properties agree to 0.2%. The transport columns differ because HSB's come
        # from older sources: HSB eta_l is up to 5% and eta_g up to 4% higher at 550 K, lambda_l within 1.2%. The tolerances are that observed spread.
        (373.15, 958.3, 0.597, 2256.7, 277.53, 12.55, 679.0, 58.91),
        (400.0, 937.5, 1.370, 2183.0, 218.9, 13.57, 685.7, 53.50),
        (430.0, 910.3, 3.020, 2092.8, 175.73, 14.716, 683.3, 47.16),
        (460.0, 879.4, 5.975, 1990.4, 147.24, 15.86, 671.3, 40.66),
        (490.0, 844.3, 10.95, 1871.5, 126.6, 17.00, 646.0, 33.90),
        (520.0, 803.8, 18.90, 1731.0, 111.05, 18.14, 618.3, 26.96),
        (550.0, 756.1, 31.52, 1562.6, 99.21, 19.33, 580.9, 19.66),
    ])
    def test_agrees_with_hsb_table_ii_1(self, t_k, rho_l, rho_v, h_fg, mu_l, mu_v, k_l, sigma):
        w = p.water_saturated(t_k - 273.15)
        assert w.rho_l == pytest.approx(rho_l, rel=0.002)
        assert w.rho_v == pytest.approx(rho_v, rel=0.002)
        assert w.h_fg / 1e3 == pytest.approx(h_fg, rel=0.002)
        assert w.mu_l * 1e6 == pytest.approx(mu_l, rel=0.05)
        assert w.mu_v * 1e6 == pytest.approx(mu_v, rel=0.04)
        assert w.k_l * 1e3 == pytest.approx(k_l, rel=0.012)
        assert w.sigma * 1e3 == pytest.approx(sigma, rel=0.007)  # worst is 0.63% at 520 K

    def test_derived_liquid_cp_agrees_with_incropera_and_hsb(self):
        # c_p = Pr k/mu from HEDH Table 8 against Incropera (app_props.tex): 0.5% (Pr and mu are rounded to 3-4 s.f.)
        for t_k in (290, 300, 320, 340, 360, 373.15, 400):
            assert p.water_saturated(t_k - 273.15).cp_l == pytest.approx(p.water_liquid(t_k - 273.15).cp, rel=0.005)
        # HSB Table II.1 at 373.15 K: 4.22 kJ/(kg K)
        assert p.water_saturated(100.0).cp_l == pytest.approx(4220, rel=0.002)

    def test_derived_vapour_cp_agrees_with_hsb(self):
        # HSB Table II.1 c_p,g at 373.15 K and 460 K: 2.03, 2.70 kJ/(kg K). 3%: it is 0.6% and 2.4% here; different compilations, as above.
        assert p.saturated_vapour("water", 100.0).cp == pytest.approx(2030, rel=0.03)
        assert p.saturated_vapour("water", 460 - 273.15).cp == pytest.approx(2700, rel=0.03)

    def test_range_enforced(self):
        p.water_saturated(0.01)
        p.water_saturated(300.0)
        with pytest.raises(ValueError):
            p.water_saturated(0.0)
        with pytest.raises(ValueError):
            p.water_saturated(300.1)

    def test_table_columns_are_smooth(self):
        # Guards transcription: properties monotonic over the table where physics requires it
        temps = [t for t in range(1, 300, 3)]
        states = [p.water_saturated(t) for t in temps]
        for a, b in zip(states, states[1:]):
            assert b.h_fg < a.h_fg
            assert b.rho_v > a.rho_v
            assert b.mu_l < a.mu_l
            assert b.mu_v > a.mu_v
            assert b.sigma < a.sigma
        ln_p = [math.log(v) for v in p._W1_P]
        reciprocal = [1 / v for v in p._W1_TK]
        slopes = [(ln_p[i + 1] - ln_p[i]) / (reciprocal[i] - reciprocal[i + 1]) for i in range(len(ln_p) - 1)]
        assert max(abs(slopes[i + 1] / slopes[i] - 1) for i in range(len(slopes) - 1)) < 0.02


class TestSaturationPressure:
    def test_water_at_atmospheric(self):
        # HEDH Table 1: 100 degC at 1.0132 bar
        assert p.saturation_temperature("water", 1.0132e5) == pytest.approx(100.0, abs=1e-6)
        assert p.saturation_pressure("water", 100.0) == pytest.approx(1.0132e5)

    @pytest.mark.parametrize("key, t_c", [("water", 45.3), ("n_pentane", 80.0), ("ammonia", 20.0), ("r22", 40.0)])
    def test_inverse(self, key, t_c):
        assert p.saturation_temperature(key, p.saturation_pressure(key, t_c)) == pytest.approx(t_c, abs=1e-8)

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            p.saturation_temperature("water", 100.0)
        with pytest.raises(ValueError):
            p.saturation_pressure("n_pentane", 0.0)
        with pytest.raises(ValueError):
            p.saturation_pressure("freon", 20.0)


class TestWaterLiquid:
    def test_node_from_app_props(self):
        # app_props.tex at 300 K
        w = p.water_liquid(300 - 273.15)
        assert isinstance(w, SinglePhaseFluidProperties)
        assert (w.rho, w.cp) == (pytest.approx(997), pytest.approx(4179))
        assert (w.mu, w.k) == (pytest.approx(855e-6), pytest.approx(0.613))

    def test_interpolated_prandtl_stays_consistent_with_table(self):
        # app_props.tex Pr at 322.5 K lies between 3.77 and 3.42
        w = p.water_liquid(322.5 - 273.15)
        assert 3.42 < w.cp * w.mu / w.k < 3.77

    def test_against_steam_table(self):
        # Incropera vs NBS/NRC: rho 0.1%; mu 1.5% and k 1.5% (the two compilations differ most at low T; 280 K is the worst)
        for t_k in (280, 300, 320, 340, 360, 380, 400):
            a, b = p.water_liquid(t_k - 273.15), p.water_saturated(t_k - 273.15)
            assert a.rho == pytest.approx(b.rho_l, rel=0.001)
            assert a.mu == pytest.approx(b.mu_l, rel=0.015)
            assert a.k == pytest.approx(b.k_l, rel=0.015)

    def test_range(self):
        p.water_liquid(1.85)
        p.water_liquid(126.85)
        with pytest.raises(ValueError):
            p.water_liquid(1.0)
        with pytest.raises(ValueError):
            p.water_liquid(127.0)


class TestAir:
    @pytest.mark.parametrize("t_k, rho, cp, mu_e7, k_e3", [
        (300, 1.1614, 1007, 184.6, 26.3),  # also in app_props.tex
        (350, 0.9950, 1009, 208.2, 30.0),
        (600, 0.5804, 1051, 305.8, 46.9),
    ])
    def test_nodes(self, t_k, rho, cp, mu_e7, k_e3):
        a = p.air(t_k - 273.15)
        assert a.rho == pytest.approx(rho)
        assert a.cp == pytest.approx(cp)
        assert a.mu == pytest.approx(mu_e7 * 1e-7)
        assert a.k == pytest.approx(k_e3 * 1e-3)

    def test_density_interpolation_is_ideal_gas(self):
        # rho T linear reproduces p/(R T): within 0.1% of 1.1614 x 300/325 at the middle of a 50 K gap
        assert p.air(325 - 273.15).rho == pytest.approx(1.1614 * 300 / 325, rel=0.001)

    def test_pressure_scales_density_only(self):
        a, b = p.air(30.0), p.air(30.0, p_pa=2 * p.ATM)
        assert b.rho == pytest.approx(2 * a.rho)
        assert b.mu == a.mu

    def test_range(self):
        with pytest.raises(ValueError):
            p.air(-74.0)
        with pytest.raises(ValueError):
            p.air(727.0)
        with pytest.raises(ValueError):
            p.air(30.0, p_pa=0.0)


# HSB Table II.1 printed Pr_l and Pr_g, used only to check the transcribed c_p, mu and k columns.
PRINTED_PR = {
    "n_pentane": ((4.29, 4.09, 3.94, 3.80, 3.67, 3.53, 3.36, 3.16, 2.84), (0.74, 0.77, 0.79, 0.81, 0.83, 0.88, 0.94, 1.03, 1.19)),
    "isopentane": ((4.47, 4.17, 4.06, 3.79, 3.58, 3.40, 3.22, 3.05, 2.99), (0.87, 0.86, 0.86, 0.87, 0.89, 0.93, 0.97, 1.07, 1.20)),
    "n_hexane": ((4.81, 4.51, 4.49, 4.36, 4.19, 4.01, 3.99, 3.83, 3.97), (0.83, 0.87, 0.91, 0.91, 0.92, 0.97, 1.02, 1.14, 1.31)),
    "n_heptane": ((5.27, 5.08, 4.87, 4.59, 4.43, 4.27, 4.28, 4.37, 4.97), (0.80, 0.84, 0.89, 0.90, 0.94, 0.97, 1.10, 1.29, 1.76)),
    "methanol": ((5.13, 4.67, 4.15, 3.61, 3.34, 2.82, 2.56, 2.42), (0.94, 0.91, 0.90, 0.92, 0.94, 0.97, 1.04, 1.35)),
    "ethanol": ((8.37, 6.88, 5.91, 5.18, 4.90, 4.64, 4.28, 4.19), (0.96, 0.95, 0.96, 0.97, 1.02, 1.20, 1.31, 1.51)),
    "ammonia": ((2.06, 1.88, 1.58, 1.39, 1.36, 1.32, 1.34, 1.41, 1.43), (1.04, 1.11, 1.17, 1.25, 1.31, 1.34, 1.49, 1.70, 1.86)),
    "r22": ((2.92, 2.88, 2.84, 2.92, 3.09, 3.45, 4.18, 6.89), (0.86, 0.89, 0.91, 0.97, 1.07, 1.22, 1.69, 2.31)),
}
# The book's Pr column carries its own citation [12] and does not always equal its c_p eta / lambda. Tolerance is 3% (columns have 2 to 4 s.f.) except where the
# book itself disagrees by more; those columns were re-read from zoomed page images and are transcribed as printed. Each exception is the worst mismatch found:
PR_TOLERANCE = {
    ("methanol", "l"): 0.095,  # systematic: columns give 3 to 9% below the printed Pr in every row, so it is not one misread cell
    ("ammonia", "l"): 0.095,  # 8.8% at 370 K, 3.1% at 270 K
    ("ammonia", "g"): 0.095,  # systematic: 4 to 9% above from 270 K up
    ("r22", "l"): 0.05,  # 4.6% at 355 K
    ("r22", "g"): 0.16,  # 15% at 340 K only
    ("n_pentane", "l"): 0.055,  # 5.3% at 440 K only
    ("isopentane", "l"): 0.035,  # 3.2% at 430 K only
}


class TestHsbTables:
    @pytest.mark.parametrize("key", sorted(PRINTED_PR))
    def test_printed_prandtl_matches_columns(self, key):
        table = p._TABLES[key]
        pr_l, pr_g = PRINTED_PR[key]
        assert len(pr_l) == len(table.t) and len(pr_g) == len(table.t)
        for i in range(len(table.t)):
            # kJ/(kg K) x uPa.s / (mW/m K) = 1
            assert table.cp_l[i] * table.mu_l[i] / table.k_l[i] == pytest.approx(pr_l[i], rel=PR_TOLERANCE.get((key, "l"), 0.03)), (key, "l", table.t[i])
            assert table.cp_v[i] * table.mu_v[i] / table.k_v[i] == pytest.approx(pr_g[i], rel=PR_TOLERANCE.get((key, "g"), 0.03)), (key, "g", table.t[i])

    def test_pentane_at_335k_matches_example_16_2(self):
        # Table II.1 node equals the data set stated in Ex. 16.2, so the two sources agree
        a = p.saturated("n_pentane", 335 - 273.15)
        b = p.pentane_hsb_example_16_2()
        for name in ("h_fg", "rho_l", "rho_v", "mu_l", "mu_v", "k_l", "cp_l", "sigma"):
            assert getattr(a, name) == pytest.approx(getattr(b, name), rel=1e-9), name

    def test_hexane_at_370k_shows_ex_16_1_used_vapour_enthalpy(self):
        # Book error documented: Ex. 16.1 states h_LG = 776.9 kJ/kg, but at 370 K Table II.1 has h_g = 776.9 and dh = 311.7.
        assert p.saturated("n_hexane", 370 - 273.15).h_fg == pytest.approx(311.7e3)

    def test_ex_16_4_properties_are_hexane_not_isopentane(self):
        # Book error documented: Ex. 16.4 names isopentane at 400 K, but its numbers are the n-hexane row at 400 K
        hexane = p.saturated("n_hexane", 400 - 273.15)
        assert (hexane.rho_l, hexane.rho_v, hexane.h_fg, hexane.mu_l, hexane.k_l) == pytest.approx((551, 14.1, 283.8e3, 128.5e-6, 82e-3))
        assert p.saturated("isopentane", 400 - 273.15).rho_l == pytest.approx(488)

    def test_interpolates_and_respects_range(self):
        mid = p.saturated("n_pentane", (335 + 350) / 2 - 273.15)
        assert 566.0 < mid.rho_l < 582.9
        with pytest.raises(ValueError):
            p.saturated("n_pentane", 309.0 - 273.15)
        with pytest.raises(ValueError):
            p.saturated("n_pentane", 441.0 - 273.15)
        with pytest.raises(ValueError):
            p.saturated("no_such_fluid", 50.0)

    def test_table_columns_have_equal_lengths(self):
        for key, table in p._TABLES.items():
            for name in ("p", "rho_l", "rho_v", "h_fg", "cp_l", "cp_v", "mu_l", "mu_v", "k_l", "k_v", "sigma"):
                assert len(getattr(table, name)) == len(table.t), (key, name)
            assert list(table.t) == sorted(table.t)
            assert list(table.p) == sorted(table.p)


class TestPresets:
    def test_every_preset_evaluates_at_both_range_limits_and_refuses_beyond(self):
        for preset in p.PRESETS.values():
            for t_c in (preset.t_min_c, preset.t_max_c):
                state = preset.at(t_c)
                assert state is not None and preset.source
            if preset.t_min_c == preset.t_max_c:
                with pytest.raises(ValueError):
                    preset.at(preset.t_min_c + 1.0)
            else:
                with pytest.raises(ValueError):
                    preset.at(preset.t_max_c + 0.5)
                with pytest.raises(ValueError):
                    preset.at(preset.t_min_c - 0.5)

    def test_kinds_match_return_types(self):
        for preset in p.PRESETS.values():
            state = preset.at(preset.t_min_c)
            expected = SinglePhaseFluidProperties if preset.kind == "single_phase" else CondensingFluidProperties
            assert isinstance(state, expected), preset.key

    def test_required_names_present(self):
        assert {"water_saturated", "water_liquid", "air", "n_pentane", "n_pentane_hsb_ex16_2"} <= set(p.PRESETS)
