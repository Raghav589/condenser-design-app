"""Benchmarks of common/condensation.py against the worked examples of HSB, read from the page images.

Expected numbers are the book's. All inputs are the book's stated data; nothing is tuned. Tolerances say why they are what they are.
Where the book is wrong, the test pins both the book's printed number and the correct one; docs/algorithms/condensation_library.md lists each.
Classes named Limit... are NOT benchmarks: they are limiting cases, hand values or trends for functions without a worked example.
"""
import math

import pytest

from condenser_app.common import condensation as c
from condenser_app.common import iteration, properties
from condenser_app.schemas import CondensingFluidProperties as Fluid


@pytest.fixture
def pentane():
    return properties.pentane_hsb_example_16_2()


def fluid(**kw):
    defaults = dict(name="test", T_sat=None, h_fg=None, rho_l=None, rho_v=None, mu_l=None, mu_v=None, k_l=None, cp_l=None, sigma=None)
    defaults.update(kw)
    return Fluid(**defaults)


@pytest.fixture
def heptane_octane():
    # HSB Example 10.1, p. 406: rho_G 18.1, rho_L 567, eta_G 7.1, eta_L 156 uN.s/m2, sigma 8.1 mN/m
    return fluid(rho_l=567.0, rho_v=18.1, mu_l=156e-6, mu_v=7.1e-6, sigma=8.1e-3)


class TestExample16_1:
    """n-hexane on a vertical plate: HSB p. 577-578. The book states h_LG = 776.9 kJ/kg (see test below)."""

    def stated(self):
        return fluid(rho_l=585.0, rho_v=7.0, mu_l=158.9e-6, k_l=91.0e-3, h_fg=776.9e3)

    def test_mean_coefficient(self):
        # Book rounds to 4 s.f.
        assert c.nusselt_mean(self.stated(), 20.0, 0.1) == pytest.approx(1482, abs=1)

    def test_condensate_and_reynolds(self):
        gamma = c.condensate_per_width(1482.0, 20.0, 0.1, 776.9e3)
        assert gamma == pytest.approx(0.003815, rel=1e-3)
        assert c.film_reynolds(gamma, 158.9e-6) == pytest.approx(96.04, rel=1e-3)

    def test_eq_16_16_form_is_the_same_result(self):
        f = self.stated()
        alpha = c.nusselt_mean(f, 20.0, 0.1)
        gamma = c.condensate_per_width(alpha, 20.0, 0.1, f.h_fg)
        assert c.nusselt_mean_from_flow(f, gamma) == pytest.approx(alpha, rel=1e-10)

    def test_book_used_vapour_enthalpy_as_latent_heat(self):
        # Table II.1, p. 987, n-hexane at 370 K: h_g = 776.9, dh_g,l = 311.7 kJ/kg. With the true latent heat the answers are 1180 and 0.00757.
        true = fluid(rho_l=585.0, rho_v=7.0, mu_l=158.9e-6, k_l=91.0e-3, h_fg=311.7e3)
        alpha = c.nusselt_mean(true, 20.0, 0.1)
        assert alpha == pytest.approx(1180, abs=1)
        assert c.condensate_per_width(alpha, 20.0, 0.1, true.h_fg) == pytest.approx(0.00757, rel=2e-3)


class TestExample16_2:
    """n-pentane, shear control, x = 0.85, m = 500, D = 0.03: HSB p. 583-587."""

    def test_friedel_terms(self, pentane):
        t = c.friedel_terms(pentane, 500.0, 0.85, 0.03)
        pg = c.phase_gradients(pentane, 500.0, 0.85, 0.03)
        assert pg.re_lo == pytest.approx(9.434e4, rel=1e-3)
        assert pg.re_go == pytest.approx(1.974e6, rel=1e-3)
        assert pg.f_lo == pytest.approx(0.00451, abs=5e-6)  # book 3 s.f.
        assert pg.f_go == pytest.approx(0.00211, abs=5e-6)
        # 128.95 uses f rounded to 0.00451; the unrounded f gives 128.89 (0.05%)
        assert t.dpdz_lo == pytest.approx(128.95, rel=1e-3)
        assert t.E == pytest.approx(31.0, abs=0.05)
        assert t.F == pytest.approx(0.559, abs=5e-4)
        assert t.H == pytest.approx(33.09, abs=0.01)
        assert t.rho_tp == pytest.approx(7.468, abs=1e-3)
        assert t.Fr == pytest.approx(15231, abs=2)
        assert t.We == pytest.approx(8.887e4, rel=1e-3)
        assert t.phi_lo2 == pytest.approx(57.08, abs=0.05)  # book carries rounded E
        assert t.dpdz == pytest.approx(7360, rel=2e-3)  # rounding of f and E, 0.1%

    def test_laminar_film(self, pentane):
        gamma = c.film_flow_in_tube(500.0, 0.85, 0.03)
        assert gamma == pytest.approx(0.5625)
        assert c.shear_stress_from_gradient(7360.0, 0.03) == pytest.approx(55.20)
        assert c.laminar_shear_coefficient(pentane, gamma, 55.20) == pytest.approx(1314, abs=1)
        delta = c.laminar_shear_film_thickness(pentane, gamma, 55.20)
        assert delta == pytest.approx(7.456e-5, rel=1e-3)
        assert delta * pentane.rho_l * 9.81 == pytest.approx(0.426, abs=1e-3)
        m_cond = math.pi * 0.03 * 1314 * 35 / 337.3e3
        assert m_cond == pytest.approx(0.01285, rel=1e-3)

    def test_turbulent_film(self, pentane):
        re_f = c.film_reynolds(c.film_flow_in_tube(500.0, 0.85, 0.03), pentane.mu_l)
        assert re_f == pytest.approx(14150, abs=1)
        delta_plus = c.kosky_staub_delta_plus(re_f)
        assert delta_plus == pytest.approx(215.9, abs=0.1)
        assert c.prandtl_liquid(pentane) == pytest.approx(4.09, abs=5e-3)
        assert c.kosky_staub_t_plus(delta_plus, 4.09) == pytest.approx(40.71, abs=0.01)
        # With the book's tau (55.20): 11,104
        assert c.kosky_staub_coefficient(pentane, re_f, 55.20) == pytest.approx(11104, rel=2e-4)
        delta = c.turbulent_film_thickness(pentane, delta_plus, 55.20)
        assert delta == pytest.approx(1.914e-4, rel=1e-3)
        # BOOK SLIP, 0.7%: it prints delta rho_L g = 1.087 N/m2, but its own 1.914e-4 x 582.9 x 9.81 = 1.0945. Harmless: still "around 2%" of tau.
        assert delta * pentane.rho_l * 9.81 == pytest.approx(1.0945, rel=2e-3)
        assert delta * pentane.rho_l * 9.81 / 55.20 == pytest.approx(0.02, abs=5e-4)

    def test_end_to_end_from_inputs_only(self, pentane):
        # Friedel, then tau, then Kosky-Staub, without using any book intermediate: 0.05% from rounding of f and E in the book
        result = c.vertical_downflow_local(pentane, 500.0, 0.85, 0.03)
        assert result.tau_i == pytest.approx(55.20, rel=2e-3)
        assert result.alpha_shear == pytest.approx(11104, rel=1e-3)
        assert result.film_thickness == pytest.approx(1.914e-4, rel=2e-3)
        # Book's note: gravity "around 2%" of tau_i
        assert result.gravity_to_shear == pytest.approx(0.02, abs=5e-4)


class TestExample16_3:
    """Boyko-Kruzhilin at the Ex. 16.2 conditions: HSB p. 588-589."""

    def test_alpha_lo(self, pentane):
        assert c.boyko_kruzhilin_alpha_lo(pentane, 500.0, 0.03) == pytest.approx(1199, abs=1)

    def test_local(self, pentane):
        # 0.06%: book multiplies by the rounded 1199
        assert c.boyko_kruzhilin_local(pentane, 500.0, 0.85, 0.03) == pytest.approx(10593, rel=1e-3)

    def test_mean_is_arithmetic_of_ends(self, pentane):
        assert c.boyko_kruzhilin_local(pentane, 500.0, 0.95, 0.03) == pytest.approx(11191, rel=1e-3)
        assert c.boyko_kruzhilin_local(pentane, 500.0, 0.05, 0.03) == pytest.approx(2820, rel=1e-3)
        assert c.boyko_kruzhilin_mean(pentane, 500.0, 0.03, 0.95, 0.05) == pytest.approx(7006, rel=1e-3)

    def test_book_comparison_with_example_16_2(self, pentane):
        # Book's note says Ex. 16.2 gives "11,179" and is "around 5% higher". Ex. 16.2 actually prints 11,104; the ratio 1.048 survives.
        ratio = 11104 / c.boyko_kruzhilin_local(pentane, 500.0, 0.85, 0.03)
        assert ratio == pytest.approx(1.05, abs=0.005)


class TestExample19_2:
    """Boyko-Kruzhilin in the air-cooled tube bundle, HSB p. 682-683 (tube-side part only)."""

    def test_values(self, pentane):
        m = 4 * 1.0 / (80 * math.pi * 0.0186**2)
        assert m == pytest.approx(46.0, abs=0.05)
        assert c.boyko_kruzhilin_alpha_lo(pentane, 46.0, 0.0186) == pytest.approx(195.7, abs=0.1)
        assert c.boyko_kruzhilin_local(pentane, 46.0, 0.85, 0.0186) == pytest.approx(1729, abs=1.5)
        assert c.boyko_kruzhilin_local(pentane, 46.0, 0.05, 0.0186) == pytest.approx(460, abs=1)
        assert c.boyko_kruzhilin_mean(pentane, 46.0, 0.0186, 0.85, 0.05) == pytest.approx(1095, abs=1)


class TestExample16_4:
    """Condensation outside a horizontal tube with iteration on wall temperature: HSB p. 590-592.

    The book says isopentane at 400 K, but its numbers are the n-hexane row of Table II.1 (see test_properties). Inputs as stated.
    """

    def stated(self):
        return fluid(k_l=82e-3, rho_l=551.0, rho_v=14.1, h_fg=283.8e3, mu_l=128.5e-6)

    def test_group_constant(self):
        # Book: alpha = 2388 / (T_sat - T_w)^(1/4)
        assert c.nusselt_horizontal_tube_mean(self.stated(), 1.0, 0.03) == pytest.approx(2388, abs=1)

    def test_iteration_converges_to_book_value(self):
        f = self.stated()
        # Water side from the book's data: 0.023 Re^0.8 Pr^(1/3) k/D reproduces its 3454 (0.1%)
        re_w = 1.0 * 996 * 0.025 / 1003e-6
        pr_w = 4183 * 1003e-6 / 598e-3
        assert re_w == pytest.approx(24825, abs=1)
        assert pr_w == pytest.approx(7.015, abs=1e-3)
        alpha_w = 0.023 * 598e-3 / 0.025 * re_w**0.8 * pr_w ** (1 / 3)
        assert alpha_w == pytest.approx(3454, rel=2e-3)
        alpha_w = 3454.0  # use the book's value below so that the comparison isolates Eq. 16.57

        def next_tw(t_w):
            alpha_o = c.nusselt_horizontal_tube_mean(f, 400.0 - t_w, 0.03)
            u = 1.0 / (1.0 / alpha_o + 0.03 / (alpha_w * 0.025))
            return 400.0 - u * (400.0 - 293.0) / alpha_o

        t_w, info = iteration.fixed_point(next_tw, 0.5 * (400.0 + 293.0), name="wall temperature", tol=1e-4)
        alpha_o = c.nusselt_horizontal_tube_mean(f, 400.0 - t_w, 0.03)
        assert info.converged
        assert t_w == pytest.approx(316.01, abs=0.02)
        assert alpha_o == pytest.approx(789, abs=0.5)
        u = 1.0 / (1.0 / alpha_o + 0.03 / (alpha_w * 0.025))
        assert u * (400 - 293) == pytest.approx(66.25e3, rel=1e-3)
        assert math.pi * 0.03 * u * 107 / 283.8e3 == pytest.approx(0.0220, abs=1e-4)


class TestExample16_5:
    """Horizontal tube, n-pentane, x = 0.075: HSB p. 595-598."""

    def test_wallis_parameter(self, pentane):
        assert c.wallis_gas_parameter(pentane, 500.0, 0.075, 0.03) == pytest.approx(1.14, abs=5e-3)

    def test_premoli_void_fraction(self, pentane):
        assert c.premoli_slip_ratio(pentane, 500.0, 0.075, 0.03) == pytest.approx(2.037, abs=1e-3)
        assert c.void_fraction_premoli(pentane, 500.0, 0.075, 0.03) == pytest.approx(0.7849, abs=2e-4)

    def test_stratifying_coefficient(self, pentane):
        assert c.stratified_coefficient(pentane, 35.0, 0.03, 0.7849) == pytest.approx(950, abs=0.5)

    def test_friedel_terms_and_the_book_arithmetic_slip(self, pentane):
        t = c.friedel_terms(pentane, 500.0, 0.075, 0.03)
        assert t.E == pytest.approx(1.097, abs=1e-3)
        assert t.F == pytest.approx(0.1301, abs=1e-4)
        assert t.rho_tp == pytest.approx(74.74, abs=0.01)
        assert t.Fr == pytest.approx(152.08, abs=0.05)
        assert t.We == pytest.approx(8880, abs=1)
        # BOOK ERROR: it prints 9.110, but its own terms 1.097 + 3.24(0.1301)(33.09)/(152.08^0.045 x 8880^0.035) give 9.19.
        book_terms = 1.097 + 3.24 * 0.1301 * 33.09 / (152.08**0.045 * 8880**0.035)
        assert book_terms == pytest.approx(9.19, abs=0.01)
        assert t.phi_lo2 == pytest.approx(book_terms, abs=0.01)
        assert t.phi_lo2 != pytest.approx(9.110, abs=0.05)

    def test_annular_coefficient_with_book_shear(self, pentane):
        re_f = c.film_reynolds(c.film_flow_in_tube(500.0, 0.075, 0.03), pentane.mu_l)
        assert re_f == pytest.approx(87264, abs=1)
        assert c.kosky_staub_delta_plus(re_f) == pytest.approx(1061, abs=2)
        assert c.kosky_staub_t_plus(1061.0, 4.09) == pytest.approx(44.69, abs=0.01)
        assert c.kosky_staub_coefficient(pentane, re_f, 8.811) == pytest.approx(4041, abs=2)

    def test_blend(self, pentane):
        assert c.palen_interpolation(1.14, 4041.0, 950.0) == pytest.approx(2928, abs=1)
        # Carrying the corrected multiplier through gives 4059 and 2945 (0.6% above the book; cause: its 9.110).
        tau = c.shear_stress_from_gradient(c.friedel_gradient(pentane, 500.0, 0.075, 0.03), 0.03)
        re_f = c.film_reynolds(c.film_flow_in_tube(500.0, 0.075, 0.03), pentane.mu_l)
        alpha_an = c.kosky_staub_coefficient(pentane, re_f, tau)
        assert alpha_an == pytest.approx(4059, abs=2)
        u = c.wallis_gas_parameter(pentane, 500.0, 0.075, 0.03)
        assert c.palen_interpolation(u, alpha_an, 950.0) == pytest.approx(2945, abs=3)


class TestExample16_6:
    """Effect of condensation on shear, Eq. 16.76: HSB p. 603-604."""

    def test_values(self, pentane):
        flux = 0.1093 / (math.pi * 0.03)
        assert flux == pytest.approx(1.160, abs=1e-3)
        assert c.mass_transfer_ratio(-1.4044) == pytest.approx(1.861, abs=1e-3)
        tau_m = c.suction_corrected_shear(55.20, 500.0, 0.85, 6.36, flux)
        assert tau_m == pytest.approx(102.73, abs=0.05)
        alpha = pentane.cp_l * math.sqrt(pentane.rho_l * 102.73) / 40.71
        assert alpha == pytest.approx(15147, abs=5)
        re_f = c.film_reynolds(c.film_flow_in_tube(500.0, 0.85, 0.03), pentane.mu_l)
        assert c.kosky_staub_coefficient(pentane, re_f, 102.73) == pytest.approx(15147, rel=1e-3)
        # Wall temperature: 335 - 390.9e3/15147
        assert 335 - flux * 337.3e3 / alpha == pytest.approx(309.2, abs=0.05)

    def test_note_on_inconsistent_rate(self):
        # Book: "0.1093, the rate calculated ... in Example 16.2"; Ex. 16.2 prints 0.1086. 0.1093 corresponds to its other stale value 11,179.
        assert math.pi * 0.03 * 11104 * 35 / 337.3e3 == pytest.approx(0.1086, abs=1e-4)
        assert math.pi * 0.03 * 11179 * 35 / 337.3e3 == pytest.approx(0.1093, abs=1e-4)


class TestExample16_13:
    """Silver-Bell-Ghaly for a hydrocarbon mixture: HSB p. 636-638, Table 16.5. Node values from Table 16.4 and the text."""

    # M_G at inlet: the text gives 192.98 and the sum of its terms is 192.99; the table prints 198.98, a misprint (the table's own mean 168.94 needs 192.9).
    # T_c at 355 K: the table gives 312.26 (heat balance agrees); the formula on p. 638 prints 312.84, a misprint.
    Q = [0.0, 16.85e6, 32.74e6, 47.45e6, 63.58e6]
    TE = [360.0, 355.0, 350.0, 345.0, 340.0]
    TC = [320.28, 312.26, 304.69, 297.68, 290.0]
    AG = [150.0, 118.0, 88.0, 60.0, 22.0]
    MG = [192.98, 144.94, 100.25, 60.58, 17.27]

    def test_increment_one_in_detail(self):
        assert c.sbg_z(168.96, 2000.0, 5.0, 16.85e6) == pytest.approx(0.1002, abs=2e-4)
        delta_a = 16.85e6 * c.sbg_integrand(0.1002, 800.0, 134.0, 357.5, 316.27)
        assert delta_a == pytest.approx(816, abs=1)

    def test_total_and_each_increment(self):
        result = c.sbg_area(self.Q, self.TE, self.TC, self.AG, self.MG, 2000.0, 800.0)
        # Book rounds each to the nearest m2 and each Z to 4 d.p.: 0.1%
        assert result.delta_a == pytest.approx([816, 721, 632, 609], rel=1.5e-3)
        assert result.z == pytest.approx([0.1002, 0.0771, 0.0547, 0.0241], abs=2e-4)
        assert result.area == pytest.approx(2778, rel=1e-3)

    def test_coolant_heat_balance(self):
        # T_c,out = 63.58e6/(500 x 4200) + 290 = 320.28
        assert 63.58e6 / (500 * 4200) + 290 == pytest.approx(320.28, abs=5e-3)


class TestExample3_3:
    """General area integral, HSB p. 165-166, Table 3.2."""

    def test_area(self):
        q = [i * 1e5 for i in range(18)]
        th = [210, 208, 206, 204, 200, 195, 190, 184, 178, 171, 163, 155, 146, 137, 126, 113, 100, 86]
        tc = [162, 154, 145, 138, 130, 125, 118, 113, 108, 105, 102, 100, 96, 93, 88, 81, 70, 59]
        u = [640, 680, 705, 855, 925, 970, 1005, 1025, 1030, 1020, 1000, 960, 925, 875, 830, 780, 730, 680]
        area = c.area_by_integration(q, u, [a - b for a, b in zip(th, tc)])
        # BOOK MISPRINT: the text says the sum is 39.355 m2, but its own column of 17 increments adds to 39.855.
        book_increments = [2.989, 2.524, 2.049, 1.658, 1.509, 1.428, 1.378, 1.381, 1.436, 1.562, 1.767, 2.028, 2.380, 2.884, 3.589, 4.286, 5.007]
        assert sum(book_increments) == pytest.approx(39.855, abs=1e-3)
        assert area == pytest.approx(39.855, abs=0.005)


class TestExample10_1:
    """Two-phase gradient by three correlations: HSB p. 406-409."""

    def test_reference_gradients(self, heptane_octane):
        pg = c.phase_gradients(heptane_octane, 407.0, 0.1, 0.0221)
        # Book rounds f to 3 s.f.: dp_GO printed 1938, unrounded 1950 (0.6%); the others 0.1%
        assert pg.dpdz_l == pytest.approx(112, abs=0.2)
        assert pg.dpdz_g == pytest.approx(34.70, abs=0.05)
        assert pg.dpdz_lo == pytest.approx(134.83, abs=0.1)
        assert pg.dpdz_go == pytest.approx(1938, rel=0.007)
        assert c.martinelli_parameter(pg.dpdz_l, pg.dpdz_g) == pytest.approx(1.796, abs=0.003)
        assert c.chisholm_parameter(pg.dpdz_go, pg.dpdz_lo) == pytest.approx(3.80, abs=0.005)

    def test_lockhart_martinelli(self, heptane_octane):
        assert c.LOCKHART_MARTINELLI_C["tt"] == 20.0
        assert c.lockhart_martinelli_phi_l2(1.796) == pytest.approx(12.44, abs=0.01)
        # Book 12.44 x 112 = 1393
        assert c.lockhart_martinelli_gradient(heptane_octane, 407.0, 0.1, 0.0221) == pytest.approx(1393, abs=2)

    def test_friedel(self, heptane_octane):
        t = c.friedel_terms(heptane_octane, 407.0, 0.1, 0.0221)
        assert t.E == pytest.approx(0.954, abs=1e-3)
        assert t.F == pytest.approx(0.162, abs=1e-3)
        assert t.H == pytest.approx(12.36, abs=0.01)
        assert t.rho_tp == pytest.approx(140.6, abs=0.05)
        assert t.Fr == pytest.approx(38.65, abs=0.01)
        assert t.We == pytest.approx(3214, abs=1)
        assert t.phi_lo2 == pytest.approx(5.101, abs=2e-3)
        assert t.dpdz == pytest.approx(687.8, abs=0.5)

    def test_chisholm(self, heptane_octane):
        # Book: Y = 3.80, B = 2.73: 1 + 13.44(0.332 + 0.018) = 5.70. Book 3 s.f.
        assert c.chisholm_phi_lo2(0.1, 3.80, 407.0) == pytest.approx(5.70, abs=0.01)
        # From inputs only: 769.0 against 768.5, 0.06% (book rounds f and B)
        assert c.chisholm_gradient(heptane_octane, 407.0, 0.1, 0.0221) == pytest.approx(768.5, rel=0.002)

    def test_chisholm_equation_needs_the_braces(self):
        # BOOK MISPRINT in Eq. 10.45: at x = 1, phi_LO^2 must equal Y^2, as only the bracketed form of the worked example gives.
        assert c.chisholm_phi_lo2(1.0, 3.8, 407.0) == pytest.approx(3.8**2)
        assert c.chisholm_phi_lo2(0.0, 3.8, 407.0) == pytest.approx(1.0)

    def test_ratio_stated_by_book(self, heptane_octane):
        # Book: 1 : 1.12 : 2.03 for Friedel, Chisholm, Lockhart-Martinelli
        f = c.friedel_gradient(heptane_octane, 407.0, 0.1, 0.0221)
        assert c.chisholm_gradient(heptane_octane, 407.0, 0.1, 0.0221) / f == pytest.approx(1.12, abs=0.01)
        assert c.lockhart_martinelli_gradient(heptane_octane, 407.0, 0.1, 0.0221) / f == pytest.approx(2.03, abs=0.01)


class TestExample10_3:
    """Void fraction and gravity: HSB p. 414-416."""

    def test_void_fractions(self, heptane_octane):
        assert c.void_fraction_homogeneous(0.1, 567.0, 18.1) == pytest.approx(0.777, abs=1e-3)
        assert c.void_fraction_lockhart_martinelli(12.44**0.5) == pytest.approx(0.716, abs=1e-3)
        assert c.premoli_slip_ratio(heptane_octane, 407.0, 0.1, 0.0221) == pytest.approx(1.682, abs=2e-3)
        assert c.void_fraction_premoli(heptane_octane, 407.0, 0.1, 0.0221) == pytest.approx(0.674, abs=5e-4)

    def test_gravity_and_total(self):
        # Eq. 10.27 is printed with a minus between the density terms; the example adds them, which is the correct form (Eq. 10.26)
        assert c.gravity_gradient(567.0, 18.1, 0.674) == pytest.approx(1933, abs=1)
        assert 687.8 + c.gravity_gradient(567.0, 18.1, 0.674) == pytest.approx(2621, abs=1)


class TestExample10_5:
    """Flooding limit: pentane reflux, D = 3 cm: HSB p. 418-420, with Eq. 18.34 on p. 674."""

    def test_critical_mass_flux(self):
        assert c.flooding_critical_mass_flux(0.03, 529.0, 19.0, C=c.FLOODING_C["hsb_reflux"]) == pytest.approx(14.6, abs=0.05)

    def test_consistency_with_wallis_definitions(self):
        m = c.flooding_critical_mass_flux(0.03, 529.0, 19.0)
        assert c.wallis_flooding_sum(m / 529.0, m / 19.0, 0.03, 529.0, 19.0) == pytest.approx(0.75)
        assert c.flooding_vapour_velocity(0.03, 529.0, 19.0, 1.0, 0.75) * 19.0 == pytest.approx(m)

    def test_constants_scale_as_c_squared(self):
        low = c.flooding_critical_mass_flux(0.03, 529.0, 19.0, C=0.75)
        assert c.flooding_critical_mass_flux(0.03, 529.0, 19.0, C=c.FLOODING_C["lecture_24"]) / low == pytest.approx((0.8 / 0.75) ** 2)
        assert c.flooding_critical_mass_flux(0.03, 529.0, 19.0, C=c.FLOODING_C["hedh_square_end"]) < low


class TestExample11_2:
    """Forced-convective part of Chen, as HSB suggests for condensation: p. 451-452. The nucleate part is not applicable."""

    def stated(self):
        return fluid(rho_l=567.0, rho_v=18.09, mu_l=156e-6, mu_v=7.11e-6, k_l=0.086, cp_l=2730.0)

    def test_values(self):
        assert c.martinelli_xtt(0.3, 567.0, 18.09, 156e-6, 7.11e-6) == pytest.approx(0.552, abs=1e-3)
        assert c.chen_f(0.552) == pytest.approx(3.949, abs=5e-3)
        assert c.chen_forced_convective(self.stated(), 678.0, 0.3, 0.0221) == pytest.approx(4388, rel=1.5e-3)


class TestEnglishBevelAgainstHedhStatement:
    def test_increases(self):
        # HEDH 3.4.3: "increases in velocity of 5, 25, and 54% for cuts of 30, 60, and 75 deg". Checks the reading of the exponent 0.32. Rounded to 1%.
        assert c.english_bevel_factor(30.0) == pytest.approx(1.05, abs=0.005)
        assert c.english_bevel_factor(60.0) == pytest.approx(1.25, abs=0.005)
        assert c.english_bevel_factor(75.0) == pytest.approx(1.54, abs=0.005)
        assert c.english_bevel_factor(0.0) == 1.0

    def test_raw_scales_as_printed_only(self):
        base = c.english_flooding_raw(25.0, 1000.0, 1.2, 72.0, 1.0, 0.0, 1.0)
        assert c.english_flooding_raw(50.0, 1000.0, 1.2, 72.0, 1.0, 0.0, 1.0) / base == pytest.approx(2**0.3)
        assert c.english_flooding_raw(25.0, 1000.0, 1.2, 72.0, 1.0, 60.0, 1.0) / base == pytest.approx(c.english_bevel_factor(60.0))
        with pytest.raises(ValueError):
            c.english_bevel_factor(90.0)


class TestFigure16_6ReadOff:
    """Coarse comparison with the curve of HSB Fig. 16.6 (and HEDH Fig. 6). Read by eye from a log plot, so 10%; not a numerical benchmark."""

    def test_points(self):
        f = fluid(k_l=1.0, rho_l=1000.0, rho_v=1.0, mu_l=1e-3, cp_l=4000.0)  # Pr = 4
        scale = c.gravity_length(f) / f.k_l
        assert c.gravity_mean(f, 10.0) * scale == pytest.approx(0.68, rel=0.1)
        assert c.gravity_mean(f, 100.0) * scale == pytest.approx(0.33, rel=0.1)
        assert c.gravity_mean(f, 1600.0) * scale == pytest.approx(0.18, rel=0.1)


# ------------------------------------------------------------------------------------------- not benchmarks from here on
class TestLimitsGravity:
    """Limiting cases and hand values. No worked example exists for these branches."""

    @pytest.fixture
    def f(self):
        return fluid(k_l=0.65, rho_l=985.0, rho_v=0.1, mu_l=500e-6, cp_l=4180.0, h_fg=2.38e6)

    def test_mean_is_four_thirds_of_local_at_the_end(self, f):
        # HSB Eqs. 16.13, 16.14
        assert c.nusselt_mean(f, 10.0, 2.0) == pytest.approx(4 / 3 * c.nusselt_local(f, 10.0, 2.0))

    def test_mean_is_the_integral_of_local(self, f):
        n = 20000
        total = sum(c.nusselt_local(f, 10.0, (i + 0.5) * 2.0 / n) for i in range(n)) * 2.0 / n / 2.0
        assert total == pytest.approx(c.nusselt_mean(f, 10.0, 2.0), rel=2e-3)

    def test_inclination(self, f):
        assert c.nusselt_mean(f, 10.0, 2.0, inclination_deg=30.0) == pytest.approx(c.nusselt_mean(f, 10.0, 2.0) * 0.5**0.25)

    def test_film_thickness_gives_local_coefficient(self, f):
        # alpha = lambda/delta (HSB Eq. 16.13)
        gamma = c.condensate_per_width(c.nusselt_mean(f, 10.0, 2.0), 10.0, 2.0, f.h_fg)
        assert f.k_l / c.nusselt_film_thickness(f, gamma) == pytest.approx(c.nusselt_local(f, 10.0, 2.0), rel=1e-10)

    def test_laminar_branches_agree_with_nusselt_forms(self, f):
        # gravity_local and gravity_mean at Re < 30 equal the z-based and Gamma-based forms
        gamma = c.condensate_per_width(c.nusselt_mean(f, 0.2, 0.01), 0.2, 0.01, f.h_fg)
        re = c.film_reynolds(gamma, f.mu_l)
        assert re < 30
        assert c.gravity_local(f, re) == pytest.approx(c.nusselt_local(f, 0.2, 0.01), rel=1e-10)
        assert c.gravity_mean(f, re) == pytest.approx(c.nusselt_mean(f, 0.2, 0.01), rel=1e-10)

    def test_regime_selection(self):
        assert c.gravity_regime(29.0) == "laminar wave-free"
        assert c.gravity_regime(30.0) == "laminar wavy"
        assert c.gravity_regime(1599.0) == "laminar wavy"
        assert c.gravity_regime(1600.0) == "turbulent"

    @pytest.mark.parametrize("pr_cp", [1000.0, 4180.0, 40000.0])
    def test_mean_is_continuous_at_30_and_1600(self, pr_cp):
        # HEDH Eq. 19 meets Eq. 12 and Eq. 21. 0.5% at 30 because HEDH rounds its constants.
        g = fluid(k_l=0.65, rho_l=985.0, rho_v=0.1, mu_l=500e-6, cp_l=pr_cp)
        assert c.gravity_mean(g, 29.999) == pytest.approx(c.gravity_mean(g, 30.0), rel=5e-3)
        assert c.gravity_mean(g, 1599.999) == pytest.approx(c.gravity_mean(g, 1600.0), rel=1e-3)

    def test_mean_equals_hedh_eq_18_integral_of_local(self, f):
        # Re_L/alpha_mean = Integral dRe/alpha over wave-free + wavy: Eq. 19 is that integral with rounded constants. 0.5%.
        re_l = 800.0
        n = 40000
        total = 0.0
        for i in range(n):
            re = (i + 0.5) * re_l / n
            total += (re_l / n) / c.gravity_local(f, re)
        assert re_l / total == pytest.approx(c.gravity_mean(f, re_l), rel=5e-3)

    def test_eq_16_27_refused_below_its_range(self, f):
        # Book calls it "overall"; evaluated at Re = 100, Pr = 4 it gives 0.043 against 0.32. Shown here with the raw formula.
        assert 100 / (8750 + 58 * 4**-0.5 * (100**0.75 - 253)) == pytest.approx(0.043, abs=1e-3)
        with pytest.raises(ValueError):
            c.butterworth_mean_turbulent(f, 100.0)

    def test_prandtl_cap(self):
        # HEDH: use Pr = 10 when Pr > 10 in the turbulent equations
        a = fluid(k_l=0.1, rho_l=800.0, rho_v=1.0, mu_l=1e-3, cp_l=1000.0)  # Pr = 10
        b = fluid(k_l=0.1, rho_l=800.0, rho_v=1.0, mu_l=1e-3, cp_l=3000.0)  # Pr = 30
        assert c.labuntsov_local(a, 5000.0) == c.labuntsov_local(b, 5000.0)

    def test_superheat_raises_the_coefficient(self):
        # Not HSB Eq. 16.24's sign; see the docstring. Hand value: zeta = 2000(50)/2e6 = 0.05.
        assert c.superheat_correction(1000.0, 2000.0, 150.0, 100.0, 2e6) == pytest.approx(1000 * 1.05**0.25)
        assert c.superheat_correction(1000.0, 2000.0, 100.0, 100.0, 2e6) == 1000.0

    def test_subcooling_parameter(self):
        assert c.subcooling_parameter(4000.0, 10.0, 2e6) == pytest.approx(0.02)


class TestLimitsShear:
    def test_kosky_staub_joins_are_continuous(self):
        # Book's fits are continuous to about 0.1%
        assert c.kosky_staub_delta_plus(50.0) == pytest.approx(c.kosky_staub_delta_plus(50.001), rel=2e-3)
        assert c.kosky_staub_delta_plus(1483.0) == pytest.approx(c.kosky_staub_delta_plus(1483.001), rel=2e-3)
        for pr in (0.7, 4.0, 10.0):
            assert c.kosky_staub_t_plus(5.0, pr) == pytest.approx(c.kosky_staub_t_plus(5.0001, pr), rel=1e-3)
            assert c.kosky_staub_t_plus(30.0, pr) == pytest.approx(c.kosky_staub_t_plus(29.9999, pr), rel=1e-3)

    def test_kosky_staub_reduces_to_laminar_shear_film(self, pentane):
        # For Re <= 50 and delta+ <= 5, Eq. 16.52 with Eqs. 16.46, 16.49 is Eq. 16.33, since 0.7071 = 1/sqrt(2)
        gamma = 0.0015
        re = c.film_reynolds(gamma, pentane.mu_l)
        assert re < 50
        assert c.kosky_staub_coefficient(pentane, re, 5.0) == pytest.approx(c.laminar_shear_coefficient(pentane, gamma, 5.0), rel=1e-3)

    def test_butterworth_laminar_is_hsb_33(self, pentane):
        gamma = 0.0015
        re = c.film_reynolds(gamma, pentane.mu_l)
        # HEDH rounds sqrt(2) to 1.41: 0.3%
        assert c.butterworth_shear_coefficient(pentane, re, 0.5) == pytest.approx(c.laminar_shear_coefficient(pentane, gamma, 0.5), rel=4e-3)

    def test_butterworth_and_kosky_staub_agree_in_turbulent_ex_16_2(self, pentane):
        # Two independent routes at Ex. 16.2 conditions; agreement is a sanity check and not a benchmark
        re = c.film_reynolds(c.film_flow_in_tube(500.0, 0.85, 0.03), pentane.mu_l)
        ratio = c.butterworth_shear_coefficient(pentane, re, 55.2) / c.kosky_staub_coefficient(pentane, re, 55.2)
        assert ratio == pytest.approx(1.0, abs=0.06)

    def test_critical_reynolds(self):
        assert c.critical_film_reynolds(0.0) == 1600.0
        assert c.critical_film_reynolds(9.04) == pytest.approx(50.0, abs=5.0)  # Eq. 22 meets Eq. 23 near 50, Fig. 7
        assert c.critical_film_reynolds(20.0) == 50.0

    def test_dimensionless_shear_hand_value(self):
        f = fluid(rho_l=1000.0, rho_v=1.0, mu_l=1e-3)
        expected = 1000 * 2.0 / (1000 * 999 * 1e-3 * 9.81) ** (2 / 3)
        assert c.dimensionless_shear(f, 2.0) == pytest.approx(expected)

    def test_henstock_hanratty_hand_value(self, pentane):
        # Worked by hand from HEDH Eqs. 31-37 at the Ex. 16.2 conditions: 98 N/m2. Not a benchmark: HEDH has no worked example. 1% (hand rounding).
        assert c.henstock_hanratty_shear(pentane, 500.0, 0.85, 0.03) == pytest.approx(98.2, rel=0.01)
        with pytest.raises(ValueError):
            c.henstock_hanratty_shear(pentane, 500.0, 0.0, 0.03)

    def test_suction_limits(self):
        assert c.suction_corrected_shear(10.0, 500.0, 0.5, 5.0, 0.0) == pytest.approx(10.0)
        assert c.suction_corrected_shear(10.0, 500.0, 0.5, 5.0, 1.0) > 10.0
        assert c.mass_transfer_ratio(1e-12) == pytest.approx(1.0)
        with pytest.raises(ValueError):
            c.suction_corrected_shear(10.0, 500.0, 0.5, 5.0, -1.0)

    def test_ackermann_limits(self):
        assert c.ackermann_factor(0.0, 2000.0, 100.0) == pytest.approx(1.0)
        assert c.ackermann_factor(0.05, 2000.0, 100.0) > 1.0  # condensation raises the coefficient (HSB Sec. 16.5.2)

    def test_chen_low_parameter_rule(self):
        assert c.chen_f(100.0) == 1.0  # 1/X_tt < 0.1

    def test_mcnaught_hand_value(self):
        assert c.mcnaught_shear(100.0, 1.0) == pytest.approx(126.0)

    def test_akers_limits(self):
        assert c.akers_equivalent_mass_flux(100.0, 0.0, 800.0, 5.0, 10.0, 20.0) == pytest.approx(100.0)
        assert c.akers_equivalent_mass_flux(100.0, 0.5, 800.0, 5.0, 10.0, 10.0) == pytest.approx(50 + 50 * 160**0.5)

    def test_tube_row_factors(self):
        # HSB Eqs. 18.2 to 18.4; one row gives 1, and for many rows the stated ordering Nusselt < Kern < Eissenberg holds
        assert c.tube_row_factor(1, "nusselt") == c.tube_row_factor(1, "kern") == pytest.approx(1.0)
        assert c.tube_row_factor(1, "eissenberg") == pytest.approx(1.02)  # the book's fit gives 1.02 at n = 1
        assert c.tube_row_factor(16, "nusselt") == pytest.approx(0.5)
        assert c.tube_row_factor(16, "nusselt") < c.tube_row_factor(16, "kern") < c.tube_row_factor(16, "eissenberg")
        with pytest.raises(ValueError):
            c.tube_row_factor(0)


class TestLimitsCombined:
    def test_eq_16_56_limits(self):
        assert c.combined_gravity_shear(3.0, 4.0) == pytest.approx(5.0)
        assert c.combined_gravity_shear(100.0, 0.0) == 100.0
        assert c.combined_gravity_shear(0.0, 100.0) == 100.0

    def test_shear_dominates_at_high_mass_flux_gravity_at_low(self, pentane):
        high = c.vertical_downflow_local(pentane, 500.0, 0.85, 0.03)
        assert high.alpha_shear > 4 * high.alpha_gravity
        assert high.alpha == pytest.approx(math.hypot(high.alpha_shear, high.alpha_gravity))
        low = c.vertical_downflow_local(pentane, 5.0, 0.05, 0.03)
        assert low.alpha_gravity > 4 * low.alpha_shear
        assert low.regime == "laminar wavy"

    def test_coefficient_falls_with_condensation_at_high_flux(self, pentane):
        # Shear falls faster than the gravity term rises at this mass flux; each step down in x gives a lower alpha
        alphas = [c.vertical_downflow_local(pentane, 500.0, x, 0.03).alpha for x in (0.95, 0.8, 0.6, 0.4, 0.2, 0.05)]
        assert alphas == sorted(alphas, reverse=True)

    def test_all_methods_and_suction(self, pentane):
        for method in c.SHEAR_METHODS:
            result = c.vertical_downflow_local(pentane, 500.0, 0.5, 0.03, shear_method=method)
            assert result.alpha >= max(result.alpha_shear, result.alpha_gravity)
        plain = c.vertical_downflow_local(pentane, 500.0, 0.85, 0.03)
        corrected = c.vertical_downflow_local(pentane, 500.0, 0.85, 0.03, condensation_flux=1.16)
        assert corrected.tau_used > plain.tau_i
        assert corrected.alpha_shear > plain.alpha_shear

    def test_bad_input(self, pentane):
        with pytest.raises(ValueError):
            c.vertical_downflow_local(pentane, 500.0, 1.0, 0.03)
        with pytest.raises(ValueError):
            c.vertical_downflow_local(pentane, 500.0, 1.2, 0.03)
        with pytest.raises(ValueError):
            c.vertical_downflow_local(pentane, 500.0, 0.5, 0.03, shear_method="nope")
        no_sigma = properties.pentane_hsb_example_16_2()
        no_sigma.sigma = None
        with pytest.raises(ValueError):
            c.vertical_downflow_local(no_sigma, 500.0, 0.5, 0.03)

    def test_palen_ends(self):
        assert c.palen_interpolation(0.4, 10.0, 2.0) == 2.0
        assert c.palen_interpolation(0.5, 10.0, 2.0) == 2.0
        assert c.palen_interpolation(1.5, 10.0, 2.0) == 10.0
        assert c.palen_interpolation(1.0, 10.0, 2.0) == pytest.approx(6.0)


class TestLimitsIntegration:
    def test_mean_over_quality_analytic(self):
        # alpha = 1 + x over 0..1: arithmetic 1.5; area-weighted 1/ln 2
        assert c.mean_coefficient_over_quality(lambda x: 1 + x, 1.0, 0.0, n=2000, weighting="arithmetic") == pytest.approx(1.5, rel=1e-6)
        assert c.mean_coefficient_over_quality(lambda x: 1 + x, 1.0, 0.0, n=2000, weighting="area") == pytest.approx(1 / math.log(2), rel=1e-6)

    def test_mean_handles_x_equals_one_at_inlet(self, pentane):
        alpha = c.mean_coefficient_over_quality(lambda x: c.vertical_downflow_local(pentane, 500.0, x, 0.03).alpha, 1.0, 0.0, n=40)
        assert 3000.0 < alpha < 13000.0

    def test_boyko_kruzhilin_end_mean_vs_true_means(self, pentane):
        # HSB's two-point mean sits below the true arithmetic mean here because alpha(x) is concave, and above the area-weighted one (spread stated in docs)
        ends = c.boyko_kruzhilin_mean(pentane, 500.0, 0.03, 0.95, 0.05)
        fn = lambda x: c.boyko_kruzhilin_local(pentane, 500.0, x, 0.03)
        arithmetic = c.mean_coefficient_over_quality(fn, 0.95, 0.05, weighting="arithmetic")
        area = c.mean_coefficient_over_quality(fn, 0.95, 0.05, weighting="area")
        assert area < arithmetic and ends < arithmetic
        assert area == pytest.approx(ends, rel=0.05)

    def test_area_integral_constant_values(self):
        assert c.area_by_integration([0.0, 1000.0], [100.0, 100.0], [10.0, 10.0]) == pytest.approx(1.0)
        with pytest.raises(ValueError):
            c.area_by_integration([0.0, 1000.0], [100.0, 100.0], [10.0, -1.0])
        with pytest.raises(ValueError):
            c.area_by_integration([0.0], [100.0], [10.0])

    def test_quality_range_checks(self):
        with pytest.raises(ValueError):
            c.mean_coefficient_over_quality(lambda x: 1.0, 0.2, 0.5)
        with pytest.raises(ValueError):
            c.mean_coefficient_over_quality(lambda x: 1.0, 0.5, 0.2, weighting="other")


class TestLimitsTwoPhase:
    def test_friedel_single_phase_limits(self, pentane):
        pg = c.phase_gradients(pentane, 500.0, 0.0, 0.03)
        assert c.friedel_gradient(pentane, 500.0, 0.0, 0.03) == pytest.approx(pg.dpdz_lo)
        assert c.friedel_gradient(pentane, 500.0, 1.0, 0.03) == pytest.approx(pg.dpdz_go)

    def test_martinelli_identity(self):
        assert c.lockhart_martinelli_phi_g2(1.7) == pytest.approx(c.lockhart_martinelli_phi_l2(1.7) * 1.7**2)

    def test_xtt_matches_blasius_gradients(self, heptane_octane):
        # X from the gradients with Blasius friction equals the closed form
        pg = c.phase_gradients(heptane_octane, 407.0, 0.1, 0.0221)
        assert c.martinelli_parameter(pg.dpdz_l, pg.dpdz_g) == pytest.approx(c.martinelli_xtt(0.1, 567.0, 18.1, 156e-6, 7.1e-6))

    def test_homogeneous_momentum_flux(self):
        # With homogeneous void, Eq. 10.27's bracket reduces to 1/rho_H of Eq. 10.22
        x, rl, rv, m = 0.3, 800.0, 5.0, 100.0
        eps = c.void_fraction_homogeneous(x, rl, rv)
        assert c.momentum_flux(m, x, rl, rv, eps) == pytest.approx(m**2 * (x / rv + (1 - x) / rl))

    def test_momentum_limits(self):
        assert c.momentum_flux(100.0, 0.0, 800.0, 5.0, 0.0) == pytest.approx(100**2 / 800)
        assert c.momentum_flux(100.0, 1.0, 800.0, 5.0, 1.0) == pytest.approx(100**2 / 5)

    def test_acceleration_is_recovery_in_condensation(self):
        assert c.acceleration_pressure_drop(100.0, 0.9, 0.99, 0.1, 0.8, 800.0, 5.0) < 0

    def test_gravity_sign(self):
        up = c.gravity_gradient(800.0, 5.0, 0.5, 1.0)
        assert c.gravity_gradient(800.0, 5.0, 0.5, -1.0) == -up
        assert c.gravity_gradient(800.0, 5.0, 0.5, 0.0) == 0.0

    def test_premoli_negative_bracket_gives_no_slip(self, pentane):
        # Source gap: Ex. 16.2 conditions at x = 0.85 make the bracket negative. S = 1 is this library's choice.
        assert c.premoli_slip_ratio(pentane, 500.0, 0.85, 0.03) == 1.0

    def test_tube_drop_hand_sum(self, pentane):
        # Uniform quality: no acceleration; friction is gradient x length; gravity is rho_TP g L
        r = c.tube_two_phase_pressure_drop(pentane, 100.0, 0.02, 3.0, 0.4, 0.4, sin_alpha=-1.0, n=4)
        eps = c.void_fraction_premoli(pentane, 100.0, 0.4, 0.02)
        assert r.dp_acceleration == pytest.approx(0.0)
        assert r.dp_friction == pytest.approx(c.friedel_gradient(pentane, 100.0, 0.4, 0.02) * 3.0)
        assert r.dp_gravity == pytest.approx(-9.81 * 3.0 * (pentane.rho_v * eps + pentane.rho_l * (1 - eps)))
        assert r.dp_total == pytest.approx(r.dp_friction + r.dp_gravity)

    def test_tube_drop_condensing_recovers_momentum_and_converges(self, pentane):
        coarse = c.tube_two_phase_pressure_drop(pentane, 100.0, 0.02, 3.0, 0.95, 0.05, n=20)
        fine = c.tube_two_phase_pressure_drop(pentane, 100.0, 0.02, 3.0, 0.95, 0.05, n=200)
        assert fine.dp_acceleration < 0
        assert fine.dp_gravity == 0.0
        assert fine.dp_friction == pytest.approx(coarse.dp_friction, rel=2e-3)
        for friction in ("chisholm", "lockhart_martinelli"):
            for void in ("homogeneous", "lockhart_martinelli"):
                assert c.tube_two_phase_pressure_drop(pentane, 100.0, 0.02, 3.0, 1.0, 0.0, friction=friction, void=void, n=10).dp_friction > 0

    def test_tube_drop_custom_profile(self, pentane):
        r = c.tube_two_phase_pressure_drop(pentane, 100.0, 0.02, 3.0, 0.95, 0.05, x_of_z=lambda s: 0.5, n=4)
        assert set(r.x) == {0.5}

    def test_tube_drop_bad_options(self, pentane):
        with pytest.raises(ValueError):
            c.tube_two_phase_pressure_drop(pentane, 100.0, 0.02, 3.0, 0.9, 0.1, friction="nope")
        with pytest.raises(ValueError):
            c.tube_two_phase_pressure_drop(pentane, 100.0, 0.02, 3.0, 0.9, 0.1, void="nope")
        with pytest.raises(ValueError):
            c.tube_two_phase_pressure_drop(pentane, 100.0, 0.02, 0.0, 0.9, 0.1)


class TestLimitsFlooding:
    def test_hedh_c_squared(self):
        assert c.flooding_c_squared_hedh(0.072) == pytest.approx(0.53)
        # Pentane sigma 11.3 mN/m with n = 0.5: C^2 = 0.53 sqrt(0.157) = 0.21, i.e. C = 0.46, far below HSB's 0.75 for the same fluid
        assert c.flooding_c_squared_hedh(0.0113) == pytest.approx(0.53 * (0.0113 / 0.072) ** 0.5)
        assert c.flooding_c_squared_hedh(0.072, f2=c.FLOODING_F2_TABLE[(60, 25)]) == pytest.approx(0.53 * 1.2)

    def test_belter_is_dimensionless_consistent_and_inverts(self):
        # Hand value of the printed expression (not a benchmark: no worked example, and see docstring on pi/4)
        v_l, d, rl, rv, s, mu, mw = 0.05, 0.0254, 1000.0, 1.2, 0.072, 1e-3, 1e-3
        delta = rl - rv
        expected = 0.286 * (d**2 * 9.81 * delta / s) ** 0.26 * (v_l * math.pi / 4 * d * (9.81 * delta**3 / s**3) ** 0.25) ** -0.22 * 2 ** -0.18
        kg = c.belter_kutateladze(v_l, d, rl, rv, s, mu, mw)
        assert kg == pytest.approx(expected)
        v_g = c.belter_vapour_velocity(kg, rl, rv, s)
        assert v_g * rv**0.5 / (9.81 * s * delta) ** 0.25 == pytest.approx(kg)

    def test_rejects_inverted_densities(self):
        with pytest.raises(ValueError):
            c.flooding_critical_mass_flux(0.03, 10.0, 20.0)
        with pytest.raises(ValueError):
            c.wallis_star(1.0, 1.0, 0.03, 10.0, 20.0)

    def test_design_fraction(self):
        assert c.REFLUX_DESIGN_FRACTION == 0.75


class TestArguments:
    def test_missing_fields_are_named(self, pentane):
        g = fluid(rho_l=500.0, rho_v=5.0)
        with pytest.raises(ValueError, match="k_l"):
            c.nusselt_mean(g, 10.0, 1.0)
        with pytest.raises(ValueError):
            c.nusselt_mean(pentane, -5.0, 1.0)
        with pytest.raises(ValueError):
            c.film_flow_in_tube(100.0, 1.5, 0.02)
        with pytest.raises(ValueError):
            c.martinelli_xtt(0.0, 500.0, 5.0, 1e-4, 1e-5)
        with pytest.raises(ValueError):
            c.gravity_local(pentane, 0.0)
