"""Preliminary mechanical sizing and vibration screening from HEDH Part 4.

BENCHMARK: HEDH 4.3.6, floating-head TEMA AJS, read from the page images. Inputs: P_shell 2000 kPa, P_tube 500 kPa, S = 121 MPa (SA-515-70), E = 0.85 (spot), corrosion 3 mm,
635 mm shell, 19.05 x 2.11 mm tube on 23.81 mm pitch, G = 678.28 mm. HEDH 4.6 has no worked example, so the vibration functions have hand values only.
"""
import math

import pytest

from condenser_app.common import mechanical as m


# ---------------------------------------------------------------------------- BENCHMARK, HEDH 4.3.6
S, E, C = 121e6, 0.85, 3e-3
P_SHELL, P_TUBE = 2000e3, 500e3


def test_hedh_4_3_6_shell_thickness_table_2():
    # R_v is the corroded radius: 635/2 + 3 = 320.5 mm. Book: shell 6.31 + 3 = 9.31, then the TEMA 9.52 minimum, 10 mm selected.
    t = m.cylinder_thickness_circumferential(P_SHELL, 320.5e-3, S, E)
    assert t * 1e3 == pytest.approx(6.31, abs=0.005)
    assert (t + C) * 1e3 == pytest.approx(9.31, abs=0.005)
    assert m.cylinder_thickness(P_SHELL, 320.5e-3, S, E) == t  # Eq. (9) can be neglected, as HEDH says


def test_hedh_4_3_6_channel_table_2():
    assert m.cylinder_thickness(P_TUBE, 320.5e-3, S, E) * 1e3 == pytest.approx(1.56, abs=0.005)


def test_hedh_4_3_6_shell_cover_differs_by_rounding_in_the_book():
    # Book prints 7.12 mm for R_v = 363 mm; its own Eq. (8) gives 2000 (363) / (121000 (0.85) - 1200) = 7.14 mm, +0.3%.
    assert m.cylinder_thickness(P_SHELL, 363e-3, S, E) * 1e3 == pytest.approx(7.12, rel=0.004)


def test_hedh_4_3_6_tube_internal_pressure():
    # R_v = 0.5 (19.05 - 2 x 2.11) = 7.415 mm, E = 1. Book 0.031 mm.
    assert m.tube_thickness_internal(P_TUBE, 7.415e-3, S) * 1e3 == pytest.approx(0.031, abs=0.0005)


def test_hedh_4_3_6_tube_external_pressure():
    # d_o/t = 9.03, B = 17 200 lbf/in2 from ASME, x 6.895 to kPa. Book 18 610 (its printed numbers give 18 581, -0.16%) and, with S = 206 MPa, 40 573 kPa.
    assert m.tube_external_allowable_a1(9.03, 17200 * 6.895) == pytest.approx(18610.0, rel=0.002)
    assert m.tube_external_allowable_a2(9.03, 206e6) == pytest.approx(40573e3, rel=0.001)


def test_hedh_4_3_6_tubesheet():
    # Book: t = (1.0 x 678.28 / 2) (2000 x 10^-3 / 121)^(1/2) = 43.6 mm. Overall 43.6 + 2 x 5 raised faces = 53.6, rounded up to 54 mm.
    assert m.tubesheet_bending(m.F_FLOATING_GASKETED, 678.28e-3, P_SHELL, S) * 1e3 == pytest.approx(43.6, abs=0.05)
    ratio, limit = m.tubesheet_shear_margin(P_SHELL, S, 19.05e-3, 23.81e-3)
    assert ratio == pytest.approx(0.0165, abs=5e-5)
    assert limit == pytest.approx(0.0639, abs=5e-4)
    assert ratio < limit  # "only bending need be considered"


def test_hedh_4_3_6_tema_minimum():
    assert m.tema_minimum_shell(0.635) == pytest.approx(9.52e-3, abs=0.03e-3)  # Table 2 prints 9.52; 4.2.6-1 prints 9.5


# --------------------------------------------------------------------------------- other rules
def test_design_stress_takes_the_smallest():
    assert m.design_stress(600e6, 255e6) == pytest.approx(170e6)  # 255/1.5 < 600/3
    assert m.design_stress(450e6, 300e6) == pytest.approx(150e6)  # 450/3 < 300/1.5
    assert m.design_stress(600e6, 300e6, sigma_r=150e6) == pytest.approx(100e6)  # 150/1.5
    assert m.design_stress(600e6, 300e6, sigma_t=55e6) == pytest.approx(50e6)  # 55/1.1


def test_longitudinal_never_governs_at_equal_efficiency():
    for p in (1e5, 1e6, 2e7):
        assert m.cylinder_thickness_longitudinal(p, 0.3, S, E) < m.cylinder_thickness_circumferential(p, 0.3, S, E)


def test_pressure_too_high_is_an_error():
    with pytest.raises(ValueError):
        m.cylinder_thickness_circumferential(200e6, 0.3, S, E)
    with pytest.raises(ValueError):
        m.cylinder_thickness(-1.0, 0.3, S, E)


@pytest.mark.parametrize("d_s, carbon, alloy", [(0.33, 9.5, 4.8), (0.737, 9.5, 4.8), (0.75, 11.1, 6.4), (0.991, 11.1, 6.4), (1.016, 12.7, 7.9), (2.032, 12.7, 7.9), (2.057, 12.7, 9.5), (2.54, 12.7, 9.5)])
def test_tema_minimum_table(d_s, carbon, alloy):
    assert m.tema_minimum_shell(d_s) == pytest.approx(carbon * 1e-3)
    assert m.tema_minimum_shell(d_s, alloy=True) == pytest.approx(alloy * 1e-3)


def test_tema_minimum_outside_the_table():
    assert m.tema_minimum_shell(0.3) is None
    assert m.tema_minimum_shell(2.6) is None


def test_tubesheet_constants_and_rough_estimate():
    assert m.F_U_TUBE_GASKETED == 1.25  # TEMA R-7.142 as quoted, HEDH 4.3.2I
    u = m.tubesheet_bending(m.F_U_TUBE_GASKETED, 0.6, 1e6, S)
    assert u == pytest.approx(1.25 * 0.6 / 2 * math.sqrt(1e6 / S))
    assert m.tubesheet_rough(0.8) == pytest.approx(0.08)
    assert m.tubesheet_shear(0.5, 19.05e-3, 23.81e-3, 2e6, S) == pytest.approx(0.31 * 0.5 / (1 - 19.05 / 23.81) * 2e6 / S)


# ------------------------------------------------------------------------------- standard sizes
def test_next_standard_shell():
    assert m.next_standard_shell(0.5) == 0.54
    assert m.next_standard_shell(0.54) == 0.54
    assert m.next_standard_shell(0.541) == 0.591
    assert m.next_standard_shell(3.1) is None


def test_table_2a_lists_have_the_table_lengths_and_rise():
    assert len(m.STANDARD_SHELLS_MM) == 36
    assert list(m.STANDARD_SHELLS_MM) == sorted(m.STANDARD_SHELLS_MM)
    for key, (full, reduced) in m.TABLE_2A_ONE_PASS.items():
        assert len(full) == len(reduced)
        assert len(full) == (24 if key[0] == 15.88 else 36)  # the 15.88 mm columns stop at 1524 mm
        assert list(full) == sorted(full) and list(reduced) == sorted(reduced)
        assert all(r < f for r, f in zip(reduced, full))


def test_table_2a_against_hedhs_own_examples():
    # HEDH 4.2.5H Example 1: 889 mm, 25.4 mm tube on 31.75 mm triangular, one pass: N_t = 671.
    assert m.table_2a_counts(25.4e-3, 31.75e-3, 0.889) == (671, 622)
    assert m.table_2a_counts(19.05e-3, 23.812e-3, 0.635) == (600, 559)
    assert m.table_2a_counts(15.88e-3, 19.84e-3, 1.6) is None  # the column ends at 1524 mm
    assert m.table_2a_counts(20.0e-3, 25.0e-3, 0.635) is None  # no such tube


def test_smallest_shell_holding_a_count():
    assert m.smallest_table_2a_shell(25.4e-3, 31.75e-3, 267) == 0.591
    assert m.smallest_table_2a_shell(25.4e-3, 31.75e-3, 268) == 0.635
    assert m.smallest_table_2a_shell(25.4e-3, 31.75e-3, 99999) is None


def test_typical_wall_settles_the_lecture_9_typo():
    assert m.typical_tube_wall(25.4e-3, True) == pytest.approx(2.77e-3)
    assert m.typical_tube_wall(25.4e-3, False) == pytest.approx(2.11e-3)
    assert m.typical_tube_wall(24.4e-3, True) is None  # the slide's misprint is not a tube
    assert m.typical_tube_wall(19.05e-3, True) == pytest.approx(2.11e-3)


# ------------------------------------------------------------------------------------- vibration
def test_metal_data_from_5_5_12():
    assert m.TUBE_METALS[m.CARBON_STEEL] == (199e9, 7900.0, True)
    assert m.TUBE_METALS[m.STAINLESS_304] == (190e9, 7900.0, False)
    assert m.TUBE_METALS[m.ADMIRALTY] == (107e9, 8550.0, False)
    assert m.TUBE_METALS[m.CUNI_90_10] == (118e9, 8900.0, False)


def test_added_mass_chart_reads():
    assert m.added_mass_coefficient(1.25, True) == 3.0  # held at the top of the chart
    assert m.added_mass_coefficient(1.40, True) == pytest.approx(2.25)
    assert m.added_mass_coefficient(1.40, False) == pytest.approx(2.17)
    assert m.added_mass_coefficient(1.50, False) == pytest.approx(1.88)
    assert m.added_mass_coefficient(2.0, False) == 1.71  # held at the end
    for r in (1.26, 1.3, 1.4, 1.5, 1.6):
        assert m.added_mass_coefficient(r, True) > m.added_mass_coefficient(r, False)  # staggered above in-line


def test_strouhal_chart_reads():
    assert m.strouhal_number(90, 1.25) == pytest.approx(0.415)
    assert m.strouhal_number(30, 1.25) == pytest.approx(0.18)
    assert m.strouhal_number(45, 1.5) == pytest.approx(0.565)


def test_frequency_constant():
    assert m.frequency_constant() == 10.0
    assert m.frequency_constant(2.0) == pytest.approx(15.97 / 2.0**1.86)
    assert m.frequency_constant(10.0) == pytest.approx(0.22, abs=0.01)  # Fig. 4.6.2(2) ends near 0.2 at lambda = 10


def test_vibration_hand_chain():
    # 19.05 x 2.11 mm carbon steel, water both sides, staggered at 1.25, span 0.6 m.
    d_o, d_i = 19.05e-3, 14.83e-3
    inertia = m.tube_moment_of_inertia(d_o, d_i)
    assert inertia == pytest.approx(math.pi * (d_o**4 - d_i**4) / 64)
    me = m.effective_mass(d_o, d_i, 7900.0, 1000.0, 1000.0, 3.0)
    assert me == pytest.approx(1000 * 1.7273e-4 + 7900 * 1.1230e-4 + 3 * 1000 * 2.8503e-4, rel=1e-3)  # about 1.91 kg/m
    fn = m.natural_frequency(10.0, 199e9, inertia, me, 0.6)
    assert fn == pytest.approx(10.0 / (2 * math.pi) * math.sqrt(199e9 * inertia / (me * 0.6**4)))
    assert 80.0 < fn < 100.0
    uc = m.critical_velocity(2.0, fn, d_o, me, 0.10, 1000.0)
    assert uc == pytest.approx(2.0 * fn * d_o * math.sqrt(me * 0.10 / (1000.0 * d_o**2)))
    assert m.vortex_shedding_frequency(0.4, 1.0, d_o) == pytest.approx(0.4 / d_o)
    assert m.buffeting_frequency(1.0, d_o, 0.8 * 23.81e-3, 23.81e-3) == pytest.approx(d_o / (0.8 * 23.81e-3**2) * (3.05 * (1 - d_o / 23.81e-3) ** 2 + 0.28))
    assert m.acoustic_frequency(2, 300.0, 1.0) == pytest.approx(300.0)


def test_span_change_of_80_percent_raises_frequency_as_hedh_says():
    # HEDH 4.6.2-3: limiting spans to 80% "increases the natural frequency by more than 50%".
    a = m.natural_frequency(10.0, 199e9, 1e-8, 1.0, 1.0)
    b = m.natural_frequency(10.0, 199e9, 1e-8, 1.0, 0.8)
    assert b / a == pytest.approx(1 / 0.8**2)
    assert b / a > 1.5


def test_interpolate_holds_the_ends():
    pts = ((1.0, 10.0), (2.0, 20.0))
    assert m.interpolate(pts, 0.5) == 10.0
    assert m.interpolate(pts, 1.5) == 15.0
    assert m.interpolate(pts, 3.0) == 20.0
