"""Preliminary mechanical sizing and flow-induced vibration screening.

Source: Heat Exchanger Design Handbook (HEDH) Part 4, read from the page images, and Table 5.5.12 for metal properties:
    4.1.1   design stress rule                        4.2.5, 4.2.6   tube, shell and tube count data (E. A. D. Saunders)
    4.3.2   index to the codes (M. Morris)            4.3.6          worked floating-head design, TEMA AJS (D. Harris)
    4.6     flow-induced vibration (J. M. Chenoweth)  3.3.5          rough tubesheet thickness, Eq. (1)
ASME Section VIII and TEMA are quoted here only as HEDH Part 4 quotes them (ASME 1977/1986 editions, TEMA 6th ed. 1978). This is preliminary sizing:
nothing in this file establishes code compliance.

Conventions: SI throughout, so pressures and stresses are both in Pa and lengths are in m. HEDH's factors of 10^-3 convert its kPa to MPa and disappear here.
``Tested`` in a docstring says whether the function reproduces HEDH 4.3.6 or has only hand values. HEDH 4.6 has no worked vibration example.
"""
from __future__ import annotations

import math
from bisect import bisect_right
from typing import Sequence


def _positive(**kwargs: float) -> None:
    for name, value in kwargs.items():
        if value <= 0:
            raise ValueError(f"{name} must be positive, got {value:g}.")


def interpolate(points: Sequence[tuple[float, float]], x: float) -> float:
    """Piecewise linear through (x, y) points in rising x. Held at the end values outside them: callers state the range of the chart."""
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    xs = [p[0] for p in points]
    k = bisect_right(xs, x) - 1
    (x0, y0), (x1, y1) = points[k], points[k + 1]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


# ----------------------------------------------------------------------------------- design stress
def design_stress(sigma_u: float, sigma_y: float, sigma_r: float | None = None, sigma_t: float | None = None) -> float:
    """S = smallest of sigma_u/3, sigma_Y/1.5, sigma_R/1.5, sigma_T/1.1. HEDH 4.1.1 Eq. (1).

    sigma_u ultimate tensile stress; sigma_Y yield or 0.2% proof stress; sigma_R stress for failure after 100 000 h at the design temperature; sigma_T stress for 1% permanent strain after 100 000 h.
    The creep terms are optional, for use above the creep range only. HEDH: the factors "are specified by the design code and may differ slightly" from these. Tested: hand value.
    """
    _positive(sigma_u=sigma_u, sigma_y=sigma_y)
    values = [sigma_u / 3.0, sigma_y / 1.5]
    if sigma_r is not None:
        values.append(sigma_r / 1.5)
    if sigma_t is not None:
        values.append(sigma_t / 1.1)
    return min(values)


# ------------------------------------------------------------------------------------- cylinders
def cylinder_thickness_circumferential(pressure: float, radius: float, stress: float, efficiency: float) -> float:
    """t = P R / (S E - 0.6 P), circumferential stress on longitudinal joints. HEDH 4.3.6 Eq. (8), quoting ASME UG-27; also HEDH 4.3.2A.

    R is the corroded inside radius and E the joint efficiency (0.85 for spot examination, ASME UW-12 as quoted by HEDH). The result excludes the corrosion allowance.
    Tested: HEDH 4.3.6 Table 2, channel 1.56 mm, shell 6.31 mm (cover: book 7.12, this gives 7.14).
    """
    _positive(pressure=pressure, radius=radius, stress=stress, efficiency=efficiency)
    denominator = stress * efficiency - 0.6 * pressure
    if denominator <= 0.0:
        raise ValueError("HEDH 4.3.6 Eq. (8): S E - 0.6 P is not positive; the pressure is too high for this allowable stress.")
    return pressure * radius / denominator


def cylinder_thickness_longitudinal(pressure: float, radius: float, stress: float, efficiency: float) -> float:
    """t = P R / (2 S E + 0.4 P), longitudinal stress on circumferential joints. HEDH 4.3.6 Eq. (9), quoting ASME UG-27. Never governs when both joints have the same E, as HEDH says. Tested: below Eq. (8) by hand."""
    _positive(pressure=pressure, radius=radius, stress=stress, efficiency=efficiency)
    return pressure * radius / (2.0 * stress * efficiency + 0.4 * pressure)


def cylinder_thickness(pressure: float, radius: float, stress: float, efficiency: float) -> float:
    """The greater of HEDH 4.3.6 Eqs. (8), (9). Excludes the corrosion allowance."""
    return max(cylinder_thickness_circumferential(pressure, radius, stress, efficiency), cylinder_thickness_longitudinal(pressure, radius, stress, efficiency))


def tube_thickness_internal(pressure: float, bore_radius: float, stress: float, efficiency: float = 1.0) -> float:
    """Tube wall for internal pressure: Eq. (8) with R the tube inside radius and E = 1 for seamless tube. HEDH 4.3.6 Sec. D, quoting ASME UG-31. Tested: HEDH 4.3.6, 0.031 mm."""
    return cylinder_thickness_circumferential(pressure, bore_radius, stress, efficiency)


def tube_external_allowable_a1(d_o_over_t: float, factor_b: float) -> float:
    """P_a1 = [2.167 / (d_o/t) - 0.0833] B. HEDH 4.3.6 Eq. (10), the ASME procedure for d_o/t < 10.

    B is a chart factor from the ASME material charts, which HEDH does not reproduce: it prints only its example's value, 17 200 lbf/in2, so the module cannot evaluate this. Pressure unit as B.
    Tested: HEDH 4.3.6 prints 18 610 kPa; its numbers give 18 581 kPa, -0.16%, a rounding in the book.
    """
    _positive(d_o_over_t=d_o_over_t, factor_b=factor_b)
    return (2.167 / d_o_over_t - 0.0833) * factor_b


def tube_external_allowable_a2(d_o_over_t: float, stress: float) -> float:
    """P_a2 = 2 S / (d_o/t) [1 - 1/(d_o/t)]. HEDH 4.3.6 Eq. (11), with S the lesser of twice the allowable and 0.9 x yield. Tested: HEDH 4.3.6, 40 573 kPa (S = 206 MPa)."""
    _positive(d_o_over_t=d_o_over_t, stress=stress)
    return 2.0 * stress / d_o_over_t * (1.0 - 1.0 / d_o_over_t)


# HEDH 4.2.6-1, TEMA class R minimum shell and head barrel thickness, as quoted. (largest D_s in the row [m], carbon and low alloy [m], alloy [m]).
TEMA_MIN_SHELL = ((0.737, 9.5e-3, 4.8e-3), (0.991, 11.1e-3, 6.4e-3), (2.032, 12.7e-3, 7.9e-3), (2.540, 12.7e-3, 9.5e-3))
TEMA_MIN_SHELL_RANGE = (0.330, 2.540)


def tema_minimum_shell(d_s: float, alloy: bool = False) -> float | None:
    """TEMA class R minimum barrel thickness as quoted in HEDH 4.2.6-1, or None outside 330 to 2540 mm, the range HEDH prints.

    HEDH's rows are 330-737, 762-991, 1016-2032, 2057-2540 mm; a diameter in a gap between rows takes the next row up. HEDH 4.3.6 Table 2 uses 9.52 mm at 635 mm, the same value before rounding.
    """
    if not TEMA_MIN_SHELL_RANGE[0] <= d_s <= TEMA_MIN_SHELL_RANGE[1]:
        return None
    for limit, carbon, other in TEMA_MIN_SHELL:
        if d_s <= limit:
            return other if alloy else carbon
    return None


# ------------------------------------------------------------------------------------ tubesheets
F_FLOATING_GASKETED = 1.0  # TEMA R-7.141 as quoted in HEDH 4.3.2I, 4.3.6-13
F_U_TUBE_GASKETED = 1.25  # TEMA R-7.142 as quoted in HEDH 4.3.2I


def tubesheet_bending(factor_f: float, gasket_diameter: float, pressure: float, stress: float) -> float:
    """T = (F G / 2) (P / S)^(1/2), effective thickness in bending. TEMA R-7.132 as quoted in HEDH 4.3.6 Eq. (17), 4.3.2I.

    F = 1.0 gasketed floating head, 1.25 gasketed U tube. G is the gasket mean diameter (U tube: gasket diameter or shell ID). P is the greater design pressure in HEDH 4.3.6, "maximum differential pressure" in 4.3.2I.
    Add corrosion and any groove. Tested: HEDH 4.3.6, 43.6 mm. Fixed tubesheets are not covered: HEDH 4.3.2I lists the clause sequence only and does not print the factors, so see tubesheet_rough.
    """
    _positive(factor_f=factor_f, gasket_diameter=gasket_diameter, pressure=pressure, stress=stress)
    return factor_f * gasket_diameter / 2.0 * math.sqrt(pressure / stress)


def tubesheet_shear_margin(pressure: float, stress: float, d_o: float, pitch: float) -> tuple[float, float]:
    """(P/S, 1.6 (1 - d_o/p)^2). Shear does not control while the first is below the second. HEDH 4.3.6-13, -14, TEMA R-7.133 as quoted. Tested: HEDH 4.3.6, 0.0165 and 0.0639."""
    _positive(pressure=pressure, stress=stress, d_o=d_o, pitch=pitch)
    return pressure / stress, 1.6 * (1.0 - d_o / pitch) ** 2


def tubesheet_shear(d_l: float, d_o: float, pitch: float, pressure: float, stress: float) -> float:
    """t = 0.31 D_L / (1 - d_o/p) (P/S). HEDH 4.3.6 Eq. (18), TEMA R-7.133 as quoted.

    HEDH does not define D_L anywhere in 4.3 and its example never evaluates this equation, so the module applies only the criterion of tubesheet_shear_margin. Tested: algebra only.
    """
    _positive(d_l=d_l, d_o=d_o, pitch=pitch, pressure=pressure, stress=stress)
    return 0.31 * d_l / (1.0 - d_o / pitch) * pressure / stress


def tubesheet_rough(d_s: float) -> float:
    """L_ts = 0.1 D_s: HEDH 3.3.5 Eq. (1), "roughly estimated (especially for preliminary designs)". Used for fixed tubesheets because the fixed-tubesheet formula is not printed in HEDH. Not a pressure calculation."""
    _positive(d_s=d_s)
    return 0.1 * d_s


# ------------------------------------------------------------------------------ standard sizes
# HEDH 4.2.5 Table 2a (fixed tubesheet), shell inside diameters [mm], read from the page images (idx 1551-1553).
STANDARD_SHELLS_MM = (203, 254, 305, 337, 387, 438, 489, 540, 591, 635, 686, 737, 787, 838, 889, 940, 991, 1067, 1143, 1219, 1295, 1372, 1448, 1524, 1600, 1676, 1753, 1829, 1981, 2134, 2286, 2438, 2591, 2743, 2896, 3048)
# One tube pass, triangular pitch: (D_o mm, pitch mm) -> (full count "0.0", reduced count "0.33"), in the order of STANDARD_SHELLS_MM. The reduced count allows for untubed areas at nozzles of D_s/3 (HEDH 4.2.5H).
TABLE_2A_ONE_PASS = {
    (15.88, 19.84): (
        (77, 125, 186, 230, 310, 401, 505, 621, 748, 869, 1019, 1180, 1354, 1539, 1737, 1946, 2167, 2520, 2901, 3308, 3742, 4203, 4690, 5205),
        (73, 119, 175, 216, 291, 376, 472, 579, 698, 810, 948, 1098, 1258, 1430, 1612, 1805, 2009, 2336, 2688, 3063, 3464, 3889, 4338, 4813),
    ),
    (19.05, 23.81): (
        (52, 85, 127, 158, 213, 276, 348, 428, 516, 600, 704, 816, 936, 1064, 1201, 1346, 1499, 1745, 2009, 2291, 2592, 2912, 3250, 3607, 3982, 4376, 4789, 5220, 6138, 7130, 8196, 9337, 10552, 11842, 13205, 14643),
        (50, 81, 120, 148, 200, 259, 325, 400, 481, 559, 655, 759, 870, 989, 1115, 1249, 1391, 1617, 1861, 2121, 2399, 2694, 3006, 3335, 3681, 4044, 4424, 4821, 5668, 6581, 7563, 8614, 9732, 10920, 12176, 13500),
    ),
    (19.05, 25.40): (
        (46, 75, 112, 138, 187, 243, 306, 376, 454, 527, 619, 717, 823, 935, 1056, 1183, 1318, 1533, 1765, 2014, 2278, 2559, 2856, 3170, 3500, 3846, 4209, 4587, 5394, 6266, 7203, 8206, 9274, 10407, 11606, 12869),
        (43, 71, 106, 130, 176, 227, 286, 351, 423, 491, 576, 667, 765, 869, 980, 1098, 1222, 1421, 1635, 1864, 2109, 2368, 2641, 2931, 3235, 3554, 3888, 4237, 4981, 5783, 6646, 7571, 8553, 9597, 10701, 11864),
    ),
    (25.40, 31.75): (
        (28, 47, 70, 87, 115, 153, 193, 238, 287, 334, 392, 455, 522, 594, 671, 752, 838, 975, 1123, 1282, 1451, 1630, 1820, 2020, 2231, 2452, 2684, 2926, 3441, 3998, 4598, 5238, 5921, 6645, 7412, 8220),
        (27, 44, 66, 82, 110, 143, 180, 220, 267, 311, 365, 423, 485, 551, 622, 697, 776, 903, 1040, 1186, 1342, 1508, 1680, 1867, 2061, 2265, 2479, 2701, 3177, 3689, 4241, 4832, 5460, 6127, 6833, 7576),
    ),
    (31.75, 39.69): (
        (17, 29, 43, 54, 74, 96, 122, 150, 181, 211, 248, 288, 331, 377, 426, 478, 532, 620, 715, 816, 924, 1038, 1159, 1287, 1422, 1563, 1711, 1866, 2195, 2551, 2934, 3343, 3779, 4243, 4732, 5249),
        (16, 27, 41, 51, 69, 90, 113, 140, 169, 196, 231, 268, 307, 350, 395, 443, 493, 574, 661, 755, 854, 960, 1071, 1189, 1313, 1443, 1580, 1722, 2026, 2353, 2706, 3083, 3484, 3911, 4362, 4837),
    ),
}
TUBE_LENGTH_RANGE = {"fixed": (3.6, 15.0), "removable": (3.6, 9.0)}  # m, HEDH 4.2.5-2, "the most common tube length range"; up to 22 m has been made


def next_standard_shell(d_s: float) -> float | None:
    """The smallest shell ID of HEDH 4.2.5 Table 2a at or above d_s [m], or None above 3048 mm."""
    for size in STANDARD_SHELLS_MM:
        if size * 1e-3 >= d_s - 1e-9:
            return size * 1e-3
    return None


def table_2a_counts(d_o: float, pitch: float, shell: float) -> tuple[int, int] | None:
    """(full, reduced) one-pass count of HEDH 4.2.5 Table 2a for a triangular layout, or None if the table has no such tube, pitch or shell. 0.05 mm tolerance on the tube and pitch."""
    index = next((k for k, size in enumerate(STANDARD_SHELLS_MM) if abs(size * 1e-3 - shell) < 1e-6), None)
    if index is None:
        return None
    for (od, p), (full, reduced) in TABLE_2A_ONE_PASS.items():
        if abs(od - d_o * 1e3) <= 0.05 and abs(p - pitch * 1e3) <= 0.05:
            return (full[index], reduced[index]) if index < len(full) else None
    return None


def smallest_table_2a_shell(d_o: float, pitch: float, n_tubes: int) -> float | None:
    """The smallest Table 2a shell [m] whose reduced count holds n_tubes, or None."""
    for size in STANDARD_SHELLS_MM:
        counts = table_2a_counts(d_o, pitch, size * 1e-3)
        if counts is not None and counts[1] >= n_tubes:
            return size * 1e-3
    return None


# Lecture 9 table, confirmed by HEDH 4.2.5-1: OD [mm] -> (carbon and low alloy, stainless and copper alloys) wall [mm]. Lecture 9 prints 24.40 and 25.04 for the third tube; HEDH prints 25.40 mm (1 in), BWG 12 and 14.
TUBE_WALL_TABLE_MM = {15.88: (1.65, 1.24), 19.05: (2.11, 1.65), 25.40: (2.77, 2.11), 31.75: (3.40, 2.77)}


def typical_tube_wall(d_o: float, carbon_or_low_alloy: bool) -> float | None:
    """Typical wall [m] for the OD, Lecture 9 and HEDH 4.2.5-1, or None if the OD is not in the table (0.05 mm tolerance)."""
    for od, (carbon, other) in TUBE_WALL_TABLE_MM.items():
        if abs(od - d_o * 1e3) <= 0.05:
            return (carbon if carbon_or_low_alloy else other) * 1e-3
    return None


# ------------------------------------------------------------------------------ flow-induced vibration
CARBON_STEEL, STAINLESS_304, ADMIRALTY, CUNI_90_10 = "carbon steel", "18Cr/8Ni stainless (S30400)", "Admiralty brass (C44300)", "copper nickel 90/10 (C70600)"
# HEDH 5.5.12 Tables 1, 2, 4, the 100 degC column: modulus of elasticity [Pa], density [kg/m3], is it carbon or low alloy (for the wall table). Read from the page images.
TUBE_METALS = {
    CARBON_STEEL: (199e9, 7900.0, True),
    STAINLESS_304: (190e9, 7900.0, False),
    ADMIRALTY: (107e9, 8550.0, False),
    CUNI_90_10: (118e9, 8900.0, False),
}

# HEDH Fig. 4.6.2(3), added mass coefficient, pixel-digitised, +-0.02. Pitch ratio -> k. Staggered above in-line. Held outside 1.25 to 1.6, where the chart ends.
ADDED_MASS_STAGGERED = ((1.255, 3.0), (1.26, 2.94), (1.28, 2.79), (1.30, 2.68), (1.33, 2.53), (1.35, 2.44), (1.40, 2.25), (1.45, 2.08), (1.50, 1.95), (1.55, 1.84), (1.60, 1.77))
ADDED_MASS_IN_LINE = ((1.25, 2.90), (1.26, 2.82), (1.28, 2.69), (1.30, 2.59), (1.33, 2.44), (1.35, 2.36), (1.40, 2.17), (1.45, 2.01), (1.50, 1.88), (1.55, 1.77), (1.60, 1.71))
# HEDH Fig. 4.6.4(2), Strouhal number by layout angle, read by eye against the gridlines, +-0.03.
STROUHAL = {
    30: ((1.25, 0.18), (1.33, 0.23), (1.50, 0.285), (1.75, 0.345), (2.00, 0.347)),
    45: ((1.25, 0.587), (1.33, 0.597), (1.50, 0.565), (1.75, 0.455), (2.00, 0.325)),
    90: ((1.25, 0.415), (1.33, 0.36), (1.50, 0.30), (1.75, 0.275), (2.00, 0.278)),
}
# Chart reads, conservative. C_n: the lowest curve of Fig. 4.6.2(2) (N > 5) on its plateau, whatever the span count; the N = 3 to 5 curves lie above it, 13 to 11 at lambda = 1.
FREQUENCY_CONSTANT_LOWEST = 10.0
# beta: the lower edge of the lowest band of Fig. 4.6.4(4), about 2.0 to 2.15 over pitch ratios 1.25 to 1.75, +-0.1. The scan does not show clearly which upper band belongs to which layout, so the lowest one is applied to every layout.
BETA_LOWEST = 2.0
LOG_DECREMENT = {"gas": 0.03, "liquid": 0.10}  # HEDH 4.6.4C: 0.03 for gases and vapours, 0.10 for liquids and two-phase mixtures


def added_mass_coefficient(pitch_ratio: float, staggered: bool) -> float:
    """k of M_s = k rho_s (pi/4) D_o^2. HEDH Fig. 4.6.2(3), from the TEMA Standards 6th ed. as quoted; digitised, +-0.02."""
    return interpolate(ADDED_MASS_STAGGERED if staggered else ADDED_MASS_IN_LINE, pitch_ratio)


def strouhal_number(layout_angle: int, pitch_ratio: float) -> float:
    """Sr = f_vs D_o / u_c. HEDH Fig. 4.6.4(2), from Chen; read by eye, +-0.03. HEDH says few data lie between pitch ratios 1.25 and 1.75."""
    return interpolate(STROUHAL[layout_angle], pitch_ratio)


def frequency_constant(span_ratio: float = 1.0) -> float:
    """C_n. For lambda > 1.2, C_n = 15.97/lambda^1.86, printed on HEDH Fig. 4.6.2(2). Below that the lowest curve, 10.0, conservative for every span count."""
    _positive(span_ratio=span_ratio)
    return 15.97 / span_ratio**1.86 if span_ratio > 1.2 else FREQUENCY_CONSTANT_LOWEST


def tube_moment_of_inertia(d_o: float, d_i: float) -> float:
    """I = pi (D_o^4 - D_i^4)/64. HEDH 4.6 nomenclature."""
    return math.pi * (d_o**4 - d_i**4) / 64.0


def effective_mass(d_o: float, d_i: float, rho_metal: float, rho_tube_fluid: float, rho_shell_fluid: float, k: float) -> float:
    """M_e = M_t + M_w + M_s [kg/m]: tube-side fluid, tube, and shell-side fluid displaced with M_s = k rho_s (pi/4) D_o^2. HEDH 4.6.2A, Fig. 3. Tested: hand value."""
    _positive(d_o=d_o, d_i=d_i, rho_metal=rho_metal, rho_tube_fluid=rho_tube_fluid, rho_shell_fluid=rho_shell_fluid, k=k)
    return rho_tube_fluid * math.pi / 4.0 * d_i**2 + rho_metal * math.pi / 4.0 * (d_o**2 - d_i**2) + k * rho_shell_fluid * math.pi / 4.0 * d_o**2


def natural_frequency(c_n: float, e: float, i: float, m_e: float, span: float) -> float:
    """f_n = (C_n/2 pi) [E I/(M_e L^4)]^0.5. HEDH 4.6.2 Eq. (1), Franklin et al. Axial load is not included (Eq. 2 needs it, and HEDH warns it can move f_n by 40%). Tested: hand value; HEDH 4.6 has no worked example."""
    _positive(c_n=c_n, e=e, i=i, m_e=m_e, span=span)
    return c_n / (2.0 * math.pi) * math.sqrt(e * i / (m_e * span**4))


def critical_velocity(beta: float, f_n: float, d_o: float, m_e: float, log_decrement: float, rho_shell: float) -> float:
    """u_crit = beta f_n D_o [M_e delta_o/(rho_s D_o^2)]^0.5. HEDH 4.6.4 Eq. (4), Connors. Tested: hand value."""
    _positive(beta=beta, f_n=f_n, d_o=d_o, m_e=m_e, log_decrement=log_decrement, rho_shell=rho_shell)
    return beta * f_n * d_o * math.sqrt(m_e * log_decrement / (rho_shell * d_o**2))


def vortex_shedding_frequency(strouhal: float, u_c: float, d_o: float) -> float:
    """f_vs = Sr u_c/D_o. HEDH 4.6.4 Eq. (2)."""
    _positive(strouhal=strouhal, d_o=d_o)
    return strouhal * u_c / d_o


def buffeting_frequency(u_c: float, d_o: float, p_l: float, p_t: float) -> float:
    """f_tb = (u_c D_o/(p_l p_t)) [3.05 (1 - D_o/p_t)^2 + 0.28]. HEDH 4.6.4 Eq. (3), Owens. HEDH: developed from gas data, "may not be applicable" to a liquid. Tested: hand value."""
    _positive(d_o=d_o, p_l=p_l, p_t=p_t)
    return u_c * d_o / (p_l * p_t) * (3.05 * (1.0 - d_o / p_t) ** 2 + 0.28)


def acoustic_frequency(mode: int, u_sound: float, d: float) -> float:
    """f_a = m u_sound/(2 d), usually with d the shell diameter. HEDH 4.6.4 Eq. (6). A vapour or gas on the shell side only. Tested: hand value."""
    _positive(mode=mode, u_sound=u_sound, d=d)
    return mode * u_sound / (2.0 * d)
