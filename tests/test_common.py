"""Shared library checks against values worked by hand or in the ME6220 tutorials."""
import math

import pytest

from condenser_app.common import iteration, singlephase as sp, thermal, units
from condenser_app.schemas import ConvergenceError


class TestUnits:
    def test_flow_and_pressure(self):
        assert units.convert(18, "kg/min", "kg/s") == pytest.approx(0.3)
        assert units.convert(3, "bar", "Pa") == pytest.approx(3e5)
        assert units.convert(250, "kJ/kg", "J/kg") == pytest.approx(250e3)

    def test_temperature(self):
        assert units.c_to_k(25) == pytest.approx(298.15)
        assert units.convert(298.15, "K", "degC") == pytest.approx(25)

    def test_incompatible_units_rejected(self):
        with pytest.raises(ValueError):
            units.convert(1, "bar", "kg/s")
        with pytest.raises(ValueError):
            units.convert(1, "furlong", "m")


class TestThermal:
    def test_lmtd_counterflow_tutorial_3_q3(self):
        assert thermal.lmtd_counterflow(357, 220, 30, 200) == pytest.approx(172.98, abs=0.01)

    def test_lmtd_equal_ends(self):
        assert thermal.lmtd(40, 40) == pytest.approx(40)

    def test_lmtd_rejects_temperature_cross(self):
        with pytest.raises(ValueError):
            thermal.lmtd(-5, 20)

    def test_isothermal_lmtd_hsb_example_3_10(self):
        # Condenser at 80 degC, water 20 -> 35 degC: HSB Example 3.10
        assert thermal.lmtd_isothermal(80, 20, 35) == pytest.approx(52.14, abs=0.01)

    def test_overall_u_hsb_example_3_1(self):
        # D_o = 25, D_i = 20 mm, k = 100, oil inside 200 + R_f 0.0008, water outside 5000 + R_f 0.0002: U = 130.2
        u = thermal.overall_u_outer(5000, 200, 0.025, 0.020, 100, r_f_o=0.0002, r_f_i=0.0008)
        assert u == pytest.approx(130.2, abs=0.1)

    def test_condensing_duty(self):
        assert thermal.condensing_duty(1.0, 337.3e3, 0.85, 0.05) == pytest.approx(269840)

    def test_effectiveness_limit(self):
        assert thermal.effectiveness_isothermal(0) == 0
        assert thermal.effectiveness_isothermal(50) == pytest.approx(1)


class TestSinglePhase:
    def test_reynolds_forms_agree(self):
        rho, u, d, mu = 820.0, 3.7266, 0.025, 0.0075
        m_dot = rho * u * math.pi * d**2 / 4
        assert sp.reynolds(rho, u, d, mu) == pytest.approx(sp.reynolds_massflow(m_dot, d, mu))

    def test_dittus_boelter_tutorial_2_q6_oil(self):
        assert sp.dittus_boelter(10185.9, 122, heating=False) == pytest.approx(156.3, abs=0.1)

    def test_friction_selection_by_re(self):
        assert sp.darcy_f(1000) == pytest.approx(0.064)
        assert sp.darcy_f(1e4) == pytest.approx(0.316e-1, abs=1e-5)
        assert sp.darcy_f(1.17e5) == pytest.approx(0.184 * 1.17e5**-0.2)

    def test_pipe_dp_hsb_example_2_12(self):
        f = 4 * 0.046 * 1.17e5**-0.2
        assert sp.pipe_dp(f, 100, 0.05, 979, 1.0) == pytest.approx(17470, rel=0.01)


class TestIteration:
    def test_fixed_point_converges(self):
        x, info = iteration.fixed_point(lambda v: math.sqrt(2 + v), 1.0, name="sqrt", tol=1e-8)
        assert x == pytest.approx(2.0)
        assert info.converged and info.iterations > 1 and info.history[0] == 1.0

    def test_fixed_point_failure_raises_with_history(self):
        with pytest.raises(ConvergenceError) as err:
            iteration.fixed_point(lambda v: v + 1, 0.0, name="diverge", max_iter=5)
        assert err.value.info.converged is False and len(err.value.info.history) == 6

    def test_bracketed_root(self):
        x, info = iteration.solve_bracketed(lambda v: v**2 - 2, 0, 2, name="root")
        assert x == pytest.approx(math.sqrt(2)) and info.converged

    def test_bracket_without_root_raises(self):
        with pytest.raises(ConvergenceError):
            iteration.solve_bracketed(lambda v: v**2 + 1, 0, 2, name="no root")
