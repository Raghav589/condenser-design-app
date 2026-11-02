"""Extreme inputs must come back as a status or a warning, never as a crash or a non-real number."""
import copy
import math

import pytest
from scipy.integrate import solve_ivp

from condenser_app.common import condensation as c
from condenser_app.modules import MODULES, direct_contact as dc
from condenser_app.schemas import Status, iter_inputs
from condenser_app.service import run_design

FACTORS = (0.0, -1.0, 0.001, 1000.0)


def _set(inputs, path, value):
    owner = inputs
    parts = path.split(".")
    for part in parts[:-1]:
        owner = getattr(owner, part)
    setattr(owner, parts[-1], value)


def _cases():
    for key, module in MODULES.items():
        for path, _owner, _name, _meta in iter_inputs(module.example_inputs()):
            if not path.startswith("solver."):  # a huge step count is slow by request, not a defect
                yield pytest.param(key, path, id=f"{key}-{path}")


def _is_real(value):
    return not isinstance(value, complex) and (not isinstance(value, float) or math.isfinite(value))


@pytest.mark.parametrize("key, path", list(_cases()))
def test_extreme_values_of_any_input_do_not_crash(key, path):
    base = MODULES[key].example_inputs()
    _owner, current = next((o, getattr(o, n)) for p, o, n, _m in iter_inputs(base) if p == path)
    if isinstance(current, bool) or not isinstance(current, (int, float)):
        return
    for factor in FACTORS:
        inputs = copy.deepcopy(base)
        _set(inputs, path, current * factor)
        result = run_design(key, inputs)
        assert result.status in Status
        for output in result.outputs:
            assert _is_real(output.value), f"{path} x {factor}: {output.key} = {output.value}"
        for step in result.steps:
            assert _is_real(step.value), f"{path} x {factor}: {step.label} = {step.value}"


def test_friedel_refuses_vapour_more_viscous_than_liquid():
    fluid = copy.deepcopy(MODULES["vertical_shell_tube"].example_inputs().condensing)
    fluid.mu_l = fluid.mu_v / 2.0
    with pytest.raises(ValueError, match="vapour viscosity below liquid viscosity"):
        c.friedel_terms(fluid, 500.0, 0.5, 0.02)


def test_vertical_module_reports_the_viscosity_problem():
    inputs = MODULES["vertical_shell_tube"].example_inputs()
    inputs.condensing.mu_l = inputs.condensing.mu_v / 2.0
    result = run_design("vertical_shell_tube", inputs)
    assert result.status is Status.FAILED
    assert "viscosity" in result.message


def test_gravity_path_with_strong_drag_matches_a_stiff_integrator():
    """The coarse RK4 passes overshoot zero velocity here. The refined answer must match Radau, an independent stiff solver."""
    v_in, t, d, rho_g, mu_g, rho_l = 5.0, 0.353, 0.001, 0.596, 0.012, 963.0

    def rhs(_t, y):
        g_net, drag = dc._drop_rates(y[1], d, rho_g, mu_g, rho_l)
        return [y[1], g_net - drag]

    reference = solve_ivp(rhs, (0.0, t), [0.0, v_in], method="Radau", rtol=1e-9, atol=1e-12)
    assert reference.success
    x_coarse, v_coarse = dc.drop_path_with_gravity(v_in, t, d, rho_g, mu_g, rho_l, 32)
    assert not isinstance(x_coarse, complex) and not isinstance(v_coarse, complex)
    x_fine, v_fine = dc.drop_path_with_gravity(v_in, t, d, rho_g, mu_g, rho_l, 8192)
    assert x_fine == pytest.approx(reference.y[0, -1], rel=1e-4)
    assert v_fine == pytest.approx(reference.y[1, -1], rel=1e-4)


def test_spray_with_strong_drag_returns_a_result():
    inputs = MODULES["direct_contact"].example_inputs()
    inputs.vapour.mu_v *= 1000.0
    result = run_design("direct_contact", inputs)
    assert result.status in (Status.OK_WITH_WARNINGS, Status.NOT_CONVERGED)
