"""Every outcome box of HSB Fig. 17.1 must be reachable by its own path through the chart."""
import pytest

from condenser_app.modules.selection import OUTCOME_TO_MODULE, OperatingData, SelectionInputs, walk_chart
from condenser_app.schemas import Status
from condenser_app.service import get_module, run_design

PATHS = {
    "Air-cooled heat exchanger": dict(q1_coolant="utility", q2_cooling_water="no"),
    "Plate heat exchanger": dict(q1_coolant="process", q3_plate_feasible="yes"),
    "Plate-fin exchanger": dict(q1_coolant="process", q3_plate_feasible="no", q4_clean_close="yes"),
    "Baffle column condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="yes", q6_multicomponent_total="no", q7_solids_countercurrent="yes"),
    "Spray condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="yes", q6_multicomponent_total="no", q7_solids_countercurrent="no"),
    "Vertical tube-side reflux condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="no", q8_reflux="yes"),
    "Vertical shell-side condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="no", q8_reflux="no", q9_horizontal="no", q10_high_p_corrosive="no", q11_low_dp="no", q12_coolant_tube_side="yes"),
    "Vertical tube-side condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="no", q8_reflux="no", q9_horizontal="no", q10_high_p_corrosive="yes"),
    "Horizontal tube-side condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="no", q8_reflux="no", q9_horizontal="yes", q10_high_p_corrosive="yes"),
    "Horizontal shell-side (E-shell) condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="no", q8_reflux="no", q9_horizontal="yes", q10_high_p_corrosive="no", q11_low_dp="no"),
    "Horizontal rod-baffle condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="no", q8_reflux="no", q9_horizontal="yes", q10_high_p_corrosive="no", q11_low_dp="yes", q13_counterflow="yes"),
    "Horizontal J-shell condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="no", q8_reflux="no", q9_horizontal="yes", q10_high_p_corrosive="no", q11_low_dp="yes", q13_counterflow="no", q14_wide_range_venting="yes"),
    "Horizontal X-shell condenser": dict(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="no", q8_reflux="no", q9_horizontal="yes", q10_high_p_corrosive="no", q11_low_dp="yes", q13_counterflow="no", q14_wide_range_venting="no"),
}


def test_every_chart_outcome_has_a_test_path():
    assert set(PATHS) == set(OUTCOME_TO_MODULE)


@pytest.mark.parametrize("outcome,answers", PATHS.items())
def test_chart_path(outcome, answers):
    assert walk_chart(SelectionInputs(**answers))[0] == outcome


def test_multicomponent_total_condensation_returns_to_shell_and_tube():
    # Question 6 yes: corrosive vapour, but direct contact is unsuitable, so shell-and-tube questions follow.
    inputs = SelectionInputs(q1_coolant="utility", q2_cooling_water="yes", q3_plate_feasible="no", q5_corrosive_dirty="yes", q6_multicomponent_total="yes", q8_reflux="yes")
    assert walk_chart(inputs)[0] == "Vertical tube-side reflux condenser"


def test_coolant_preference_decides_vertical_when_low_dp_also_answered_yes():
    inputs = SelectionInputs(**{**PATHS["Vertical shell-side condenser"], "q11_low_dp": "yes", "q12_coolant_tube_side": "no"})
    assert walk_chart(inputs)[0] == "Vertical tube-side condenser"


def test_process_coolant_skips_cooling_water_question():
    outcome, path = walk_chart(SelectionInputs(**PATHS["Plate heat exchanger"]))
    assert "q2_cooling_water" not in [name for name, _ in path]


def test_missing_required_answer_is_invalid_input():
    result = run_design("selection", SelectionInputs(q1_coolant="utility"))
    assert result.status is Status.INVALID_INPUT
    assert "cooling water" in result.input_issues[0].message


def test_answer_needed_only_on_this_path_is_an_input_issue():
    # q2 is optional in general, but this path needs it: the user is told which question to answer.
    result = run_design("selection", SelectionInputs(q1_coolant="utility", q3_plate_feasible="no", q5_corrosive_dirty="no"))
    assert result.status is Status.INVALID_INPUT
    assert "cooling water" in result.input_issues[0].message


def test_invalid_choice_rejected():
    result = run_design("selection", SelectionInputs(q1_coolant="maybe", q3_plate_feasible="no", q5_corrosive_dirty="no"))
    assert result.status is Status.INVALID_INPUT


def test_result_traces_the_path_and_maps_to_module():
    result = run_design("selection", SelectionInputs(**PATHS["Vertical tube-side condenser"]))
    assert result.get("recommended_type") == "Vertical tube-side condenser"
    assert result.get("design_module") == "vertical_shell_tube"
    assert len(result.steps) == 7


def test_outcome_without_a_module_warns():
    result = run_design("selection", SelectionInputs(**PATHS["Horizontal J-shell condenser"]))
    assert result.get("design_module") == "none"
    assert any("no design module" in w for w in result.warnings)


def test_plate_limits_warning():
    inputs = SelectionInputs(**PATHS["Plate heat exchanger"], data=OperatingData(pressure_vapour=40e5, temperature_max=250.0, flammable_or_toxic="yes"))
    result = run_design("selection", inputs)
    assert result.get("recommended_type") == "Plate heat exchanger"
    assert any("30 bar" in w and "gasket" in w for w in result.warnings)


def test_allocation_freezing_goes_to_shell():
    inputs = SelectionInputs(**PATHS["Horizontal shell-side (E-shell) condenser"], data=OperatingData(freezing_risk="yes"))
    result = run_design("selection", inputs)
    advice = [q.value for q in result.outputs if q.key.startswith("allocation_")]
    assert any("freezing" in a.lower() and "shell" in a for a in advice)


def test_air_cooled_path_needs_no_later_answers():
    result = run_design("selection", SelectionInputs(q1_coolant="utility", q2_cooling_water="no"))
    assert result.get("recommended_type") == "Air-cooled heat exchanger"
    assert result.get("design_module") == "air_cooled"


def test_form_shows_only_questions_on_the_path():
    module = get_module("selection")
    start = module.relevant_fields(SelectionInputs())
    assert "q1_coolant" in start and "q2_cooling_water" not in start
    utility = module.relevant_fields(SelectionInputs(q1_coolant="utility"))
    assert {"q1_coolant", "q2_cooling_water"} <= utility and "q3_plate_feasible" not in utility
    process = module.relevant_fields(SelectionInputs(q1_coolant="process"))
    assert {"q1_coolant", "q3_plate_feasible"} <= process and "q2_cooling_water" not in process
    deep = module.relevant_fields(SelectionInputs(q1_coolant="process", q3_plate_feasible="no", q4_clean_close="no", q5_corrosive_dirty="no", q8_reflux="no", q9_horizontal="yes", q10_high_p_corrosive="no"))
    assert {"q9_horizontal", "q10_high_p_corrosive", "q11_low_dp"} <= deep and "q12_coolant_tube_side" not in deep


def test_result_states_its_design_basis():
    result = run_design("selection", SelectionInputs(**PATHS["Spray condenser"]))
    assert any(a.startswith("Design basis: HSB Fig. 17.1") for a in result.design_basis)
    text = " ".join(result.design_basis + result.scope_limits).lower()
    assert "supplied" not in text and "project drive" not in text and "not available" not in text


def test_stale_answer_off_the_path_is_hidden():
    # q2 answered earlier, then q1 changed to process: q2 is no longer on the path.
    module = get_module("selection")
    assert "q2_cooling_water" not in module.relevant_fields(SelectionInputs(q1_coolant="process", q2_cooling_water="no"))
