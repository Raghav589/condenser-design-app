"""Drives the real Streamlit app headlessly: every module must work through the form, not just through its API."""
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

from condenser_app.schemas import Status, iter_inputs
from condenser_app.service import get_module, available_modules, run_design

APP = str(Path(__file__).resolve().parent.parent / "app.py")
KEYS = list(available_modules())


def open_page(module_key: str) -> AppTest:
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert not at.exception
    at.sidebar.radio[0].set_value(get_module(module_key).title).run()
    assert not at.exception
    return at


def press(at: AppTest, key: str) -> AppTest:
    at.button(key=key).click().run()
    assert not at.exception
    return at


def test_home_lists_every_module():
    at = AppTest.from_file(APP, default_timeout=60).run()
    assert not at.exception
    labels = [b.label for b in at.button]
    for title in available_modules().values():
        assert title in labels


@pytest.mark.parametrize("module_key", KEYS)
def test_example_through_the_form_matches_the_direct_call(module_key):
    """Load example -> Calculate on the form gives the same result as the module on the same inputs."""
    module = get_module(module_key)
    at = open_page(module_key)
    press(at, f"example:{module_key}")
    press(at, f"calculate:{module_key}")

    form_result, form_inputs = at.session_state[f"result:{module_key}"]
    direct = run_design(module_key, module.example_inputs())
    assert form_result.status is direct.status
    assert form_result.status not in (Status.INVALID_INPUT, Status.FAILED), form_result.message or form_result.input_issues
    for form_q, direct_q in zip(form_result.outputs, direct.outputs, strict=True):
        assert form_q.key == direct_q.key
        if isinstance(direct_q.value, float):
            assert form_q.value == pytest.approx(direct_q.value, rel=1e-9), form_q.key
        else:
            assert form_q.value == direct_q.value, form_q.key


@pytest.mark.parametrize("module_key", KEYS)
def test_empty_form_reports_issues_without_crashing(module_key):
    at = open_page(module_key)
    press(at, f"calculate:{module_key}")
    result, _ = at.session_state[f"result:{module_key}"]
    assert result.status is Status.INVALID_INPUT
    assert result.input_issues


@pytest.mark.parametrize("module_key", KEYS)
def test_declared_defaults_survive_the_form(module_key):
    """A field with a declared default must reach the module unless the user changes it."""
    module = get_module(module_key)
    at = open_page(module_key)
    shown = module.relevant_fields(module.new_inputs())
    press(at, f"calculate:{module_key}")
    _, form_inputs = at.session_state[f"result:{module_key}"]
    for (path, owner, name, _), (_, form_owner, _, _) in zip(iter_inputs(module.new_inputs()), iter_inputs(form_inputs)):
        default = getattr(owner, name)
        if default in (None, "") or (shown is not None and path not in shown):  # blank text means "not set"
            continue
        expected = pytest.approx(default) if isinstance(default, float) else default
        assert getattr(form_owner, name) == expected, path


def test_reset_clears_the_form():
    at = open_page("selection")
    press(at, "example:selection")
    assert at.session_state["selection:q1_coolant"] == "utility"
    press(at, "reset:selection")
    assert at.selectbox(key="selection:q1_coolant").value == "(not set)"


def test_selection_only_asks_questions_on_the_path():
    at = open_page("selection")
    labels = [s.label for s in at.selectbox]
    assert any(label.startswith("1.") for label in labels) and not any(label.startswith("2.") for label in labels)
    at.selectbox(key="selection:q1_coolant").select("utility").run()
    labels = [s.label for s in at.selectbox]
    assert any(label.startswith("2.") for label in labels) and not any(label.startswith("3.") for label in labels)


def test_selection_links_to_the_design_module():
    at = open_page("selection")
    at.selectbox(key="selection:q1_coolant").select("process").run()
    at.selectbox(key="selection:q3_plate_feasible").select("yes").run()
    press(at, "calculate:selection")
    link = [b for b in at.button if b.label.startswith("Design this condenser")]
    assert len(link) == 1
    link[0].click().run()
    assert at.sidebar.radio[0].value == get_module("plate").title


@pytest.mark.parametrize("answers,module_key,field,value", [
    ({"q1_coolant": "utility", "q2_cooling_water": "yes", "q3_plate_feasible": "no", "q5_corrosive_dirty": "yes", "q6_multicomponent_total": "no", "q7_solids_countercurrent": "no"},
     "direct_contact", "condenser_type", "spray"),
    ({"q1_coolant": "utility", "q2_cooling_water": "yes", "q3_plate_feasible": "no", "q5_corrosive_dirty": "no", "q8_reflux": "yes"},
     "vertical_shell_tube", "variant", "reflux_upflow"),
])
def test_chart_outcome_preselects_the_variant(answers, module_key, field, value):
    at = open_page("selection")
    for name, answer in answers.items():
        at.selectbox(key=f"selection:{name}").select(answer).run()
    press(at, "calculate:selection")
    next(b for b in at.button if b.label.startswith("Design this condenser")).click().run()
    assert not at.exception
    assert at.sidebar.radio[0].value == get_module(module_key).title
    assert at.selectbox(key=f"{module_key}:{field}").value == value


def test_direct_contact_form_shows_only_the_chosen_type():
    at = open_page("direct_contact")
    labels = lambda: [n.label for n in at.number_input]
    assert not any("drop diameter" in x or "Packing size" in x for x in labels())
    at.selectbox(key="direct_contact:condenser_type").select("spray").run()
    assert any("drop diameter" in x for x in labels()) and not any("Packing size" in x for x in labels())
    at.selectbox(key="direct_contact:condenser_type").select("packed_bed").run()
    assert any("Packing size" in x for x in labels()) and not any("drop diameter" in x for x in labels())


def test_display_units_convert_on_entry():
    # 2257 kJ/kg on the form must reach the module as 2.257e6 J/kg.
    at = open_page("direct_contact")
    press(at, "example:direct_contact")
    assert at.session_state["direct_contact:vapour.h_fg"] == pytest.approx(2257.0)
    press(at, "calculate:direct_contact")
    _, inputs = at.session_state["result:direct_contact"]
    assert inputs.vapour.h_fg == pytest.approx(2.257e6)
