"""Reports must build for every state a result can be in, with symbols the PDF font cannot draw."""
import io

import openpyxl
import pytest

from condenser_app.modules.selection import OperatingData, SelectionInputs
from condenser_app.reporting import _safe, input_rows, to_excel, to_pdf
from condenser_app.modules import MODULES
from condenser_app.schemas import CalcStep, Status
from condenser_app.service import run_design

GOOD = SelectionInputs(q1_coolant="utility", q2_cooling_water="no", data=OperatingData(pressure_vapour=2e5))


def test_excel_contains_results_and_notes():
    result = run_design("selection", GOOD)
    book = openpyxl.load_workbook(io.BytesIO(to_excel(result, GOOD)))
    assert {"Summary", "Inputs", "Results", "Calculation steps", "Checks", "Notes"} <= set(book.sheetnames)
    cells = [c.value for row in book["Results"].iter_rows() for c in row]
    assert "Air-cooled heat exchanger" in cells


def test_inputs_are_reported_in_display_units():
    pressure = next(r for r in input_rows(GOOD) if r["Field"] == "data.pressure_vapour")
    assert pressure["Value"] == pytest.approx(2.0) and pressure["Unit"] == "bar"


def test_pdf_builds():
    result = run_design("selection", GOOD)
    assert to_pdf(result, GOOD).startswith(b"%PDF")


def test_reports_build_for_invalid_input():
    bad = SelectionInputs()
    result = run_design("selection", bad)
    assert result.status is Status.INVALID_INPUT
    assert to_pdf(result, bad).startswith(b"%PDF")
    assert to_excel(result, bad)


def test_safe_text_handles_symbols_and_markup():
    assert _safe("ΔT ≥ 5 & ρ < 3") == "DeltaT &gt;= 5 &amp; rho &lt; 3"
    assert _safe(1234.56789) == "1234.6"


@pytest.mark.parametrize("module_key", list(MODULES))
def test_every_step_reads_correctly(module_key):
    """No repeated left-hand side such as 'x = x = ...', and every step ends in its result."""
    result = run_design(module_key, MODULES[module_key].example_inputs())
    assert result.steps
    for step in result.steps:
        chain = step.chain()
        assert chain, step.label
        if step.symbol:
            assert not chain.startswith(f"{step.symbol} = {step.symbol} ="), chain
            if step.value is not None:
                assert chain.splitlines()[-1].startswith(f"{step.symbol} = "), chain


def test_chain_when_the_formula_states_its_own_left_side():
    step = CalcStep("Overall", "U", "1/U = 1/h_o + 1/h_i", "h_o = 100, h_i = 100", 50.0, "W/m2.K")
    assert step.chain() == "1/U = 1/h_o + 1/h_i\n  with h_o = 100, h_i = 100\nU = 50 W/m2.K"


def test_chain_when_the_formula_is_only_the_right_side():
    step = CalcStep("Reynolds", "Re", "4 m / (pi D mu)", "4(0.2)/(pi 0.025 5.5e-4)", 18519.0)
    assert step.chain() == "Re = 4 m / (pi D mu)\n  = 4(0.2)/(pi 0.025 5.5e-4)\nRe = 18519"
