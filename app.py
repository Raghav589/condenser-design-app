"""ME6220 condenser design application. Run with:  streamlit run app.py"""
from __future__ import annotations

import copy
from dataclasses import fields, is_dataclass
from pathlib import Path

import pandas as pd
import streamlit as st

from condenser_app.common.units import convert
from condenser_app.modules.selection import OUTCOME_PRESET
from condenser_app.reporting import DISCLAIMER, STATUS_TEXT, to_excel, to_pdf
from condenser_app.schemas import Category, CheckState, DesignResult, Status, iter_inputs
from condenser_app.service import available_modules, get_module, run_design

st.set_page_config(page_title="Condenser Design", layout="wide")

HOME = "Home"
ABOUT = "Sources and limits"
NOT_SET = "(not set)"


# ------------------------------------------------------------------ form state
def _key(module_key: str, path: str) -> str:
    return f"{module_key}:{path}"


def _is_number(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _to_shown(value, meta):
    """Internal value -> what the widget holds. Number widgets hold floats, in the display unit."""
    if not _is_number(value):
        return value
    return float(convert(value, meta["unit"], meta["display_unit"]) if meta["unit"] != meta["display_unit"] else value)


def _load_example(module_key: str) -> None:
    _reset(module_key)
    for path, owner, name, meta in iter_inputs(get_module(module_key).example_inputs()):
        value = getattr(owner, name)
        if value is not None:
            st.session_state[_key(module_key, path)] = _to_shown(value, meta)
    st.session_state.pop(f"result:{module_key}", None)


def _reset(module_key: str) -> None:
    for key in [k for k in st.session_state if k.startswith(f"{module_key}:") or k == f"result:{module_key}"]:
        del st.session_state[key]


def _read_form(module_key: str):
    """Build the input dataclass from what is on the form, converting display units to internal units."""
    inputs = get_module(module_key).new_inputs()
    for path, owner, name, meta in iter_inputs(inputs):
        widget_key = _key(module_key, path)
        if widget_key not in st.session_state:
            continue
        value = st.session_state[widget_key]
        if value == NOT_SET or value == "":
            value = None
        elif _is_number(value):
            if meta["unit"] != meta["display_unit"]:
                value = convert(value, meta["display_unit"], meta["unit"])
            if "int" in str(type(owner).__dataclass_fields__[name].type):
                value = int(round(value))
        setattr(owner, name, value)
    return inputs


def _widget(module_key: str, path: str, meta: dict, default, is_text: bool, column) -> None:
    key = _key(module_key, path)
    if key not in st.session_state and default is not None:
        st.session_state[key] = _to_shown(default, meta)  # seed, so the field's declared default is not lost
    unit = meta["display_unit"]
    label = f"{meta['label']} [{unit}]" if unit else meta["label"]
    if meta["required"]:
        label += " *"
    help_text = " ".join(x for x in (meta["help"], f"Source: {meta['source']}" if meta["source"] else "") if x) or None
    if meta["choices"]:
        column.selectbox(label, (NOT_SET,) + tuple(meta["choices"]), key=key, help=help_text)
    elif is_text:
        column.text_input(label, key=key, help=help_text)
    else:
        column.number_input(label, value=None, format="%g", step=None, key=key, help=help_text)


def _render_group(module_key: str, obj, prefix: str, shown: set[str] | None) -> None:
    """Draw one dataclass: its own fields in two columns, then each nested group under a heading."""
    leaves = [f for f in fields(obj) if "label" in f.metadata and (shown is None or f"{prefix}{f.name}" in shown)]
    columns = st.columns(2)
    for i, f in enumerate(leaves):
        _widget(module_key, f"{prefix}{f.name}", dict(f.metadata), getattr(obj, f.name), "str" in str(f.type), columns[i % 2])
    for f in fields(obj):
        value = getattr(obj, f.name)
        if is_dataclass(value) and not isinstance(value, type):
            sub = f"{prefix}{f.name}."
            if shown is None or any(p.startswith(sub) for p in shown):
                st.markdown(f"**{f.metadata.get('group', f.name)}**")
                _render_group(module_key, value, sub, shown)


# ------------------------------------------------------------------ results
def _status_banner(result: DesignResult) -> None:
    text = STATUS_TEXT[result.status]
    if result.status is Status.OK:
        st.success(text)
    elif result.status is Status.OK_WITH_WARNINGS:
        st.warning(text)
    else:
        st.error(text)
    if result.message:
        st.write(result.message)


def _fmt(value) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.5g}"
    return str(value)


def _show_results(module_key: str, result: DesignResult, inputs) -> None:
    st.divider()
    st.subheader("Results")
    _status_banner(result)
    for issue in result.input_issues:
        st.error(issue.message)
    if result.status is Status.INVALID_INPUT:
        return

    for w in result.warnings:
        st.warning(w)

    tab_names = ["Summary", "Calculation steps", f"Checks ({len(result.checks)})", "Convergence", "Notes and sources", "Export"]
    summary, steps, checks, convergence, notes, export = st.tabs(tab_names)

    with summary:
        for category in Category:
            rows = [q for q in result.outputs if q.category is category]
            if rows:
                st.markdown(f"**{category.value.capitalize()}**")
                st.dataframe(
                    pd.DataFrame({"Result": [q.label for q in rows], "Value": [_fmt(q.value) for q in rows], "Unit": [q.unit for q in rows], "Source": [q.source for q in rows]}),
                    hide_index=True, width="stretch",
                )
        if not result.outputs:
            st.info("No results were produced.")
        _design_link(result)

    with steps:
        st.caption("Each step shows the formula, then the numbers, then the result.")
        for category in Category:
            rows = [s for s in result.steps if s.category is category]
            if not rows:
                continue
            st.markdown(f"**{category.value.capitalize()}**")
            for s in rows:
                st.markdown(f"{s.label}" + (f"  \n<small>{s.source}</small>" if s.source else ""), unsafe_allow_html=True)
                st.code(s.chain(), language=None, wrap_lines=True)

    with checks:
        marks = {CheckState.PASS: "PASS", CheckState.FAIL: "FAIL", CheckState.NOT_ASSESSED: "NOT ASSESSED"}
        if result.checks:
            st.dataframe(
                pd.DataFrame([{"State": marks[c.state], "Check": c.name, "Detail": c.detail, "Category": c.category.value, "Source": c.source} for c in result.checks]),
                hide_index=True, width="stretch",
            )
        else:
            st.info("This module recorded no engineering checks.")

    with convergence:
        if not result.convergence:
            st.info("This calculation did not need an iteration.")
        for c in result.convergence:
            st.markdown(f"**{c.name}**: {'converged' if c.converged else 'NOT converged'} in {c.iterations} passes; tolerance {c.tolerance:g}, final residual {c.residual:.3g}")
            if c.rows:
                st.dataframe(pd.DataFrame(c.rows, columns=c.columns or None), hide_index=True, width="stretch")
            elif len(c.history) > 1:
                st.line_chart(pd.DataFrame({"value": c.history}))

    with notes:
        for title, items in (("Assumptions", result.assumptions), ("Source gaps", result.source_gaps), ("References", result.references)):
            if items:
                st.markdown(f"**{title}**")
                for item in items:
                    st.markdown(f"- {item}")
        st.caption(DISCLAIMER)

    with export:
        left, right = st.columns(2)
        left.download_button("Download Excel report", to_excel(result, inputs), file_name=f"{module_key}_report.xlsx",
                             mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        right.download_button("Download PDF report", to_pdf(result, inputs), file_name=f"{module_key}_report.pdf", mime="application/pdf")


def _go_to(page: str) -> None:
    st.session_state["page"] = page


def _go_to_design(module_key: str, outcome: str) -> None:
    """Open the design module for the chart's outcome, with the matching type or variant already chosen."""
    if outcome in OUTCOME_PRESET:
        field, value = OUTCOME_PRESET[outcome]
        st.session_state[_key(module_key, field)] = value
    _go_to(get_module(module_key).title)


def _design_link(result: DesignResult) -> None:
    """After selection: one click to the design module for the recommended type."""
    if result.module != "selection":
        return
    target = next((q.value for q in result.outputs if q.key == "design_module"), "none")
    modules = available_modules()
    if target in modules:
        outcome = next(q.value for q in result.outputs if q.key == "recommended_type")
        st.button(f"Design this condenser: {modules[target]}", type="primary", on_click=_go_to_design, args=(target, outcome))
    elif target != "none":
        st.info(f"The design module '{target}' is not available yet.")


# ------------------------------------------------------------------ pages
def page_home() -> None:
    st.title("Condenser Design")
    st.write(
        "Enter the operating conditions, let the selection chart recommend a condenser type or choose one yourself, "
        "and get a preliminary design you can trace line by line."
    )
    st.markdown("**Input data → Validation → Condenser selection → Design calculations → Engineering checks → Results → Report export**")
    modules = available_modules()
    for key, title in modules.items():
        st.button(title, key=f"home:{key}", on_click=_go_to, args=(title,), help=get_module(key).description)
    st.warning(
        "The project's supplied algorithms, data-entry sheets, numerical examples and standards were not available when this was built. "
        "The modules follow the ME6220 lectures, Hewitt, Shires & Bott and the Heat Exchanger Design Handbook. See 'Sources and limits'."
    )
    st.caption(DISCLAIMER)


def page_module(module_key: str) -> None:
    module = get_module(module_key)
    st.title(module.title)
    st.write(module.description)

    left, right, _ = st.columns([2, 1, 4])
    left.button("Load example", on_click=_load_example, args=(module_key,), key=f"example:{module_key}", help="Fill the form with this module's cited example.")
    right.button("Reset", on_click=_reset, args=(module_key,), key=f"reset:{module_key}")

    # Build what is on the form so far first, so the module can decide which fields apply.
    shown = module.relevant_fields(_read_form(module_key))
    _render_group(module_key, module.new_inputs(), "", shown)
    st.caption("Fields marked with an asterisk are required.")

    if st.button("Calculate", type="primary", key=f"calculate:{module_key}"):
        inputs = _read_form(module_key)
        st.session_state[f"result:{module_key}"] = (run_design(module_key, inputs), copy.deepcopy(inputs))
    stored = st.session_state.get(f"result:{module_key}")
    if stored:
        if stored[1] != _read_form(module_key):
            st.info("The inputs have changed since this result was calculated. Press Calculate to update it.")
        _show_results(module_key, *stored)


def page_about() -> None:
    st.title("Sources and limits")
    readme = Path(__file__).parent / "docs" / "gaps_and_assumptions.md"
    if readme.exists():
        st.markdown(readme.read_text(encoding="utf-8"))
    else:
        st.info("docs/gaps_and_assumptions.md is not present.")


# ------------------------------------------------------------------ navigation
modules = available_modules()
title_to_key = {title: key for key, title in modules.items()}
pages = [HOME, *modules.values(), ABOUT]
if st.session_state.get("page") not in pages:
    st.session_state["page"] = HOME
page = st.sidebar.radio("Page", pages, key="page")

if page == HOME:
    page_home()
elif page == ABOUT:
    page_about()
else:
    page_module(title_to_key[page])
