# Condenser Design

One application for preliminary condenser design, built for ME6220 Heat Exchanger Design at IIT Madras. Enter the operating conditions, let the selection chart recommend a condenser type or choose one yourself, and get a design you can trace line by line, with its checks and warnings.

**Workflow:** input data → validation → condenser selection → design calculations → engineering checks → results → report export.

| Module | What it does |
|---|---|
| Selection and fluid allocation | Walks HSB's selection chart, recommends a type, and advises which fluid goes where |
| Direct-contact | Coolant flow, mixed outlet state and height for spray, jet, sheet, fan, baffle column, packed bed and bubble types; balances only for pool |
| Plate | Duty, plate count, coefficients, U and both pressure drops |
| Vertical shell-and-tube | In-tube downflow, unbaffled shell-side and reflux variants; zones, tube and shell geometry, Bell-Delaware coolant side, pressure drops |
| Vertical downflow air-cooled | Air flow, finned-tube bundle, tube length, pressure drops and fan power |

## Basis and limits

The modules follow the ME6220 lectures, Hewitt, Shires & Bott and the Heat Exchanger Design Handbook, including tables and charts digitised from those books.

Not every path has a worked example behind it: no source works a complete vertical shell-and-tube, condensing plate or vertical downflow air-cooled example. See [`docs/benchmarks.md`](docs/benchmarks.md) and [`docs/sources_and_scope.md`](docs/sources_and_scope.md).

This is a preliminary tool. Mechanical sizing is a screening calculation from HEDH Part 4 and does not establish code compliance.

## Setup

Tested on Python 3.14 on Windows. On Windows, clone into a short path such as `C:\projects\condenser`: a deeply nested folder can exceed the 260-character path limit and make `pip install` fail partway through Streamlit.

```bash
git clone https://github.com/Raghav589/condenser-design-app.git
cd condenser-design-app
python -m venv .venv
```

Activate the environment. On Windows:

```bash
.venv\Scripts\activate
```

On macOS or Linux:

```bash
source .venv/bin/activate
```

Then install the dependencies and run:

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Using it

1. Open a module from the sidebar. To see it working, press **Load example**.
2. Fill in the form, then press **Calculate**.
3. Read the result tabs:
   - **Summary:** results by category
   - **Calculation steps:** each step's formula, then the numbers, then the result
   - **Checks:** pass, fail or not assessed, with the criterion
   - **Convergence:** every iteration, pass by pass
   - **Notes and sources:** design basis, scope limits, references
   - **Export:** Excel or PDF

From the selection page, **Design this condenser** opens the right module with the type already chosen.

The examples are real textbook cases where one exists. Some show failed checks on purpose, because the book's case falls outside its own stated ranges. The air-cooled example is one of these.

## Tests

```bash
python -m pytest -q
```

## Documentation

| Document | Contents |
|---|---|
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | Team workflow and Git commands |
| [`docs/architecture.md`](docs/architecture.md) | Structure, interfaces and data model |
| [`docs/reference_inventory.md`](docs/reference_inventory.md) | Sources used |
| [`docs/algorithms/`](docs/algorithms/) | Each module's steps, equations, variables and scope |
| [`docs/benchmarks.md`](docs/benchmarks.md) | What is and isn't verified |
| [`docs/variable_map.md`](docs/variable_map.md) | Every input field (generated) |
| [`docs/standards_map.md`](docs/standards_map.md) | Which standard governs what |
| [`docs/sources_and_scope.md`](docs/sources_and_scope.md) | Design basis, scope, errors found in the sources |

Textbooks, standards and slides must not be committed. `.gitignore` blocks the usual formats.
