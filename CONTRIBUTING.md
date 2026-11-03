# Contributing

Three people, one repository, one application. `main` always holds the integrated, working app.

## Who owns what

This split keeps people in separate files.

| Member | Owns | Files |
|---|---|---|
| A (integration lead) | Selection, shared schemas and library, service, app, reports | `condenser_app/schemas.py`, `common/`, `service.py`, `reporting.py`, `modules/selection.py`, `app.py` |
| B | Direct-contact and plate modules | `modules/direct_contact.py`, `modules/plate.py`, their tests and docs |
| C | Vertical shell-and-tube and air-cooled modules | `modules/vertical_shell_tube.py`, `modules/air_cooled.py`, their tests and docs |

## Rules

1. Never push to `main`. Work on a feature branch and open a pull request.
2. Stay in your own files where you can. Modules register themselves, so adding one never touches a shared list.
3. Every module follows `modules/base.py` and returns a `DesignResult`.
4. **Do not change a shared signature alone.** `schemas.py`, `modules/base.py`, `service.py` and existing functions in `common/` are used by everyone. Agree the change with the team first. Adding a new function is fine.
5. Commit messages say what changed and why.
6. All tests must pass before a merge: `python -m pytest -q`.
7. Resolve merge conflicts carefully, and never overwrite someone else's work without asking them.
8. When engineering logic changes, update its test and `docs/algorithms/<module>.md` in the same pull request.

## Engineering rules

- Every equation, constant and table value names its source (book and section or equation, or lecture number). Nothing from memory.
- If a source has no procedure for a step, record it in `result.scope_limits` or require the value as an input. Do not fill the gap with a guess.
- Benchmark answers belong in tests, never in engineering code. Do not tune a calculation to match a book; find out why it differs and document it.
- Report tests as passing only after actually running them.
- Do not claim code compliance. A check with no usable procedure in the books is recorded as "not assessed".
- Do not commit textbooks, standards or slides. `.gitignore` blocks common formats; keep them outside the repository or in `references/`.

## Showing the working

Results are meant to be followed, so:

- Every `result.step` is a full chain: `formula` holds the symbols and `substitution` the numbers. Never name a correlation without writing it.
- Record each validity check where the correlation is used, with its criterion: `result.check("Dittus-Boelter range", ok, "Re = 18500 > 1e4; 0.7 <= Pr = 3.58 <= 160; L/D = 200 > 60")`.
- When choosing between alternatives, say why, in the step label or `result.design_basis`.
- Run iterations through `common.iteration`, and state why it iterates, the first guess, the update rule and the stopping rule.
- Read equations from page images, because OCR garbles them. Check a book's arithmetic before trusting a discrepancy either way.
- After changing inputs, run `python tools/generate_variable_map.py`.

## Git commands

```bash
git clone https://github.com/Raghav589/condenser-design-app.git
cd condenser-design-app
python -m venv .venv
```

Activate the environment and install, on Windows:

```bash
.venv\Scripts\activate
pip install -r requirements.txt
```

On macOS or Linux:

```bash
source .venv/bin/activate
pip install -r requirements.txt
```

Start a piece of work from an up-to-date `main`:

```bash
git checkout main
git pull origin main
git checkout -b feature/plate-pressure-drop
```

Test, commit and push:

```bash
python -m pytest -q
git add condenser_app/modules/plate.py tests/test_plate.py docs/algorithms/plate.md
git commit -m "Add port pressure drop to plate condenser"
git push -u origin feature/plate-pressure-drop
```

Then open a pull request on GitHub into `main`, ask a teammate to review it, and merge once tests pass.

If `main` moved while you worked:

```bash
git fetch origin
git merge origin/main
python -m pytest -q
git push
```

Read the `git status` output before committing and confirm no PDFs or slides are listed.
