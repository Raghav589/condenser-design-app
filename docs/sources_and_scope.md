# Sources and scope

## Sources

The modules follow the ME6220 lectures, Hewitt, Shires & Bott (HSB) and the Heat Exchanger Design Handbook (HEDH). `docs/reference_inventory.md` lists them. Each result names the source of every step.

This is a preliminary design tool for coursework. Mechanical sizing is a screening calculation from HEDH Part 4. It does not establish code compliance, and checks that would need the code itself show as **not assessed**.

## What each module rests on

| Module | Basis, and how far it has been tested |
|---|---|
| Selection | HSB Fig. 17.1. The chart is broad guidance, and the book itself says not to follow it slavishly |
| Direct-contact | Spray, sheet, packed-bed and bubble paths reproduce HSB's examples. Baffle-column stage count is built from sourced parts, because no stage procedure is published. Pool condensers: balances only, because no sizing procedure is published |
| Plate | Book correlations built in (HSB Ex. 8.1-8.2, HEDH intermating), or entered coefficients. HEDH's approximate-design curves are digitised and reported as a cross-check. The condensing path uses tube correlations on plates, as HSB reports from Kumar. No source works a complete condensing plate example |
| Vertical shell-and-tube | Bell-Delaware coolant side, zoned condensing side, standard shell and tube-count tables, HEDH Part 4 mechanical and vibration screening. No source works a complete vertical condenser example, so the whole path is tested by hand values, closure and trends |
| Air-cooled | HSB Ch. 9 and 19, HEDH 3.8. The regime rule of HEDH 3.4.6 B is reported at every station. Tube-side methods differ by up to about 50% in length, and the spread is reported. No source works a complete vertical downflow example |

## Design basis

- **Vertical in-tube coefficient:** HSB Eq. 16.56 with Kosky-Staub shear. Boyko-Kruzhilin is reported as a cross-check.
- **Suction correction:** off by default, because HSB's text is ambivalent. It is an input, and the other setting is always reported.
- **Friedel in downflow:** the printed form is used, as HSB does in Ex. 16.2, because neither book prints the downflow form. HSB prints 3.24 and (1-x)^0.24; HEDH 2.3.2 prints 3.23 and (1-x)^0.224. The code uses HSB's, because its examples do.
- **Spray height:** sized with gravity included (2.15 m in HSB Ex. 22.1), because the no-gravity height (1.62 m) is too short for that drop. The book's value is also reported. The maximum drop size is the design drop, as HEDH 2.6.8 advises for preliminary design.
- **Plate Eq. 20.6:** evaluated at a common mass flux. HSB's wording is "local conditions", but with each phase's actual flow the bracket gives absurd values and the wrong limit at x = 0.
- **Plate F factor:** a closed form derived from Lecture 15's channel equations. It matches Lecture 15's chart within 0.003.
- **Reflux flooding:** designed to 75% of the lower of HSB's and HEDH's limits, which differ by 2.7 times for pentane.
- **Shell-to-baffle clearance:** HEDH Eq. (19), the line it captions on Fig. 13.
- **Vibration screen:** the lowest beta band of HEDH Fig. 4.6.4(4), because the scan does not separate the upper bands. C_n is the lowest curve of Fig. 4.6.2(2).
- **Non-condensables in the vertical module:** HEDH 3.4.9B's factor 1 - 0.18C, below 5 wt%. HEDH has no simple method between 5 and 80 wt%, so higher values are refused.

## Chart data read from the books

Digitised values are stored in the module that uses them, each labelled with its figure and read uncertainty.

| Chart | Used for | Read uncertainty |
|---|---|---|
| HEDH 3.7.10 Figs 1-3 | Plate approximate design | About 1.5% |
| HEDH 3.3.5 Fig. 4 | Baffle cut against spacing | About 1% B_c |
| HEDH 3.3.5 Figs 12, 14 | Tube-to-baffle and bundle-to-shell clearance | 0.01 and 0.5 mm |
| HEDH 4.6 Figs 2-4 and the Strouhal chart | Vibration screen | See `docs/algorithms/vertical_shell_tube.md` |

## Scope

| Area | Status |
|---|---|
| Pressure boundary | Shell and tube thickness, TEMA minimum shell and tubesheet sizing from HEDH Part 4. Fixed tubesheets use HEDH's rough 0.1 D_s, because the fixed-tubesheet formula is not printed |
| Flanges, nozzle reinforcement, expansion joints, plate stress, frames | Not assessed. No usable procedure in the books |
| Vent sizing and fogging in tubular units | Not calculated. HEDH says there are very few rules. Its guidance is reported in the checks |
| Mixtures | Vertical module only, with an entered condensation curve |
| Plate condensers: several vapour passes, subcooling, desuperheating, mixtures | No procedure in HSB or HEDH |
| Vertical: several tube passes, cocurrent coolant, baffled shell-side condensation, flooded subcooling | No usable procedure in HSB or HEDH |
| Air-cooled: superheat, reflux section, multiple tube passes, in-line high fins | Outside this module |
| Horizontal condensers and plate-fin | The selection chart can recommend them; no design module exists |
| Standard shell | The next standard shell and its tube count are reported. The design is not re-run on it |
| Cost estimation | Outside this application |

Not in any of the books: the fraction of terminal velocity that limits vapour velocity, drop-distribution parameters (HEDH and HSB say nozzle makers provide them), packing coefficient tables, the range of HSB Eq. 22.24, tail-pipe friction and seal allowances, and a plate stress method.

## Errors found in the sources

These were found by comparing a book's numbers with its own equations. Each is listed with its evidence in `docs/algorithms/`. A sample was re-checked independently against page images, including HSB Ex. 22.4's c_p and Ex. 16.1's latent heat.

| Where | Issue |
|---|---|
| HSB Ex. 22.4 | c_p printed as 4709; 4209 is correct, so Ja′ is 24.1, not 27.0 |
| HSB Ex. 16.1 | Uses 776.9 kJ/kg as n-hexane's latent heat; its Table II.1 gives 311.7 |
| HSB Ex. 16.4 | Says isopentane, but the data are n-hexane |
| HSB Eq. 16.24, 22.6, 22.18, 22.20, 22.42, 10.27, 10.45, 16.54, 20.7 | Sign or printing errors |
| HSB Eq. 19.3 | Prints +10 log 2πR²; HEDH's minus is the one its example uses |
| HSB Ex. 3.3 | Prints 39.355; its column adds to 39.855 |
| HSB Ex. 19.2 | Air balance does not close |
| HSB Eq. 2.229 | Diagonal-plane test is backwards |
| HEDH 4.3.6 | Prints 7.12 mm where its numbers give 7.14 mm, and 18 610 kPa where they give 18 581 kPa |
| HEDH 2.6.8 Eq. (81), Fig. 25 | ν^(1/2) fails a dimensional check; the axis prints Pe where Eq. (87) has Pe^(1/2) |
| Lecture 6 | Coupled-exchanger bracket; three-point method has a stray "3" |
| Lecture 8 | "1−P1g" should be 1−P1j; connecting list repeats one line and omits another |
| Lecture 9 | 24.40 and 25.04 mm should be 25.40 mm (HEDH 4.2.5-1 settles this); "1542 mm" should be 1524 mm |
| Lecture 23 | Baffle column captioned as a bubble column |
| Lecture 24 | Silver-Bell-Ghaly dated "1872" |

## More detail

- `docs/reference_inventory.md`: the sources
- `docs/benchmarks.md`: what is and isn't verified
- `docs/standards_map.md`: which standard governs what
- `docs/algorithms/`: each module's steps, equations, variables and scope
- `docs/variable_map.md`: every input field
