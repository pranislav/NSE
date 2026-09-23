# MMS report TODO

## Context / resume reference

- Chat ID needed to resume this context window: `codex resume 01a0b55d-94b9-78f2-bc7f-8d13b4b80af5`
- Relevant template: [`mms_verification_report_template.tex`](../mms_verification_report_template.tex).

## Proposed report scope

- Write a dissertation-ready MMS verification report/chapter of roughly 12--16
  pages of main text, with full convergence tables and generated source terms
  in appendices.
- Keep the scientific sequence: verification scope -> mathematical model ->
  manufactured problem -> numerical procedure -> flow results -> temperature
  results -> spatial diagnostics -> limitations -> conclusions.
- Treat commit history as provenance, not as the report's narrative structure.

## Data status

### Flow verification

- [x] Results exist for all core flow combinations:
  - [x] `Re = 100`, degrees 1, 2, and 3.
  - [x] `Re = 7500`, degrees 1, 2, and 3.
- [x] Each of these has raw solution output, global error tables, fitted
  convergence rates, and individual norm-wise plots.
- [x] Available flow norms: velocity `L2`, pressure `L2`, and velocity `H1`.

### Coupled temperature verification

- [x] `Re = 100`, degrees 1, 2, and 3 have final coupled result sets in
  `solns/mms_re100_deg{1,2,3}_finer_integral_error_temperature_tighter/`.
- [x] `Re = 7500`, degrees 1, 2, and 3: coupled temperature MMS results are
  missing. Generate these if the report claims coupled flow--temperature
  verification at both Reynolds numbers.
- [x] Degree-1 temperature convergence at `Re = 100` is stable:
  fitted `L2(T)` rate = `2.9989 +/- 0.0003`.
- [ ] Degree-2 and degree-3 temperature data reach a fine-mesh numerical
  floor. Do not present the current full-range fitted rates as asymptotic
  discretisation rates without further analysis.

## Analysis TODO

- [ ] Select one authoritative result set per reported parameter combination;
  mark older `mms_try*`, legacy, and non-tightened-temperature directories as
  exploratory/archival.
- [ ] Produce a final expected-versus-observed flow-rate table for `Re=100`
  and `Re=7500`, including fitted-slope standard errors and stated fit ranges.
- [ ] Reanalyse temperature convergence for degrees 2 and 3:
  - determine the first floor-dominated refinement level;
  - fit only defensible pre-floor points, or rerun with a demonstrated
    tolerance criterion;
  - report raw errors and local rates regardless of the chosen fit.
- [ ] Make a final composite temperature-convergence figure from the existing
  `error-global-q*_L2_temperature.png` outputs.
- [ ] Decide whether to add temperature to the per-cell integral-error plot.
  The solver writes `temperature_L2_cell_error`, but
  `scripts/visualize_mms_integral_errors.py` currently visualises only
  velocity `L2`, velocity `H1`, and pressure `L2` cell errors.
- [x] If required by the final claim, run coupled MMS at `Re=7500` for degrees
  1--3 and perform the same convergence analysis.

## Suggested figures and tables

- [ ] Exact fields and source terms: `mms_visualization.png`.
- [ ] One compact flow convergence-rate comparison table.
- [ ] Flow log--log convergence plots, grouped cleanly by degree or Reynolds
  number; use consistent axes and labels.
- [ ] Temperature `L2` convergence figure that visibly identifies excluded or
  floor-dominated levels.
- [ ] A cropped/rebuilt version of `mms_integral_errors.png`, rather than the
  dense full 3x4 overview, for the main text.
- [ ] Put complete per-cycle generated tables in an appendix.

## Report-writing safeguards

- [ ] Distinguish code verification from validation against experiments.
- [ ] State that the current MMS is smooth, steady, two-dimensional, and on a
  fluid unit square.
- [ ] Do not claim verification of fluid--solid conjugate interface coupling;
  that would require a separate manufactured fluid--solid problem.
- [ ] Explain pressure mean-value normalisation before pressure-error results.
- [ ] Keep exact fields, C++ functions, and LaTex equations synchronised via
  `scripts/expressions_mms.py` and `scripts/generate_mms.py`.

## Key sources

- Exact MMS definition: `scripts/expressions_mms.py`
- Generated solver functions: `src/mms_generated.h`
- deal.II MMS wrappers: `src/mms.h`
- Error computation/output: `src/conjugate_heat_transfer_solver.cc`
- Rate fitting: `scripts/convergence_rate_analysis.py`
- Report template: `mms_verification_report_template.tex`
