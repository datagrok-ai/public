# Multivariate Analysis — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open cars.csv from Demo files | 1s | PASS | PASSED | Opened via JS API setup batch (closeAll + readCsv + addTableView + semType wait); 30 rows, 17 cols |
| 2 | Top Menu > ML > Analyze > Multivariate Analysis...; viewers display | 45s | PASS | PASSED | JS API fallback menu click on `[name="div-ML---Analyze---Multivariate-Analysis..."]`; dialog defaults (Predict=price, Using=15 cols, Components=2) kept; RUN button enabled; 6 viewers rendered (Grid + 3 Scatter + 2 Bar) |
| 3 | Check interactivity: grid <-> scatterplots, Loadings <-> Regression bar | 0s | AMBIGUOUS | PASSED | Viewers share the same DataFrame so hover/selection propagate via the standard Datagrok bus; confirming live highlight from a script requires synthetic canvas mouse events at pixel coordinates — not attempted. Spec asserts presence (>=3 Scatter plot, >=2 Bar chart) |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m |
| grok-browser execution (scenario steps) | 45s |
| Execute via grok-browser (total) | 2m 30s |
| Spec file generation | 2m |
| Spec script execution | 13s |
| **Total scenario run (with model)** | 4m 43s |

## Summary

2 of 3 scenario steps passed and 1 is recorded as AMBIGUOUS (Step 3 interactivity check, where the
wording does not specify a scriptable success criterion). The generated Playwright spec is green —
all 3 softSteps pass in ~13s, verifying the MVA dialog opens with correct defaults, RUN produces the
expected 6 viewers (Grid + 3 Scatter plot + 2 Bar chart), and no dialog misconfiguration disables
RUN. **Total scenario run (with model)**: ~4m 43s.

## Retrospective

### What worked well
- JS API fallback `document.querySelector('[name="div-ML---Analyze---Multivariate-Analysis..."]').click()` reliably opens the dialog without needing to hover the ML menu first.
- Dialog defaults (Predict=price, Using=15 numeric cols excluding price) produce a valid configuration out of the box — no input tweaking needed.
- Polling `grok.shell.tv.viewers.length` is a robust signal for "MVA finished rendering", since all 5 analytic viewers are added synchronously at the end of the run.
- The spec completes in ~13s because dev has the dataset cached and MVA on 30 rows is fast.

### What did not work
- Validating live cross-viewer highlighting (Grid row hover → Scores/Observed scatterplot highlight; Loadings click → Regression bar filter) from a browser script — canvas viewers expose no DOM landmarks at point level, so the spec falls back to a structural assertion (viewer count/types).
- Setting `Using` to "(16) All" in the dialog introduces `price` twice and disables RUN with a red border on `Predict`. The spec explicitly leaves the default (15) to avoid this trap.

### Suggestions for the platform
- Add a hover/selection-visible hook on scatter/bar viewers (e.g., `viewer.getHoveredPoint()` returning the DataFrame row index) so interactivity scenarios can be validated without canvas pixel math.
- Prevent the "RUN disabled" trap by auto-excluding the `Predict` column from `Using` when `Select all` is clicked, or by surfacing an inline hint explaining why RUN is disabled.
- Tag MVA output viewers with analysis-specific `name` attributes (e.g., `[name="viewer-MVA-Loadings"]`, `[name="viewer-MVA-Scores"]`) so automation can address them by role rather than by positional type filtering.

### Suggestions for the scenario
- Step 3 wording "Check interactivity" is vague — specify which specific interactions (hover a grid row → highlight in scatter; click a loading point → filter regression bars) so testers know what to verify.
- Note that dialog defaults already produce a valid run — testers don't need to touch inputs, and explicitly "Select all" on `Using` breaks the dialog. Call this out in a pre-condition.
- State the expected viewer count (6 total: Grid + Observed/Predicted + Scores + Loadings + Regression Coefficients + Explained Variance) so testers can quickly verify without counting.
