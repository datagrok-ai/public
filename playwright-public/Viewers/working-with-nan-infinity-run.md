# Working with NaN and Infinity — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | SPGI_v2_infinity.csv — open dataset, verify infinity values present | 25s | PASS | PASSED | 3624 rows, 88 cols; `Chemical Space X` and `Chemical Space Y` confirmed to have Infinity values |
| 2 | Add Scatter plot (Chemical Space X / Chemical Space Y) | 18s | PASS | PASSED | Viewer renders; clusters visible; Infinity rows silently excluded from plot |
| 3 | Add Histogram on Chemical Space X (infinity column) | 12s | PASS | PASSED | Distribution renders finite range only; infinity excluded from binning |
| 4 | Save layout for SPGI dataset | 8s | PASS | PASSED | Layout saved via `grok.dapi.layouts.save` |
| 5 | Run console script — demog, NaN/Infinity set, scatter plot x=height y=weight size=age color=race | 20s | PASS | PASSED | All 4 column assignments confirmed; NaN stored as internal missing encoding; Infinity stored as null; regression=true, markerType=square |
| 6 | Scatter plot renders correctly — NaN/Infinity rows excluded | 15s | PASS | PASSED | Regression lines per race category shown; canvas visible without crash |
| 7 | Add Histogram on height column (NaN column) | 12s | PASS | PASSED | Histogram renders clean distribution; NaN row excluded from bins |
| 8 | Add Histogram on weight column (Infinity column) | 12s | PASS | PASSED | Histogram renders clean distribution; Infinity (stored as null) excluded from bins |
| 9 | Save layout for demog dataset | 6s | PASS | PASSED | Layout saved; both layouts cleaned up after run |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~4 min |
| grok-browser execution (scenario steps) | ~2 min |
| Execute via grok-browser (total) | ~6 min |
| Spec file generation | ~3 min |
| Spec script execution | 1m 24s |
| **Total scenario run (with model)** | ~11 min |

## Summary

All 9 spec steps PASSED in 1m 24s. NaN and Infinity values in numeric columns are handled gracefully across Scatter Plot and Histogram viewers — affected rows are silently excluded from rendering without errors or visual artifacts. The total scenario run took approximately 11 minutes including model thinking and spec execution.

## Retrospective

### What worked well
- Scatter plot excluded NaN/Infinity rows cleanly; regression lines computed on finite data only
- Histogram binning naturally excluded infinite/missing values with no special configuration
- `grok.data.demo.demog()` + manual `col.set(0, NaN/Infinity)` reproduced test conditions reliably
- All 4 scatter plot column assignments (x=height, y=weight, size=age, color=race) verified

### What did not work
- `t.col('height').set(0, NaN)` — NaN is not stored as IEEE NaN; getter returns a tiny float (~2.7e-34), Datagrok's internal missing-value encoding for `double` columns. Treated as missing by all viewers.
- `t.col('weight').set(0, Infinity)` — Infinity stored as `null` (missing). Datagrok does not store actual Infinity in double columns.
- Step 6 initially failed (`showRegressionLine=false`) due to timing: `setOptions()` is async; fixed by adding 2s wait.
- Step 7 initially failed (strict mode violation) because the SPGI tab's Histogram element stays in the DOM when switching to the demog tab; fixed by counting viewers in `grok.shell.tv` before/after click.

### Suggestions for the platform
- Document (or warn at `set()` call time) that `NaN` and `Infinity` are silently coerced to missing values — surprising behavior for developers testing edge cases.
- A `console.warn` when `set()` receives `NaN` or `Infinity` would improve debuggability.

### Suggestions for the scenario
- "Add viewers+layouts" is ambiguous — specifying "Scatter Plot + Histogram on height and weight" would make the test deterministic.
- Each step could state the expected outcome: "verify that rows with NaN/Infinity are excluded from the plot and do not cause errors."
- The scenario could note that `Infinity` in a column becomes a missing value — testers may be surprised to find no Infinity markers in the viewer.
