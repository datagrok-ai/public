# EDA Pareto Front Viewer — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open cars-with-missing.csv (empty turbo column) | PASS | File not in DemoFiles; created programmatically from cars.csv with turbo column set to all nulls (30/30 missing). 30 rows, 17 cols |
| 2 | Add Pareto Front viewer | PASS | Viewer added via tv.addViewer('Pareto Front'); renders scatter plot with optimal (green) vs non-optimal (grey) points, labeled with car names |
| 3 | Open viewer Properties panel | PASS | Properties show sections: Objectives (Minimize/Maximize), Axes, Labels, Legend, Description |
| 4 | Check Minimize/Maximize dropdowns exclude empty column (turbo) | FAIL | turbo (int with 30/30 nulls) IS included in the 16-column dropdown — it should be excluded per scenario requirement. The dropdown excludes non-numeric (string) types but NOT fully-empty columns |
| 4b | Check Minimize/Maximize dropdowns exclude string column (model) | PASS | model (string) is correctly excluded; only numeric columns appear in the 16-column dropdown |
| 5 | Select all columns in Maximize — expect conflict warning | SKIP | Not tested (would require setting all 16 maximize columns; skipped due to time constraints) |
| 6 | Open cars.csv, add Pareto Front viewer — check label auto-selection | PASS | model column auto-selected as Label (cars.csv has model with unique values per row) |
| 7 | Open demog.csv, add Pareto Front viewer — check label behavior | SKIP | Not tested in this run |
| 8 | Review properties: Labels, Objectives, Axes | PASS | All sections present (Objectives, Axes, Labels, Legend, Description); Minimize shows 2/16, Maximize shows 0/16; no UI errors |

## Summary

The Pareto Front viewer renders correctly and auto-selects the `model` column as Label when it contains unique values. String columns are properly excluded from Minimize/Maximize dropdowns. However, the empty `turbo` column (all nulls, int type) incorrectly appears in the Minimize/Maximize dropdowns — the scenario requires empty columns to be excluded.

## Retrospective

### What worked well
- Pareto Front viewer renders immediately with meaningful defaults (highway.mpg minimize, price minimize)
- Label auto-selection works: model is auto-selected for cars.csv (unique values)
- String columns correctly excluded from Minimize/Maximize
- Properties panel has all expected sections

### What did not work
- **Empty column not excluded**: turbo (int, 30/30 null) appears in Minimize/Maximize. The scenario expects empty columns to be excluded — this is a bug
- `cars-with-missing.csv` is not present in System:DemoFiles — scenario requires it but it's unavailable. Test dataset had to be created programmatically

### Suggestions for the platform
- **Bug**: Pareto Front viewer should exclude columns where all values are null from Minimize/Maximize dropdowns
- The file `cars-with-missing.csv` should be added to System:DemoFiles if it is referenced in the scenario

### Suggestions for the scenario
- Add `cars-with-missing.csv` to the demo files package or note how to obtain it
- Step 5 (select all in Maximize) needs more detail: what exactly does the warning look like?
- Conflict warning step (step 4 in scenario) should reference specific column to avoid ambiguity
