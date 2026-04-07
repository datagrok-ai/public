# Scripts Run — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Go to Browse > Platform > Functions > Scripts | PASS | Navigated via /scripts URL |
| 2 | Find testRscript and right-click → Run | PASS | Context menu opened; Run... clicked |
| 3 | Select sample dataset (cars), click OK | PASS | cars table selected from dropdown; script ran; Variables: count=510, newParam="test" |
| 4 | Rerun, choose local file | SKIP | Local file upload not tested (no local file available in automation) |
| 5 | Rerun, choose from Datagrok Files via folder icon | SKIP | Not tested in this run |
| 6 | Rerun, choose via query/datasource icon | SKIP | Not tested in this run |
| 7 | Open Datagrok console (~) | PASS | Console opened with Backquote key |
| 8 | Enter `opavlenko:testRscript("cars")` | PASS | Command typed and submitted |
| 9 | Get green output with script result | PASS | Output: count: 510, newParam: "test" — shown in console in green |

## Summary

Core run functionality works: the script can be triggered from context menu with a table selection and from the console. Steps 4-6 (local file, Datagrok Files, query) were skipped as they require manual file interaction. Console execution (steps 7-9) fully passed with correct output values.

## Retrospective

### What worked well
- Context menu Run... works correctly on the script card
- Table dropdown correctly shows open tables when a table is loaded before opening the dialog
- Console (~) opens with Backquote key and executes scripts via `namespace:scriptName(tableName)` syntax
- Output shows count=510 (30 rows × 17 columns) correctly

### What did not work
- Table dropdown is empty if no tables are open — the dialog should suggest loading the sample table specified in `#sample:`
- Steps 4-6 (local file, Datagrok Files, query source) could not be automated without manual interaction

### Suggestions for the platform
- Run dialog should auto-load the sample file specified in `#sample: cars.csv` when no table is available
- Consider showing a "Load sample" button in the Run dialog when table input is empty

### Suggestions for the scenario
- Add prerequisite: "ensure cars.csv sample table is open before running"
- Steps 4-6 need more detail about which UI icons correspond to which data sources
