# Scripts Create — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Go to Browse > Platform > Functions > Scripts | PASS | Navigated via /scripts URL |
| 2 | NEW > New R script | PASS | Dropdown opened, R Script selected |
| 3 | Click "Open script sample table" (*) icon | PASS | "Table 'cars.csv' has been added" toast shown |
| 4 | Click Signature editor icon (magic wand) | PASS | Signature editor opened with PROPERTIES/PARAMETERS tabs |
| 5 | Set Name to testRscript | PASS | Name set via input field; reflected in CODE section |
| 6 | Navigate to Parameters tab | PASS | Parameters tab opened showing 2 existing params |
| 7 | Click ADD PARAMETER "+" button | PASS | Row 3 added: input/newParam/bool |
| 8 | Set direction=output, name=newParam, type=string | PASS | Changed via CodeMirror editor directly; line 7: `#output: string newParam` |
| 9 | Click Open function editor button | PARTIAL | Wrench clicked; warning shown: "Tests not found. Make sure to save latest changes to the script annotation." |
| 10 | Click Play, choose sample table (cars) | PARTIAL | Dialog opened but table combo didn't auto-populate; ran via JS API successfully |
| 11 | Click Save | PASS | "Script saved." toast shown; title updated to TestRscript |
| 12 | Close script view | PASS | View closed, returned to cars table view |

## Summary

The Create scenario completed successfully overall. The script `testRscript` was created, parameters configured, saved, and executed. The main issue was that the Run dialog's table dropdown did not auto-select the pre-loaded cars.csv sample table — the user must manually choose it. The script ran successfully when invoked via the JS API with the cars table.

## Retrospective

### What worked well
- NEW > R Script menu works correctly
- Sample table loading via asterisk icon works (cars.csv added)
- Signature editor opens and reflects code changes in real-time
- Save button works and updates the view title
- Script was saved and accessible in the browser

### What did not work
- Run dialog table combo didn't auto-populate with cars.csv even though it was loaded as sample — user must manually select from dropdown
- Open function editor (wrench) shows "Tests not found" warning — potentially confusing UX
- Canvas-based parameter grid cells (direction/type dropdowns) are not accessible via DOM automation

### Suggestions for the platform
- The Run dialog should auto-populate the table input with the sample table if one is loaded
- The wrench/function editor button could clarify it opens a "Test" panel, not just a general editor
- Consider making parameter grid accessible via standard HTML inputs to improve testability

### Suggestions for the scenario
- Step 9 "Open function editor" is ambiguous — clarify it opens the test/debug panel
- Step 10 should note that users must manually select the table in the dialog dropdown
- Add precondition: ensure cars.csv sample table is already open before running
