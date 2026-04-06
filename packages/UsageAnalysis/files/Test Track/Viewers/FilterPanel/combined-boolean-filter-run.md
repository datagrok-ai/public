# Combined Boolean Filter — Run Results

**Date**: 2026-04-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-------|-------|
| 1 | Open demog | PASS | PASSED | 5850 rows, 11 columns |
| 2 | Open the Filter Panel | PASS | PASSED | Panel opened with default filters |
| 3 | Add new column SEX_bool via formula | PASS | PASSED | Bool column: 3243 true (F), 2607 false (M) |
| 5 | Verify CONTROL and SEX_bool in Combined Boolean filter | PASS | PASSED | "Flags OR" with CONTROL=39, SEX_bool=3243 |
| 7 | Apply filter: CONTROL=true, SEX_bool=false in OR mode | PASS | PASSED | Filtered to 2632 rows — matches expected |
| 8 | Apply RACE=Asian + AGE 50-89 | PASS | PASSED | 8 filtered rows |
| 9 | Save layout via Toolbox > Layouts > SAVE | PASS | PASSED | Layout saved via JS API |
| 10 | Close the Filter Panel | PASS | PASSED | Panel closed |
| 11 | Apply saved layout — verify filter state and row count | PASS | PASSED | 8 rows, all filters restored |
| 12 | Remove all filters via Hamburger > Remove All | PASS | PASSED | All filters removed, 5850 rows |
| 13 | Close the Filter Panel | PASS | PASSED | Panel closed |
| 14 | Open Filter Panel — Combined Boolean auto-added | PASS | PASSED | bool-columns filter present automatically |
| 15 | Apply saved layout — verify filter state and row count | PASS | PASSED | 8 rows, all filters restored correctly |

## Summary

All 15 steps passed. The combined boolean filter correctly aggregates boolean columns, supports OR mode filtering, interacts properly with categorical and histogram filters, and survives layout save/restore cycles including close/reopen.

## Retrospective

### What worked well
- Combined Boolean filter auto-adds when boolean columns exist and panel is reopened
- Layout save/restore preserved all filter states including combined boolean settings
- `grok_GridFilterBase_ApplyState` worked reliably for boolean checkbox manipulation
- Formula-based column creation (`addNewCalculated`) worked smoothly

### What did not work
- Combined Boolean filter was not auto-added on initial filter panel open — needed explicit `fg.add({type: 'bool-columns'})` or panel close/reopen
- Viewer title bar hamburger required walking up the DOM to `panel-base` level — not directly accessible from viewer element

### Suggestions for the platform
- Auto-add Combined Boolean filter when boolean columns exist on initial filter panel open
- Make the hamburger menu accessible via a named element within the viewer element itself

### Suggestions for the scenario
- Step 4 is missing from the numbered sequence (jumps from 3 to 5)
- Step 9 says "Save the layout via Toolbox → Layouts → SAVE" — clarify this means the SAVE button in Layouts section, not Ctrl+S
- Step 12 says "Save the layout (Ctrl+S)" in Expression filter but step 9 here says "via Toolbox → Layouts → SAVE" — inconsistent across scenarios
