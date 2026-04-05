# Cloned Views — Run Results

**Date**: 2026-04-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-------|-------|
| 1 | Open spgi-100 (SPGI.csv) | PASS | PASSED | 3624 rows, 88 columns, Chem dataset |
| 2 | Close Context Panel and Context Help | PASS | PASSED | Via JS API |
| 3 | Open Filter Panel | PASS | PASSED | 43 filter cards created |
| 4 | Navigate to Competition assay filter | PASS | PASSED | Filter found |
| 5 | Filter out missing values for Competition assay | PASS | PASSED | Via fg.updateOrAdd with filterOutMissingValues |
| 6 | Set Stereo Category to S_ACHIR | PASS | PASSED | 705 rows after both filters |
| 7 | Sketch c1ccccc1 in Structure filter | PASS | PASSED | 391 rows, used sketcher dialog |
| 8 | Clone View via View > Layout > Clone View | PASS | PASSED | "Table copy" created |
| 9 | Verify cloned view matches original | PASS | PASSED | 391 rows, panel open, structure filter present |
| 10 | Turn all filters off (global checkbox) | PASS | PASSED | 3624 rows visible |
| 11 | Turn filters back on | PASS | PASSED | 391 rows restored |
| 12 | Clear Structure, set to C1CCCCC1 | PASS | PASSED | 10 rows (cyclohexane substructure) |
| 13 | Remove Structure filter from clone | PASS | PASSED | Count stayed at 10, original still has Structure filter |
| 14 | Save layout | PASS | PASSED | Layout saved via grok.dapi.layouts.save |
| 15 | Close Filter Panel | PASS | PASSED | Closed via fg.close() |
| 16 | Apply saved layout | PASS | PASSED | Panel opens, no Structure filter, Filtered: 10 |

## Summary

All 16 steps passed. Cloned views correctly inherit the full filter panel state including categorical, histogram, and structure filters. The global filter toggle works on cloned views. Layout save/restore preserves the filter configuration — the removed Structure filter stays absent after layout restore. Both views share the same DataFrame, so filter value changes in one view affect the other.

## Retrospective

### What worked well
- Clone View correctly duplicates all filter panel state (3 active filters)
- Global filter toggle (checkbox) enables/disables all filters properly
- Structure filter sketcher dialog works reliably for both setting and clearing
- Layout save/restore correctly captures filter panel state without the removed Structure filter
- Filter panel close via `fg.close()` works when title bar X icon is not accessible

### What did not work
- Closing filter panel via title bar X icon failed — the `icon-times` selector only matched per-filter remove icons, not the viewer-level close. Had to fall back to `fg.close()` JS API

### Suggestions for the platform
- The filter panel viewer should have a clearly distinct close icon selector (e.g., on a `.d4-viewer-title-bar`)

### Suggestions for the scenario
- Step 13 wording "filtered state should not change" is ambiguous — the count stays the same because views share a DataFrame, not because removing a filter has no effect. Suggest clarifying: "the filtered row count should remain the same since both views share the same DataFrame"
- Step 2 "Close the Context Panel and the Context Help" is only achievable via JS API — consider noting this
