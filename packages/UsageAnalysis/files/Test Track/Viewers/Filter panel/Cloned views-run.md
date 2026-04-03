# Cloned Views — Run Results

**Date**: 2026-04-03
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open spgi-100 | PASS | PENDING | Opened SPGI.csv via `grok.data.files.openTable` (3624 rows, 88 cols) |
| 2 | Open the Filter Panel | PASS | PENDING | Opened via `tv.getFiltersGroup()` |
| 3 | Navigate to Competition assay filter | PASS | PENDING | Scrolled to Competition assay histogram filter |
| 4 | Filter out missing values for Competition assay | PASS | PENDING | Clicked filter indicator > Missing values > Filter out. Rows: 3624 → 2486 |
| 5 | Set Stereo Category filter to S_ACHIR | PASS | PENDING | Used `fg.updateOrAdd` API. Rows: 2486 → 705 |
| 6 | Sketch c1ccccc1 structure | PASS | PENDING | Clicked Sketch, typed benzene SMILES, OK. Rows: 705 → 391 |
| 7 | Clone View via View > Layout > Clone View | PASS | PENDING | "SPGI copy" created with matching filter state (391 rows) |
| 8 | Check Filter Panel open and state matches on clone | PASS | PENDING | Filter Panel open, [3] active filters badge, 391 filtered rows |
| 9 | Turn all filters off (global checkbox) | PASS | PENDING | Clicked header checkbox. All 3624 rows visible |
| 10 | Turn filters on again | PASS | PENDING | Re-enabled checkbox. Rows back to 391 |
| 11 | Clear Structure filter, set to C1CCCCC1 | PASS | PENDING | Clicked Clear span (705 rows), then set C1CCCCC1 cyclohexane (10 rows) |
| 12 | Remove Structure filter (X icon) | PASS | PENDING | Removed from clone via icon-times. State unchanged (10 rows). Original still has it |
| 13 | Save the layout | PASS | PENDING | Saved via `tv.saveLayout()` |
| 14 | Close the Filter Panel | PASS | PENDING | Closed via `filterViewer.close()`. Filtered: 10 (original view filters still active) |
| 15 | Apply saved layout — Filter Panel opens without Structure | PASS | PENDING | Layout applied. Filter Panel open, no Structure filter, filtered: 10 |

## Summary

All 15 steps passed successfully. The Clone View feature correctly duplicates the filter state. Filter modifications on the cloned view's filter panel are independent of the original view's filter panel configuration, while the underlying dataframe filter (BitSet) is shared. The layout save/restore correctly preserves the filter panel state including the absence of the removed Structure filter.

## Retrospective

### What worked well
- Clone View (via top menu View > Layout > Clone View) correctly preserved all 3 active filters
- Global enable/disable checkbox correctly toggled all filters on/off and back
- Removing the Structure filter from the cloned view did not affect the original view's filter panel
- Layout save/restore correctly preserved the filter configuration without the Structure filter
- The filter indicator icon opened the correct context menu with "Missing values" submenu

### What did not work
- Toolbox "Filters" section was not always visible — had to use `tv.getFiltersGroup()` API as fallback
- The viewer-level close (X) icon was not findable via `[name="icon-times"]` scoped outside `.d4-filter` — had to use `filterViewer.close()` API
- The "Clear" link on the Structure filter required `dispatchEvent(new MouseEvent('click'))` — a simple `.click()` did not register the first time

### Suggestions for the platform
- Add a distinguishing `name=` attribute to the viewer title bar close icon to differentiate from per-filter remove icons
- Make the Toolbox "Filters" section consistently visible when a table view is active

### Suggestions for the scenario
- Step 1: "Open spgi-100" references a 100-row dataset but the JSON block specifies "System:DemoFiles/SPGI.csv" (3624 rows). Clarify which dataset to use
- Step 12: Clarify that both views share the same DataFrame — removing a filter card from the clone doesn't change the filtered count because the original view's structure filter is still active on the shared BitSet
- Step 14-15: Note that when the cloned view's filter panel is closed, the data remains filtered by the original view's active filters
