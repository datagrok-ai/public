# Cloned Views — Run Results

**Date**: 2026-04-06
**URL**: https://dev.datagrok.ai
**Dataset**: System:DemoFiles/SPGI.csv (spgi-100.csv not found on dev; SPGI.csv has same columns, 3624 rows)
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open spgi-100 | PASS | PASSED | Used SPGI.csv (3624 rows, 88+2 cols). spgi-100.csv missing on dev server (502/not found) |
| 3 | Open the Filter Panel | PASS | PASSED | Opened via `tv.getFiltersGroup()`, 43 filter cards created |
| 4 | Navigate to Competition assay filter | PASS | PASSED | Filter card present as Histogram type |
| 5 | Filter out missing values for Competition assay | PASS | PASSED | Used `fg.updateOrAdd({filterOutMissingValues: true})` — UI submenu unreliable per reference |
| 6 | Set Stereo Category to S_ACHIR | PASS | PASSED | Used `fg.updateOrAdd({type: CATEGORICAL, selected: ['S_ACHIR']})`. Result: 705 filtered rows |
| 7 | Sketch c1ccccc1 structure | PASS | PASSED | Clicked `.sketch-link`, typed SMILES via keyboard, pressed Enter, clicked OK. Result: 391 filtered rows |
| 8 | Clone View via View > Layout > Clone View | PASS | PASSED | Used dispatchEvent for menu navigation. "Table copy" opened |
| 9 | Verify cloned view filter state matches original | PASS | PASSED | Filter panel open, 391 rows, structure canvas present, S_ACHIR and Competition assay active |
| 10 | Turn all filters off | PASS | PASSED | Clicked `.d4-filter-group-header input[type="checkbox"]`. All 3624 rows visible |
| 11 | Turn filters back on | PASS | PASSED | Re-clicked checkbox. Back to 391 filtered rows |
| 12 | Clear Structure, set to C1CCCCC1 | PASS | PASSED | Clicked `.chem-clear-sketcher-button` (705 rows), then set C1CCCCC1 (10 rows) |
| 13 | Remove Structure filter — state should not change | PASS | PASSED | Clicked X icon on Structure card. Count stayed at 10 (original view still filtering) |
| 14 | Save the layout | PASS | PASSED | Saved via `grok.dapi.layouts.save()` |
| 15 | Close the Filter Panel | PASS | PASSED | Closed via `fg.close()`. Title bar X icon not accessible via DOM |
| 16 | Apply saved layout — panel opens without Structure | PASS | PASSED | Filter panel restored, no Structure filter, 10 filtered rows preserved |

## Summary

All 14 scenario steps passed on dev server. The cloned view correctly inherits filter state from the original view, including the structure substructure filter. Removing the Structure filter from the cloned view does not affect filtering because the original view's Structure filter remains active on the shared DataFrame. Layout save/restore correctly preserves the filter configuration without the removed Structure filter.

## Retrospective

### What worked well
- JS API for filter manipulation (`fg.updateOrAdd`, `fg.close()`) is reliable and fast
- Clone View via menu `dispatchEvent` works correctly
- Layout save/restore preserves filter state accurately
- `.chem-clear-sketcher-button` is always visible (no hover needed) for clearing structure filters
- Global filter toggle checkbox works as expected on cloned views

### What did not work
- `spgi-100.csv` does not exist on dev server — had to use `SPGI.csv` instead. The scenario JSON should reference the correct file
- Filter panel close via title bar X icon: the `.grok-font-icon-close` element is not accessible via standard DOM selectors or a11y tree. Had to use `fg.close()` instead
- Sketcher dialog SMILES input requires real keyboard events — `input.value = '...'` + dispatchEvent doesn't trigger Dart change listeners reliably

### Suggestions for the platform
- Add a `name=` attribute to the filter panel's close icon in the title bar for easier automation
- Consider adding a JS API for setting structure filter SMILES directly (e.g., `fg.updateOrAdd({type: 'Chem:substructureFilter', column: 'Structure', smiles: 'c1ccccc1'})`)

### Suggestions for the scenario
- Update dataset reference from `spgi-100.csv` to `SPGI.csv` — the file doesn't exist on dev/public servers
- Step numbering jumps from 1 to 3 (step 2 is missing)
- Clarify step 13: "filtered state should not change" — specify this means the row count stays the same because the original view's Structure filter is still active on the shared DataFrame
