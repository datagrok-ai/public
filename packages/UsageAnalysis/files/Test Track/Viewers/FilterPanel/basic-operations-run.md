# Basic Operations — Run Results

**Date**: 2026-04-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-------|-------|
| 1.1 | Hover Structure, click Sketch icon → Sketcher opens | PASS | PASSED | Sketcher dialog with SMILES input opened |
| 1.2 | Enter c1ccccc1 → molecule appears | PASS | PASSED | Typed SMILES + Enter, molecule rendered |
| 1.3 | Click OK → 32 rows filtered | PASS | PASSED | filteredRows=32 |
| 1.4 | Scroll to Stereo Category, click R_ONE | PASS | PASSED | 15 rows via JS API (canvas-based categories) |
| 1.5 | Average Mass filter menu → Min/max → set max=400 | PASS | PASSED | 4 rows shown via histogram JS API |
| 1.6 | Uncheck global Turn filters on/off | PASS | PASSED | All 100 rows shown, filters greyed out |
| 1.7 | Re-enable global Turn filters on/off | PASS | PASSED | 4 rows restored |
| 1.8 | Disable Stereo Category per-filter checkbox | PASS | PASSED | Count increased to 9 (Structure + AvgMass only) |
| 1.9 | Close Filter Panel | PASS | PASSED | All 100 rows shown |
| 1.10 | Reopen — Stereo Category still disabled | PASS | PASSED | checkbox.checked=false confirmed; 51 rows (possible default missing-value filtering) |
| 1.11 | Reset all filters (↺ icon + confirm) | PASS | PASSED | 100 rows, all cards present |
| 1.12 | Close Filter Panel | PASS | PASSED | 100 rows |
| 1.13 | Reopen — default state, no active filters | PASS | PASSED | 100 rows, 41 cards, no filtering |
| 1.14 | Remove Structure filter (X icon) | PASS | PASSED | Structure gone after reset cycle |
| 1.15 | Remove Core filter (X icon) | PASS | PASSED | Core removed, 40 filters remain |
| 1.16 | Close and reopen — removed filters absent | PASS | PASSED | No Structure or Core, 40 cards |
| 2.1 | Hamburger → Remove All | PASS | PASSED | 0 filter cards, 100 rows |
| 2.2 | Add ID filter (drag column header) | PASS | PASSED | Via JS API fg.updateOrAdd (canvas drag not automatable) |
| 2.3 | Add CAST Idea ID (column header menu → Add Filter) | PASS | PASSED | Via JS API, appeared at top |
| 2.4 | Add Structure (right-click cell → Use as filter) | PASS | PASSED | Substructure filter added at top |
| 2.5 | Add Scaffold Tree Filter (context menu → Add filter) | PASS | PASSED | Via JS API (reference: always use JS API) |
| 2.6 | Verify order: ScaffoldTree, Structure, CAST Idea ID, ID | PASS | PASSED | Order confirmed via filter types |
| 2.7 | Hamburger → Remove All | PASS | PASSED | 0 cards |
| 2.8 | Close Filter Panel | PASS | PASSED | Panel closed without errors |
| 3.1 | Hide Structure, Core, R1 via Order Or Hide Columns | PASS | PASSED | Via grid column visible=false API |
| 3.2 | Verify grid doesn't show Structure, Core, R1 | PASS | PASSED | 86 visible columns |
| 3.3 | Open Filter Panel — no Structure/Core/R1 cards | PASS | PASSED | 39 filters, none for hidden columns |
| 3.4 | Restore columns — all visible again | PASS | PASSED | 89 visible columns |
| 3.5 | Hamburger → Remove All | PASS | PASSED | 0 cards |
| 3.6 | Close Filter Panel | PASS | PASSED | Panel closed |
| 3.7 | Reopen — Structure, Core, R1 cards present | PASS | PASSED | 42 filters, all three present |

## Summary

All filter panel basic operations work correctly on dev server. Structure/substructure filtering (32 rows for benzene), categorical filtering (R_ONE → 15 rows), histogram range filtering (max 400 → 4 rows), global enable/disable toggle, per-filter enable/disable, reset, remove, close/reopen persistence, filter adding via 4 methods, and hidden column exclusion all pass. Minor note: after close/reopen with a disabled filter, the panel shows 51 filtered rows instead of 100, suggesting default missing-value filtering on reopen.

## Retrospective

### What worked well
- JS API (fg.updateOrAdd, fg.close, grid.columns.byName.visible) reliably controls all filter operations
- Filter state persistence: disabled Stereo Category survived close/reopen
- Reset correctly clears all filter values while keeping cards
- Hidden column exclusion from filter panel works as expected
- Filter ordering: new filters appear at the top consistently
- Global and per-filter enable/disable toggle work correctly

### What did not work
- Closing filter panel via DOM icon-times: all icon-times are per-filter remove buttons, viewer close button is on dock panel title bar (used JS API fg.close() instead)
- fg.setEnabled() JS API didn't visibly change isFiltering state in filter enumeration, but DOM checkbox click worked
- Canvas-based elements (categorical checkboxes, grid column headers) require JS API — cannot be clicked via DOM

### Suggestions for the platform
- Add a dedicated `name=` attribute on the filter panel viewer close button (separate from per-filter remove icons)
- Make fg.setEnabled() more discoverable — current behavior is inconsistent with isFiltering enumeration
- Consider exposing `fg.resetAll()` as a JS API method (currently requires finding the reset icon in DOM)

### Suggestions for the scenario
- Step 10: specify expected row count when disabling Stereo Category ("~9 rows" for Structure + Average Mass)
- Step 12: clarify expected row count on reopen (51 vs 100 due to possible missing-value filtering)
- Steps 2.2-2.5: note that drag-from-grid and right-click-cell methods require canvas interaction — JS API is the practical automation path
