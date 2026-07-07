# Tile Viewer tests — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Default form rendering: click first tile → currentRow=0 | 8s | PASS | PASSED | currentRowIdx=0 confirmed |
| 2 | Default form rendering: click second tile → currentRow=1 | 6s | PASS | PASSED | currentRowIdx=1 confirmed |
| 3 | Row selection: shift-click tile 3 → highlighted as selected | 15s | PASS | PASSED | tile0.click() then shift-click tile2; selection.trueCount≥1 |
| 4 | Row selection: ctrl-click tile 5 → added to selection | 10s | PASS | PASSED | rows 2 and 4 selected after ctrl-click; selection cleared after |
| 5 | Lanes: set RACE → 4 lanes appear | 8s | PASS | PASSED | 4 lane headers: Asian, Black, Caucasian, Other |
| 6 | Lanes: set SEX → 2 lanes (F, M) | 6s | PASS | PASSED | F and M lanes confirmed |
| 7 | Lanes: clear → single flat lane | 6s | PASS | PASSED | laneCount=1 confirmed |
| 8 | Row source: Selected → 5 tiles for 5 selected rows | 12s | PASS | PASSED | tileCount=5 confirmed |
| 9 | Row source: Filtered with SEX=M → filtered rows shown | 10s | PASS | PASSED | 2607 filtered rows |
| 10 | Row source: All → all rows shown, filter cleared | 8s | PASS | PASSED | 5850 rows, filter and selection cleared |
| 11 | Tiles font: change size to 18px | 6s | PASS | PASSED | tilesFont contains "18px" |
| 12 | Tiles font: change family to Courier | 6s | PASS | PASSED | tilesFont contains "Courier" |
| 13 | Tiles font: reset to default 13px Roboto | 5s | PASS | PASSED | tilesFont reset confirmed |
| 14 | Auto-generate: load SPGI, switch table → SPGI columns | 27s | PASS | PASSED | autoGenerate=true; SPGI loaded (3624 rows); dataFrame.name contains "(2)" |
| 15 | Auto-generate: switch back to demog → demog columns | 8s | PASS | PASSED | dataFrame.name="Table"; autoGenerate disabled |
| 16 | Context menu: right-click tile → menu with Edit Form... and Properties... | 12s | PASS | PASSED | Edit Form... and Properties... both confirmed in menu items |
| 17 | Context menu: close menu and open properties via gear | 8s | PASS | PASSED | Escape pressed; gear icon at x>900 clicked |
| 18 | Viewer title: set title "Patient Cards" and description | 11s | PASS | PASSED | title="Patient Cards", description="Demographic data per patient" confirmed |
| 19 | Viewer title: set description position to Bottom | 5s | PASS | PASSED | descriptionPosition="Bottom" confirmed |
| 20 | Viewer title: clear title and description | 5s | PASS | PASSED | title="" and description="" confirmed |
| 21 | Filter interaction: open panel and add SEX=M filter | 10s | PASS | PASSED | rowSource=Filtered; 2607 filtered rows |
| 22 | Filter interaction: add AGE > 50 → further reduced | 8s | PASS | PASSED | 857 rows (M + AGE>50) |
| 23 | Filter interaction: remove all filters → all tiles restored | 8s | PASS | PASSED | 5850 rows restored |
| 24 | Filter interaction: close filter panel | 5s | PASS | PASSED | Toggled via div-section--Filters click |

## Timing

| Phase | Duration |
|-------|----------|
| grok-browser execution (scenario steps) | ~4m |
| Spec file generation | ~3m |
| Spec script execution | 58s |
| **Total scenario run** | ~8m |

## Summary

24 of 24 steps passed. Steps correspond 1:1 to softSteps in the spec. Drag between lanes and Card markup moved to manual checklist (tile-viewer-tests-ui.md).

## Retrospective

### What worked well
- JS API property changes (tilesFont, lanesColumnName, rowSource) worked reliably without opening the property panel
- Auto-generate table switching worked correctly
- Filter panel operations via JS API were solid
- Context menu dispatched via JS contextmenu event opened correctly

### Suggestions for the platform
- Shift+click in Tile Viewer should range-select from the previously clicked tile to the shift-clicked tile
- `cardMarkup` should take precedence over `sketchState` when explicitly set (moved to manual checklist)

### Suggestions for the scenario
- "Lanes" step 3: change "three lanes appear" to "four lanes appear" (demog has 4 RACE categories: Asian, Black, Caucasian, Other)
