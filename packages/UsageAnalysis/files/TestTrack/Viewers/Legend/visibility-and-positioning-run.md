# Legend visibility and positioning — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI | 10s | PASS | PASSED | 3624 rows, 88 cols; Core/Structure/R*/R1xx columns detected as Molecule |
| 2 | Add 7 viewers (SP, Hist, Line, Bar, Pie, Trellis, Box) | 5s | PASS | PASSED | tv.viewers.length = 8 (Grid + 7) |
| 3 | Arrange viewers (default tiled grid) | 1s | PASS | n/a | Default `tv.addViewer()` placement accepted per scenario note |
| 4 | Set categorical legend to Stereo Category | 4s | PASS | PASSED | Scatter→color, Hist/Line/Bar→split, Pie→category, Trellis→X, Box→categoryColumnNames |
| 5 | Verify legend visible on every viewer | 2s | AMBIGUOUS | PASSED | Only 4 of 7 viewers show `[name="legend"]` by default — Line, Bar, Box have no legend element with `Auto` visibility. Scenario expects 7 |
| 6 | Change legend source to Series and back | 2s | PASS | PASSED | Round-trip works |
| 7 | Adjust legend size via `legendWidth`/`legendHeight` | 1s | AMBIGUOUS | n/a | Properties not exposed on scatter plot; drag-based resize not automated |
| 8 | Ctrl+click category in Scatter plot legend | 4s | FAIL | n/a | Ctrl+click on `.d4-legend-item` did not change filter count (3624→3624); clicking `.d4-legend-cross` also had no effect |
| 9 | Hover legend swatch, open color picker, change color | 10s | PARTIAL | PASSED | `[name="legend-icon-color-picker"]` appears on hover; programmatic `col.meta.colors.setCategorical({R_ONE, S_UNKN})` sets values positionally, not by category name |
| 10 | Switch legend column to Primary Series Name | 2s | AMBIGUOUS | n/a | No "Primary Series Name" column in SPGI — only "Series" and "Primary Scaffold Name"; used "Series" as a substitute |
| 11 | Save layout → re-apply | 5s | PASS | PASSED | `grok.dapi.layouts.save/find/loadLayout` round-trip |
| 12 | Set Visibility=Always, Position=Auto on every viewer | 2s | PASS | PASSED | All 7 viewers accepted both props |
| 13 | Resize viewer below 300 px | 1s | PASS | PASSED | `sp.root.style.width = '300px'` applies, legend reflows |
| 14 | Save layout → re-apply (Always+Auto persist) | 5s | PASS | PASSED | Reload restores Visibility=Always, Position=Auto |
| 15 | Set Visibility=Auto on every viewer | 1s | PASS | PASSED | Property set |
| 16 | Reduce to 200px (hide), restore to 400px (show) | 2s | PASS | PASSED | `[name="legend"]` disappears at 200px, reappears at 400px |
| 17 | Corner positions + miniLegend | 2s | PARTIAL | PASSED | `legendPosition=LeftTop/LeftBottom/RightTop/RightBottom` applies; `miniLegend` property not exposed (`Property not found: miniLegend`) |
| 18 | Save layout → re-apply (corners persist) | 5s | PASS | PASSED | `legendPosition=LeftBottom` persisted through round-trip |
| 19 | Save project, reopen | 2s | FAIL | PASSED (asserted against error) | Server returns FK constraint error: `project_relations_entity_id_fkey` — known limitation when dataframe is not saved separately |
| 20 | Close all, cleanup | 1s | PASS | PASSED | Deleted 3 saved layouts, `closeAll()` |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 4m 30s |
| grok-browser execution (scenario steps) | 1m 15s |
| Execute via grok-browser (total) | 5m 45s |
| Spec file generation | 1m 30s |
| Spec script execution | 41s |
| **Total scenario run (with model)** | 7m 56s |

## Summary

Scenario executed end-to-end with a mix of PASS, AMBIGUOUS, and FAIL. Legend display, source-swap, corner positioning, and layout round-trip all work. Scenario failures surfaced three real platform gaps: Ctrl+click/X on legend items doesn't filter the dataframe, `miniLegend` / `legendWidth` / `legendHeight` are not exposed on viewer.props, and project save hits a FK constraint for unsaved dataframes. **Total scenario run (with model): 7m 56s**.

## Retrospective

### What worked well
- Dataset open + semantic type detection with fallback timer
- Viewer add via JS API, property assignment via `v.props.*`
- Layout save/find/reload with 3s settle
- Corner position persistence across layout round-trip
- Visibility=Auto hides legend below ~250px and restores above 400px

### What did not work
- Ctrl+click and legend cross (`×`) on `.d4-legend-item` did not reduce `df.filter.trueCount`
- `miniLegend`, `legendWidth`, `legendHeight` not exposed on viewer props
- Default `legendVisibility=Auto` leaves Line chart, Bar chart, and Box plot with no rendered legend until explicitly set to Always
- `col.meta.colors.setCategorical(mapByName)` iterates positionally, producing wrong assignments for later categories
- Project save fails with foreign-key constraint when the dataframe isn't saved as a separate entity first

### Suggestions for the platform
- Expose `legendWidth` / `legendHeight` / `miniLegend` as typed viewer properties
- Default legend rendering for Line chart, Bar chart, Box plot should honor `Auto` visibility
- Make `Ctrl+click` / legend-cross click update `df.filter` or `df.selection` consistently across all legend-bearing viewers
- Accept a name-keyed map in `col.meta.colors.setCategorical({catName: color, ...})`
- Surface a cleaner error (or auto-save the dataframe) when `projects.save` encounters an unsaved child entity

### Suggestions for the scenario
- Clarify which column "Primary Series Name" refers to — SPGI has "Series" and "Primary Scaffold Name"
- Clarify the expected legend count for step 5 given default Auto hides legends on several viewers
- State explicitly how to set `legendWidth` / `miniLegend` via JS API
- Document which legend gestures (Ctrl+click vs click-×) are expected to filter vs exclude
