# Box plot viewer — Run Results

**Date**: 2026-04-06
**URL**: http://localhost:8888
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open SPGI, SPGI-linked1, SPGI-linked2 | PASS | 5s | PASSED | All 3 datasets opened (3624, 3624, 224 rows). Tables renamed. |
| 2 | Go to Tables > SPGI | PASS | 1s | PASSED | Switched to SPGI view by matching 3624 rows, 88 cols |
| 3 | Click Box plot in Viewers tab | PASS | 2s | PASSED | Box plot added via `[name="icon-box-plot"]` click |
| 4 | Click Gear icon — Property Pane opens | PASS | 1s | PASSED | Gear icon found in `.panel-titlebar` ancestor, not inside `[name="viewer-Box-plot"]` |
| 5a | Tables switching (SPGI, linked2, linked1) | PASS | 3s | PASSED | Successfully switched `boxPlot.dataFrame` across all 3 tables |
| 5b | Set Filter to ${Average Mass} > 225 | PASS | 2s | PASSED | Filter applied via `boxPlot.props.filter`, 120 filtered rows shown |
| 5c | Set Color to Link Column 1 | PASS | 1s | PASSED | `markerColorColumnName` set to `link column 1` from linked table |
| 5d | Change Bin Color Aggr Type | PASS | 1s | PASSED | Changed from `avg` to `med`, set `binColorColumnName` to `Average Mass` |
| 5e | Save to Layout. Check | PASS | 5s | PASSED | All Data properties preserved: filter, markerColor, binColorAggr |
| 6a | Change axes | PASS | 1s | PASSED | Value → Average Mass, Category → Series |
| 6b | Invert Y Axis | PASS | 1s | PASSED | `invertYAxis` set to true, Y axis values go top-to-bottom |
| 6c | Save to Layout. Check | PASS | 5s | PASSED | Axes and invertY preserved in layout |
| 7a | Right-click box plot — context menu Tooltip | PASS | 2s | PASSED | Context menu has Tooltip section with Edit, Color Scheme, Use as Group Tooltip |
| 7b | Property Pane > Tooltip. Save to Layout | PASS | 5s | PASSED | `showTooltip = "show custom tooltip"` preserved in layout |
| 7c | Reset tooltip, Save to Layout | PASS | 5s | PASSED | `showTooltip = "inherit from table"` preserved in layout |
| 8a | Set Color to Link Column 1 | PASS | 1s | PASSED | Same as 5c |
| 8b | Change Bin Color Aggr Type | PASS | 1s | PASSED | Changed to `max` |
| 8c | Save to Layout. Check. Close all. | PASS | 5s | PASSED | Coloring preserved, all views closed |
| 9a | Open SPGI_v2, add Box Plot, adjust viewport | PASS | 4s | PASSED | Set valueMin=0, valueMax=10 |
| 9b | Save project, reopen, verify viewport preserved | SKIP | - | N/A | Project save via API threw `ApiException`; tested via layout instead |
| 9c | Save layout, reload, verify viewport NOT preserved | AMBIGUOUS | 6s | PASSED | Viewport IS preserved in layout (valueMin=0, valueMax=10 restored). Scenario says it should NOT be. |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~75s |
| Spec file generation | ~5s |
| Spec script execution | ~104s |

## Summary

Most Box Plot viewer functionality works correctly. Data properties (table switching, filter, color, bin color aggr type), axes (change columns, invert Y), tooltips (context menu, property pane), and coloring all pass and are correctly preserved in layouts. Step 9 (visualization zoom) could not fully test project save/restore due to an API exception, and the layout save behavior differs from the scenario expectation — viewport zoom IS preserved in layouts, while the scenario says it should NOT be.

## Retrospective

### What worked well
- All property-based viewer configuration works correctly
- Layout save/restore reliably preserves all Box Plot settings
- Table switching between linked tables works seamlessly
- Context menu is comprehensive with all expected options
- Linked table columns (e.g., `link column 1`) are accessible for color settings

### What did not work
- **Project save via API** — `grok.dapi.projects.save()` threw `ApiException`. Could not test project-based viewport preservation.
- **Range slider interaction** — Mouse event dispatch on range slider handles didn't trigger viewport zoom; had to use `props.valueMin/valueMax` directly.
- **Gear icon location** — The gear icon is not inside `[name="viewer-Box-plot"]` but in the `.panel-titlebar` of the `.panel-base` ancestor. This differs from the reference docs.

### Suggestions for the platform
- Expose project save/restore in a simpler API for test automation
- Consider making viewer title bar icons accessible within the `[name="viewer-*"]` container for simpler selector scoping
- Range slider zoom events should be dispatchable programmatically for test automation

### Suggestions for the scenario
- Step 1 says "Open SPGI, SPGI-linked1, SPGI-linked2" but doesn't mention linking the tables — add explicit linking step
- Column name "Link Column 1" should match actual column name `link column 1` (case mismatch)
- Step 9 references "SPGI_v2 dataset" which doesn't exist as a separate file — clarify if this is SPGI renamed
- Step 9 says viewport should NOT be preserved in layout — verify this is still the expected behavior, as current implementation preserves it
- Step 9 project save/restore instructions are vague — specify the API or UI method to use
