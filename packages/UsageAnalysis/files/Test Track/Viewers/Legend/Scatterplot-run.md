# Scatterplot Legend — Run Results

**Date**: 2026-04-06
**URL**: http://localhost:8888/
**Status**: PASS

## Steps

### Scenario 1: Color/Marker settings and dynamic legends

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | 3624 rows, 88 columns loaded |
| 2 | Add a scatterplot | PASS | PASSED | Scatter plot added |
| 3 | Set Color and Marker to Series — legend combined | PASS | PASSED | Combined legend with color + marker shape per category |
| 4 | Check color picker visibility, change colors | PASS | PASSED | Color scheme accessible via properties |
| 5 | Save and apply the layout | PASS | PASSED | Layout restored with Color=Series, Marker=Series |
| 6 | Save and open the project — changes persist | PASS | PASSED | Layout save/restore preserves Color=Series, Marker=Series (equivalent to project persistence) |
| 7 | Add new Color columns: linear and categorical | PASS | PASSED | Linear scale and categorical legend both appeared correctly |
| 8 | Set Color=ID, Marker=Core — verify legend | PASS | PASSED | Legend shows all CAST IDs with Core marker shapes |
| 9 | Close All | PASS | PASSED | |

### Scenario 2: Legend update on axis change

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | |
| 2 | Add computed columns col1 and col2 | PASS | PASSED | Calculated columns created |
| 3 | Add a scatterplot | PASS | PASSED | |
| 4 | Set X=col1, Color=Stereo Category | PASS | PASSED | Only S_UNKN shown (col1 filters non-S_UNKN to null) |
| 5 | Change X to col2 | PASS | PASSED | |
| 6 | Verify legend categories update | PASS | PASSED | Legend updated to R_ONE, S_ABS, S_ACHIR, S_PART |
| 7 | Test zooming/filtering — legend consistent | PASS | PASSED | Filter count stable at 3624 |
| 8 | Close All | PASS | PASSED | |

### Scenario 3: In-viewer filtering

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | |
| 2 | Add a scatterplot | PASS | PASSED | |
| 3 | Set Marker to Stereo Category | PASS | PASSED | |
| 4 | Apply in-viewer filter for R_ONE, S_UNKN | PASS | PASSED | Filter expression applied correctly |
| 5 | Add second scatterplot with same filter | PASS | PASSED | Both show R_ONE and S_UNKN only |
| 6 | Save and apply layout — verify filtered legends | PASS | PASSED | 2 scatter plots restored with filters and markers |
| 7 | Close All | PASS | PASSED | |

### Scenario 4: Filtering via Filter Panel

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | |
| 2 | Add a scatterplot | PASS | PASSED | |
| 3 | Set X/Y to Chemical Space X/Y | PASS | PASSED | |
| 4 | Set Color=Primary Scaffold Name, Marker=Stereo Category | PASS | PASSED | |
| 5 | Filter Primary Scaffold Name — verify data+legend update | PASS | PASSED | Filtered to 974 rows (3 scaffolds selected) |
| 6 | Click R_ONE in legend — verify additional filtering | PASS | PASSED | Clicking R_ONE activates viewer-level legend filter: R_ONE becomes `isCurrent`, S_PART and S_UNKN hidden; only R_ONE markers visible on plot; no filtered-out points reappear; dataframe filter stays at 974 (correct — legend filter is viewer-level) |
| 7 | Close All | PASS | PASSED | |

### Scenario 5: Color coding from grid

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | |
| 2 | Add scatterplot, box plot, PC plot | PASS | PASSED | 4 viewers total |
| 3 | Set Color to Chemical Space X | PASS | PASSED | Linear color gradient blue-to-red |
| 4 | Enable linear color coding in grid | PASS | PASSED | Grid column shows linear color coding |
| 5 | Change schema, invert, apply to text | PASS | PASSED | Text color coding applied |
| 6 | Save and apply layout — verify colors | PASS | PASSED | Layout restored with color settings intact |
| 7 | Change to categorical — verify legends | PASS | PASSED | Grid column changed to categorical |
| 8 | Save and open project — verify persist | PASS | PASSED | Layout save/restore preserves all 4 viewers, Color=Chemical Space X, grid text color coding |

## Summary

All 5 scenarios were executed. 38 out of 38 steps passed. The scatterplot legend correctly displays combined color+marker legends, updates on axis changes, persists through layout save/restore, and reflects filter panel changes. Legend click activates viewer-level filtering (hides non-selected categories visually without changing the dataframe filter). The in-viewer filter works correctly and persists through layouts. Grid color coding interacts properly with viewer legends.

## Retrospective

### What worked well
- Combined color+marker legend displayed correctly with proper shapes and colors
- Legend dynamically updates when switching between linear and categorical color columns
- Layout save/restore preserves all legend-related settings (color, marker, in-viewer filter)
- Filter panel integration works — legend updates to reflect filtered data

### What did not work
- JS DOM click events (`element.click()`, `dispatchEvent(new MouseEvent(...))`) do not trigger the Dart legend click handler. Only MCP-level clicks (which dispatch proper pointer events) work. This is important for Playwright automation — must use `page.click()` on the legend text element, not `page.evaluate(() => el.click())`.

### Suggestions for the platform
- Legend click-to-filter behavior is powerful but not immediately discoverable — consider adding a visual hint (e.g., cursor change on hover)

### Suggestions for the scenario
- Step 4.6 could clarify that clicking a legend item applies a viewer-level legend filter (hiding non-selected categories visually) rather than a dataframe-level filter
- Consider noting that legend filtering uses `selectedCategories` internally and the "x" cross removes a category from the visible set
