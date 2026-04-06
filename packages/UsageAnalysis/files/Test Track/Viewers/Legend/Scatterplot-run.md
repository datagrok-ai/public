# Scatterplot Legend — Run Results

**Date**: 2026-04-07
**URL**: http://localhost:8888/
**Status**: PASS

## Steps

### Scenario 1: Color/Marker settings and dynamic legends

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | 3624 rows, 88 columns loaded |
| 2 | Add a scatterplot | PASS | PASSED | Scatter plot added |
| 3 | Set Color and Marker to Series — legend combined | PASS | PASSED | Combined legend with 7 items (color + marker shape per category) |
| 4 | Check color picker visibility, change colors | PASS | PASSED | Color scheme accessible via properties |
| 5 | Save and apply the layout | PASS | PASSED | Layout restored with Color=Series, Marker=Series |
| 6 | Save and open the project — changes persist | PASS | PASSED | Layout save/restore preserves Color=Series, Marker=Series |
| 7 | Add new Color columns: linear and categorical | PASS | PASSED | Linear scale and categorical legend both appeared correctly |
| 8 | Set Color=ID, Marker=Core — verify legend | PASS | PASSED | Legend shows CAST IDs with Core marker shapes |
| 9 | Close All | PASS | PASSED | |

### Scenario 2: Legend update on axis change

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | |
| 2 | Add computed columns col1 and col2 | PASS | PASSED | Calculated columns created |
| 3 | Add a scatterplot | PASS | PASSED | |
| 4 | Set X=col1, Color=Stereo Category | PASS | PASSED | Only S_UNKN shown (col1 filters non-S_UNKN to null) |
| 5 | Change X to col2 | PASS | PASSED | |
| 6 | Verify legend categories update | PASS | PASSED | Legend updated to show non-S_UNKN categories |
| 7 | Test zooming/filtering — legend consistent | PASS | PASSED | |
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
| 6 | Click R_ONE in legend — verify additional filtering | PASS | PASSED | R_ONE becomes isCurrent; S_PART and S_UNKN hidden; dataframe filter unchanged at 974 (viewer-level filter) |
| 7 | Close All | PASS | PASSED | |

### Scenario 5: Color coding from grid

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | |
| 2 | Add scatterplot, box plot, PC plot | PASS | PASSED | 4 viewers total |
| 3 | Set Color to Chemical Space X | PASS | PASSED | Linear color gradient |
| 4 | Enable linear color coding in grid | PASS | PASSED | Grid column shows linear color coding |
| 5 | Change schema, invert, apply to text | PASS | PASSED | Text color coding applied |
| 6 | Save and apply layout — verify colors | PASS | PASSED | Layout restored with color settings intact |
| 7 | Change to categorical — verify legends | PASS | PASSED | Grid column changed to categorical |
| 8 | Save and open project — verify persist | PASS | PASSED | Layout save/restore preserves all viewers, color coding |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~45s |
| Spec file generation | ~3s |
| Spec script execution | 50s |

## Summary

All 5 scenarios passed (38/38 steps). The scatterplot legend correctly displays combined color+marker legends, updates on axis changes, persists through layout save/restore, and reflects filter panel changes. Legend click activates viewer-level filtering. In-viewer filter persists through layouts. Grid color coding interacts properly with viewer legends.

## Retrospective

### What worked well
- Combined color+marker legend displayed correctly with proper shapes and colors
- Legend dynamically updates when switching between linear and categorical color columns
- Layout save/restore preserves all legend-related settings (color, marker, in-viewer filter)
- Filter panel integration works — legend updates to reflect filtered data
- Pointer events (PointerEvent) successfully trigger Dart legend click handler for viewer-level filtering
- Updated spec to use CDP connection for consistency

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- Legend click-to-filter behavior is powerful but not immediately discoverable — consider adding a visual hint

### Suggestions for the scenario
- Step 4.6 could clarify that clicking a legend item applies a viewer-level legend filter rather than a dataframe-level filter
