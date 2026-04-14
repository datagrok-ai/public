# Histogram tests (Playwright) -- Run Results

**Date**: 2026-04-08
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-------|-------|
| **Bins configuration** | | | | |
| 1 | Set Value to AGE | PASS | PASSED | |
| 3 | Set Bins to 5 | PASS | PASSED | |
| 4 | Set Bins to 100 | PASS | PASSED | |
| 5 | Set Bins to 1 | PASS | PASSED | |
| 6 | Reset Bins to 20 | PASS | PASSED | |
| 7 | Bin Width Ratio = 1.0 | PASS | PASSED | |
| 8 | Bin Width Ratio = 0.3 | PASS | PASSED | |
| 9 | Bin Width Ratio = 0.8 | PASS | PASSED | |
| **Split column** | | | | |
| 1 | Split = SEX | PASS | PASSED | |
| 2 | Normalize Values = true | PASS | PASSED | |
| 3 | Normalize Values = false | PASS | PASSED | |
| 4 | Show Markers = false | PASS | PASSED | |
| 5 | Spline Tension = 5 | PASS | PASSED | |
| 6 | Split = RACE | PASS | PASSED | |
| 7 | Split = None | PASS | PASSED | |
| 8 | Split = SEX | PASS | PASSED | |
| 9 | Split Stack = true | PASS | PASSED | |
| 10 | Show Values = true | PASS | PASSED | |
| 11 | Split Stack = false | PASS | PASSED | |
| 12 | Show Distribution Lines = true | PASS | PASSED | |
| 13 | Show Distribution Lines = false | PASS | PASSED | |
| **Color coding** | | | | |
| 1 | Split = None | PASS | PASSED | |
| 2 | Color Column = WEIGHT | PASS | PASSED | |
| 3 | Color Aggr Type = min | PASS | PASSED | |
| 4 | Color Aggr Type = max | PASS | PASSED | |
| 5 | Invert Color Scheme = true | PASS | PASSED | |
| 6 | Invert Color Scheme = false | PASS | PASSED | |
| 7 | Color Column = None | PASS | PASSED | |
| **Value range** | | | | |
| 1 | Set Value to AGE | PASS | PASSED | |
| 2 | Value Min = 30 | PASS | PASSED | |
| 3 | Value Max = 60 | PASS | PASSED | |
| 4 | Clear Min and Max | PASS | PASSED | |
| 5 | Show Range Inputs = true | PASS | PASSED | |
| 6 | Type in range inputs | AMBIGUOUS | SKIPPED | Inputs appear on hover |
| 7 | Show Range Inputs = false | PASS | PASSED | |
| **Spline mode** | | | | |
| 1 | Spline = true | PASS | PASSED | |
| 2 | Fill Spline = true | PASS | PASSED | |
| 3 | Fill Spline = false | PASS | PASSED | |
| 4 | Spline = false | PASS | PASSED | |
| **Appearance** | | | | |
| 1 | Show X Axis = true | PASS | PASSED | |
| 2 | Show Y Axis = true | PASS | PASSED | |
| 3 | Disable both | PASS | PASSED | |
| 4 | X Axis Height = 30 | PASS | PASSED | |
| 5 | Show Column Selector = false | PASS | PASSED | |
| 6 | Show Bin Selector = false | PASS | PASSED | |
| 7 | Show Split Selector = false | PASS | PASSED | |
| 8 | Show Range Slider = false | PASS | PASSED | |
| 9 | Re-enable all four | PASS | PASSED | |
| **Labels** | | | | |
| 1 | Split = SEX | PASS | PASSED | |
| 2 | Legend Visibility = Never | PASS | PASSED | |
| 3 | Legend Visibility = Always | PASS | PASSED | |
| 4 | Legend Position = RightTop | PASS | PASSED | |
| 5 | Split = None | PASS | PASSED | |
| 6 | Show Title = true | PASS | PASSED | |
| 7 | Title = "Age Distribution" | PASS | PASSED | |
| 8 | Description = "Shows distribution..." | PASS | PASSED | |
| 9 | Description Visibility = Always | PASS | PASSED | |
| 10 | Description Position = Bottom | PASS | PASSED | |
| **Bin selection** | | | | |
| 2-6 | Click bins, Ctrl+click, Shift+drag | AMBIGUOUS | SKIPPED | Canvas-based |
| 7 | Split = SEX | PASS | PASSED | |
| 8 | Split Stack = true | PASS | PASSED | |
| 9 | Click stacked segment | AMBIGUOUS | SKIPPED | Canvas-based |
| 10 | Split Stack = false | PASS | PASSED | |
| 11 | Click spline line | AMBIGUOUS | SKIPPED | Canvas-based |
| **Filtering** | | | | |
| 1 | Filter SEX = M | PASS | PASSED | 2607 rows |
| 2 | Show Filtered Out Rows = false | PASS | PASSED | |
| 3 | Show Filtered Out Rows = true | PASS | PASSED | |
| 4 | Change Filtered Out Color | PASS | PASSED | |
| 5 | Remove filter | PASS | PASSED | 5850 rows |
| 6 | Set Value to AGE | PASS | PASSED | |
| 7 | Filtering Enabled = true | PASS | PASSED | |
| 8 | Drag range slider | AMBIGUOUS | SKIPPED | Canvas interaction |
| 9 | Normalize To Filter = false | PASS | PASSED | |
| 10 | Normalize To Filter = true | PASS | PASSED | |
| 11 | Zoom To Range = false | PASS | PASSED | |
| 12 | Zoom To Range = true | PASS | PASSED | |
| 13 | Bin To Range = true | PASS | PASSED | |
| 14 | Filtering Enabled = false | PASS | PASSED | 5850 rows restored |
| 15 | Filter RACE = Asian+Other | PASS | PASSED | 426 rows |
| 16 | Reset RACE filter | PASS | PASSED | 5850 rows |
| **Context menu** | | | | |
| 1-4 | Histogram context menu | AMBIGUOUS | SKIPPED | Canvas right-click shows view menu |
| **Layout persistence** | | | | |
| 1 | Set Value = WEIGHT | PASS | PASSED | |
| 2 | Set Bins = 15 | PASS | PASSED | |
| 3 | Set Split = RACE | PASS | PASSED | |
| 4 | Enable Split Stack | PASS | PASSED | |
| 5 | Save layout | PASS | PASSED | |
| 6 | Close histogram | PASS | PASSED | |
| 7-8 | Apply and verify layout | PASS | PASSED | All settings restored |
| **Data properties** | | | | |
| 1 | Select 5 rows | PASS | PASSED | |
| 2 | Row Source = Selected | PASS | PASSED | |
| 3 | Row Source = All | PASS | PASSED | |
| 4 | Filter formula = ${AGE} > 40 | PASS | PASSED | |
| 5 | Clear filter formula | PASS | PASSED | |
| 6-9 | Table switching | SKIP | SKIPPED | Complex table prop |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~120s |
| Spec file generation | ~30s |
| Spec script execution | 50s (PASSED) |

## Summary

Most histogram property-based tests passed successfully. All property setters (bins, split, color,
spline, appearance, labels, filtering, layout) worked correctly via JS API. Canvas-based interactions
(bin clicks, range slider drag, context menu on histogram area) could not be automated. Layout
save/restore preserved all histogram settings correctly.

## Retrospective

### What worked well
- All histogram properties accessible and modifiable via JS API
- Filter integration works correctly (categorical + histogram filtering)
- Layout persistence preserves all histogram settings including split, stack, bins
- Row Source switching works for Selected/All views
- Viewer filter formula works correctly

### What did not work
- Histogram context menu (Show Filtered Out Rows, Selection group) not reachable via DOM right-click dispatch -- always shows the view-level menu instead
- Canvas-based bin selection (click, Ctrl+click, Shift+drag) not automatable
- Range slider drag interaction not automatable

### Suggestions for the platform
- Consider making histogram's own context menu accessible via right-click on the viewer container, not just on specific canvas regions
- Add `name=` attributes to histogram bin elements if feasible

### Suggestions for the scenario
- Separate canvas-dependent bin selection steps from property-based tests
- Add expected row counts for bin selection verification
- Clarify which context menu items should appear at which click position
