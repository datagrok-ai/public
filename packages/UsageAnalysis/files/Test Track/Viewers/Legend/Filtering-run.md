# Filtering — Run Results

**Date**: 2026-04-06
**URL**: http://localhost:8888/
**Status**: PARTIAL

## Steps

### Section 1: Filtering

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | 3624 rows, 88 columns loaded |
| 2 | Add viewers and set legend to Stereo Category | PASS | PASSED | 7 viewers added: scatter, histogram, line, bar, pie, trellis, box |
| 3 | Set legend to Stereo Category for each viewer | PASS | PASSED | All viewers show 5-category legend |
| 4 | Open Filter Panel and apply filters — check legend | PASS | PASSED | Categorical (R_ONE, S_ABS) + numerical (CAST Idea ID range) → 352 rows; legends show only R_ONE, S_ABS |
| 5 | Save and apply layout — check legend | PASS | PASSED | Layout restored with 352 filtered rows, all viewers and legends intact |
| 6 | Reset filters | PASS | PASSED | All 3624 rows restored after closing/reopening filter panel |
| 7 | Set in-viewer Filter to R_ONE/S_UNKN — check legend | PARTIAL | PASSED | Histogram, line chart, bar chart, box plot correctly show R_ONE/S_UNKN only; scatter plot, pie chart, trellis plot still show all 5 categories |
| 8 | Apply additional Filter Panel filters with in-viewer filter | PASS | PASSED | Combined filtering: panel (R_ONE, S_ABS, S_UNKN + range) + in-viewer (R_ONE, S_UNKN) → 630 rows; viewers show intersection correctly |
| 9 | Filter via viewers (categorical filter) — legend shows only visible values | PASS | PASSED | Filtered to R_ONE, S_PART → all legends show only R_ONE, S_PART |
| 10 | Save and apply layout — check legend | PASS | PASSED | Layout restored with 1832 rows, legends preserved |
| 11 | Set different Row Source values — check legend | PASS | PASSED | Scatter plot (Selected, 100 rows): legend shows R_ONE, S_ABS; Histogram (FilteredSelected): correct; others unchanged |

### Section 2: Bar chart edge case

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Configure bar chart (Value=CAST Idea ID, Category=Stereo Category, Stack=Primary scaffold name) | PASS | PASSED | Bar chart shows 5 Stereo Category bars stacked by scaffold |
| 2 | Uncheck Include nulls | PASS | PASSED | No null bar visible |
| 3 | Filter Primary scaffold name — legend shows all displayed categories | PASS | PASSED | Deselected UNSUBSTDIAZA, HYDROXY_AMINOPIPERIDINE; legend correctly shows only AMINOPYRROLIDINE, EXTENDED_AMINOPYRROLIDINE, TRISUBSTITUTED |

## Summary

Most filtering scenarios work correctly. Legends update properly when filters are applied via the Filter Panel, and layout save/restore preserves both filter state and legend configuration. The in-viewer `filter` property works for histogram, line chart, bar chart, and box plot but does not affect scatter plot, pie chart, or trellis plot legends. The bar chart edge case with stacked bars and filtered categories behaves correctly.

## Retrospective

### What worked well
- Legend updates immediately when filters change across all viewer types
- Layout save/restore reliably preserves filter state and legend configuration
- Bar chart stacked legend correctly reflects only visible scaffold categories after filtering
- Row Source changes correctly affect individual viewer legends

### What did not work
- In-viewer `filter` property has no effect on scatter plot, pie chart, and trellis plot — the property may not be supported for these viewer types, or the legend does not respond to it

### Suggestions for the platform
- The in-viewer `filter` property should be supported consistently across all viewer types, or documentation should clarify which viewers support it
- Consider adding a `resetAll()` method to the filter group API for easier programmatic filter reset

### Suggestions for the scenario
- Step 7 ("Set in-viewer Filter") should specify which viewers support the `filter` property, as behavior differs
- Step 9 ("Filter table via viewers") could clarify how to trigger zoom-based filtering on scatter plot programmatically
- Steps could be more specific about expected legend content (e.g., "legend should show exactly R_ONE and S_UNKN")
