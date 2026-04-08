# Legend Filtering — Run Results

**Date**: 2026-04-08
**URL**: http://localhost:8888/
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1.1 | Open SPGI | PASS | 8s | PASSED | 3624 rows, 88 columns |
| 1.2 | Add viewers | PASS | 7s | PASSED | 8 viewers total (Grid + 7 added via toolbox) |
| 1.3 | Set legend to Stereo Category | PASS | 3s | PASSED | Color/Split/Stack set on all viewers |
| 1.4 | Open Filter Panel, apply filter, check legend | PASS | 4s | PASSED | Filtered to R_ONE+S_ACHIR (2247 rows); legends show only 2 categories |
| 1.5 | Save and apply layout -- check legend | PASS | 5s | PASSED | Layout restored; programmatic filter reset as expected; legend shows all 5 |
| 1.6 | Reset filters | PASS | 1s | PASSED | All 3624 rows visible |
| 1.7 | Set in-viewer Filter to R_ONE, S_UNKN | PASS | 2s | PASSED | Scatter Plot legend shows R_ONE+S_UNKN; dataframe still 3624 rows |
| 1.8 | Apply additional Filter Panel filters | PASS | 3s | PASSED | SP legend: R_ONE+S_UNKN (intersection); Histogram: R_ONE+S_ABS+S_UNKN |
| 1.9 | Filter via viewers (legend click) | PASS | 2s | PASSED | Legend click: R_ONE current, others not; dataframe unchanged |
| 1.10 | Save and apply layout -- check legend | PASS | 5s | PASSED | 9 viewers (incl Filters) restored; legend shows all 5 categories |
| 1.11 | Set different Row Source values | PASS | 3s | PASSED | All: 5 cats; Filtered: 2; Selected: 1 |
| 2.1 | Bar chart: set Value, Category, Stack | PASS | 3s | PASSED | Value=CAST Idea ID, Split=Stereo Category, Stack=Primary scaffold name |
| 2.2 | Uncheck Include nulls | PASS | 1s | PASSED | includeNulls=false |
| 2.3 | Filter scaffold names, check bar chart legend | PASS | 3s | PASSED | 3 filtered scaffold categories shown in bar chart legend |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~55s |
| Spec file generation | ~5s |
| Spec script execution | 35.9s |

## Summary

All legend filtering behaviors work correctly. Legends dynamically update when dataframe filters, in-viewer filters, or row source changes are applied. The intersection of in-viewer and dataframe filters is correctly reflected in legends. Different Row Source values (All, Filtered, Selected) show appropriate legend categories. The bar chart edge case with stacked categories and Include nulls=false correctly shows only filtered scaffold categories in the legend.

## Retrospective

### What worked well
- Legend dynamically updates to reflect filtered categories across all viewer types
- In-viewer filter and dataframe filter interact correctly (intersection shown in legend)
- Row Source property correctly controls which categories appear in legend
- Bar chart stacked legend correctly reflects filtered scaffold categories

### What did not work
- Programmatic filter (`df.filter.init()`) is not preserved through layout save/restore -- expected behavior but worth noting
- `legendVisibility: 'Auto'` hides legends when viewers are small (Filter Panel open); had to set `'Always'` to verify legend content
- `FilterGroup.filters[i].column` property is undefined; must use `filterColumnName` accessor

### Suggestions for the platform
- Consider preserving Filter Panel state in layout save/restore so filtered legends are preserved
- Document `filterColumnName` as the correct accessor for filter column name (not `.column.name`)

### Suggestions for the scenario
- Step 4 "apply filters (structure, numerical, categorical)" is vague -- specify which filters to apply
- Step 9 "zoom for scatterplot" is not clearly testable via automation (zoom is a canvas interaction)
- Step 11 "different Row Source values" -- specify which values to test (All, Filtered, Selected, etc.)
- Consider separating Section 1 and Section 2 into separate scenario files for clarity
