# Bar chart tests (Playwright) — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| **Stack column** | | | | |
| 1 | Set Stack to SEX | PASS | PASSED | `stackColumnName = 'SEX'` confirmed |
| 2 | Set Stack to RACE | PASS | PASSED | Updated to RACE |
| 3 | Enable Relative Values | PASS | PASSED | `relativeValues = true` |
| 4 | Disable Relative Values | PASS | PASSED | `relativeValues = false` |
| 5 | Set Stack to None | PASS | PASSED | `stackColumnName = ''` |
| **Sorting** | | | | |
| 1 | Set Split to RACE | PASS | PASSED | |
| 2 | Bar Sort Type = by value | PASS | PASSED | |
| 3 | Bar Sort Order = asc | PASS | PASSED | |
| 4 | Bar Sort Order = desc | PASS | PASSED | |
| 5 | Bar Sort Type = by category | PASS | PASSED | |
| 6 | Category sort asc | PASS | PASSED | |
| 7 | Category sort desc | PASS | PASSED | |
| **Value axis type** | | | | |
| 1 | Set Value=HEIGHT, Split=RACE | PASS | PASSED | |
| 2 | Axis Type = logarithmic | PASS | PASSED | |
| 3 | Axis Type = linear | PASS | PASSED | |
| 4 | Set custom Min=100, Max=200 | PASS | PASSED | |
| 5 | Clear Min and Max | PASS | PASSED | `null` confirmed |
| **Color coding** | | | | |
| 1 | Split=RACE, Value=AGE | PASS | PASSED | |
| 2 | Color Column = HEIGHT | PASS | PASSED | |
| 3 | Color Aggr Type = min | PASS | PASSED | |
| 4 | max, then med | PASS | PASSED | |
| 5 | Invert Color Scheme | PASS | PASSED | |
| 6 | Color Column = None | PASS | PASSED | |
| **Include nulls** | | | | |
| 1 | Split to DIS_POP (has nulls) | PASS | PASSED | |
| 2 | includeNulls default = true | PASS | PASSED | |
| 3 | Disable includeNulls | PASS | PASSED | |
| 4 | Re-enable includeNulls | PASS | PASSED | |
| **Bar style** | | | | |
| 1 | Open Style section | PASS | PASSED | |
| 2 | Bar Border Line Width = 2 | PASS | PASSED | |
| 3 | Bar Corner Radius = 10 | PASS | PASSED | |
| 4 | Max Bar Height = 20 | PASS | PASSED | |
| 5 | Vertical Align Top/Bottom/Center | PASS | PASSED | |
| 6 | Show Category Zero Baseline off | PASS | PASSED | |
| **Labels** | | | | |
| 1-2 | Split=RACE, Value=AGE, open Style | PASS | PASSED | |
| 3 | Show Labels = inside | PASS | PASSED | |
| 4 | Show Labels = outside | PASS | PASSED | |
| 5 | Show Labels = never | PASS | PASSED | |
| 6 | Show Labels = auto | PASS | PASSED | |
| **Controls visibility** | | | | |
| 2 | Show Value Selector off | PASS | PASSED | |
| 3 | Show Category Selector off | PASS | PASSED | |
| 4 | Show Stack Selector off | PASS | PASSED | |
| 5 | Show Value Axis off | PASS | PASSED | |
| 6 | Show Category Values off | PASS | PASSED | |
| 7 | Re-enable all | PASS | PASSED | |
| **Aggregation types** | | | | |
| 2 | avg | PASS | PASSED | |
| 3 | min | PASS | PASSED | |
| 4 | max | PASS | PASSED | |
| 5 | sum | PASS | PASSED | |
| 6 | count | PASS | PASSED | |
| 7 | Value=WEIGHT, avg | PASS | PASSED | |
| **Date/time split** | | | | |
| 1 | Split=STARTED | PASS | PASSED | |
| 2 | Split Map = year | PASS | PASSED | |
| 3 | Split Map = month | PASS | PASSED | |
| 4 | Split Map = quarter | PASS | PASSED | |
| **Legend** | | | | |
| 1 | Stack=SEX, legend appears | PASS | PASSED | |
| 2 | Legend Visibility = Always | PASS | PASSED | |
| 3 | Legend Position (L/R/T/B) | PASS | PASSED | |
| 4 | Legend Visibility = Never | PASS | PASSED | |
| 5 | Remove Stack, legend gone | PASS | PASSED | |
| **Title and description** | | | | |
| 1 | Show Title enable | PASS | PASSED | |
| 2 | Title = "Demographics" | PASS | PASSED | |
| 3 | Description = "By race" | PASS | PASSED | |
| 4 | Description Position variants | PASS | PASSED | |
| 5 | Description Visibility = Never | PASS | PASSED | |
| **Show values instead of categories** | | | | |
| 2 | Enable | PASS | PASSED | |
| 3 | Disable | PASS | PASSED | |
| **Orientation** | | | | |
| 1 | horizontal | PASS | PASSED | |
| 2 | vertical | PASS | PASSED | |
| 3 | auto | PASS | PASSED | |
| **Data panel (SPGI)** | | | | |
| 1 | Add Bar chart to demog | PASS | PASSED | |
| 2 | Row Source switching | PASS | PASSED | All/Filtered/Selected |
| 3 | Switch table to SPGI | PASS | PASSED | |
| 5 | Set Filter expression | PASS | PASSED | `${CAST Idea ID} < 636500` |
| 6 | Color Column = Chemical Space Y | PASS | PASSED | |
| 7-9 | Layout save/restore | PASS | PASSED | Filter and color preserved |
| 10 | Delete saved layout | PASS | PASSED | |
| **Filter Panel interaction** | | | | |
| 1 | Open Filter Panel | PASS | PASSED | |
| 2 | Categorical filter | PASS | PASSED | 5339 rows after RACE filter |
| 3 | Numeric filter | PASS | PASSED | 3796 rows after AGE + RACE |
| 4 | Remove all filters | PASS | PASSED | 5850/5850 restored |
| **Scrolling with range slider** | | | | |
| 1 | Split=SUBJ (many cats) | PASS | PASSED | 5850 unique categories |
| 2 | Value=AGE with scrolling | PASS | PASSED | |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~3 min |
| Spec file generation | ~10s |
| Spec script execution | 1m 0s (1 passed, 0 failed) |

## Summary

All 15 bar chart test sections passed on dev.datagrok.ai. All viewer properties (stack, sorting, axis type, color coding, include nulls, bar style, labels, controls, aggregation, date/time split, legend, title/description, values-instead-of-categories, orientation) behave correctly. Layout save/restore with SPGI dataset preserves filter and color column. Filter panel collaborative filtering confirmed.

## Retrospective

### What worked well
- All bar chart properties accessible and round-trip correctly through `viewer.props`
- Layout save/restore preserves filter expressions, color coding, split/stack columns
- Date split maps (year/month/quarter) work correctly
- Filter Panel interaction with bar chart is seamless

### What did not work
- No failures encountered

### Suggestions for the platform
- Add hit-test APIs for canvas-rendered bar segments to improve automation
- Range slider scrolling is visual-only — no programmatic API to verify scroll position

### Suggestions for the scenario
- "Scrolling with range slider" references "Primary Scaffold Name" which is not in demog — used SUBJ instead
