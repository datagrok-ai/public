# Box plot tests (Playwright) — Run Results

**Date**: 2026-04-10
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Plot style: box vs violin — set Value=AGE, Cat1=RACE, switch violin, bins 50/500, IQR width 10, back to box | PASS | All property changes confirmed |
| 2 | Two-level categories — Cat1=SEX, Cat2=RACE, toggle minor/all categories, clear Cat2 | PASS | All toggles work correctly |
| 3 | Statistics display — verify defaults, enable all stats, disable showStatistics | PASS | showStatistics=true by default, all individual stats toggleable |
| 4 | Markers — markerType=square, size=10, opacity=80, markersCol=RACE, sizeCol=WEIGHT, clear both | PASS | All marker properties set and cleared |
| 5 | Marker and Bin color coding — markerColor AGE/RACE, colorAxisType log, invertColorScheme, binColor WEIGHT, aggrType min/max/med | PASS | All color properties work |
| 6 | Value axis configuration — axisType log, invertYAxis, valueMin/Max 20/60, clear, reset to linear | PASS | All axis properties confirmed |
| 7 | Zoom by filter with filter panel — zoomValuesByFilter default/toggle, AGE filter 40-70, reset | PASS | Filter reduced to 3697 rows, reset to 5850 |
| 8 | Show empty categories — toggle off/on | PASS | Property toggles correctly |
| 9 | Box plot components — toggle off all 6 (mean/median/upper/lower/inside/outside), re-enable all | PASS | All six toggles work |
| 10 | Controls visibility — disable all 6 selectors/axes, re-enable all | PASS | All six controls toggleable |
| 11 | Title and description — showTitle, title text, description text, visibility mode Always/Never, position Bottom | PASS | All description properties work |
| 12 | Date category mapping — Cat1=STARTED, map=month/quarter, back to RACE | PASS | Date categorization works |
| 13 | Style customization — whiskerLineWidth=4, whiskerWidthRatio 1.0/0.3, autoLayout off/on, axisUseColumnFormat off/on | PASS | All style properties confirmed |
| 14 | Viewer filter formula — set `${AGE} > 40`, clear | PASS | Filter formula set and cleared |
| 15 | P-value (t-test) — Cat1=SEX, toggle showPValue off/on, T key toggle off/on, Cat1=RACE for Alexander-Govern | PASS | T key correctly toggles p-value, all modes work |
| 16 | Legend — markerColor=RACE, legendVisibility Never/Always, legendPosition RightTop/LeftBottom, clear color | PASS | Legend properties work correctly |
| 17 | Layout save and restore — configure (WEIGHT/RACE/SEX/totalCount/violin), save, reopen, load, verify, delete | PASS | All properties preserved after layout load |
| 18 | Visualization zoom (project and layout) — zoom via valueMin/Max, reset, project save, layout save/load | AMBIGUOUS | Zoom and reset work. Project save fails with "Unable to add entity" error on dev server. Layout DOES preserve zoom (valueMin=30, valueMax=60 after restore), but scenario expects zoom NOT preserved in layout — potential behavior discrepancy |
| 19 | Data properties and table switching (spgi-100) — open both, switch table, set Value=Average Mass, Cat1=Series, filter, binColor=TPSA | PASS | Table switching and all properties work on spgi-100 |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~180s |
| Spec file generation | existing spec reused |
| Spec script execution | 47s (PASSED) |

## Summary

17 of 19 sections passed. Section 18 (Visualization zoom) is AMBIGUOUS: project save fails on dev server, and layout restore preserves zoom state contrary to scenario expectation. All box plot property manipulations, filter interactions, T-key toggle, layout save/restore, and table switching work correctly.

## Retrospective

### What worked well
- All box plot properties are fully accessible and modifiable via JS API
- Layout save/restore correctly preserves all viewer configuration
- T key shortcut for p-value toggle works reliably
- Table switching between demog and spgi-100 works seamlessly

### What did not work
- Project save via `grok.dapi.projects.save()` fails with "Unable to add entity" error — likely a server-side or permissions issue on dev
- Layout zoom preservation behavior differs from scenario expectation (layout DOES preserve valueMin/valueMax)

### Suggestions for the platform
- Investigate project save API error on dev server
- Clarify whether layout should preserve viewport zoom (valueMin/valueMax) — current behavior preserves it

### Suggestions for the scenario
- Section 18: Clarify expected behavior for zoom preservation in layouts vs projects
- Section 18: Consider splitting zoom/project/layout tests into separate sections for easier debugging
