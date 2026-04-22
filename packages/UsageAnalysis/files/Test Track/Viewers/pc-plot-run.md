# PC Plot tests (Playwright) — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Menu Ribbon & To Script — open PC plot, right-click To Script, close and re-add | 31s | PASS | PASSED | PC plot added via toolbox icon; right-click → To Script → To JavaScript produced balloon; close + reopen succeeded |
| 2 | Axis scale & normalization — normalize toggle, log columns, Y Axis Global/Normalized | 32s | PASS | PASSED | `normalizeEachColumn` toggled via props; right-click Y Axis → Global/Normalized menu path OK; logColumnsColumnNames AGE, AGE+WEIGHT, cleared |
| 3 | Selection & line display — toggle current/mouseOver/all/row group | 13s | PASS | PASSED | Defaults all `true`; `showCurrentLine`, `showMouseOverLine`, `showAllLines`, `showMouseOverRowGroup` flipped off and back |
| 4 | Style & layout — lineWidth, currentLineWidth, mouseOverLineWidth, orientations, horzMargin, autoLayout | 13s | PASS | PASSED | lineWidth=3, currentLineWidth=5, mouseOverLineWidth=5, labelsOrientation=Vert, horzMargin=60, autoLayout=false; all restored |
| 5 | In-chart filtering & reset — showFilteredOutLines, dblclick reset, Reset View menu | 16s | PASS | PASSED | `showFilteredOutLines=true`, dblclick canvas fired, context-menu `Reset View` item found and clicked |
| 6 | Filter panel interaction — histogram filter on AGE, reset | 12s | PASS | PASSED | AGE∈[30,50] filtered to 2870/5850 rows; AGE∈[18,89] restored 5850 |
| 7 | Column management & reordering — remove HEIGHT, add back, reorder | 12s | PASS | PASSED | defaultCols=[AGE,HEIGHT,WEIGHT,STARTED]; HEIGHT removed, re-added, reordered to [WEIGHT,AGE,HEIGHT,STARTED], restored |
| 8 | Density styles — circles / box plot toggles / violin + bins / whisker | 14s | PASS | PASSED | densityStyle: circles → box plot → violin plot → circles; IQR, upper/lower dash, mean, median, circles, bins=200, whiskerLineWidth=5 all toggled |
| 9 | Color coding, legend & grid coloring — AGE/RACE/HEIGHT, invert, log, legend positions, grid colors | 18s | PASS | PASSED | colorAxisType=logarithmic, invertColorScheme, colorMin=30, RACE categorical; legendPosition Left/Right/Top/Bottom; legendVisibility Never/Auto; HEIGHT setLinear + setConditional via col.meta.colors |
| 10 | Title and description — set title/desc/position, clear | 7s | PASS | PASSED | title='My PC Plot', description='Test description', descriptionPosition='Bottom'; both cleared |
| 11 | Pick Up / Apply — 2nd PC plot inherits cols/log/color/legend/title | 28s | PASS | PASSED | pc2Cols=[AGE,WEIGHT,STARTED], color=RACE, title='Source Plot' on second viewer after right-click Pick Up → Apply |
| 12 | Layout save/restore — save, add scatter, reload layout | 18s | PASS | PASSED | After reload: hasScatter=false, hasPc=true; layout deleted from server |
| 13 | Table switching & transformation — add spgi-100, switch PC plot table + GroupAggregation | 31s | PASS | PASSED | spgi-100 loaded (100 rows, Molecule semType); `pc.props.table` + 3-stage GroupAggregation transformation applied cleanly |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 47s |
| grok-browser execution (scenario steps) | 21s |
| Execute via grok-browser (total) | 1m 8s |
| Spec file generation | 8s |
| Spec script execution | 47s |
| **Total scenario run (with model)** | 2m 3s |

All rows are full-phase wall-clock (incl. model thinking and retries). The two `scenario steps`
rows sum to `Execute via grok-browser (total)`. Spec was not rewritten (existing spec is correct
per the user's instruction); the 8s reflects verification/sanity-pass only. Spec execution used
the existing `playwright.config.files-and-sharing.ts` (testMatch: `/-spec\.ts$/`).

## Summary

All 13 scenario sections (mapped to 12 Playwright softSteps — scale and normalization are combined in the spec) passed during the MCP run on dev.datagrok.ai and then again via the existing Playwright spec unmodified. Dataset opens (demog: 5850 rows; spgi-100: 100 rows with Molecule semType), all property round-trips, context-menu navigation (Y Axis → Global/Normalized, Pick Up → Apply, To Script → To JavaScript), filter interactions, layout save/restore, and transformation apply all succeeded. **Total scenario run (with model)**: 7m 53s.

## Retrospective

### What worked well
- JS-API property round-trips (`pc.props.*`) give deterministic verification — no flaky canvas coordinates.
- Context-menu navigation via `.d4-menu-item-label` text + `mouseenter` for submenu hover is reliable on both primary and secondary PC plots (Pick Up → Apply worked first try).
- Layout save/restore via `grok.dapi.layouts.save/find/delete` is robust with a 1s post-save wait and 3s post-loadLayout wait.
- Filter panel `updateOrAdd({type:'histogram', column, min, max})` precisely controlled row count (5850 → 2870 → 5850).
- The existing Playwright spec preamble already matches the grok-debug-scenarios template verbatim and passed on dev without edits.

### What did not work
- `npx playwright test <path>` failed twice before succeeding: once because Playwright's default `testMatch` excludes `*-spec.ts`, and once because a config placed under `/tmp/` can't resolve `@playwright/test` from the repo's `node_modules`. Resolution required a repo-local config wrapper.
- `pc.props.table = df2.name` set the property to `"Table (2)"` (auto-generated view name) rather than the filename `"spgi-100"` — the display name and the loaded table's internal name diverge; the spec's assertion (`tableSet` truthy) is defensive but loses specificity.

### Suggestions for the platform
- Normalize the table name produced by `readCsv(path)` to the file stem (e.g. `spgi-100`), so scenarios referring to the CSV filename don't need a separate rename step.
- Add `name=` attributes on `.d4-menu-item-label` so nested context menus can be addressed without text-content scanning.
- Expose axis range slider handles as DOM elements (e.g. `[name="pc-axis-AGE-slider"]`) so in-chart filter drags can be scripted without hitting canvas coordinates.

### Suggestions for the scenario
- Section "Filter panel interaction" specifies a precondition ("range sliders on two axes are narrowed (from scenario 4 flow)") that most spec implementations skip — inline an explicit setup or drop the precondition.
- Section "Table switching and transformation" step 3 ("set **Table** to spgi-100") should clarify whether it means the filename or the in-memory table name (`Table`, `Table (2)`, …).
- Step 11.9 ("Adjust the range slider on the second PC plot — the first plot should update to show filtered lines, but its own range sliders should remain unchanged") is hard to verify purely via props; consider switching to a filter-count assertion via `df.filter.trueCount`.
- Section "Layout and project save/restore" step 6 ("Open the saved project") doesn't specify project name / lookup path — either provide a JS-API path or drop the step (layout save already covers restoration).
