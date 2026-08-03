# Viewers: uniform default tooltip — Run Results

**Date**: 2026-05-12
**URL**: https://dev.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open a scatter plot and a box plot | 25s | PASS | PASSED | Clicked `[name="icon-scatter-plot"]` then `[name="icon-box-plot"]` in the Toolbox. `grok.shell.tv.viewers` reports `["Grid","Scatter plot","Box plot"]`. SP auto-picked `x=WEIGHT, y=HEIGHT`, BP auto-picked `value=HEIGHT, cat1=DIS_POP` — both viewers `showTooltip='inherit from table'`, `rowTooltip=''` |
| 2 | Hover viewers (incl. grid): same column set in same order | 90s | FAIL | FAILED | All three tooltips have the **same set** of 7 columns (`USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT`) but **different order**: Scatter plot prepends Y/X → `HEIGHT, WEIGHT, USUBJID, AGE, SEX, RACE, DIS_POP`; Box plot keeps natural order → `USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT`; Grid (default `showVisibleColumnsInTooltip=false`) tooltips only the columns that get clipped off-screen — content differs depending on layout |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~11m |
| grok-browser execution (scenario steps) | ~45s |
| Execute via grok-browser (total) | ~12m |
| Spec file generation | ~3m |
| Spec script execution | 24s |
| **Total scenario run (with model)** | ~20m |

## Summary

Step 1 PASS, Step 2 FAIL. Scatter plot and Box plot both populate their tooltip from
`tooltip.dart#getRowTooltip` with the same 7 natural columns, but Scatter plot's
`dataValues='Merge'` causes its data columns (Y, X) to be prepended via `getExpandedColumns`,
while Box plot has no data-column merge — yielding a different column order across viewers.
The Grid uses an entirely different mechanism (`showVisibleColumnsInTooltip`, default `false`)
that lists only off-screen columns, so its tooltip content is layout-dependent and rarely
matches either chart's tooltip. The scenario's "same set in the same order" expectation is
not met by the platform's default tooltip behavior. Total scenario run (with model): ~20m.

## Retrospective

### What worked well
- Subscribing to `viewer.onTooltipCreated` from JS gave a reliable snapshot of each viewer's
  tooltip content without depending on real pointer behavior.
- Reading `core/client/d4/lib/src/widgets/tooltip/tooltip.dart#getRowTooltip` (lines 480-570)
  identified the exact divergence point: `getExpandedColumns` prepends `getDataColumns()` for
  viewers that opt in via `dataColumnsFilter`/`dataValuesMergeOption`, but only Scatter Plot
  uses `Merge` by default; Box Plot leaves the natural column order intact.
- The Playwright spec correctly fails at `expect(bpTip.cols).toEqual(spTip.cols)`, which is
  the actual scenario assertion — confirming the codified test would catch a future regression
  toward (or fix of) uniform ordering.
- Using `page.mouse.move()` with a deliberate "move away → settle → move in → wait 1.1s past
  the tooltip debounce" cadence is the only reliable way to trigger row tooltips on canvas
  viewers from Playwright; pure `dispatchEvent('mousemove')` chains keep resetting the
  debounce and never fire.

### What did not work
- Initial spec attempts (3 of 4 runs) used synthesized `dispatchEvent('mousemove')` to drive
  scatter plot tooltips. The tooltip's 200 ms debounce never elapsed because subsequent
  dispatches reset it. The fix was to switch to `page.mouse.move()` with explicit `waitForTimeout`
  longer than the debounce.
- The scatter plot listens on the **data canvas** (first `<canvas>` child); the box plot
  listens on the **overlay canvas** (second). My first spec uniformly used the overlay for
  both, which silently never fired for the scatter plot. Verified against
  `core/client/d4/lib/src/viewers/box_plot/box_plot_core.dart` (`overlay.onMouseMove.listen`)
  and the scatter plot reference (`scatterplot.md` "Dual canvas layers" note).
- Grid hover in Playwright produced a tooltip *only* because the side-by-side layout pushed
  trailing columns off-screen — `showVisibleColumnsInTooltip=false` then surfaces them.
  In the MCP run with a different grid width, the grid produced no tooltip at all. Both
  outcomes already fail the scenario's "same set" expectation, just for different reasons.

### Suggestions for the platform
- Either (a) make Box Plot honor the same `dataValues='Merge'` semantics by exposing the
  value column (and category columns, if non-null) so they prepend the natural tooltip, or
  (b) standardize all chart viewers on the natural-order tooltip so axis-merge is opt-in
  rather than only-Scatter-Plot. Today every other XY-style viewer ignores `dataValues`,
  silently making the Scatter Plot the odd one out for users comparing tooltips.
- Either make the Grid default `showVisibleColumnsInTooltip=true` (so its tooltip lists the
  same columns as the charts) or expose a `dataframe.tooltipColumns` configuration that all
  viewers (Grid included) share. Today the Grid's "default tooltip" is fundamentally a
  different feature from the chart row tooltip — they happen to overlap only when columns
  clip off-screen.
- `tooltip.dart#getRowTooltip` is the single decision point for tooltip composition. A small
  Dart unit test verifying "ordered column list is identical across SP/BP/Line/Bar for a
  table with no explicit `.row-tooltip` tag" would catch regressions automatically.
- Add a JS API such as `viewer.getRowTooltipColumns()` (returns `Column[]`) — the Scatter
  Plot already exposes `getRowTooltip(row)` but the Box Plot does not. Parity would let
  test code verify tooltip contents without driving real pointer events.

### Suggestions for the scenario
- Specify the dataset to use (e.g., `System:DemoFiles/demog.csv`). Without one the scenario
  silently depends on whichever dataset the tester happens to have open, and "same set of
  columns" only makes sense relative to a specific schema.
- Clarify what "same set of columns in the same order" means precisely: is column order
  intended to follow the table's natural order, the viewer's data-column order (Y, X first),
  or `.tooltip` / `.row-tooltip` tag order? The current wording is satisfied by three
  mutually exclusive interpretations.
- Explicitly include the **Grid** rule: state whether the grid is expected to show all
  columns or only off-screen ones, since the grid uses `showVisibleColumnsInTooltip` rather
  than the chart-style tooltip pipeline.
- Note that if `.row-tooltip` is set on the dataframe, scatter plot's `dataValues='Merge'`
  still prepends Y, X — testers should reset the dataframe tooltip configuration before
  running this scenario to get a clean comparison.
