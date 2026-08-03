# Grid: include visible columns in tooltip — Run Results

**Date**: 2026-05-12
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open linked datasets (energy_uk.csv) | 8s | PASS | PASSED | 68 rows, 3 cols (value double, source string, target string) |
| 2 | Open grid properties and find Show Visible Columns In Tooltip | 6s | PASS | PASSED | Settings gear opened the property panel; `[name="prop-show-visible-columns-in-tooltip"]` row found in Tooltip category, default `false` |
| 3 | Unchecked + all cols visible: no tooltip on cell hover | 6s | PASS | PASSED | Tooltip element stays `display:none` with empty text — matches scenario |
| 4 | Unchecked + last column clipped: tooltip shows clipped column | 7s | PASS | PASSED | After widening `value` to 1100px so `target` falls off-screen, hover renders `.d4-tooltip` with text `target … Losses` — and `source` (still visible) is NOT in the tooltip. Required setting `showTooltip = 'inherit from table'` first; see "What did not work" below |
| 5 | Enable Show Visible Columns In Tooltip | 4s | PASS | PASSED | Clicking the checkbox toggled `grid.props.showVisibleColumnsInTooltip` to `true` |
| 6 | Checked: tooltip lists all visible columns in both states | 10s | PASS | PASSED | All three columns (`value`, `source`, `target`) appear in the tooltip whether the grid is fully visible or `target` is clipped — content identical across the two states |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~14m |
| grok-browser execution (scenario steps) | ~45s |
| Execute via grok-browser (total) | ~15m |
| Spec file generation | ~5m |
| Spec script execution | 22s |
| **Total scenario run (with model)** | ~21m |

## Summary

All six scenario steps PASS. The grid's `Show Visible Columns In Tooltip` property behaves exactly as documented: when unchecked the cell tooltip lists only columns that are clipped off-screen (or nothing if all columns fit); when checked it always lists every visible column, and the listed set is identical regardless of whether the last column is clipped. The diagnostic loop uncovered that the grid's default `Show Tooltip = 'show custom tooltip'` short-circuits the entire `showVisibleColumnsInTooltip` filter when `Row Tooltip` is empty — the scenario implicitly assumes `Show Tooltip = 'inherit from table'`, which the spec now sets explicitly. Total scenario run (with model): ~21m.

## Retrospective

### What worked well
- The accordion-like property-grid table (`[name="prop-category-tooltip"]` + `[name="prop-show-visible-columns-in-tooltip"]`) is reliable once the category is expanded; `state:'attached'` waits avoid the `display:table-row` hidden-row trap.
- `page.mouse.move(0,0)` → `page.mouse.move(box.x - 50, …)` → `page.mouse.move(box.x + box.width/2, …)` is enough movement to fire `mouseLeave` → `mouseEnter` and clear the 200 ms tooltip debounce.
- Reading `core/client/d4/lib/src/viewers/grid/grid_core.dart` (`showCellTooltip`) and `core/client/d4/lib/src/widgets/tooltip/tooltip.dart` (`getRowTooltip`, lines 525-537) clarified that the cell-tooltip filter (which is what `showVisibleColumnsInTooltip` flows into) is only applied under specific `Show Tooltip` modes — without that read the diagnosis loop would have stayed AMBIGUOUS.
- Using `ui.tooltip.show(...)` from JS to confirm `.d4-tooltip` is the right element to watch ruled out "wrong selector" as a cause.

### What did not work
- The grid's default `Show Tooltip = 'show custom tooltip'` plus default empty `Row Tooltip` causes `getRowTooltip` (tooltip.dart:525-537) to return `null` from `tooltipColumns.length == 0 && preText == null`. The cell-tooltip filter that `showCellTooltip` builds is never consulted, so `showVisibleColumnsInTooltip` has no observable effect in this mode. This is what made the first spec run land on AMBIGUOUS for steps 4 and 6 — the spec was honoring the user-facing default but the platform behavior never triggers there.
- `grok test` without the local `public/tools/` build fails with `UsageAnalysis:test returned null instead of a DataFrame` because the globally installed `datagrok-tools@6.1.10` predates `--skip-puppeteer`. Worked around by `npm install` + `npm run build` in `public/tools/` and running `node public/tools/bin/grok.js test …` directly.

### Suggestions for the platform
- Either honor `showVisibleColumnsInTooltip` in `'show custom tooltip'` mode (when `rowTooltip` is empty, apply the same column-filter logic instead of the early `return null`), or update the property's UI hint/dependsOn to make it clear the flag only affects `'inherit from table'`. Today the property silently does nothing under the default grid configuration, which is exactly the case end users start in.
- Update the globally-distributed `datagrok-tools` to a 6.2+ build so `--skip-puppeteer` is available without a manual local build dance.
- Consider adding a JS-API entry point such as `grid.showCellTooltipAt(row, col)` so test code can verify tooltip content without relying on real pointer motion — useful for headless CI even when CDP-driven hover is fine.

### Suggestions for the scenario
- Either set `Show Tooltip = 'inherit from table'` explicitly as a precondition, or note that the default mode currently suppresses the behavior under test (so future readers don't mistake the "no tooltip" outcome for a regression).
- Step 4 wording ("extend a column's width to push the last column(s) out of sight") is ambiguous about whether partial clipping counts. Reword as "until the last column is fully off-screen and a horizontal scrollbar appears."
- Step 6 says the tooltip "remains the same" in both states — clarify what "the same" means (same column set, same row values). The spec now asserts that all three column names are present in both states.
- Add an explicit pre-condition that the **Tooltip** category in the property panel may be collapsed by default; instruct the tester to click the header to expand before searching for the property.
