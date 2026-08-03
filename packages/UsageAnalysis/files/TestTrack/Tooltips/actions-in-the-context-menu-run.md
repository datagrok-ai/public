# Viewers: actions in the context menu — Run Results

**Date**: 2026-05-12
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Close all, open demog dataset | 6s | PASS | PASSED | 5850 rows, 11 cols (USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY) |
| 2 | Add Histogram, Line Chart, Bar Chart, Trellis plot | 5s | PASS | PASSED | Clicking the four toolbox icons added all viewers; `tv.viewers` reports `["Grid","Histogram","Line chart","Bar chart","Trellis plot"]` |
| 3 | Grid: Tooltip section has Hide, Edit..., Use as Group Tooltip, Remove Group Tooltip | 4s | PASS | PASSED | All four labels present; Grid has one extra `Current Column` item below them — a grid-specific addition, not a regression |
| 4 | Histogram: Tooltip section has the four expected items | 3s | PASS | PASSED | Top-level groups: General, Reset View, ... Tooltip, Properties...; Tooltip submenu labels match exactly |
| 5 | Line chart: Tooltip section has the four expected items | 3s | PASS | PASSED | Top-level groups include General, Reset View, Tools, Data, ..., Tooltip; Tooltip submenu labels match exactly |
| 6 | Bar chart: Tooltip section has the four expected items | 3s | PASS | PASSED | Top-level groups include General, On Click, Orientation, ..., Tooltip; Tooltip submenu labels match exactly |
| 7 | Trellis plot: Tooltip section has the four expected items | 3s | PASS | PASSED | Top-level groups: General, On Click, Pick Up / Apply, Tooltip, Properties..., To Script; Tooltip submenu labels match exactly |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~2m |
| grok-browser execution (scenario steps) | ~30s |
| Execute via grok-browser (total) | ~2m 30s |
| Spec file generation | ~3m |
| Spec script execution | 19s |
| **Total scenario run (with model)** | ~9m |

## Summary

All five viewers (Grid, Histogram, Line chart, Bar chart, Trellis plot) expose a `Tooltip` group in their right-click context menu, and that submenu contains the four expected items: `Hide`, `Edit...`, `Use as Group Tooltip`, `Remove Group Tooltip`. Grid adds one extra child item, `Current Column`, which is a grid-specific affordance and does not contradict the scenario. The Playwright spec passed on the second attempt; the first run failed because synthetic `dispatchEvent` right-clicks didn't trigger Dart's submenu DOM creation — switching to `page.mouse.click(..., {button: 'right'})` plus an explicit `mouseenter`/`mouseover`/`mousemove` on the Tooltip group resolved it. Total scenario run (with model): ~9m.

## Retrospective

### What worked well
- The same `.d4-menu-item.d4-menu-group` + label-text lookup pattern from `references/widgets/menu.md` worked across all five viewer types — no per-viewer special casing was needed once the menu was open.
- Submenu items remain populated in the DOM even when the container has `display:none`, so reading labels after a brief wait did not require driving a real cursor over the `Tooltip` group during the MCP run.
- `softStep` lets every viewer's check run independently and aggregates failures — useful here because the first spec run would otherwise have stopped after the first viewer.

### What did not work
- In the spec's first run, synthetic `target.dispatchEvent(new MouseEvent('contextmenu', …))` opened the popup but Dart never created the `.d4-menu-item-container` submenu DOM under the `Tooltip` group. The labels array came back empty for every viewer. In the MCP browser the same `dispatchEvent` path did populate the submenu — the difference appears to be Playwright's stricter `event.isTrusted` semantics or the Dart event-loop scheduling under headed Playwright. Using `page.mouse.click(x, y, {button: 'right'})` plus an explicit hover on the Tooltip group fixed it.
- The first `grok test` invocation still falls back to `UsageAnalysis:test returned null instead of a DataFrame` when run via the globally installed `datagrok-tools@6.1.10` (no `--skip-puppeteer` support). Worked around by invoking `node public/tools/bin/grok.js test …` directly — same workaround as noted in `default-tooltip-run.md`.

### Suggestions for the platform
- Either honor `event.isTrusted = false` synthetic right-clicks in the Dart context-menu code path, or document that test code must use real-mouse Playwright events for context menus. The current behaviour is silent: the popup opens, top-level group labels are present, but submenu children never materialise, which is a very misleading partial-state from a test author's perspective.
- Ship `datagrok-tools` 6.2+ (with `--skip-puppeteer`) as the global `npm i -g datagrok-tools` build so contributors don't need to keep a locally-built `public/tools/bin/grok.js` on PATH for Test Track replays.

### Suggestions for the scenario
- The scenario lists four expected items but is silent on whether *additional* items are allowed. Grid's `Current Column` extra item passed under the current "contains all four" reading; consider clarifying whether the section must contain exactly the four items or at least the four items.
- Add a sentence noting that "right-click on each viewer" means on the chart area (canvas), not on the title bar — clicking on the title bar opens the hamburger menu rather than the context menu, and produces a different set of items.
- Mention which viewer types are explicitly in scope. The current list (Grid, Histogram, Line Chart, Bar Chart, Trellis plot) is enumerated in the steps but could be promoted to a pre-condition table so the reader knows to add exactly those viewers before testing.
