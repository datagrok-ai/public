# Viewers: default tooltip visibility — Run Results

**Date**: 2026-05-12
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Close all. Open demog dataset | 5s | PASS | PASSED | Setup via JS API (closeAll + readCsv + addTableView); 5850 rows, 11 columns |
| 2 | Open grid properties; enable Show Visible Columns In Tooltip | 8s | PASS | PASSED | Clicked grid gear; ticked `[name="prop-show-visible-columns-in-tooltip"]` checkbox; verified `grid.props.showVisibleColumnsInTooltip === true` |
| 3 | Open Scatter, Box, Histogram, Line, Bar, Trellis | 7s | PASS | PASSED | Clicked each Toolbox icon; final `tv.viewers` has all 7 types (Grid + 6) |
| 4 | Right-click viewer → Tooltip > Hide | 4m | PASS | PASSED | Submenu opens reliably only with `mouseover+mouseenter` then incremental `mousemove` slope sequence; click set df tag `.tooltip-visibility="false"` (affects row-based tooltips: grid/scatter/box) |
| 5 | Right-click viewer → Tooltip > Show Custom (replaces Hide) | 8s | PASS | PASSED | Submenu now showed `Show Custom` (Hide gone); click set `.tooltip-visibility="true"` |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 4m 10s |
| grok-browser execution (scenario steps) | 25s |
| Execute via grok-browser (total) | 4m 35s |
| Spec file generation | 2m |
| Spec script execution | 22s |
| **Total scenario run (with model)** | 7m |

## Summary

All five scenario steps passed end-to-end against dev.datagrok.ai and the generated Playwright spec
runs in 18.8s (1 test passed). The `Tooltip > Hide` and `Tooltip > Show Custom` context-menu actions
toggle a dataframe-level tag (`.tooltip-visibility`), which is the observable side-effect the spec
asserts. Most of the model time was spent figuring out the reliable way to open the Tooltip submenu.
**Total scenario run (with model): ~7m.**

## Retrospective

### What worked well
- `[name="prop-show-visible-columns-in-tooltip"]` is a stable selector for the grid's tooltip toggle.
- `[name="div-Tooltip"]` + `[name="div-Tooltip---Hide"]` / `[name="div-Tooltip---Show-Custom"]` are stable selectors once the submenu is rendered.
- The `Show Custom` ↔ `Hide` label swap is a clean DOM-visible signal that lets the spec assert the action took effect without sampling tooltips on canvas.
- Asserting on the dataframe tag `.tooltip-visibility` is much more reliable than hovering over a canvas and waiting for a `.d4-tooltip` DOM element.

### What did not work
- A single `mousemove` (or even `mouseover+mouseenter+mousemove`) on `[name="div-Tooltip"]` did **not** open the submenu — Dart's hover-slope tracking ignored the synthetic events.
- `tooltipGroup.click()` on the menu **closes** the popup entirely instead of opening the submenu — definitely not what a user gets.
- The viewer's own `showTooltip` property does not change when `Tooltip > Hide` is invoked, which initially misled me into thinking the click had no effect. The actual side-effect lives on `df.tags['.tooltip-visibility']`.

### Suggestions for the platform
- Document (or expose a JS API helper for) the `.tooltip-visibility` dataframe tag — it's the single source of truth for the `Hide` / `Show Custom` toggle but isn't discoverable from a viewer's property panel.
- Consider adding a `name=` attribute on the dynamic submenu container or emitting a `d4-submenu-show` debug-attribute so automated tests can wait for `[data-submenu-open="Tooltip"]` instead of probing element visibility in a polling loop.
- Synthetic `mousemove` should open submenus the same way real movement does — the slope-tracking logic in `menu.dart` is currently the single biggest source of flakiness for automated UI tests.

### Suggestions for the scenario
- Clarify in step 4 which viewer to right-click on (scenario reads as if it could be any one viewer); current wording made me unsure whether to apply `Hide` once or once per affected viewer.
- Add a one-line note that `Tooltip > Hide` is a dataframe-level toggle (not per-viewer) — otherwise the expected-result list of "grid, scatter, box plot only" reads as a per-viewer override rather than as the consequence of which viewer types render row tooltips.
- Pre-condition could state explicitly that the dataframe must have `.tooltip-visibility` unset on entry, so that re-runs from a dirty session don't silently start from `false`.
