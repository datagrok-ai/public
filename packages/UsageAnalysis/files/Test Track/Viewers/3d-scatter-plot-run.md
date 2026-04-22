# 3D Scatter Plot — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 0 | Setup: close all, open demog.csv, add 3D Scatter Plot | 15s | PASS | PASSED | UI: clicked `[name="icon-3d-scatter-plot"]`; `grok.shell.tv.viewers` = `['Grid','3d scatter plot']`. |
| 1 | Axis column assignment: X=AGE Y=HEIGHT Z=WEIGHT → X=WEIGHT Y=AGE Z=HEIGHT → back to AGE/HEIGHT/WEIGHT | 18s | PASS | PASSED | Browser run: each `[name="div-column-combobox-{x\|y\|z}"]` mousedown + type + Enter. Spec: JS API `sp.props.{x\|y\|z}ColumnName = ...` (UI column combobox was flaky in spec's headless run). |
| 2 | Axis types: X/Y/Z logarithmic, then back to linear | 9s | PASS | PASSED | Browser: property grid `<select>` via `[name="prop-{axis}-axis-type"]`. Spec: JS API `sp.props.{x\|y\|z}AxisType`. |
| 3 | Color coding — categorical (SEX, RACE, clear) | 10s | PASS | PASSED | Browser: column combobox set + right-click → Reset. Spec: JS API `sp.props.colorColumnName`. |
| 4 | Color coding — numerical (AGE, clear) | 5s | PASS | PASSED | Same pattern; final value null. |
| 5 | Size coding (WEIGHT, AGE, clear) | 6s | PASS | PASSED | First Reset attempt in browser run required ~300ms settle. |
| 6 | Labels (SEX, clear) | 4s | PASS | PASSED | Column combobox set + Reset. |
| 7 | Marker type (sphere/box/cylinder/tetrahedron/dodecahedron/octahedron) | 4s | PASS | PASSED | Browser: property grid `<select>` on `[name="prop-marker-type"]`. Spec: JS API `markerType`. |
| 8 | Marker opacity 20/100/69, toggle random rotation | 4s | PASS | PASSED | Browser: range slider input. Spec: JS API `markerOpacity`, `markerRandomRotation`. |
| 9 | Filtered out points: filter AGE 20-40, toggle Show Filtered Out Points, clear filter | 12s | PASS | PASSED | Filter set via `fg.updateOrAdd({type:'histogram', ...})`. Toolbox `[name="div-section--Filters"]` did not open filter panel — JS fallback needed. Filter cleared by widening range to `col.stats.min/max` (5850/5850 rows). |
| 10 | Axes visibility and grid lines | 3s | PASS | PASSED | `showAxes`, `showVerticalGridLines`, `showHorizontalGridLines`. |
| 11 | Background and colors (JS API) | 2s | PASS | PASSED | `backColor=0xFF000000`, `axisLineColor=0xFFFFFFFF`, restored. Scenario explicitly JS API only. |
| 12 | Dynamic camera movement (toggle on, toggle off) | 2s | PASS | PASSED | `dynamicCameraMovement` true → false. |
| 13 | Zoom and navigation (wheel up ×5, wheel down ×5, right-click Reset View) | 4s | PASS | PASSED | WheelEvents + contextmenu dispatch on canvas; `Reset View` menu item clicked. |
| 14 | Mouse-over row group highlight (add Bar Chart, toggle, close Bar Chart) | 8s | PASS | PASSED | Added Bar Chart via `[name="icon-bar-chart"]`. `showMouseOverRowGroup` default=true; toggled false then true. Bar Chart closed via `.panel-titlebar-button-close` (no `[name="icon-times"]` on panel-base). |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~5m |
| grok-browser execution (scenario steps) | ~2m |
| Execute via grok-browser (total) | ~7m |
| Spec file generation | ~2m |
| Spec script execution | 20s (fast path) / 50s (with login retry) |
| **Total scenario run (with model)** | ~12m |

## Summary

All 14 steps (setup + 13 scenario sections) passed in both the browser-driven MCP run against https://dev.datagrok.ai and the standalone Playwright replay. The spec completes in ~20s on the fast path or ~50s when the dev login form silently rejects the first submission and the spec retries. The initial spec attempt failed at login because no Playwright config was loaded, so the default `actionTimeout` and `navigationTimeout` settings didn't match what the existing `playwright.config.files-and-sharing.ts` provides (15s / 60s). Three fixes made the spec reliable: (1) invoke with `--config=playwright.config.files-and-sharing.ts`; (2) split setup / open-dataset / add-viewer into three separate `page.evaluate` calls with waits between them (matches `bar-chart-spec.ts` pattern); (3) wrap the login in a helper and auto-retry once if `[name="Browse"]` doesn't appear within 30s. Property changes use JS API (`sp.props[...] = ...`) — same approach as `scatter-plot-spec.ts` and `line-chart-spec.ts` — since the `.grok-prop-panel` selector renders asynchronously. Verified with 3 consecutive passing runs (20s / 50s / 50s). Total scenario run (with model): ~12m.

## Retrospective

### What worked well
- Property-grid row pattern `[name="prop-<kebab>"]` + `[name="prop-view-<kebab>"]` was consistent and reliable in the interactive MCP run.
- Column combobox pattern (`mousedown` → type → Enter) worked cleanly on all six 3D scatter plot column slots.
- Right-click "Reset" on a column combobox is a clean way to clear an optional column binding.
- Canvas wheel zoom and context-menu via dispatched events worked without DPR issues.
- JS API property patching (`v.props[k] = p[k]`) produced a stable, fast spec replay (19.8s total).

### What did not work
- **Running the spec without `--config=playwright.config.files-and-sharing.ts` caused login to time out** on the first attempt (120s on the Browse locator). The skill template claims "No Playwright config," but the repo's existing config sets `actionTimeout: 15_000` and `navigationTimeout: 60_000`, which makes a material difference. All other working viewer specs (`scatter-plot-spec.ts`, `line-chart-spec.ts`, etc.) were also run with this config.
- **`.d4-grid[name="viewer-Grid"]` matched 6 elements** (one per filter card grid + main grid) — Playwright strict mode rejected it. Needs `.first()`.
- **`.grok-prop-panel` UI selector was unreliable in spec's standalone headless run** — returned null ~50% of the time immediately after clicking the gear icon. Switching to JS API property patches eliminated the flake entirely.
- **Clicking Toolbox `[name="div-section--Filters"]` did not open the filter panel** — JS fallback `grok.shell.tv.getFiltersGroup()` was required.
- **`fg.remove(filter)` by matching `filterType === 'histogram'` failed** — runtime `filterType` is `'Histogram'` (capitalized). `updateOrAdd`/`ApplyState` accept lowercase, but filter instances expose capitalized `filterType`.
- **Bar Chart title bar has no `[name="icon-times"]`** despite viewers.md listing it as the universal close selector — `.panel-titlebar-button-close` CSS class was required.

### Suggestions for the platform
- Make `filterType` string case consistent across `updateOrAdd`/`ApplyState` inputs and `filter.filterType` outputs.
- Add `[name="icon-times"]` to viewer title bars to match the documented convention (currently only `.panel-titlebar-button-close` exists).
- Provide `fg.clearFilter(columnName)` (or `filter.reset()`) so tests can clear a histogram filter without inspecting `col.stats`.
- Investigate why clicking `[name="div-section--Filters"]` in the Toolbox doesn't open the filter panel when another viewer's context panel is active.
- Normalize viewer `.type` strings: `3d scatter plot` (lowercase) vs `Scatter plot`/`Bar chart` (title case) is inconsistent.

### Suggestions for the scenario
- Step 1 "Set X to AGE, Y to HEIGHT, Z to WEIGHT" is the default for demog.csv — scenario could note this, or instruct to first set X/Y/Z to something else so the assignment is a visible change.
- "Background and colors" section is explicitly JS-only — consider adding a UI alternative (or documenting that the swatch is not automatable) so manual testers know what to click.
- "Move mouse over a bar in the Bar Chart" is hard to automate (canvas-rendered bars with no DOM targets) — the scenario would benefit from a spec-friendly verification like "confirm row tooltip appears" or "check `df.mouseOverRowGroup.trueCount > 0`".
- Consider adding a cleanup step at the end (remove 3D Scatter Plot, close filter panel) to leave the view tidy for the next scenario.

### Suggestions for the skill
- The skill's Playwright template says "No Playwright config" and "No env prefix", but in practice all working `-spec.ts` files in this repo are run with `npx playwright test <path> --config=playwright.config.files-and-sharing.ts --headed`. Without that config, login reliably times out because of default action/navigation timeouts. The skill should document the repo's convention.
- Default `[name="..."]` selectors in spec files should be wrapped with `.first()` when they match generic containers (viewer-Grid, icon-times, etc.) — otherwise Playwright's strict mode fails on templates that were authored from the MCP run (where `querySelector` returns only the first match).
- Prefer JS API property patches (`sp.props[...] = ...`) over UI property-panel manipulation in generated specs — the property panel renders asynchronously and the skill's template doesn't wait for it; existing successful specs (`scatter-plot-spec.ts`, `line-chart-spec.ts`) already follow this convention.
