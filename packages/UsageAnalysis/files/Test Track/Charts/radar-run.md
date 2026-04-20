# Radar viewer (Charts package) — Run Results

**Date**: 2026-04-20
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open earthquakes.csv + add Radar viewer | 35s | PASS | PASSED | Add viewer ribbon icon (`i.svg-add-viewer`) opens dialog; clicked Radar card. Dataset path corrected to `System:DemoFiles/geo/earthquakes.csv`. |
| 2 | Open demog.csv + add Radar viewer | 18s | PASS | PASSED | Same flow; Radar viewer present in current TableView. |
| 3 | Open Radar settings (gear icon) | 6s | PASS | PASSED | Gear lives on `panel-base` ancestor (Radar is a JsViewer in a panel). Property categories: Data, Selection, Color, Value, Misc, Style, Legend, Description. |
| 3a | Switching tables (Table dropdown) | 5s | PASS | PASSED | `setOptions({table: 'Table'})` rebound viewer; values columns auto-changed to earthquakes columns. Restored to demog. |
| 3b | Selection check-boxes | 3s | PASS | PASSED | Toggled `showCurrentRow` / `showMouseOverRow` / `showMouseOverRowGroup` off→on; values reflected. |
| 3c | Increase/decrease Values | 4s | PASS | PASSED | Reduced to 2 columns, increased back to 4. Viewer redrew with corresponding axes. |
| 3d | Style (color) changes | 3s | PASS | PASSED | Background min/max → light pink/red, line color → blue. Verified visually in screenshot. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~3m |
| grok-browser execution (scenario steps) | ~1m 14s |
| Execute via grok-browser (total) | ~4m 14s |
| Spec file generation | ~50s |
| Spec script execution | 10.8s |
| **Total scenario run (with model)** | ~5m 15s |

## Summary

The Radar viewer (Charts package) reproduced cleanly on dev for both earthquakes.csv and demog.csv. Switching tables, selection check-boxes, value-column count changes, and color/style changes all reflected on the viewer. Generated Playwright spec passed end-to-end in 10.8s. **Total scenario run (with model)**: ~5m 15s.

## Retrospective

### What worked well
- `i.svg-add-viewer` ribbon icon reliably opens the Add Viewer gallery, and the gallery card click works via synthetic `MouseEvent` chain (mousedown → mouseup → click → dblclick).
- `setOptions(...)` + `getOptions().look` is a solid, fast verification path for JsViewer property changes.
- Reusing the live MCP Chrome profile via `chromium.connectOverCDP('http://127.0.0.1:9222')` + `browser.contexts()[0]` keeps auth and avoids the login screen.

### What did not work
- Initial `readCsv('System:DemoFiles/earthquakes.csv')` failed with HTTP 502 — actual path is `System:DemoFiles/geo/earthquakes.csv`. Scenario JSON only listed `demog.csv`, so the earthquakes path had to be discovered via `grok.dapi.files.list` with recursive search.
- Radar (a JsViewer from Charts) lacks the standard d4 viewer title bar; the gear icon is one level up on the wrapping `.panel-base`. Standard `[name="viewer-Radar"] [name="icon-font-icon-settings"]` selector does not match — need to use `radar.closest('.panel-base')` first.
- `tv.viewers.toList()` is not a function on dev — must iterate with `for...of`.

### Suggestions for the platform
- JsViewers should expose the same in-viewer title-bar icons (gear / close / hamburger / help) as native viewers when `body.selenium` is set, so a single selector pattern works across viewer types.
- 502 from a missing file path is misleading — should be 404 with a clear message.

### Suggestions for the scenario
- Add `System:DemoFiles/geo/earthquakes.csv` to the `datasets` JSON block (currently only `demog.csv`).
- Clarify Step 1/2: the "Add viewer" icon is on the **top toolbar** (ribbon), not the Toolbox Viewers section.
- Step 3 says "Check all properties" — split into explicit sub-steps mirroring the four bullet points, each with a measurable expected outcome.

---
{
  "order": 28,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/geo/earthquakes.csv"]
}
