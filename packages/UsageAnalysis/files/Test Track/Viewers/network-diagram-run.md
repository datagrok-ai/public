# Network Diagram — Run Results

**Date**: 2026-04-14
**URL**: https://dev.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open demog dataset | PASS | ~3s | PASSED | 5850 rows via `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')` |
| 2 | Click Network diagram in Toolbox, auto-pick columns | PASS | ~2s | PASSED | Auto-picked `SEX` / `CONTROL` (fewest categories). Canvas drawn, no console errors. Selector names are `div-column-combobox-node1` / `node2` (lowercase), not `Node1ColumnName` as reference claims |
| 3 | Switch Node 1 = RACE, Node 2 = DEMOG | PASS | ~1s | PASSED | Column combobox popup on `mousedown` rendered as an empty Grid backdrop (not the expected ColumnGrid) → fell back to `viewer.props.node1/2ColumnName`. Graph rebuilt |
| 4 | Click a node — rows for incident edges selected | SKIP | — | SKIPPED | Canvas-based; synthetic PointerEvents do not reach vis.js/Hammer.js. Only verified `selectRowsOnClick = true` default |
| 5 | Shift+click, Ctrl+click nodes | SKIP | — | SKIPPED | Same limitation as step 4 |
| 6 | Click an edge — backing rows selected | SKIP | — | SKIPPED | Same limitation. `selectEdgesOnClick = true` verified |
| 7 | Double-click empty canvas — clear + fit | SKIP | — | SKIPPED | Same limitation |
| 8 | Click gear on viewer title bar | FAIL→JS | — | JS API | No `[name="icon-font-icon-settings"]` icon in the viewer DOM (title bar not rendered for this viewer). Fell back to JS API to set props |
| 9 | Data: edgeColor=AGE(avg), edgeWidth=WEIGHT(avg), node1Size=AGE, node1Color=SEX | PASS | ~1.5s | PASSED | Gradient edges & coloured Node 1 circles visible in screenshot |
| 10a | Style: showColumnSelectors off → on | PASS | ~1s | PASSED | `display:flex` → `display:none` on both combo boxes; toggled back |
| 10b | Style: showArrows = to | PASS | <1s | PASSED | Property set |
| 10c | Style: suspendSimulation = true | PASS | <1s | PASSED | Property set |
| 11 | Filter AGE > 40, toggle showFilteredOutNodes | PASS | ~1.4s | PASSED | filter trueCount 3715 / 5850; `showFilteredOutNodes` flipped false → true |
| 12 | Close viewer via × icon | FAIL→JS | ~0.5s | JS API | No `icon-times` in viewer DOM. Fell back to `viewer.close()`; 0 Network diagrams remain; no warnings |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~405 s |
| Spec file generation | ~3 s |
| Spec script execution | 7.4 s (1 passed) |

## Summary

Viewer core functionality (auto-pick, column rebinding, property-driven styling, filter integration,
clean close) all work. Four steps relating to node/edge clicks on the canvas could not be
reproduced because vis.js listens through Hammer.js, which ignores synthetic
`PointerEvent` / `MouseEvent` dispatches. Title-bar icons (gear, ×) are not rendered in the DOM of
this viewer even with `body.selenium`, forcing JS API fallback for open-settings and close.

## Retrospective

### What worked well
- JS API fallback for property setting is fast and reliable (edgeColor, edgeWidth, node1Size, etc.)
- Filter integration via `df.filter.init` immediately propagates to the viewer
- Column combo boxes expose `name=div-column-combobox-node1/node2` attributes

### What did not work
- **Node/edge clicks on the canvas** — vis.js uses Hammer.js for pointer events; synthetic
  `PointerEvent`/`MouseEvent` dispatched on the canvas do not trigger `_onClick`. Result:
  `df.selection.trueCount` stayed 0.
- **Column selector popup** — `mousedown` opened a backdrop, but the popup rendered as an
  empty Grid root instead of the expected ColumnGrid. Typing/searching inside was not possible.
- **Gear icon / × icon not in DOM** — the viewer has no rendered title bar in the DOM (only
  the canvas region). `body.selenium` does not help here.

### Suggestions for the platform
- Expose the underlying `vis.Network` instance on the Dart viewer object
  (e.g. `nd.dart.visNetwork`) or add helpers `nd.clickNode(value)` / `nd.clickEdge(from, to)` so
  automated tests can exercise click behaviour without canvas hit-testing.
- Always render viewer title bar elements (gear, close) in the DOM with `body.selenium`,
  even for canvas-only viewers — they are required for scripted automation.
- Align column combo box `name=` attributes with the reference doc
  (`div-column-combobox-Node1ColumnName-`) or update the reference to match the actual
  `div-column-combobox-node1` attribute.
- Investigate why the Network Diagram column-selector popup renders as an empty Grid backdrop
  instead of the ColumnGrid populated with the dataframe's columns.

### Spec runnability notes
- Initial spec used `page.goto(baseUrl)` (box-plot pattern), which produced a fresh context
  without a logged-in session and timed out waiting for `grok` to initialise.
- Fixed by switching to the `calendar-spec.ts` pattern: `chromium.connectOverCDP(...)` and
  reuse the already-open Datagrok tab from the live Chrome session.
- Also added `.first()` to the `.d4-grid[name="viewer-Grid"]` locator to avoid Playwright
  strict-mode violations from filter cells that share the same `name=` attribute.

### Suggestions for the scenario
- Pre-condition: mention that steps 4–7 require real user interaction with a canvas — they are
  not reliably scriptable. Suggest a JS-API-only variant for CI runs.
- Step 2 says auto-picked columns are e.g. `sex` / `race`; on demog the heuristic picks
  `SEX` / `CONTROL`. Update the example.
- Step 8 mentions the **Gear** icon; for canvas-only viewers like Network Diagram, the gear
  is not in the DOM — document the F4 or sidebar-Settings alternative for scripted runs.
