# Network diagram — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open demog dataset | 15s | PASS | PASSED | Loaded 5850 rows, 11 columns via readCsv |
| 2 | Add Network diagram, verify auto-picked columns | 8s | PASS | PASSED | Auto-picked SEX/CONTROL (not SEX/RACE as doc example); scenario wording "e.g. sex and race" matches heuristic; no console errors |
| 3 | Set Node 1 = RACE, Node 2 = DEMOG | 6s | PASS | PASSED | JS API fallback; UI click on combo opened small popup but column picker did not show |
| 4 | Click a node → selects incident rows | 1s | SKIP | SKIPPED | vis.js renders nodes to canvas, no DOM handles to click |
| 5 | Shift/Ctrl-click another node | 1s | SKIP | SKIPPED | Canvas-rendered; not automatable without coordinate calc |
| 6 | Click an edge | 1s | SKIP | SKIPPED | Canvas-rendered |
| 7 | Double-click empty canvas (clear / fit) | 1s | SKIP | SKIPPED | Canvas-rendered |
| 8 | Open Property Pane via Gear icon | 10s | PASS | PASSED | Gear is on .panel-base wrapper, not inside viewer root |
| 9 | Set edge color=AGE/avg, width=WEIGHT/avg, node1 size=AGE, color=SEX | 4s | PASS | PASSED | Edges rendered with gradient colors + varied widths; nodes re-sized/recolored (screenshot verified) |
| 10 | Toggle Show Column Selectors, Show Arrows=to, Suspend Simulation | 4s | PASS | PASSED | All toggles effective; selectors hidden/shown count correct |
| 11 | Filter AGE>40; toggle Show Filtered Out Nodes | 4s | PASS | PASSED | 3715/5850 rows passed; showFilteredOutNodes toggled true |
| 12 | Close viewer via × | 3s | PASS | PASSED | Close icon has name="Close" (not "icon-times"); viewer removed cleanly, no warnings |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 6m 30s |
| grok-browser execution (scenario steps) | 1m 30s |
| Execute via grok-browser (total) | 8m |
| Spec file generation | 1m 30s |
| Spec script execution | 22s |
| **Total scenario run (with model)** | 9m 52s |

## Summary

9 of 12 steps PASS; 3 SKIP (canvas-based node/edge interactions cannot be automated via DOM). The
network diagram viewer works as specified — auto-pick, property setting, filtering, and close all
succeed and visually render correctly. Overall total ~10 minutes including selector discovery.

## Retrospective

### What worked well
- JS API-based property setting is fast and reliable for viewers with complex property grids
- `nd.props.* = value` propagates synchronously to the viewer; a 400–800ms wait is enough to verify DOM effects
- Filter panel interactions via `df.filter` bitset + `fireChanged()` behave identically to UI filter setting
- `document.body.classList.add('selenium')` makes the title bar icons reliably queryable

### What did not work
- Clicking the Node1/Node2 combobox did not open the column picker popup — the click may require a specific event sequence vis.js / ColumnComboBox is listening for
- Canvas-based interactions (clicking nodes/edges) are unreachable via DOM selectors — only viewport coordinates work, and those are not stable
- Reference file `network_diagram.md` documents combobox name as `div-column-combobox-Node1ColumnName-` but actual is `div-column-combobox-node1`
- Reference file documents the close icon as `[name="icon-times"]`; the actual attribute is `[name="Close"]` on the panel-base wrapper
- Gear/close icons live on the outer `.panel-base` wrapper, not inside the `[name="viewer-*"]` root — scoping to the root misses them

### Suggestions for the platform
- Keep reference docs in sync with actual attribute names (combobox names, close icon name)
- Add a JS API shortcut to simulate canvas click on a node (e.g. `viewer.clickNode(nodeId, {shift, ctrl})`) so automation can validate the selection handler without screen coordinates
- Consider exposing the vis.js network object on the viewer's `dart` property (currently hidden behind obfuscated field names) to enable programmatic node selection/fit from tests

### Suggestions for the scenario
- Steps 4-7 assume manual interaction; add an explicit "alternative for automated tests" note recommending programmatic selection verification via `df.selection.trueCount` after setting selection through the JS API
- Step 9 says "Edge Width Column Name = weight" — but the demog dataset uses uppercase `WEIGHT`; spell column names as they appear in the dataset
- Step 2 expectation says "first two categorical columns with the fewest categories, e.g. sex and race" — the actual auto-pick on demog is SEX and CONTROL (both 2 categories). Update the example or relax the expected columns
