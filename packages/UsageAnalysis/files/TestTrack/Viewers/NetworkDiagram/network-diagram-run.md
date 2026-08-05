# Network diagram (Playwright) — Run Results

**Date**: 2026-08-05
**URL**: https://dev.datagrok.ai (1.28.0)
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add Network diagram from the Viewers toolbox | PASS | Node columns auto-picked (SEX / CONTROL); canvas carries content |
| 2 | Suspend simulation freezes the layout | PASS | Idle canvas delta after freezing is 0 |
| 3 | Switch Node 1 to RACE | PASS | Selector text changes, ~110k pixels redrawn |
| 4 | Colour nodes by SEX and size them by AGE | PASS | Property grid shows both, diagram repaints |
| 5 | Colour edges by AGE | PASS | Edge colours change |
| 6 | Clicking a node selects the rows behind it | PASS | 3243 rows selected by a real click |
| 7 | Click-selection switches off stop the clicks from selecting | PASS | Both **Select Rows On Click** and **Select Edges On Click** off → 0 selected |
| 8 | Show Column Selectors hides the on-viewer selectors | PASS | DOM visibility round-trip |
| 9 | Show Arrows draws directions on the edges | PASS | Arrow heads repaint immediately; setting survives a rebuild |
| 10 | Show Filtered Out Nodes | PASS | Pixel count grows when the option goes on — the filtered-away nodes are drawn again |
| 11 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 39–41s (two consecutive runs) |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `waitForViewerRepaint`,
`waitForPropertyValue`, `expect.poll`), so the run costs what the UI actually
takes and a missing repaint fails with a message naming what was expected.

* This viewer draws into a plain 2d canvas, so repaints are proved with real pixel
  histograms (`snapshotCanvasColors` / `diffCanvasColors`) rather than screenshot
  hashes. **Suspend Simulation** is switched on early precisely so those diffs
  mean something — while the physics runs, the canvas changes on its own and any
  "it repainted" assertion would pass for free.
* Nodes have no DOM handles and the underlying vis.js network is not reachable
  from the viewer object, so node positions are found from the pixels: saturated
  colour blobs bucketed into 40px cells, largest first.
* `Edge Width` exists **twice** in the property grid — a column selector under
  Data and a number under Misc, both rendered as `tr[name="prop-edge-width"]`.
  Property-grid helpers therefore take an optional category to disambiguate;
  without it they silently drive the first match.

## Filed bugs

| Property | Ticket | State on dev 1.28.0 (2026-08-05) |
|---|---|---|
| **Show Arrows** | GROK-20617 | **Fixed** — `none` → `to` repaints straight away |
| **Show Filtered Out Nodes** | GROK-20618 | **Fixed** — with `SEX = F` applied the pixel count goes 1786 → 2662 when the option is switched on |
| **Edge Width** | GROK-17125 | Still open (covers Shapes and Edge Width) |

Both fixes were caught by the `knownOpenBug()` wrapper: on 2026-08-05 the run
threw `[KNOWN_BUG_FIXED:GROK-20617]`, and the wrappers were replaced with plain
hard assertions. In Jira, GROK-20618 was in *Ready for QA* and GROK-20617 was
still *Open* at that point — the run is what proved both.

The Show Filtered Out Nodes step also had to be repaired, and the reason is worth
keeping: it filtered `SEX = F` while the graph was built on **RACE**, so no node
was ever filtered away and the canvas legitimately did not change. The step
"reproduced" a bug that its own state guaranteed, and it hid the fix. It now
forces **Node 1** to the filtered column first, and asserts that the pixel count
GROWS with the option on instead of merely asserting "something repainted".

## Not automated

* Double-click node expansion, node dragging, node images, colour-scheme pickers,
  and the remaining physics toggles — kept as a manual list in `network-diagram.md`.
