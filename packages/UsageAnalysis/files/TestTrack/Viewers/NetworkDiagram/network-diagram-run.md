# Network diagram (Playwright) — Run Results

**Date**: 2026-08-04
**URL**: https://dev.datagrok.ai
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
| 9 | Show Arrows draws directions on the edges | PASS | `KNOWN_BUG_REPRODUCED:GROK-20617`; arrows drawn after the graph is rebuilt |
| 10 | Show Filtered Out Nodes | PASS | `KNOWN_BUG_REPRODUCED:GROK-20618` |
| 11 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 47s |

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

Three Misc properties commit their value but leave the diagram untouched until the
graph is rebuilt (e.g. by changing a node column):

| Property | Observed | Ticket |
|---|---|---|
| **Show Arrows** | `none` → `to` gives a canvas delta of 0; the arrow heads appear only after a rebuild | GROK-20617 |
| **Show Filtered Out Nodes** | with `SEX = F` applied, canvas pixel count is 2483 both off and on | GROK-20618 |
| **Edge Width** | `4` → `14` gives a canvas delta of 0, and the value is reset back to `4` by the next rebuild | GROK-17125 (already open, covers Shapes and Edge Width) |

Measured on dev 2026-08-04, with and without the simulation suspended, and with a
forced repaint (pointer move over the canvas) in between. Data-category properties
on the same viewer repaint immediately, so this is specific to these settings.

Both steps keep the **real** assertion — that the canvas repaints as soon as the
property is set — wrapped in `knownOpenBug()` from `helpers/known-open-bug.ts`:

* while the bug reproduces the run stays green and logs
  `[KNOWN_BUG_REPRODUCED:GROK-2061x]`;
* the moment the fix lands the assertion starts passing and the wrapper throws
  `[KNOWN_BUG_FIXED:GROK-2061x]`, so the spec tells us to set
  `related_bugs[].status: fixed` and replace the wrapper with a plain `expect`.

The Show Arrows step additionally keeps a hard assertion on what a user gets
today: after a rebuild the arrows are drawn and the setting survives.

## Not automated

* Double-click node expansion, node dragging, node images, colour-scheme pickers,
  and the remaining physics toggles — kept as a manual list in `network-diagram.md`.
