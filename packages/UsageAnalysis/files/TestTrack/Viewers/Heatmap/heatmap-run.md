# Heat map (Playwright) — Run Results

**Date**: 2026-08-04
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Add Heat map from the Viewers toolbox | PASS | Content canvas carries ~537k coloured pixels; **Is Heatmap** is on |
| 2 | Heatmap Colors off falls back to plain cells | PASS | `KNOWN_BUG_REPRODUCED:GROK-20619` |
| 3 | Global Color Scaling recolours the cells | PASS | Canvas repaints on both flips |
| 4 | Col Labels Orientation rotates the header | PASS | Auto → Vert → Horz → Auto |
| 5 | Max Heatmap Columns limits the columns on screen | PASS | 100 → 3 redraws the map, then back to 100 |
| 6 | Show Heatmap Scrollbars hides the range sliders | PASS | `x-slider` hidden, then visible again |
| 7 | Is Heatmap off turns the viewer into a plain grid | PASS | Canvas repaints in both directions |
| 8 | Row Height is not offered in heatmap mode | PASS | `KNOWN_BUG_REPRODUCED:GROK-20619` |
| 9 | Clicking a cell makes its row current | PASS | Current row moves to the clicked cell's row |
| 10 | Alt+drag zooms into an area and the slider follows | PASS | Slider window shrinks; double-click resets to the full extent |
| 11 | Filtering the table redraws the heat map | PASS | `SEX = M` |
| 12 | Close the viewer from its title bar | PASS | Viewer removed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec run (full) | 30s |

## Notes

No fixed `waitForTimeout` pauses: every wait blocks on the condition the step
cares about (`waitForCanvasChange`, `waitForCanvasQuiet`, `waitForViewerRepaint`,
`waitForPropertyValue`, `expect.poll`), so the run costs what the UI actually
takes and a missing repaint fails with a message naming what was expected.

* The viewer hosts three canvases — a 10px scrollbar strip, the content canvas
  (`canvas[name="canvas"]`) and an overlay. Pixel reads target the content one
  explicitly; the shared canvas helpers take a selector for exactly this reason.
* Alt+drag zoom is asserted through the range slider's window size rather than a
  screenshot: zooming shrinks the span, and double-clicking the slider resets it
  to the **full** extent — wider than where the test started, because the viewer
  opens already scrolled on a 5850-row table.

## Filed bugs — GROK-20619

| Setting | Observed |
|---|---|
| **Heatmap Colors** | Unchecking it leaves the cells fully coloured; the content canvas moves by ~10 pixels out of 537k, and no cell text appears |
| **Row Height** | `28` → `20` commits the value but leaves the content canvas byte-identical |

**Row Height is a grid-only option** — its tooltip says so, and it is not meant to
apply in heatmap mode. The agreed fix is to hide it there rather than make it work,
so the bug is that the property is offered at all and silently accepts a value.
Prior art: GROK-11515 (*Viewers | Heatmap: the 'Row Height' property doesn't work*,
closed in 2022) treated it as a defect in the property itself.

Not related to GROK-19087 (*Heatmap: some Styles properties are broken*), which
covers Max font size, Col labels orientation and allow col resizing — different
properties, and under Style rather than Misc and Rows.

Both steps keep the **desired** assertion wrapped in `knownOpenBug()` from
`helpers/known-open-bug.ts`:

* Heatmap Colors — the cells must stop being colour-filled (a repaint of more than
  1000 px, well above the ~10 px the current-cell outline accounts for);
* Row Height — the property must not have a row in the property grid at all, since
  the agreed fix is to hide this grid-only option in heatmap mode.

While the bug reproduces both log `[KNOWN_BUG_REPRODUCED:GROK-20619]` and the run
is green. When the fix lands they throw `[KNOWN_BUG_FIXED:GROK-20619]` — at which
point the Heatmap Colors wrapper becomes a plain `expect` and the Row Height step
is simply deleted.

## Not automated

* Right-click content panning and custom sorting on a categorical column (needs
  SPGI) — kept in `heatmap-ui.md`.
