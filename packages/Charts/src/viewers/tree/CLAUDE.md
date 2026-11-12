# Tree

An echarts `tree` series that draws a categorical hierarchy: a root labelled **`All`**, one level
per hierarchy column, one node per distinct value under its parent. Registered by Charts as `Tree`,
so a feature reaches it as **`tree viewer`**.

| File | Purpose |
|---|---|
| `tree-viewer.ts` | The viewer: hierarchy pick, the tree it feeds echarts, click-to-select/filter, branch painting, the automation signals |
| `tree-status.ts` | The areas and readings `getWidgetStatus` reports |

## What a click does

Two different things, depending on what is under the pointer:

- **the line into a node** — the zrender click handler (`getTargetPath`) takes the element created
  just before the clicked line, which is the node it leads into, and selects (or filters) that
  node's rows; Shift/Ctrl add to the selection.
- **the node's symbol** — echarts folds or unfolds the node. `getTargetPath` finds nothing for a
  symbol, so it selects nothing.

Nodes deeper than `initialTreeDepth` start folded. A folded node's descendants have no layout and
no graphic element: they are not drawn and the status does not report them.

## Automation surface (`getWidgetStatus`)

`parts`: `root`, and `canvas` — the echarts canvas, absent while the message is shown. Areas are in
the canvas's pixel space. A branch is addressed by the node names from the root down, root label
included, joined with ` | `: `branch All | false | F | Asian`. An empty value is `(empty)`.

| Area | What |
|---|---|
| `view` | the chart canvas |
| `branch <path>` | a 4 px square on the line into that node — the click target that selects its rows. The point is taken along the line (middle first) and kept only where zrender's own `findHover` hits that line, not a label or another line; a branch with no such point has no area. The root has no line into it and no area |

| Reading | What |
|---|---|
| `branches` | nodes drawn, the root excluded |
| `rows of branch <path>` | the row count the node represents (the series value) |
| `rows shown` | `filter.trueCount` — the viewer's own combined filter |

## `isRenderPending` and `onRendered`

Both come from `RenderSignals` (`utils/utils.ts`), capped at 120 frames here. `render` is queued;
`_render` re-creates the chart and ends at `setOption(..., lazyUpdate: true)`, and every branch
repaint (`applySelectionFilterChange`) is another lazy `setOption`. The tree then animates its
nodes and lines into place for `animationDurationUpdate` (500 ms). The flag covers:

- the queued renders, settling on a failed one too — so one failure no longer stalls every later render;
- the lazy frame and the animation (`zr.animation.isFinished()`);
- the 50 ms data debounce and the 10 ms reset-filter debounce, from the first event of a burst;
- a resize, from the size event through the frame it is applied on;
- a fold or unfold: echarts applies it before the viewer's click handler runs, and the folded
  branches animate out — they are already gone from `branches` and the areas.
