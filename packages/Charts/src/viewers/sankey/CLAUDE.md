# Sankey

A d3-sankey diagram in svg: one node per distinct source or target value, one link per source →
target pair (per pair and color with Link Color set). Registered by Charts as `Sankey`, so a feature reaches it as
**`sankey viewer`**.

| File | Purpose |
|---|---|
| `sankey.ts` | The viewer: column pick, the graph, the svg, selection by click, the automation signals |
| `sankey-status.ts` | The areas and readings `getWidgetStatus` reports |

## What it draws

`prepareData` walks the viewer's filter and draws one link per source → target pair (one per pair
and color with Link Color set), its Value summed over the pair's rows: d3-sankey's layout grows
faster than the links, and a link per row froze a table of 10,000 rows for seconds. A plain string
column is read as its category indexes and an int or double Value as its raw array; any other column
(virtual, a changed type) through `getString`/`get`. Hovering or clicking a link covers its pair's
rows of its color, a node every row that has it as source or target. A graph with a cycle, or a table without two
string columns of at most 50 categories and a numeric one, shows a message instead, and so does a
Value with a negative value or with only zeros among the links (d3-sankey has no meaningful layout
for either; only zeros gives NaN coordinates). A row whose
Value is missing weighs nothing (read raw, the int null would sink both of its nodes). A graph with
no links — no rows, or none with both a source and a target — draws nothing: d3-sankey sizes the
nodes by their links.

## Automation surface (`getWidgetStatus`)

There is no canvas: `parts.root` is the svg (the viewer root when there is none), and every area is
in px of the svg's box. Everything is read from the svg the last render built, so an empty or
message state reports `0`/`''` and no area of an earlier frame.

| Area | What |
|---|---|
| `node <name>` | the node's rect |
| `link <source> -> <target>` | a 4 px square on the topmost link of that pair, at a point along it (middle first) where `elementFromPoint` finds that pair's link, not a crossing one |

| Reading | What |
|---|---|
| `nodes` | node rects drawn |
| `node names` | their names, comma-separated, in drawn order |
| `links` | link paths drawn — one per pair, or per pair and color with Link Color set |
| `rows shown` | `filter.trueCount` — the viewer's own combined filter |

## `isRenderPending` and `onRendered`

Both come from `RenderSignals` (`utils/utils.ts`). `render` rebuilds the svg synchronously and fires
`onRendered` on its way out, on every path. The selection and size subscriptions are debounced 50
ms; each holds the flag from the first event of a burst until its render, whatever renders between.
