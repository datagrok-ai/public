# Sunburst

An echarts `sunburst` series that draws a categorical hierarchy as nested rings: one ring per
hierarchy column, one sector per distinct value under its parent. Registered by Charts as
`Sunburst`, so a feature reaches it as **`sunburst viewer`**.

| File | Purpose |
|---|---|
| `sunburst-viewer.ts` | The viewer: hierarchy pick, the tree it feeds echarts, click-to-select/filter, molecule labels, the three automation signals |
| `sunburst-status.ts` | The areas and readings `getWidgetStatus` reports |

## What it draws

`_render` keeps only the hierarchy columns that are usable — when the viewer's own filter passes
`CATEGORIES_NUMBER` (500) rows, a column with more than 500 categories is dropped — then takes at
most `MAXIMUM_COLUMN_NUMBER` (20) of them. That surviving list is `eligibleHierarchyNames`, and it
is what the status reports as `hierarchy columns`: the property is the request, this is the answer.
With none left the viewer shows a message instead of a chart.

`render` is token-queued: each call takes the next token and chains onto `renderQueue`, and
`_renderWithToken` returns without drawing if a newer token has since been issued. So a burst of
property changes collapses to one frame.

## A sector is addressed by its path

The click handler builds a sector's identity from `params.treePathInfo` minus the root, and the
selection is keyed by that path joined with `|`. The status reports the same names joined with
` | ` instead, for readability: `segment Cancer | Male`. A one-column sunburst therefore names its
sectors by that column alone (`segment Cancer`).

## Automation surface (`getWidgetStatus`)

Built on demand by `sunburst-status.ts`. Every `hitArea` is in the chart canvas's own pixel space.

A sector's area is **not** its bounding box — for a wide sector the box's centre is the hole, not
the ring. It is a 6 px square centred on the sector's mid-radius, mid-angle point, and that point is
offered to the sector's own `contain()` before it is reported. zrender builds a sector with
`cy + r*sin(angle)`, so that candidate goes first and the opposite sign is the fallback; a sector
that accepts neither is left unplaced rather than reported wrongly. This replaces the
radius-by-direction sweep a spec previously needed to land on a ring.

| Area | What |
|---|---|
| `view` | the chart canvas |
| `segment <path>` | a click target inside that sector — `segment Cancer \| Male` |

| Reading | What |
|---|---|
| `hierarchy columns` | the columns that survived the category and count limits, in ring order |
| `on click` | `Select` or `Filter` |
| `include nulls` | the `includeNulls` property |
| `rows shown` | `filter.trueCount` — the viewer's own combined filter |
| `segments` | sectors the layout actually placed |
| `segment names` | their paths, comma-separated |
| `rows of segment <path>` | that sector's value, as the last layout computed it |

`isRenderPending` covers both the queued renders and the frame echarts lays the sectors out on:
`_render` ends at `setOption(..., lazyUpdate: true)`, so the series is built on the next zrender
frame and neither the promise nor the canvas says the sunburst is there. `renderFinished` re-checks
`laidOutSegmentCount` across a few frames (capped, and short-circuited by `renderError`) before
`onRendered` fires — the same shape as the word cloud's. `renderError` reads the message the viewer
put in place of the chart, and while it is set the status reports no areas from the previous frame.

The status reads echarts internals (`getData()`, `getItemGraphicEl`, `tree.getNodeByDataIndex`), so
every step is optional-chained. A sector it cannot name or cannot place contributes nothing, which
is the honest answer — never a guessed rectangle.
