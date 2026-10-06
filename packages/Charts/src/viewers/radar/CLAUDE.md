# Radar

An echarts `radar` chart: one axis per value column, one line per row. Registered by Charts as
`Radar`, so a feature reaches it as **`radar viewer`**.

| File | Purpose |
|---|---|
| `radar-viewer.ts` | The viewer: value columns, the row lines, percentile areas, legend, the automation signals |
| `radar-status.ts` | The areas and readings `getWidgetStatus` reports |

## What it draws

Three series: the min and max percentile areas (0 and 1) and the row lines (2). `createSeriesData`
walks the viewer's filter, skips a row whose color category the legend has deselected, and stops
at `MAXIMUM_ROW_NUMBER` (1000) **lines** — the cap applies after the legend, so on a table with more
than 1000 rows of the picked category the chart still draws 1000 lines. The current row and the
mouse-over row are pushed again on top as highlight lines. When the filter passes more than 1000
rows the viewer adds the `Only first 1000 shown` notice above the chart; with no usable value column
it covers the chart with `The Radar viewer requires a minimum of 1 numerical column.`

## Automation surface (`getWidgetStatus`)

`parts`: `root`, and `canvas` — the echarts canvas, absent while the error covers the chart.

| Area | What |
|---|---|
| `view` | the chart canvas |

| Reading | What |
|---|---|
| `axes` | the indicator axes the radar coordinate system was last laid out with, comma-separated, in drawn order |
| `rows shown` | row lines drawn — the rows the legend let through, up to the 1000-line cap. The current and mouse-over highlight lines (drawn on top, with a symbol, legend or not) are not counted |
| `message` | the error covering the chart, else the `Only first 1000 shown` notice, else `''` |

Under the error `axes` is `''` and `rows shown` is `0`: nothing of the previous frame is reported.

## `isRenderPending` and `onRendered`

Both come from the `RenderSignals` `EChartViewer` holds (`utils/utils.ts`). `render` ends at
`setOption(..., lazyUpdate: true)`: the model takes the option at once, the axes and lines are built
on the next zrender frame. The flag is held while `render` runs and then until echarts owes no frame
(`echartsFramePending`: no lazy `setOption` left unapplied, nothing left unpainted), on every path,
the error and a thrown render included. The size subscription is debounced 50 ms and renders on the
frame after, so the flag is held from the first size event of a burst through that render.
