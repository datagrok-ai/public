# Timelines

An echarts `custom` series that draws one lane per split-by value and, in it, one bar per
start–end interval (or a marker for a point event). Registered by Charts as `Timelines`, so a
feature reaches it as **`timelines viewer`**.

| File | Purpose |
|---|---|
| `timelines-viewer.ts` | The viewer: column pick, the intervals, legend filter, zoom state, the automation signals |
| `timelines-status.ts` | The areas and readings `getWidgetStatus` reports |
| `echarts-options.ts` | The initial option: a category y axis, four dataZooms, the custom series |

## What it draws

`getSeriesData` walks the viewer's filter and merges rows with the same lane, event and interval
into one item; `render` then drops the items whose color category the legend has deselected. The
x axis is zoomed by an inside and a slider dataZoom in `weakFilter` mode (an item is dropped only
when it lies wholly outside the window), the y axis by two in `filter` mode. The mouse wheel over
the plot drives the inside zooms. With no split-by or time column the viewer hides the chart and
shows its message in `titleDiv`.

## Automation surface (`getWidgetStatus`)

`parts`: `root`, and `canvas` — the echarts canvas, absent while the message is shown. Areas are in
the canvas's pixel space.

| Area | What |
|---|---|
| `view` | the plot area — the grid rect, axis labels and legend excluded — so a wheel at its centre zooms |
| `lane <value>` | the box of a lane label the y axis drew; the category axis skips labels that would overlap, and those lanes have no area. A click selects the lane's rows. An empty value is `(empty)` |

| Reading | What |
|---|---|
| `lanes` | lanes inside the y axis' current zoom window |
| `intervals` | intervals and markers the series drew: after the legend filter and both zoom windows |
| `rows shown` | `filter.trueCount` — the rows the viewer draws from |

Lane labels are found as the axis view's `Text` elements AxisBuilder tags `label_<ordinal>`, and the
ordinal is looked up in the axis' categories, so a truncated label still names its full lane.

## `isRenderPending` and `onRendered`

Both come from the `RenderSignals` `EChartViewer` holds (`utils/utils.ts`). `setOption` is not lazy
here, so `render` has laid the chart out by the time it returns; the flag is held through it and
until zrender has painted. The base class's 50 ms selection and data debounces hold it from the
first event of a burst, and its resize until the repaint. A mouse wheel is outside `render`
altogether: the inside zoom dispatches its action throttled (20 ms), so a `mousewheel` holds the
flag until the `dataZoom` event lands, or 100 ms for a wheel that zoomed nothing.
