# Band chart

An echarts line chart that shades the area between a lower and an upper column along X — a prediction
or confidence band — with an optional line through it (a median, a fit), one band per category of the
split column. Registered by Charts as `Band chart`, so a feature reaches it as **`band chart viewer`**,
and a u2 view as `viewers.byType('Band chart', …)` or the spec tag `u2-viewer-band-chart`.

| File | Purpose |
|---|---|
| `band-chart-viewer.ts` | The viewer: columns, the stacked band series, the line, tooltip, the automation signals |
| `band-chart-status.ts` | The areas and readings `getWidgetStatus` reports |

## What it draws

`collectBands` walks the viewer's filter, skips a row with X, lower or upper missing, groups the rest by
the split category (one band, named `''`, without a split column) and sorts each group by X. echarts
has no area between two series, and a `stack` on two value axes stacks X, not Y; so each band is a
`custom` series that draws one polygon — the upper edge forward, the lower edge back — through
`api.coord`, which also makes `yAxisType: 'logarithmic'` right. The polygon is drawn on the first item
only; every item carries `[x, lower, upper]` (`encode.y: [1, 2]`) so the band sets the axis extent. The
line column, when set, is its own line series in the band's color. A band's series share its name,
so a legend click hides the band and its line together. The tooltip reads the values from `bands`. A band's colour is the split column's colour coding at one of its rows
(`DG.Color.getRowColor`, as the platform's viewers read a coded column), else the categorical palette in
band order.

Axis titles: `xAxisTitle` / `yAxisTitle`, or when empty the X column name and the line column name (else
`lower – upper`); the column's `units` tag (`col.meta.units`: the X column's, the line's or else the
upper's for Y) follows in brackets (`axisTitle`) and in the tooltip; on `yAxisType: 'logarithmic'` the Y title adds
`log scale` and the axis labels its minor ticks too, so the scale is read from the axis itself. Reference lines are the table's
own formula lines (`df.meta.formulaLines`, as the scatter plot and line chart read them): `refLines`
keeps the visible constant ones (`DG.FormulaLinesHelper.getMeta` names the column, not the constant) — `${X column} = 19` vertical, `${any other column} = 1500` horizontal,
dashed unless `style: 'solid'`, labelled with the line's title — and skips other formulas (they need the
platform's evaluator). They ride an empty scatter series as echarts `markLine`s; that series' radius-0
points keep a horizontal line inside the Y axis. The viewer re-renders on `onMetadataChanged`, since
formula lines and units are tags.

With no X / lower / upper declared, `onTableAttached` takes the first three numerical columns, so the
viewer shows something on any table (the `Viewers` category runs it on demog). Without three usable
columns it covers the chart with `The Band chart needs an X, a lower and an upper numerical column.`

## Automation surface (`getWidgetStatus`)

`parts`: `root`, and `canvas` — the echarts canvas, absent while the error covers the chart.

| Area | What |
|---|---|
| `view` | the chart canvas |

| Reading | What |
|---|---|
| `bands` | the band names drawn, comma-separated, in drawn order (`''` for the single band without a split) |
| `points` | X values drawn, over all bands |
| `reference lines` | the formula lines drawn, `x = 19, y = 1500`, comma-separated |
| `message` | the error covering the chart, else `''` |

## `isRenderPending` and `onRendered`

From the `RenderSignals` `EChartViewer` holds, as in the other echarts viewers: `render` runs `draw`
under `signals.render`, and the data, selection and filter subscriptions are debounced 50 ms through
`signals.debounce`, so the flag is held from the first change of a burst.
