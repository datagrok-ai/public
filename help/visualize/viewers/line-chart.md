---
title: "Line chart"
---

Line chart shows data points as connected line segments. It is commonly used to track trends, changes over time, and compare multiple data series. 

> Developers: To add the viewer from the console, use:
`grok.shell.tv.addViewer('Line chart');`

General:

|                       |                 |
|-----------------------|-----------------|
| Right click           | Context menu    |
| Alt+Mouse Drag        | Zoom            |
| Alt+F                 | Show in full screen |
| Shift+Mouse Drag      | Select          |
| Up, Down, Left, Right | Scroll          |
| Ctrl -                | Zoom out        |
| Ctrl +                | Zoom in         |
| Shift+drag axes       | Select segments |

![Line Chart](../../uploads/gifs/line-chart.gif "Line chart")

## Custom aggregated tooltip

The line chart supports tooltips for both single and aggregated values. The default [single-value tooltip](../table-view-1.md#tooltips) is inherited from the grid. The default aggregated tooltip displays axis values and can be customized to show additional summary statistics.

To configure a custom aggregated tooltip, go to **Context menu > Tooltip > Edit**. This opens the aggregated tooltip dialog, where you can select the columns and corresponding aggregation functions. Available aggregation functions depend on the column type.

<details> 
   <summary>Supported aggregation functions</summary>

| Column type   | Aggregation functions                 |
|---------------|--------------------------------------|
| Numerical   | `min`, `max`, `sum`, `avg`, `med`, `geomean`,<br></br> `count`, `values`, `unique`, `nulls`,<br></br>`stdev`, `variance`, `skew`, `kurt`,<br></br>`q1`, `q2`, `q3`, `first`|
| Categorical     | `first`, `count`, `values`, `unique,nulls`, `values or unique count`,<br></br> `concat all`, `concat unique`,<br></br>  l`ongest`, `shortest`,<br></br> `most frequent`, `concat counts`  |
| Date          | `first`, `count`, `values`, `unique`, `nulls`,<br></br>`min`, `max`, `avg`,<br></br>  `range`                    |
</details> 

![](img/line-chart-aggregated-tooltip.gif)

## Statistical Process Control

Line chart supports Statistical Process Control (SPC) features out of the box, including contol limits and adjustable rules.
![SPC Chart](../../uploads/gifs/line-chart-spc.png "SPC Chart")

## Videos

[![Line Chart](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=934s)

See also:

* [Viewers](../viewers/viewers.md)
* [Table View](../table-view-1.md)
* [JS API: Line chart](https://public.datagrok.ai/js/samples/ui/viewers/types/line-chart)
* [Community: Visualization-related updates](https://community.datagrok.ai/t/visualization-related-updates/521)
