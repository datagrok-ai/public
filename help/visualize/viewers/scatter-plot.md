---
title: "Scatterplot"
format: mdx
---

A scatterplot displays data points on the horizontal (X) and vertical (Y) axes
to show the relationship between two variables. By using marker color, shape,
and size, you can show up to three additional data dimensions. Use the
scatterplot to explore patterns and relationships between variables in your
data.

A scatterplot is a [chemically-aware viewer](../../datagrok/solutions/domains/chem/chemically-aware-viewers#scatter-plot)
 and can be used to explore a chemical space.

:::tip

To show the relationship between three variables, use a [3D Scatterplot](3d-scatter-plot.md).

::: 

## Controls

|                        |                      |
|------------------------|----------------------|
| Context menu           | Right-click         |
| Zoom                   | Alt+Mouse Drag       |
| Zoom in                | Mouse Wheel Up or Plus|
| Zoom out               | Mouse Wheel Down or Minus|
| Double-click           | Reset view             |
| Select                 | Shift+Mouse Drag       |
| Invert selected        | Ctrl+Mouse Click       |
| Scroll                 |  Up, Down, Left, Right |
| Toggle lasso tool      | L        |
| Toggle regression line | R       |
| Show in full screen    | Alt+F        |

![Scatterplot](../../uploads/gifs/scatter-plot.gif)

## Adding and configuring a scatterplot

To add a scatterplot, click the **Scatterplot** icon on the **Toolbox**.

Use the viewer controls to select columns for each axis and marker color and
size. For additional configurations, click the **Gear** icon on top of the
viewer and set your preferences in the **Context Panel**. Here, you can change
the background, adjust the legend position, add or remove labels, show drop
lines, and more. You can also access key settings from the context menu by
right-clicking.

<img alt="Molecular structures rendering" src={require('./img/Axes.gif').default}
width="800px"/>

:::note developers

To add a scatterplot from the **Console**, use
`grok.shell.tv.addViewer('Scatterplot');`

:::

### Calculations and trends

A scatterplot can show reference lines that represent formulas or equations.
These lines are used to emphasize specific areas on the chart or data. Common
examples include a regression line, value bands, and so on. 

To toggle a regression line, press the **R** key. 

To show a custom formula line, right-click a scatterplot, then choose **Tools**
> **Formula Lines...** This action opens a **Formula Lines** dialog. Here, enter
your formula and configure the line settings. Your formula should refer to the
columns on the **X** and **Y** axes. The syntax for the formula is similar to
that used to [Add New Column](../../transform/add-new-column.md).

<img alt="Formula lines" src={require('./img/formula-lines.gif').default}
width="800px"/>

:::note developers

You can [add formula lines programmatically](https://datagrok.ai/help/develop/how-to/show-formula-lines).

:::

### Filtering and selection

To adjust the data display settings, use the **Data** info pane from the
**Context Panel**. Here, you can specify which rows to show on the scatterplot,
define zoom and filtering actions, and so on.

For manual selection settings, use the **Misc** info pane. Here, you can choose
the selector you want and define actions for the mouse drag. 

### Tooltip

By default, a scatterplot inherits the tooltip from the grid. However, you can
customize the scatterplot's tooltip to show the data you want using the
**Tooltip** info panel or via the context menu.

In addition, a scatterplot itself can be used as a 
[group tooltip](../../datagrok/navigation/views/table-view.md#group-tooltips), which may be especially useful when 
dealing with grouped or clustered data or when the screen space is limited.

![Group Tooltip](../../uploads/viewers/viewer-group-tooltip.png "Group Tooltip")

To use this feature: 

1. Add and configure a scatterplot.
1. From the context menu, select **Tooltip** > **Use as Group Tooltip**.
1. Optional. Close the scatterplot.

### Connecting lines

You can set a column that defines order in which points are connected. 
Below, we see the (gdp, life expectancy) trajectory of different countries over time.

![](img/scatter-plot-lines.png)

## Videos

[![ScatterPlot](../../uploads/youtube/visualizations2.png "Open on
Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=214s)

See also:

* [Column selectors](column-selectors.md)
* [Viewers](viewers.md)
* [Table view](../../datagrok/navigation/views/table-view.md)
* [JS API:
  Scatterplot](https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot)
