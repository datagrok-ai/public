---
title: "Scatter plot"
---

Scatter plot displays the relationship between two variables using Cartesian
coordinates. Data is shown as points on the plot, positioned according to the
values of the two variables being studied. To show the relationship between
three variables, use [3D Scatter plot](3d-scatter-plot.md).

You can make a scatter plot more informative by using colors, sizes, and marker
shapes. In the example below, you see additional insights about `AGE`,
`DIS_POP`, and `SEX`.

![Scatter plot](scatter-plot.png)

## Configuring a scatter plot

To configure a scatter plot, use the selectors on the viewer or the settings on
the **Context Pane**. To do that, click the **Gear** icon on top of the viewer
and manage the scatter plot’s settings.

### Filtering

A scatter plot can respond to filters and apply filters to other viewers (**Context Pane** > **Data**):

* Use **Row Source** to define the rows to display on the scatter plot. Choices are: `All`, `Filtered`, `Selected`, `SelectedOrCurrent`, `FilteredSelected`, `MouseOverGroup`, `CurrentRow`, `MouseOverRow`.
* Apply **Zoom and Filter** to define how the scatter plot responds during zooming. Сhoices include: `zoom by filter`, `filter by zoom`, `pack and zoom by filter`.

<!--<img alt="Filtering" src={require('./filtering.gif').default} width="800px"/>-->
![Filtering](./filtering.gif "Filtering")

### Axes

You can customize scatter plot axes using settings from the **X** and **Y** info panels: 

* Apply **X** (**Y**) **Axis Type** to switch between `linear` and `logarithmic`
  scales. Log scale is aplicable for positive numeric values only.
* Set **Min** and **Max** values for the axes.
* Enable **Invert X** (**Y**) **Axis**  to reverse the direction of the axis on
  a viewer. 

>Note: The key settings from the **Context Panel** are replicated in the context
>menu for scatter plot axes.

#### Molecular rendering

A scatter plot renders molecular structures on the axes. When you hover over a structure, a tick mark shows its position on the axis, and its label appears above the others on the axis.

<!--<img alt="Molecular structures rendering" src={require('./rendering.gif').default} width="800px"/>-->
![Molecular structures rendering](./rendering.gif "Molecular structures rendering")

### Selection

A scatter plot offers two data selection modes: **Rectangle Marquee** (by default)
and **Lasso** (**L** or **Misc** > **Lasso tool**).  By default, to select points on a
scatter plot, you need to use Shift+Mouse Drag. To enable selecting an area
using only mouse drag, switch the **Mouse Drag** setting to `select`. 

<!--<img alt="Selection" src={require('./selection.gif').default}width="800px"/>-->
![Selection](./selection.gif "Selection")

### Tooltip

By default, a scatter plot inherits a tooltip from the grid. You can set a
custom tooltip for the scatter plot  from the **Tooltip** info panel:

1. Set **Show tooltip** to `show custom tooltip`.
1. In **Row Tooltip**, select the desired columns.
1. Select the suitable choice for **Data Values** to combine values from the
   axes and **Data Values** in the tooltip.

>Note: The key settings from the **Context Panel** are replicated in the context
>menu for the scatter plot (**Context Menu > Tooltip**).

<!--<img alt="Tooltip" src={require('./tooltip.gif').default} width="800px"/>-->
![Tooltip](./tooltip.gif "Tooltip")

Learn also about [Goup tooltip](https://datagrok.ai/help/visualize/viewers/#group-tooltips)

## Advanced features

### Formula lines

Scatter plot supports formula lines (COntext Menu > Tools > Formula Lines…). These lines enable you to visually represent mathematical formulas or equations on the plot. This can help illustrate trends, relationships, or patterns between variables in the data. 

>Developers: Learn more about how to [Show formula lines](https://datagrok.ai/help/develop/how-to/show-formula-lines)

<!--<img alt="Formula lines" src={require('./formula-lines.gif').default} width="800px"/>-->
![Formula lines](./formula-lines.gif "Formula lines")

### Regression line

Hit 'R' to toggle the regression line.

![Scatter plot](../../uploads/gifs/scatter-plot.gif "scatter plot")

## Viewer controls

| Action                 | Control              |
|------------------------|----------------------|
| Zoom                   | Alt+Mouse Drag       |
| Zoom in  | Mouse Wheel Up or Plus |
| Zoom out  | Mouse Wheel Down or Minus  |
| Select | Shift+Mouse Drag       |
| Invert selected | Ctrl+Mouse Click        |
| Scroll  |  Up, Down, Left, Right       |
| Toggle lasso tool | L        |
| Toggle regression line | R       |
| Show in full screen | Alt+F        |

## Videos

[![Scatter Plot](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=214s)

See also:

* [Column selectors](column-selectors.md)
* [Viewers](../viewers/viewers.md)
* [Table view](../../datagrok/navigation/table-view.md)
* [JS API: Scatter plot](https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot)
