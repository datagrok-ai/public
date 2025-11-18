---
title: "Pareto front viewer"
---

The Pareto front viewer visualizes a graphical representation of a [Pareto front](https://en.wikipedia.org/wiki/Pareto_front) for [multi-objective optimization](https://en.wikipedia.org/wiki/Multi-objective_optimization). Each point represents a solution evaluated across several objectives, highlighting the trade-offs and non-dominated set of optimal candidates.

![Pareto Front](img/pareto-front.png "Pareto Front")

>Note: To use the Pareto front viewer, install the package
[EDA](https://github.com/datagrok-ai/public/tree/master/packages/EDA).

## Add a Pareto front viewer

1. On the menu ribbon, click the **Add viewer** icon. A dialog opens.
2. In the dialog, select **Pareto Front**.

> Developers: To add the viewer from the console, use:
 `grok.shell.tv.addViewer('Pareto front');`

Initially, the viewer selects two numerical columns as objectives to be minimized. The resulting Pareto-optimal points are displayed on the scatter plot in green. By default, the viewer uses:

* The columns corresponding to the optimized features for the axes
* Category columns with unique values for the labels

## Configuring the Pareto front viewer

To configure the Pareto front, click the **Gear** icon on top of the viewer and use
the info panels on the **Context Panel**. For example, you can:

* **Specify the columns with objectives to minimize** using the `Minimize` property.
* **Select the columns representing objectives to maximize** using the `Maximize` property.
* **Customize the chart axes** using the `X Axis` and `Y Axis` properties.
* **Customize the label columns** to display next to the markers using the `Label Columns` property.

## Interactivity

The Pareto front viewer responds to row selection and data filtering. When you change the target columns, it performs multi-objective optimization and displays the resulting points on the chart.

![Word Cloud](img/pareto-front-viewer.gif)

## See also

* [Viewers](../viewers/viewers.md)
* [Scatter plot](../viewers/scatter-plot.md)
