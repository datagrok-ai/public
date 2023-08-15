---
title: "Sunburst"
---

A [sunburst diagram](https://en.wikipedia.org/wiki/Pie_chart#Ring),
also referred to as a Ring Chart, Sunburst Chart, or Multi-level Pie Chart,
is a visual representation used to illustrate hierarchical categories as a series of rings.

The central circle represents the root node,
and the categories within the hierarchy extend outward in rings.

Each ring corresponds to a specific attribute or column.
For example, in the picture below, the inner circle represents the "Sex" column with categories "Male" and "Female".
The next ring represents the boolean "Control" status, indicating whether an individual is in the control group or not.
Finally, the outer circle represents the "Race" column.

![Sunburst viewer](sunburst-viewer.png)

The rings in a sunburst diagram are divided and sliced
based on their hierarchical relationship to the parent slice.
Each sector on the outer rings represents a combination of categories from the inner rings as well.

The size of the sector indicates the relative representation
of this specific combination of categories within the dataset.

In the provided picture, the large orange sector on the outer ring represents the
proportion of individuals who belong to the "Caucasian" race,
have a "False" control status, and are classified as "Male" in terms of sex.

## Hierarchical data exploration

Use the Sunburst viewer to explore the hierarchical nature of your dataset:

* To **zoom** to a category, click it.
* To **zoom out**, click the blue circle in the middle.
* To **reset view**, press Ctrl+Shift+A.

![Sunburst interactive data exploration](sunburst-interactive.gif)

In addition to drilling down into categories, you can select a category by using Ctrl+Click.
When you select a category in the Sunburst viewer,
other viewers within the application will respond to the data selection
and dynamically update their representation.

![Sunburst categories selection](sunburst-categories-selection.gif)

## Creating a Sunburst viewer

To create a **Sunburst** viewer, navigate to the **Main Menu**
and select **Add > Javascript Viewers > Charts > Sunburst**.

When you add a sunburst viewer in Datagrok, the platform selects all columns
with categorical data.

## Configuring a Sunburst viewer

You can set the used columns,
and customize some minor visualization options.
To do that, click the **Gear** icon on top of the viewer and use the **Data**
and **Misc** info pane
on the **Context Panel** to manage the viewer’s settings.

For example, you can:

* Select data source for the viewer via **table** control.
* Set the categorical columns to display via **Hierarchy** control.
* Change the margins via **Top**, **Left**, **Bottom**, and **Right** controls.
* Set the animation speed via **Animation Duration** option.

## Interaction with other viewers

The **Sunburst** viewer responds to data filters and instantly changes the visualization.

Selecting a category in the viewer selects
the corresponding rows in the grid, scatter-plot, and other viewers.

However, selecting rows in the grid or other viewers is reflected on the **Sunburst** viewer
only in case if **_all_** data from the particular subcategory are selected.

The highlighting of selected data similar to the **Pie chart viewer** will be available in further releases.

## Viewer controls

| Action                                   | Control                                       |
|------------------------------------------|-----------------------------------------------|
| Select the category and drill down to it | Click the category                            |
| Return to one level up                   | Click the blue circle in the center of viewer |
| Clear the category selection             | Ctrl+Shift+A                                  |

## See also

* [Viewers](../viewers/viewers.md)
* [Pie Chart](pie-chart.md)
