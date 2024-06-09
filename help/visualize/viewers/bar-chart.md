---
title: "Bar chart"
---

A bar chart presents grouped data as rectangular bars with lengths proportional to the values that they represent.
Unlike histograms which you can apply to display the distribution of numerical data, bar charts are primarily designed
for categorical values.

> Developers: To add the viewer from the console, use:
 `grok.shell.tv.addViewer('Bar chart');`

|                |                                   |
|----------------|-----------------------------------|
| Click on a bar | [Select or filter](../viewers/viewers.md) |
| Right click    | Context menu                      |
| Double-click   | Reset View                        |
| Alt+drag       | Zoom                              |
| Alt+F          | Show in full screen               |

## Features

### Stacked bars and relative values

Use the 'Relative Values' property in combination with the 'Stack' property to analyze the distribution of the stacked
values:

![Relative values in a bar chart](img/bar-chart-relative-values.gif "Relative values in a bar chart")


### Dates and years, quarters, months

You can categorize DateTime columns using special functions, such as 'Year', 'Month', 'Quarter', '
Year - Month' and 'Year - Quarter':

![Dates in a bar chart](img/bar-chart-dates.gif "Dates in a bar chart")

## Videos

[![Bar Chart](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=684s)

See also:

* [Column selectors](column-selectors.md)
* [Table View](../../datagrok/navigation/views/table-view.md)
* [Viewers](../viewers/viewers.md)
* [JS API: Bar Chart](https://public.datagrok.ai/js/samples/ui/viewers/types/bar-chart)
