---
title: "Matrix plot"
---

Use Matrix Plot to assess the relationship among many pairs of columns at the same time.

> Developers: To add the viewer from the console, use:
`grok.shell.tv.addViewer('Matrix plot');`

General:

|             |                     |
|-------------|---------------------|
| Right click | Context menu        |
| Alt+F       | Show in full screen |

![Matrix Plot](../../uploads/viewers/matrix-plot.png "Matrix Plot")

## Videos

[![Matrix Plot](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=1653s)


## Properties

| Property | Type | Description |
|----------|------|-------------|
| **General** | | |
| X Column Names | list | Columns to use on the X axis |
| Y Column Names | list | Column to use on the Y axis |
| Font | string |  |
| Cell Plot Type | string |  |
| Show X Axes | boolean |  |
| Show Y Axes | boolean |  |
| Back Color | number |  |
| Inner Viewer Look | lookandfeel |  |
| Row Source | string | Determines the rows shown on the plot. |
| Allow Dynamic Menus | boolean |  |
| Show Context Menu | boolean | Properties common for all viewers todo: use code generation |
| Title | string |  |
| Description | string | Viewer description that gets shown at the *Descriptor Position*. Markup is supported. |
| Help | string | Help to be shown when user clicks on the ''?'' icon on top. Could either be in markdown, or a URL (starting with ''/'' or ''http''). |
| Description Position | flexposition |  |
| Description Visibility Mode | visibilitymode |  |
| **Style** | | |
| Auto Layout | boolean |  |
| **Data** | | |
| Filter | string | Formula that filters out rows to show. Examples: `${AGE}` > 20 or `${WEIGHT / 2)}` > 100, `${SEVERITY}` == ''Medium'', `${RACE}`.endsWith(''sian'') |
| Table | string |  |
| **Description** | | |
| Show Title | boolean |  |


See also:

* [Viewers](../viewers/viewers.md)
* [Table view](../table-view-1.md)
* [JS API: Matrix plot](https://public.datagrok.ai/js/samples/ui/viewers/types/matrix-plot)
* [Community: Visualization-related updates](https://community.datagrok.ai/t/visualization-related-updates/521)