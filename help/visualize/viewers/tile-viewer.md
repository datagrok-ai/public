---
title: "Tile viewer"
---

Visualizes rows as a collection of forms that are positioned as tiles.

To edit the form, right-click and select "Edit form...".

![Tile viewer](../../uploads/gifs/tile-viewer.gif "Tile Viewer")

## Videos

[![Tile Viewer](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=3199s)


## Properties

| Property | Type | Description |
|----------|------|-------------|
| **General** | | |
| Column Name | string |  |
| Header | boolean |  |
| Visible | boolean |  |
| Show Column Name | boolean |  |
| Pos | number |  |
| Lanes Column Name | string |  |
| Card Markup | string |  |
| Allow Drag Between Lanes | boolean |  |
| Auto Generate | boolean | Whether the form auto-generates whenever columns change |
| Lanes | list |  |
| Row Source | string | Determines the rows shown on the plot. |
| Allow Dynamic Menus | boolean |  |
| Show Context Menu | boolean | Properties common for all viewers todo: use code generation |
| Title | string |  |
| Description | string | Viewer description that gets shown at the *Descriptor Position*. Markup is supported. |
| Help | string | Help to be shown when user clicks on the ''?'' icon on top. Could either be in markdown, or a URL (starting with ''/'' or ''http''). |
| Description Position | flexposition |  |
| Description Visibility Mode | visibilitymode |  |
| **Style** | | |
| Tiles Font | string |  |
| **Data** | | |
| Filter | string | Formula that filters out rows to show. Examples: `${AGE}` > 20 or `${WEIGHT / 2)}` > 100, `${SEVERITY}` == ''Medium'', `${RACE}`.endsWith(''sian'') |
| Table | string |  |
| **Description** | | |
| Show Title | boolean |  |

See also:

* [Form](form.md)
* [JS API: Tile viewer](https://public.datagrok.ai/js/samples/ui/viewers/types/tile-viewer)
* [Community: Visualization-related updates](https://community.datagrok.ai/t/visualization-related-updates/521)