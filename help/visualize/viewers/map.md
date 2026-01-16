---
title: "Map viewer"
---

Map viewer shows geospatial data on a map as either markers, or a heat map.

## Adding

1. Go to **Tables** and expand the **Viewers** panel.
1. Locate the **Map Viewer** icon and click it.

> Developers: To add the viewer from the console, use:
`grok.shell.tv.addViewer('Map');`

When you add a map viewer, it tries to automatically detect columns that contain
longitude and latitude values.

## Settings

To configure a map, click the **Gear** icon on top of the viewer and use the
info panels on the **Context Pane**. For example, you can:

* *Color-code* points using the `Color` property
* *Size-code* points using the `Size` property
* *Control point size* by using `Marker Min Size` and `Marker Max Size` properties

## Interactivity

* Shows only filtered rows
* Shows selected rows in orange
* Synchronizes current record upon clicking on the point

![Map viewer](img/map-viewer.gif)

## Custom file viewers

Сustom file viewers give user the opportunity to view files of various geographic extensions.
To see it, open the file manager, hover your mouse over the file and click on it, as shown below.

![Custom file viewer](img/map-custom-file-viewer.gif)

The system provides interaction with files of the following extensions:

* .geojson
* .topojson
* .kmz
* .kml

## Controls

|Action              |        Control                |
|------------------------|----------------------|
| Zoom in                                            | Mouse Wheel Up or Plus          |
| Zoom out                                         | Mouse Wheel Down or Minus  |
| Add a point to selection                | Shift+Click the point                   |
| Select multiple points                    | Ctrl+Mouse Drag                       |

## See also

* [Viewers](../viewers/viewers.md)
* [Globe](globe.md)
* [Community: Visualization-related updates](https://community.datagrok.ai/t/visualization-related-updates/521)

## Properties

| Property | Type | Description |
|----------|------|-------------|
| **General** | | |
| Region Column Name | string |  |
| Color Column Name | string |  |
| Color Aggr Type | string |  |
| Show Color Scale | boolean |  |
| Show Color Selector | boolean |  |
| Allow Pan Zoom | boolean |  |
| Linear Color Scheme | list |  |
| Categorical Color Scheme | list |  |
| Row Source | string | Determines the rows shown on the plot. |
| Allow Dynamic Menus | boolean |  |
| Show Context Menu | boolean | Properties common for all viewers todo: use code generation |
| Title | string |  |
| Description | string | Viewer description that gets shown at the *Descriptor Position*. Markup is supported. |
| Help | string | Help to be shown when user clicks on the ''?'' icon on top. Could either be in markdown, or a URL (starting with ''/'' or ''http''). |
| Description Position | flexposition |  |
| Description Visibility Mode | visibilitymode |  |
| **Data** | | |
| Filter | string | Formula that filters out rows to show. Examples: `${AGE}` > 20 or `${WEIGHT / 2)}` > 100, `${SEVERITY}` == ''Medium'', `${RACE}`.endsWith(''sian'') |
| Table | string |  |
| **Style** | | |
| Controls Font | string | Viewer controls elements font. |
| **Description** | | |
| Show Title | boolean |  |

## See also

* [Viewers](../viewers/viewers.md)
* [Globe](globe.md)
* [Community: Visualization-related updates](https://community.datagrok.ai/t/visualization-related-updates/521)