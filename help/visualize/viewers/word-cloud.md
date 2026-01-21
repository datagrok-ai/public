---
title: "Word cloud"
---

Word cloud (a tag cloud) shows the frequency of individual words using font size
and color. Use it to see the most popular concepts, highlight important textual
data points, compare data, etc.

![Word Cloud](img/word-cloud.png "Word Cloud")

## Adding

1. Go to **Tables** and expand the **Viewers** panel.
1. Locate the **Word Cloud** icon and click it.

Another variant:

1. Go to **Add viewer**
1. Locate the **Word cloud** viewer and click it.

Initially, the viewer picks up the first string column in the corresponding
table and builds a word cloud.

> Developers: To add the viewer from the console, use:
 `grok.shell.tv.addViewer('Word cloud');`

## Settings

To configure a word cloud, click the **Gear** icon on top of the viewer and use
the info panels on the **Context Panel**. For example, you can:

* **Select the word column** using the `Column` property
* **Select the text size** using the `Min Text Size` and `Max Text Size` properties
* **Select the rotation degree** using the `Min Rotation Degree` and `Max Rotation Degree`
properties
* **Limit the drawing out of bound** using the **Draw Out Of Bound** setting.

## Interactivity

A word cloud viewer doesn't respond to the row selection and data filtering. You
can use it to filter other viewers.

![Word Cloud](img/word-cloud.gif)

## See also

* [Viewers](../viewers/viewers.md)
* [Table View](../table-view-1.md)
* [Community: Visualization-related updates](https://community.datagrok.ai/t/visualization-related-updates/521)
