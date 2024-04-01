---
title: "Markup viewer"
---

Use this viewer to host any text, arbitrary HTML content, or [Markdown-formatted text](https://en.wikipedia.org/wiki/Markdown). In
most casees, the viewer will auto-detect content type. Use the "mode" property to explicitly specify it.

Properties:

|                     |         |
|---------------------|---------|
| Content             |     |
| Mode                | Text, Html, Markdown, or Auto |
| Markup Enabled      | When true, the rendered HTML is processed by the [Markup](../../datagrok/navigation/markup.md) engine |

Context menu:

|                       |                 |
|-----------------------|-----------------|
| Edit content...       | Opens a dialog for editing viewer's content.   |

![Markup Viewer](img/markup-viewer.png "Markup Viewer")

Here is how to embed iframes:

![Markup Viewer](img/markup-iframe-embedding.gif "iframe embedding")

## Videos

[![Markup Viewer](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=3052s)

See also:

* [Viewers](../viewers/viewers.md)
* [Table view](../../datagrok/navigation/table-view.md)
* [Flex view](../../datagrok/navigation/flex-view.md)
* [JS API: Markup](https://public.datagrok.ai/js/samples/ui/viewers/types/markup)
