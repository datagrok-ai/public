<!-- TITLE: Markup Viewer -->
<!-- SUBTITLE: -->

# Markup Viewer

Use this viewer to host any text, arbitrary HTML content, or [markdown-formatted text](../features/markdown.md). In most casees,
the viewer will auto-detect content type. Use the "mode" property to explicitly specify it.

Properties:

|                     |         |
|---------------------|---------|
| Content             |     |
| Mode                | Text, Html, Markdown, or Auto |
| Markup Enabled      | When true, the rendered HTML is processed by the [Markup](../features/markup.md) engine |

Context menu:

|                       |                 |
|-----------------------|-----------------|
| Edit content...       | Opens a dialog for editing viewer's content.   |


![Markup Viewer](markup-viewer.png "Markup Viewer") 

Here is how to embed iframes:

![Markup Viewer](markup-iframe-embedding.gif "iframe embedding") 

See also: 
  
  * [Viewers](../viewers/viewers.md)
  * [Table View](../views/table-view.md)
  * [Flex View](../views/flex-view.md)
  * [JS API: Markup](https://public.datagrok.ai/js/samples/ui/viewers/markup)
