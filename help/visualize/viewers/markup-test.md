<!-- TITLE: Tests: Markup viewer -->
<!-- SUBTITLE: -->

# Tests: Markup

Use this viewer to host any text, arbitrary HTML content, or [markdown-formatted text](../../datagrok/markdown.md). In
most casees, the viewer will auto-detect content type. Use the "mode" property to explicitly specify it.

## Testing scenario

1. Open "demog" dataset

1. Add [Markup viewer](markup.md) (from "Viewers" section
   on [Toolbox](../../datagrok/navigation.md#toolbox), "Add" menu or list on Toolbar)

* [Markup viewer](markup.md) is added on layout
* Help switched to [Markup viewer](markup.md) page
* Sample text present on [markup](markup.md)

1. Open "Edit content" dialog from [Markup viewer](markup.md) menu

* "Edit" dialog is open
* Sample text is present in the dialog input field

1. Enter [Markup](markup.md) text to the input field *

1. Click "OK" button

* [Markup](markup.md) appeared viewer from an external source
* All external viewer functionality is available inside the platform

1. Add "Markup view" from "+" menu

* New tab with Markup view added
* There is sample text on view
* Markup text shown on [Property Panel](../../datagrok/navigation.md#properties)

1. Enter markup text with the external viewer to [Property Panel](../../datagrok/navigation.md#properties) field *

* Viewer from an external source displayed on markup view
* All external viewer functionality is available inside the platform1. Return to view
  with [Source Tree](../../access/data-source.md)

(*):

```html
<iframe fremeborder="0" id="iframe_opkomst" src="https://dirkmjk.nl/files/articles/2016/opkomst/en.html"
width="100%" height="100%">
</iframe>
```

See also:

* [Viewers](../viewers.md)
* [Viewers test](viewers-test.md)
