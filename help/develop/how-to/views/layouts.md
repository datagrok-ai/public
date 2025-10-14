---
title: "Use layouts"
---

[Layouts](../../../visualize/view-layout.md) define the way visualizations are positioned in a table view, allowing to
reuse them for different datasets. Layouts contain viewer settings, positions, and relevant metadata that determine
where a layout can be further suggested.

A layout can be created either manually by adding viewers, setting their properties, and docking them,
or [programmatically](../viewers/manipulate-viewers.md).

Table of contents:

- [Creating layouts](#creating-layouts)
- [Saving and searching](#saving-and-searching)
- [Applying layouts to new data](#applying-layouts-to-new-data)
- [Storing metadata](#storing-metadata)
- [Project layouts](#project-layouts)
- [Saving data with layouts using the REST API](#saving-data-with-layouts-using-the-rest-api)

## Creating layouts

Layouts are created from views, and views, in turn, can be restored from layouts. When you save a layout to a server
repository (by choosing `View | Layout | Save to Gallery` in the top menu or hitting _Ctrl + S_) or download it as a
file (`View | Layout | Download`), internally you are working with objects of `ViewLayout` class. There are several ways
to obtain its instance: saving existing view layouts and constructing them from `JSON` strings. Let's start by getting a
layout of the currently open view:

```js
let layout = grok.shell.v.saveLayout();
```

This is quite explicit, there is just one caveat: this method can only be applied to
[table views](../../../datagrok/navigation/views/table-view.md). The same holds for its counterpart `loadLayout`
method that applies a previously saved layout to the given view. Here is an example:

```js
let view = grok.shell.addTableView(grok.data.demo.demog());
let layout = view.saveLayout();
view.addViewer('Histogram', { value: 'age' });
view.grid.columns.rowHeader.width = 100;
view.loadLayout(layout);
```

In the above code snippet, we modify the view after saving its initial layout, so
`loadLayout` in the last line rolls back these changes.

In addition, there is a way to create a layout from `JSON`. The `JSON` describing it may come from a
[user data storage](https://public.datagrok.ai/js/samples/ui/views/layouts), a file containing a downloaded layout, or
directly from serializing an instance via `layout.toJson()`.

```js
let layoutJson = ''; // a JSON string
let tableId = '';    // e.g. JSON.parse(layoutJson).viewStateMap.tableId

grok.data.openTable(tableId).then(t => {
  let view = grok.shell.addTableView(t);
  view.loadLayout(DG.ViewLayout.fromJson(layoutJson));
});
```

The `JSON` representation of a layout has the following main fields:

```json
{
  "#type": "ViewLayout",
  "viewStateMap": {
    ...
  },
  "columns": [
    ...
  ],
  "name": "LayoutName",
  "friendlyName": "Layout Name",
  "author": {
    ...
  },
  "createdOn": "YYYY-MM-DDTHH:mm:ss.sssZ"
}
```

Some of these attributes are common to all [entities](../../../datagrok/concepts/objects.md), so we will mention only the
layout-specific ones. The `viewStateMap` field contains the essential part describing in which positions viewers are
docked and which viewer properties are set. The metadata needed for further layout application is stored in `columns`.
However, a bare view state is enough to reuse the layout:

```js
let layout = view.saveLayout();
// Reapply the layout after a series of changes
view.loadLayout(DG.ViewLayout.fromViewState(layout.viewState));
```

## Saving and searching

The `grok.dapi.layouts` endpoint provides common functionality inherited from
[HttpDataSource](https://datagrok.ai/api/js/api/dg/classes/HttpDataSource) that is responsible for handling collections of
entities stored on the server. Developers can save layouts, find them by id, filter the list of entities according
to [certain criteria](../../../datagrok/navigation/views/table-view.md#search), and so on.

```js
grok.dapi.layouts.list().then(layouts => grok.shell.info(`Total: ${layouts.length}`));
```

## Applying layouts to new data

Since layouts are designed to be reusable, it is essential to determine whether they can be applied to a new dataset.
There is a special method that finds a list of appropriate layouts for a given table:

```js
let df = grok.data.demo.demog();
let view = grok.shell.addTableView(df);
grok.dapi.layouts.getApplicable(df).then(layouts => view.loadLayout(layouts[0]));
```

This method checks whether all the columns the layout was originally applied to can be mapped to the columns of the
specified table. The matching mechanism consists of the following steps:

1. Column names and column types match
2. Both columns have the same [layout-id](../../../govern/catalog/tags.md#layout-id)
3. Both columns have the same [semantic type](../../../govern/catalog/tags.md#quality)

## Storing metadata

Layouts remember the columns they were constructed from (the field `columns`) and the metadata associated with those
columns so that they can be applied to the original or similar data. In particular, metadata includes column names,
types, and tags (such as `layout-id` or `quality`). This is what the mapping rules described above rely on. Apart from
that, layouts as
[entities](../../../datagrok/concepts/objects.md) are capable of storing metadata in a form of properties.

See also:

- [Upload data with layouts using the server API](../data/upload-data.md#layout)
- [View layout](../../../visualize/view-layout.md)
- [Table view](../../../datagrok/navigation/views/table-view.md)
- [User data storage](../data/user-settings-storage.md)
- [JavaScript API Samples: Layout permissions and metadata](https://public.datagrok.ai/js/samples/dapi/layouts-and-permissions)
- [JavaScript API Samples: Saving layouts to user data storage](https://public.datagrok.ai/js/samples/ui/views/layouts)

## Project layouts

Project layout is a combination of layouts for views stored in the project. It can be accessed manually from right-click in **Browse** > **_Path to table_** > **Export layout**. 

## Saving data with layouts using the REST API

To create a dashboard consisting of a dataset that resides externally, and a pre-created layout
(common case for visualizing a dataset created as a result of a data pipeline), use
the [data upload API](../../packages/rest-api.md).
