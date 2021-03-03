<!-- TITLE: Use Layouts -->

# Layouts

[Layouts](../../visualize/view-layout.md) capture how visualizations are displayed relative to each
other in a table view. This allows reusing visual templates across different data. Such UI-first
approach, when you sketch out the desired view and save the state to come back to it later, forms a
contrast to the programmatic step-by-step construction of the view in the code (see JavaScript API
for [viewers](manipulate-viewers.md) and [grid](customize-grid.md)). Layouts contain viewer settings,
positions, and relevant metadata that determine where a layout can be further suggested.

## Creating Layouts

Layouts are created from views, and views, in turn, can be restored from layouts. When you save a
layout to a server repository (by choosing `View | Layout | Save to Gallery` in the top menu or
hitting _Ctrl + S_) or download it as a file (`View | Layout | Download`), internally you are
working with objects of `ViewLayout` class. There are several ways to obtain its instance: saving
existing views and constructing layouts from `JSON`. Let's start by getting a layout of the
currently open view:

```js
let layout = grok.shell.v.saveLayout();
```

This is quite explicit, there is just one caveat: this method can only be applied to
[table views](../../overview/table-view.md). The same holds for its counterpart `loadLayout`
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
[user data storage](https://public.datagrok.ai/js/samples/ui/views/layouts), a file containing a
downloaded layout, or from directly serializing an instance via `layout.toJson()`.

```js
let tableId = '';
let layoutJson = '';

grok.data.openTable(tableId).then(t => {
  let view = grok.shell.addTableView(t);
  view.loadLayout(DG.ViewLayout.fromJson(layoutJson));
});
```
<!-- TODO: grok.dapi.layouts entrypoint, applying layouts to new data, storing metadata in layouts -->

See also:
  - [View Layout](../../visualize/view-layout.md)
  - [Table View](../../overview/table-view.md)
  - [User Data Storage](../user-data-storage.md)
  - [JavaScript API Samples: Layout Permissions and Metadata](https://public.datagrok.ai/js/samples/dapi/layouts-and-permissions)
  - [JavaScript API Samples: Saving Layouts to User Data Storage](https://public.datagrok.ai/js/samples/ui/views/layouts)
