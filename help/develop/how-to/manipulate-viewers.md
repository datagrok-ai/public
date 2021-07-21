<!-- TITLE: Manipulate Viewers -->
<!-- SUBTITLE: -->

# API for Manipulating Viewers

[Viewers](../../visualize/viewers.md) are main visual components of the Datagrok platform. Our JavaScript API exposes functionality for manipulating native viewers, such as [scatter plot](../../visualize/viewers/scatter-plot.md) or [histogram](../../visualize/viewers/histogram.md), as well as for [developing custom viewers](develop-custom-viewer.md).

Table of contents:
  - [Adding Viewers to Table Views](#adding-viewers-to-table-views)
  - [Adding Viewers to Views](#adding-viewers-to-views)
  - [Working with Properties](#working-with-properties)
  - [Docking Viewers](#docking-viewers)

## Adding Viewers to Table Views

There are several options for attaching a viewer instance to your layout. First of all, you can add it directly to a table view:

```javascript
let data = grok.data.demo.demog();
let view = grok.shell.addTableView(data);
view.addViewer('Histogram', { value: 'age' });
```

To avoid hardcoding the name of a viewer, you can reach it from the `DG` namespace, e.g., `DG.VIEWER.HISTOGRAM` or `DG.VIEWER.PIE_CHART`. For native viewers, Datagrok's API provides a bunch of handy wrapper methods, so this code snippet is equivalent to the one given above:

```javascript
let data = grok.data.demo.demog();
let view = grok.shell.addTableView(data);
view.histogram({ value: 'age' });
```

In both cases, the `options` parameter is not required, besides, it is possible to specify them later on with the viewer's [setOptions](https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot) method. The `addViewer` method differs in that it can be used with all viewers, both native and custom-built ones:

```javascript
view.addViewer('Leaflet').setOptions({
    latitudeColumnName: 'lat',
    longitudeColumnName: 'lon',
    renderType: 'heat map'
});
```

Another way to add a viewer starts with creating it by viewer type:

```javascript
let data = grok.data.demo.demog();
let view = grok.shell.addTableView(data);
let viewer = DG.Viewer.fromType(DG.VIEWER.LINE_CHART, data);
view.addViewer(viewer);
```

Notice that the first parameter to `addViewer` may be either a string with the corresponding name or an instance of the `Viewer` class.

A dataframe on which the view is built has handy plotting methods of its own:

```javascript
let data = grok.data.demo.demog();
let view = grok.shell.addTableView(data);
let viewer = data.plot.scatter({ x: 'weight', y: 'height', color: 'disease' });
view.addViewer(viewer);
```

There is a slight difference between native and custom viewers. `DG.Viewer.fromType` is reserved for
the standard viewer types. When creating an instance of the custom one, use `DataFrame.plot.fromType`
instead. Instances of `JsViewer`'s subclass are obtained asynchronously on the client, as their package
should be initialized first. As a rule of thumb, you should wait until a JsViewer object is constructed.
However, if your plans are limited to adding it to a table view without further adjustments,
a synchronous call of `addViewer` will do.

Examples:

<ul style="column-count: 2;">
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/bar-chart" target="_blank">Bar Chart</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/markup" target="_blank">Markup</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/box-plot" target="_blank">Box Plot</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/matrix-plot" target="_blank">Matrix Plot</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/calendar" target="_blank">Calendar</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/network-diagram" target="_blank">Network Diagram</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/corr-plot" target="_blank">Correlation Plot</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/pc-plot" target="_blank">PC Plot</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/density-plot" target="_blank">Density Plot</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/pie-chart" target="_blank">Pie Chart</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/filters" target="_blank">Filters</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot" target="_blank">Scatter Plot</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/form" target="_blank">Form Viewer</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot-3d" target="_blank">3D Scatter Plot</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/globe" target="_blank">Globe Viewer</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/shape-map" target="_blank">Shape Map</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/google-map" target="_blank">Google Map</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/statistics" target="_blank">Statistics</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/grid" target="_blank">Grid</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/tile-viewer" target="_blank">Tile Viewer</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/heat-map" target="_blank">Heat Map</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/tree-map" target="_blank">Tree Map</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/histogram" target="_blank">Histogram</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/trellis-plot" target="_blank">Trellis Plot</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/line-chart" target="_blank">Line Chart</a>
  </li>
  <li>
    <a href="https://public.datagrok.ai/js/samples/ui/viewers/types/word-cloud" target="_blank">Word Cloud</a>
  </li>
</ul>

## Adding Viewers to Views

You can add a viewer to any container, such as a [main application view](/develop/how-to/build-an-app.md#the-main-view) or a [dialog window](https://public.datagrok.ai/js/samples/ui/dialogs/dialogs).

Here is how to create a new [view](develop/how-to/custom-views.md) and add a plot to it:

```javascript
let v = DG.Viewer.scatterPlot(grok.data.demo.demog());
grok.shell.newView('foo').append(v.root);
```

Here is an example of adding controls and viewers to a dialog box:

```javascript
ui.dialog()
  .add(ui.div([
    "Here's a plot in the dialog",
    DG.Viewer.scatterPlot(table)
  ]))
.show();
```

Check the [ChaRPy](https://github.com/datagrok-ai/public/blob/master/packages/ChaRPy/src/package.js) package for a similar example.

## Docking Viewers

Just like other visual components that occupy a window of their own, viewers can be docked to a particular position:

```javascript
let data = grok.data.demo.demog();
let view = grok.shell.addTableView(data);
let viewer = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, data);
view.addViewer(viewer);
view.dockManager.dock(viewer, 'right');
```

The list of positions consists of the following options: `left | right | top | down | fill`. You can refer to them from the `DG` namespace, e.g., `DG.DOCK_TYPE.RIGHT`. Notice that here the viewer will be placed independently of the table view:

```javascript
grok.shell.dockManager.dock(viewer, DG.DOCK_TYPE.RIGHT);
```

Examples:
  * [Docking Viewers](https://public.datagrok.ai/js/samples/ui/docking/docking-table-view)

See also:
  * [Viewers](../../visualize/viewers.md)
  * [How to Develop Custom Viewer](develop-custom-viewer.md)
