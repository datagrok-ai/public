<!-- TITLE: Manipulate Viewers -->
<!-- SUBTITLE: -->

# API for Manipulating Viewers

[Viewers](../../visualize/viewers.md) are main visual components of the Datagrok platform. Our JavaScript API exposes functionality for manipulating native viewers, such as [scatter plot](../../visualize/viewers/scatter-plot.md) or [histogram](../../visualize/viewers/histogram.md), as well as for [developing custom viewers](develop-custom-viewer.md).

## Adding Viewers

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

Examples:
  * [Bar Chart](https://public.datagrok.ai/js/samples/ui/viewers/types/bar-chart)
  * [Box Plot](https://public.datagrok.ai/js/samples/ui/viewers/types/box-plot)
  * [Calendar](https://public.datagrok.ai/js/samples/ui/viewers/types/calendar)
  * [Correlation Plot](https://public.datagrok.ai/js/samples/ui/viewers/types/corr-plot)
  * [Density Plot](https://public.datagrok.ai/js/samples/ui/viewers/types/density-plot)
  * [Filters](https://public.datagrok.ai/js/samples/ui/viewers/types/filters)
  * [Form Viewer](https://public.datagrok.ai/js/samples/ui/viewers/types/form)
  * [Globe Viewer](https://public.datagrok.ai/js/samples/ui/viewers/types/globe)
  * [Google Map](https://public.datagrok.ai/js/samples/ui/viewers/types/google-map)
  * [Grid](https://public.datagrok.ai/js/samples/ui/viewers/types/grid)
  * [Heat Map](https://public.datagrok.ai/js/samples/ui/viewers/types/heat-map)
  * [Histogram](https://public.datagrok.ai/js/samples/ui/viewers/types/histogram)
  * [Line Chart](https://public.datagrok.ai/js/samples/ui/viewers/types/line-chart)
  * [Markup](https://public.datagrok.ai/js/samples/ui/viewers/types/markup)
  * [Matrix Plot](https://public.datagrok.ai/js/samples/ui/viewers/types/matrix-plot)
  * [Network Diagram](https://public.datagrok.ai/js/samples/ui/viewers/types/network-diagram)
  * [PC Plot](https://public.datagrok.ai/js/samples/ui/viewers/types/pc-plot)
  * [Pie Chart](https://public.datagrok.ai/js/samples/ui/viewers/types/pie-chart)
  * [Scatter Plot](https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot)
  * [3D Scatter Plot](https://public.datagrok.ai/js/samples/ui/viewers/types/scatter-plot-3d)
  * [Shape Map](https://public.datagrok.ai/js/samples/ui/viewers/types/shape-map)
  * [Statistics](https://public.datagrok.ai/js/samples/ui/viewers/types/statistics)
  * [Tile Viewer](https://public.datagrok.ai/js/samples/ui/viewers/types/tile-viewer)
  * [Tree Map](https://public.datagrok.ai/js/samples/ui/viewers/types/tree-map)
  * [Trellis Plot](https://public.datagrok.ai/js/samples/ui/viewers/types/trellis-plot)
  * [Word Cloud](https://public.datagrok.ai/js/samples/ui/viewers/types/word-cloud)

## Docking Viewers

Just like other visual components that occupy a window of their own, viewers can be docked to a particular position:

```javascript
let data = grok.data.demo.demog();
let view = grok.shell.addTableView(data);
let viewer = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, data);
view.addViewer(viewer);
view.dockManager.dock(viewer, 'right');
```

The list of positions consists of the following options: `left | right | top | bottom | fill`. Notice that here the viewer will be placed independently of the table view:

```javascript
grok.shell.dockManager.dock(viewer, 'right');
```

Examples:
  * [Docking Viewers](https://public.datagrok.ai/js/samples/ui/docking/docking-table-view)

See also:
  * [Viewers](../../visualize/viewers.md)
  * [How to Develop Custom Viewer](develop-custom-viewer.md)
