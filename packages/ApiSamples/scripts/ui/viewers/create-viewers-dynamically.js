//api: DG.Viewer.fromType, DG.DataFramePlotHelper.bar, DG.ScatterPlotViewer.zoom
// Creating viewers dynamically: a known DG.VIEWER type returns that viewer's class
// (here ScatterPlotViewer, so `zoom` completes); a type known only at runtime returns DG.Viewer.

let t = grok.data.demo.demog();

let scatterPlot = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, t, {x: 'height', y: 'weight'});
scatterPlot.zoom(140, 40, 200, 120);

let typeName = 'Bar chart';
let byName = DG.Viewer.fromType(typeName, t);

// fluent API
let barChart = t.plot.bar();

ui.divV([scatterPlot.root, byName.root, barChart.root]);
