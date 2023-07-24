let showStats = true;  // external storage in this case
let showStatsProperty = DG.Property.bool(
  'showStatistics',
  (plot) => showStats,
  (plot, value) => showStats = value)
showStatsProperty.category = 'XXX';
DG.Property.registerAttachedProperty('ScatterPlotLook', showStatsProperty);

grok.events.onViewerAdded.subscribe((data) => {
  console.log(data.args.viewer);
});

// this property is now editable
let view = grok.shell.addTableView(grok.data.demo.demog());
let plot = view.scatterPlot();
// plot.onAfterDrawScene.subscribe((_) => {
//   if (showStats) {
//     plot.
//   }
// });
