const df = grok.data.demo.demog();
const tv = grok.shell.addTableView(df);

tv.addViewer(DG.VIEWER.TRELLIS_PLOT, {
  viewerType: 'Summary',
  innerViewerLook: {
    columnNames: ['weight', 'height'],
    aggregations: [DG.STATS.MIN, DG.STATS.MED],
    visualizationType: 'bars',
    colorColumnName: 'age',
    colorAggrType: 'stdev',
    invertColorScheme: true
  }
});