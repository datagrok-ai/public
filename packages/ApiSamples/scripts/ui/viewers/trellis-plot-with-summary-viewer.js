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
    invertColorScheme: true,
    colorSchemes: [DG.Color.categoricalPalette, [DG.Color.darkGray, DG.Color.blue, DG.Color.green, DG.Color.darkGreen,
      DG.Color.yellow, DG.Color.olive, DG.Color.orange]],
  }
});