// Multiple time series plotted on one line chart

let df = grok.data.demo.biosensor(100);
let unpivoted = df.unpivot(['time'], ['x', 'y', 'z']);
let lineChart = DG.Viewer.lineChart(unpivoted, {
  xColumnName: 'time',
  yColumnNames: ['Value'],
  split: 'Category'
});
grok.shell.newView('Lines').append(lineChart.root);