// A line chart with several variables plotted on one chart

let df = DG.DataFrame.fromColumns([
  DG.Column.fromType(DG.TYPE.INT, 'T', 0),
  DG.Column.fromType(DG.TYPE.FLOAT, 'Speed', 0),
  DG.Column.fromType(DG.TYPE.FLOAT, 'Velocity', 0),
  DG.Column.fromType(DG.TYPE.STRING, 'Object', 0)
]);

df.rows.addNew([1, 0.1, 0.04, 'Object1']); df.rows.addNew([1, 0.3, 0.05, 'Object2']); df.rows.addNew([1, 0.4, 0.02, 'Object3']);
df.rows.addNew([2, 0.2, 0.10, 'Object1']); df.rows.addNew([2, 0.1, 0.08, 'Object2']); df.rows.addNew([2, 0.3, 0.06, 'Object3']);
df.rows.addNew([3, 0.3, 0.05, 'Object1']); df.rows.addNew([3, 0.5, 0.03, 'Object2']); df.rows.addNew([3, 0.7, 0.09, 'Object3']);

let lineChart = DG.Viewer.fromType(DG.VIEWER.LINE_CHART, df);
lineChart.setOptions({
  xColumnName: 'T',
  yColumnNames: ['Speed'],
  // For two combined charts:
  // yColumnNames: ['Speed', 'Velocity'],
  split: 'Object'
});

grok.shell.newView('Lines').append(lineChart.root);
