/**
 * Create viewers and use them anywhere via static methods of the DG.Viewer class
 * See also: adding viewers to the TableView
 */

let t = grok.data.demo.demog();
grok.shell.newView('foo').append(ui.divH([
  DG.Viewer.scatterPlot(t, {x: 'height', y: 'weight'}),
  DG.Viewer.histogram(t, {value: 'height'}),
  DG.Viewer.barChart(t),
]));