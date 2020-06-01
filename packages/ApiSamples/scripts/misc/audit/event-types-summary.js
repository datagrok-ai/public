let view = grok.shell.newView('Usage');

grok.data.query('UsageAnalysis:EventsSummaryOnDate', {'date': 'today'})
  .then(t => view.append(DG.Viewer.barChart(t)));
