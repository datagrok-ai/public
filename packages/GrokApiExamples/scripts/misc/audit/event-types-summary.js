let view = grok.newView('Usage');

grok.query('UsageAnalysis:EventsSummaryOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.barChart(t)));
