let view = grok.shell.newView('Usage');

grok.data.query('UsageAnalysis:EventsOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.scatterPlot(t)));
