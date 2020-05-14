let view = grok.newView('Usage');

grok.query('UsageAnalysis:EventsOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.scatterPlot(t)));
