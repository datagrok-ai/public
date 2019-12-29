let view = grok.newView('Usage');

grok.query('EventsSummaryOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.barChart(t)));
