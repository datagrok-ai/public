let view = grok.newView('Usage');

grok.query('UsageAnalysis:ErrorsOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.grid(t)));
