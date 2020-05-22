let view = grok.shell.newView('Usage');

grok.data.query('UsageAnalysis:ErrorsOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.grid(t)));
