let view = grok.newView('Usage');

grok.query('ManualActivityByDate', {'date': 'today'})
  .then(t => view.append(Viewer.grid(t)));
