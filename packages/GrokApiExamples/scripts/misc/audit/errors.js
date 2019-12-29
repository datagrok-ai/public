let view = grok.newView('Usage');

grok.query('ErrorsOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.grid(t)));
