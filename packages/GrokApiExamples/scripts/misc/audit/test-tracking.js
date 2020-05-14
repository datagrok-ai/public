let view = grok.newView('Usage');

grok.query('TestTrack:ManualActivityByDate', {'date': 'today'})
  .then(t => view.append(Viewer.grid(t)));
