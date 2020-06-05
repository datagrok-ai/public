let view = grok.shell.newView('Usage');

grok.data.query('TestTrack:ManualActivityByDate', {'date': 'today'})
  .then(t => view.append(DG.Viewer.grid(t).root));
