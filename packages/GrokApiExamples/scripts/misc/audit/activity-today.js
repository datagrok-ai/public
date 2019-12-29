let view = grok.newView('Usage');

grok.query('EventsOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.scatterPlot(t)));
