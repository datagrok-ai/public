let view = gr.newView('Usage');

gr.query('EventsOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.scatterPlot(t)));
