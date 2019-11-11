let view = gr.newView('Usage');

gr.query('EventsSummaryOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.barChart(t)));
