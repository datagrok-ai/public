let view = gr.newView('Usage');

gr.query('ManualActivityByDate', {'date': 'today'})
  .then(t => view.append(Viewer.grid(t)));
