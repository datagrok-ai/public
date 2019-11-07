let view = gr.newView('Usage');

gr.query('ErrorsOnDate', {'date': 'today'})
  .then(t => view.append(Viewer.grid(t)));
