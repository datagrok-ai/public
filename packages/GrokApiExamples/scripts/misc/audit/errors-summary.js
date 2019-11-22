let view = gr.newView('Usage');

gr.query('ErrorsSummaryOnDate', {'date': 'this year'})
    .then(t => gr.addTableView(t));