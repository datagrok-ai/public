let view = grok.newView('Usage');

grok.query('ErrorsSummaryOnDate', {'date': 'this year'})
    .then(t => grok.addTableView(t));