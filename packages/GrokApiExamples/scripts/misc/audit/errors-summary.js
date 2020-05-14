let view = grok.newView('Usage');

grok.query('UsageAnalysis:ErrorsSummaryOnDate', {'date': 'this year'})
    .then(t => grok.addTableView(t));