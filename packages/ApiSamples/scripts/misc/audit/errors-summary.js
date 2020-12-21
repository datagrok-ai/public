let view = grok.shell.newView('Usage');

grok.data.query('UsageAnalysis:ErrorsSummaryOnDate', {'date': 'this year'})
  .then(t => grok.shell.addTableView(t));
