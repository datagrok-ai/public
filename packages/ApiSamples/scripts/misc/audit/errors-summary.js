grok.data.query('UsageAnalysis:TopErrors', {date: 'this week'})
  .then((t) => grok.shell.addTableView(t));