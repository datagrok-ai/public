grok.data.query('UsageAnalysis:TopPackageErrors', {date: 'this week', users: ['all']})
  .then((t) => grok.shell.addTableView(t));