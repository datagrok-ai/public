let errors = await grok.data.query('UsageAnalysis:TopErrors', {date: 'this week'})
grok.shell.addTableView(errors)