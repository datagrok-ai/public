const df = await grok.data.query('UsageAnalysis:Tests');
grok.shell.addTableView(df);