const df = await grok.data.query('UsageAnalysis:TestsDashboard');
grok.shell.addTableView(df);