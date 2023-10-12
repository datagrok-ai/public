const df = await grok.data.query('UsageAnalysis:TestsToday');
grok.shell.addTableView(df);