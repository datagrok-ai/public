const date = new Date().toLocaleDateString('en-US', {month: '2-digit', day: '2-digit', year: 'numeric'});
const df = await grok.data.query('UsageAnalysis:TestsToday', {date: date});
grok.shell.addTableView(df);