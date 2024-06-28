// Demonstrates markdown-like link in grid

const df = DG.DataFrame.fromColumns([
  DG.Column.fromList(DG.TYPE.STRING, 'markdown link', ['To learn more about [Datagrok](https://datagrok.ai/) open the [platform help](https://datagrok.ai/help/home)']),
]);
grok.shell.addTableView(df);