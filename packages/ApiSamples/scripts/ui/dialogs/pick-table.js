// Platform table-picker dialogs: 'Select a file' over the file shares,
// 'Select a database query' over the database connections.
// Each promise resolves with the picked DataFrame, or null when the dialog is closed without a pick.

const df = await ui.pickTableFromFiles(); // or: await ui.pickTableFromQuery()
if (df !== null)
  grok.shell.addTableView(df);
else
  grok.shell.info('No table picked');
