// The Space row of the Save project dialog in a dialog of your own: where a project goes (a space, or your own
// projects when empty), and Move or Link when a saved project changes space. Null when no package provides it.
// The last argument words the empty box and its tooltip for what the dialog saves.

const row = await ui.spaceSaveRow(null, 'Kinase panel', null,
  {placeholder: 'Your projects', tooltip: 'Empty: kept in your own projects. Pick a space to share it there.'});
if (row === null)
  grok.shell.warning('No package provides the Space row');
else {
  ui.dialog('Save to space')
    .add(row)
    .onOK(async () => grok.shell.info(JSON.stringify(await row.resolveTarget())))
    .show();
}
