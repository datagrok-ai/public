// DG.DomainGrid — an editable grid over a domain table, straight from the
// PLATFORM: no package installed, no library imported, available in every
// session (the `@datagrok-libraries/domain-ui` classes are these same ones).
//
// Free: the table's registered handler decorates the grid, in-grid edits are
// tracked / validated / highlighted (amber, red, conflict orange), the save bar
// writes the whole batch as ONE transaction with the 409 flows, and everything
// is permission-gated down to per-column writability.

const grid = await DG.DomainGrid.create(grok.dapi.domains.table('grit.issue'),
  {query: {limit: 50}});
grok.shell.newView('Issues', [grid.root]);

// The editing state lives in the FRAME, in three service columns whose single
// writer is DG.DomainFrameEditor — and which never reach an export, an upload
// or a saved project.
console.log(DG.DomainFrameEditor.SERVICE_COLUMNS.join(', '));

if (grid.editor.dataFrame.rowCount > 0) {
  grid.editor.setValue(0, 'title', 'Renamed from a sample');
  console.log(`row 0 is now "${grid.editor.stateOf(0)}", ` +
    `${grid.editor.changeCount} unsaved change(s)`);
  grid.editor.discard();                      // await grid.editor.save() would write it
}

// A grid's machine surface: its status, its dataFrame, and REAL platform Funcs
// from the shared widget-action vocabulary (DomainUi:Save, Discard, AddRow, ...).
grok.shell.info(`${grid.getWidgetStatus().description}\n` +
  `Actions: ${grid.getFunctions().map((f) => f.name).join(', ')}`);
