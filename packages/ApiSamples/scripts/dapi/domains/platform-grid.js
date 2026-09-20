// An editable grid over a domain table, straight from the PLATFORM: any DG.Grid
// the app owns hosts a DG.DomainFrameEditor through Grid.attachEditor — no
// package installed, no library imported, available in every session.
//
// Free: the table's registered handler decorates the grid
// (DomainObjectHandler.decorateGrid), in-grid edits are tracked / validated /
// highlighted (amber, red, conflict orange), the whole batch saves as ONE
// transaction with the 409 flows, and everything is permission-gated down to
// per-column writability (a table the caller cannot edit degrades to read-only).

const issues = grok.dapi.domains.table('grit.issue');
const df = await issues.queryDf({limit: 50});
const editor = await DG.DomainFrameEditor.attach(df, issues, {query: {limit: 50}});

const grid = DG.Grid.create(df);
DG.DomainObjectHandler.decorateGrid(grid, 'grit.issue', df);
grid.attachEditor(editor);
grok.shell.newView('Issues', [grid.root]);

// An edit of a cell the caller may not write is refused ONCE per attempt, in the
// table's own words ('Title is read-only', 'Issue #4 is read-only for you'): the
// grid balloons it, and onRefused hands it to a host with a status line of its own.
editor.onRefused.subscribe((r) => console.log(`refused at row ${r.row}: ${r.message}`));

// The editing state lives in the FRAME, in three service columns whose single
// writer is DG.DomainFrameEditor — and which never reach an export, an upload
// or a saved project. The grid hides every `~` column and locks it.
console.log(DG.DomainFrameEditor.SERVICE_COLUMNS.join(', '));

if (df.rowCount > 0) {
  editor.setValue(0, 'title', 'Renamed from a sample');
  console.log(`row 0 is now "${editor.stateOf(0)}", ${editor.changeCount} unsaved change(s)`);
  editor.discard();                           // await editor.save() would write it
}

// A row that does not exist yet carries a draft id ('~new:…') another row may
// reference; the save (one /transaction) resolves it and reports the server id.
const row = editor.addRow({title: 'Draft from a sample'});
console.log(`draft id: ${df.get('id', row)}`);
editor.discard();

// Detaching the grid drops the editor wiring with the viewer's own subscriptions.
grid.detachEditor();
