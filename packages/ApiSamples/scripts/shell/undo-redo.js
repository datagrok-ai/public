// Making your own operations undoable (Ctrl+Z / Ctrl+Shift+Z).
//
// Records live in memory only, are capped by count and by size, and are dropped when the
// table or view they belong to is closed. Capture values, not row indexes of live objects.

let demog = grok.data.demo.demog(100);
let view = grok.shell.addTableView(demog);
let column = demog.col('age');

let before = column.get(0);
column.set(0, 100);

DG.UndoService.push('Set age', () => column.set(0, before), {
  redo: () => column.set(0, 100),
  context: demog
});

grok.shell.info(`Ctrl+Z to undo "${DG.UndoService.undoName}"`);

// the same thing, with the operation performed by the service
DG.UndoService.run('Rename table', () => demog.name = 'patients', () => demog.name = 'demog',
  { context: demog });

grok.shell.onUndo.subscribe((args) => grok.shell.info(`Undone: ${args.name}`));

// grok.shell.undo();  // same as pressing Ctrl+Z
// grok.shell.canUndo, grok.shell.canRedo — enablement
