//api: grok.data.linkTables, DG.TableLink.detach
// Removing a link between tables (see also link-tables.js).
//
// "linkTables" returns the link it made. Here, selecting rows in the master table filters the
// detail table until "detach" removes the link; detaching also shows every detail row again.

const master = DG.DataFrame.fromColumns([DG.Column.fromStrings('key', ['a', 'b', 'c'])]);
const detail = DG.DataFrame.fromColumns([DG.Column.fromStrings('key', ['a', 'a', 'b', 'c', 'c', 'c'])]);

const link = grok.data.linkTables(master, detail, ['key'], ['key'], [DG.SYNC_TYPE.SELECTION_TO_FILTER]);
master.selection.set(0, true);
grok.shell.info(`Linked: ${detail.filter.trueCount} of ${detail.rowCount} detail rows shown`);

link.detach();
master.selection.set(2, true);
grok.shell.info(`Detached: ${detail.filter.trueCount} of ${detail.rowCount} detail rows shown`);
