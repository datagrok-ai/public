//api: DG.U2.Scope.root, DG.U2.Scope.build, DG.U2.Scope.task, DG.U2.Scope.subscribe
// The owner verbs. Scope.root starts a tree on purpose; build runs code with the owner ambient; task does async
// work and gets back into the owner after an await (t.build skips once the owner is gone); subscribe owns a stream.
const session = DG.U2.Scope.root('demog loader');
const host = ui.divV([]);
const view = grok.shell.newView('u2 owner verbs', [host]);
host.append(session.build(() => new DG.U2.Control(ui.h2('Demog'))).root);

const load = () => session.task(async (t) => {
  const df = await grok.data.getDemoTable('demog.csv');
  t.build(() => host.append(ui.divText(`${df.rowCount} rows`)));
}, {serial: 'load'});
load();
load();
load();   // one load in flight: the second waits, the third replaces it — two loads in all

session.subscribe(grok.events.onViewRemoved, (v) => {
  if (v.id === view.id)
    session.dispose();   // the loads still running skip their t.build
});
