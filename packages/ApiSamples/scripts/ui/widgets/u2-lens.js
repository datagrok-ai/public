//api: DG.U2.lens, DG.U2.afterFlush
// lens: a writable signal over a part of another signal. afterFlush: run something once the batch that changed
// signals has ended and the effects it set off have run.
const person = DG.U2.signal({name: 'Ann', age: 41});
const name = DG.U2.lens(person, (p) => p.name, (p, value) => ({...p, name: value}));
name.value = 'Bo';
grok.shell.info(`person is ${JSON.stringify(person.value)}`);

const rows = DG.U2.signal(['first']);
const owner = DG.U2.Scope.root('lens sample');
const list = ui.divV([]);
owner.effect(() => list.replaceChildren(...rows.value.map((r) => ui.divText(r))));
DG.U2.batch(() => {
  rows.value = [...rows.value, 'second'];
  DG.U2.afterFlush(() => grok.shell.info(`after the batch: ${list.children.length} rows`));
  grok.shell.info(`inside the batch: ${list.children.length} row`);
});
owner.dispose();
