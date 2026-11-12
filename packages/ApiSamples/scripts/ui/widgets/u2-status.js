//api: DG.U2.status, DG.U2.state, DG.U2.items, DG.U2.revealItem
// What a test reads off a u2 control: its state block, its status, its items, and an item found by label.
const list = new DG.U2.Control();
list.root.dataset.u2 = 'sample-list';
list.root.setAttribute('aria-label', 'Stages');
for (const stage of ['Discovery', 'Preclinical', 'Phase I']) {
  const row = ui.divText(stage);
  row.dataset.u2Item = stage.toLowerCase();
  list.root.append(row);
}
grok.shell.newView('u2 status', [list.root]);

const state = DG.U2.state(list.root);
grok.shell.info(`${state.kind} "${state.label}": enabled ${state.enabled}, visible ${state.visible}`);
grok.shell.info(`status parts: ${Object.keys(DG.U2.status(list.root).parts).length}, ` +
  `items: ${DG.U2.items(list.root)?.count ?? 'a plain control lists none'}`);
DG.U2.revealItem(list.root, 'phase i').then((row) => row.style.fontWeight = 'bold');
