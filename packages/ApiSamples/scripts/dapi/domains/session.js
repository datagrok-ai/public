// DG.DomainSession — the unit of work: several DomainFrameEditors save as ONE
// /transaction, so a child row can point at a parent row that does not exist yet.
// The parent's draft id ('~new:<uuid>', stamped by addRow) travels as the insert's
// `ref`; the server resolves it and the session hands every editor the real id back
// through DomainSaveResult.assigned. Any failure rolls the whole batch back.

const items = grok.dapi.domains.table('apitests.item');
const events = grok.dapi.domains.table('apitests.item_event');
const sku = `SKU-${Date.now()}`;

// Two editors over two tables. Each owns its frame and its pending batch.
const parent = await DG.DomainFrameEditor.create(items, {query: {filter: `sku = "${sku}"`}});
const child = await DG.DomainFrameEditor.create(events, {query: {filter: `kind = "${sku}"`}});

// The parent does not exist yet: its id is a draft the child refers to as a plain value.
const row = parent.addRow({sku: sku, name: 'Widget', quantity: 1});
const draft = parent.dataFrame.get('id', row);
child.addRow({item_id: draft, kind: sku, amount: 1});

const session = new DG.DomainSession([parent, child]);
grok.shell.info(`${session.changeCount} pending changes, dirty: ${session.isDirty}`);

session.onSaved.subscribe((r) => grok.shell.info(
  `inserted ${r.inserted}, updated ${r.updated}, deleted ${r.deleted}; ${draft} → ${r.assigned[draft]}`));

// ONE transaction: both inserts share a tx_id, and the reference is resolved server-side.
if (await session.save()) {
  const itemId = parent.dataFrame.get('id', 0);
  const event = await events.get(child.dataFrame.get('id', 0));
  grok.shell.info(`item ${itemId}; the event points at ${event.item_id}`);
}
session.dispose();

await items.deleteWhere({property: 'sku', operator: '=', value: sku});
parent.detach();
child.detach();
