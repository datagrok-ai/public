// Soft delete is reversible: `deleted` picks the scope of a read ('exclude' — the
// default, 'include', 'only' — the trash), and restore is the Delete right spent
// backwards. Off 'exclude' the rows carry ~is_deleted (a service column: hidden by
// Grid.attachEditor, never exported), and the restore is audited as 'undelete'.

const items = grok.dapi.domains.table('apitests.item');
const [row] = await items.insert({sku: `TRASH-${Date.now()}`, name: 'Restorable'});
await items.delete(row.id);

// The trash of this table: only the deleted rows, each saying so. count/exists take the
// same scope, so a trash list's total agrees with its rows.
const [deleted] = await items.query({filter: `id = "${row.id}"`, deleted: 'only'});
grok.shell.info(`in the trash: ${deleted.name} (~is_deleted: ${deleted['~is_deleted']}), ` +
  `${await items.count(null, {deleted: 'only'})} deleted rows in all`);

// One trashed row is addressable too — `get` takes the same scope, so a trash entity
// page can open a deleted row (read-only until it is restored).
grok.shell.info(`get() in the default scope: ${await items.get(row.id)}; ` +
  `in the trash: ${(await items.get(row.id, {deleted: 'only'})).name}`);

// A frame over both scopes — a service column (DG.DOMAIN_SERVICE_COLUMNS) is one the
// grid hides and no export carries.
const df = await items.queryDf({filter: `id = "${row.id}"`, deleted: 'include'});
const service = DG.DOMAIN_SERVICE_COLUMNS.filter((c) => df.col(c) !== null);
grok.shell.info(`csv without ${service.join(', ')}: ${df.toCsv()}`);

// Restore: the row comes back with the next version, audited as 'undelete'.
// A row whose reference points at a still-deleted parent is refused (DomainRestrictError
// naming that column) — restore the parent first.
const restored = await items.restore(row.id);
grok.shell.info(`restored ${restored.id} at version ${restored.version}; ` +
  `ops: ${(await items.audit(row.id)).map((a) => a.op).join(', ')}`);

// The same restore as one op of a transaction — the very code path the route runs,
// so a parent and its child come back together (the server orders the parent first).
await items.delete(row.id);
const [viaTransaction] = await grok.dapi.domains.transaction('apitests',
  [{op: 'restore', table: 'item', id: row.id}]);
grok.shell.info(`restored in a transaction at version ${viaTransaction.version}`);

// …and in an editor a restore is a PENDING change like any other: markRestored stages
// it, save() lands it together with everything else the batch holds, discard() takes
// it back. Only a row the SERVER holds as deleted (DG.DOMAIN_DELETED_COLUMN) can be.
await items.delete(row.id);
const editor = await DG.DomainFrameEditor.create(items,
  {query: {filter: `id = "${row.id}"`, deleted: 'only'}});
editor.markRestored(0);
grok.shell.info(`staged: state "${editor.stateOf(0)}", ${editor.changeCount} pending change`);
await editor.save();
grok.shell.info(`saved: ${DG.DOMAIN_DELETED_COLUMN} is now ` +
  `${editor.dataFrame.get(DG.DOMAIN_DELETED_COLUMN, 0)}`);
editor.detach();

await items.delete(row.id);
