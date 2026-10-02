// Table-level access of the CURRENT user (DomainTableClient.access): one
// server-composed answer — {can: {view, insert, edit, delete, share}, fields:
// {<column>: 'editable' | 'readonly'}} — from the predicates the writes apply
// (grants on the securing entity, column security, relation travel), cached per
// table until grok.dapi.domains.invalidateUiCaches() after a grant change. Gate
// UI on it to spare users the common 403s; the server still enforces every write.
// A column the caller may not SEE is absent from `fields` altogether.

if (!(await grok.dapi.domains.schemas.list()).some((s) => s.name === 'apitests'))
  return grok.shell.info('Deploy the ApiTests package first (it declares the apitests domain schema)');

const items = grok.dapi.domains.table('apitests.item');
const access = await items.access();
// {can, fields, travelableRelations, securingTable, securityMode, audit, hasBusinessKey, support}
grok.shell.info(`${access.securingTable} governs access (${access.securityMode} mode); ` +
  `expandable relations: ${access.travelableRelations.join(', ') || 'none'}`);

// `support` is the other half, and it is NOT about this caller: what the TABLE can do at all,
// computed once on the server from the registration. `can` is permission (a grant flips it),
// `support` is storage (no grant flips it) — so a client installs an optional affordance only
// where support declares it, and never probes the table's shape to guess:
//   writes      insert / update / batch / updateWhere (false on a read-only registration)
//   deleted     `deleted: 'include' | 'only'` reads and the `~is_deleted` column
//   restore     restore(id) — writes AND deleted
//   audit       audit(id) / auditLog(), and row watch, which needs the trail
//   ancestors   pathTo(id) and the `under` filter term (a hierarchy table only)
//   probe       updated_on is there, so a live list can poll for the newest write
//   watch       watch(id) subscribes to change notifications
//   systemColumns  the system columns this table physically has (a Core registration lists fewer)
const s = access.support;
grok.shell.info(`This table ${s.writes ? 'accepts' : 'refuses'} writes, ` +
  `${s.deleted ? 'has' : 'has no'} trash, ${s.ancestors ? 'is' : 'is not'} a tree; ` +
  `system columns: ${s.systemColumns.join(', ')}`);

// Gate the affordance, do not try and catch: an op the table cannot do fails with a
// DomainUnsupportedError (422) that no permission would have fixed.
const history = s.audit ? await items.auditLog({limit: 5}) : [];
grok.shell.info(s.audit ? `Last ${history.length} table-wide audit entries` : 'This table keeps no history');

// A form over the columns this caller may WRITE — a value in any other column
// fails the whole insert ('forbidden-column').
const inputs = {
  sku: ui.input.string('Sku', {value: `ACC-${Date.now()}`}),
  name: ui.input.string('Name', {value: 'Access sample'}),
  quantity: ui.input.int('Quantity', {value: 1}),
};
const editable = Object.keys(inputs).filter((c) => access.fields[c] === 'editable');

// The create button is offered only to callers the server would let insert.
const create = ui.bigButton('CREATE', async () => {
  const values = {};
  for (const c of editable)
    values[c] = inputs[c].value;
  const [row] = await items.insert(values);
  // Per-row truth rides the rows: `withAccess` adds ~can_edit / ~can_delete /
  // ~can_share (DG.DOMAIN_ACCESS_COLUMNS) where the table-level flags are false
  // negatives (a promoted row of a row-mode table, a master row's children).
  const fresh = await items.get(row.id, {withAccess: true});
  grok.shell.info(`Inserted ${row.id} (can edit: ${fresh['~can_edit']}, can delete: ${fresh['~can_delete']}) — deleting it again`);
  await items.delete(row.id);
});

ui.dialog('apitests.item')
  .add(ui.inputs(editable.map((c) => inputs[c])))
  .add(access.can.insert ? create : ui.divText('You cannot add rows to this table.'))
  .show();
