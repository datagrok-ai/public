// Table-level capabilities of the CURRENT user (DomainTableClient.capabilities):
// one server-composed answer — the predicates the writes apply (View/Edit/Delete/
// Share on the securing entity, column security, relation travel) — cached per
// table until grok.dapi.domains.invalidateUiCaches() after a grant change. Gate
// UI on it to spare users the common 403s; the server still enforces every write.

if (!(await grok.dapi.domains.schemas.list()).some((s) => s.name === 'apitests'))
  return grok.shell.info('Deploy the ApiTests package first (it declares the apitests domain schema)');

const items = grok.dapi.domains.table('apitests.item');
const caps = await items.capabilities();
// {canView, canInsert, canEdit, canDelete, canShareTable, writableColumns,
//  travelableRelations, securingTable, securityMode, audit, hasBusinessKey}
grok.shell.info(`${caps.securingTable} governs access (${caps.securityMode} mode); ` +
  `expandable relations: ${caps.travelableRelations.join(', ') || 'none'}`);

// A form over the columns this caller may WRITE — a value in any other column
// fails the whole insert ('forbidden-column').
const inputs = {
  sku: ui.input.string('Sku', {value: `CAP-${Date.now()}`}),
  name: ui.input.string('Name', {value: 'Capabilities sample'}),
  quantity: ui.input.int('Quantity', {value: 1}),
};
const writable = Object.keys(inputs).filter((c) => caps.writableColumns.includes(c));

// The create button is offered only to callers the server would let insert.
const create = ui.bigButton('CREATE', async () => {
  const values = {};
  for (const c of writable)
    values[c] = inputs[c].value;
  const [row] = await items.insert(values);
  grok.shell.info(`Inserted ${row.id} — deleting it again`);
  await items.delete(row.id);
});

ui.dialog('apitests.item')
  .add(ui.inputs(writable.map((c) => inputs[c])))
  .add(caps.canInsert ? create : ui.divText('You cannot add rows to this table.'))
  .show();
