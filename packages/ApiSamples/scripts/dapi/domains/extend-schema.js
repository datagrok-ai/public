// Extending a plugin's domain schema with your own tables and columns.
// A plugin opts in from its databases/<schema>/schema.json — "extensible": {"tables": true}
// at the root lets users add tables, "extensible": true on a table lets them add columns.
// Users then need the 'Extend' permission on the SCHEMA entity.
// Extensions live only in the server registry, so `grok api` codegen never sees them —
// address them through the generic client, by their LOGICAL name.

const schema = grok.dapi.domains.schema('apitests');

// The registry entity carries the opt-in flag and the extension counter;
// schema.manifest() shows the same flags in manifest vocabulary.
const info = (await grok.dapi.domains.schemas.list()).find((s) => s.name === 'apitests');

if (!info?.extensibleTables) {
  grok.shell.warning("The 'apitests' schema does not accept user extensions " +
    '(add "extensible": {"tables": true} to its schema.json)');
} else {
  // 1. Delegate extending. 'Extend' is grantable on schema entities only —
  //    on a table or column schema it is a 400. Requires Share on the schema.
  const analysts = (await grok.dapi.groups.filter('Analysts').first());
  if (analysts != null)
    await schema.grant(analysts.id, 'Extend');

  // 2. Add a table of your own, and a column of your own to a plugin table, in one apply.
  //    'extend' is full state for YOUR columns of that table: omitting one you added
  //    proposes its drop. Plugin objects are immutable here. ifVersion tracks ext_version
  //    (NOT the plugin's version) — pass it to lose cleanly against a parallel editor.
  const plan = await schema.apply({
    tables: {
      customers: {columns: {name: {type: 'string', isName: true}}},
    },
    extend: {
      item: {columns: {customer_id: {type: 'ref', ref: 'customers', onDelete: 'setnull'}}},
    },
    ifVersion: `${info.extVersion ?? 0}`,
  }, {dryRun: true});
  grok.shell.info(`plan: creates ${plan.creates.tables.join(', ') || '—'}, ` +
    `columns ${JSON.stringify(plan.creates.columns)}, extVersion ${plan.extVersion}`);

  // Columns you add to a PLUGIN table stay nullable, non-unique, never the display-name
  // column, and never cascade deletes — restrict/setnull only. Your own TABLES keep the
  // full manifest vocabulary.
  await schema.apply({
    tables: {customers: {columns: {name: {type: 'string', isName: true}}}},
    extend: {item: {columns: {customer_id: {type: 'ref', ref: 'customers', onDelete: 'setnull'}}}},
  });

  // 3. Use the new objects like any other. The extension column is physically prefixed
  //    server-side, but every API surface — insert, query, filter, patch, batch, d42 —
  //    keeps using the logical name.
  const customers = grok.dapi.domains.table('apitests.customers');
  const items = grok.dapi.domains.table('apitests.item');

  const [acme] = await customers.insert({name: 'Acme'});
  const [row] = await items.insert({sku: `SKU-${Date.now()}`, name: 'Widget', customer_id: acme.id});
  const mine = await items.query().where('customer_id', '=', acme.id).top(10);
  grok.shell.info(`${mine.length} item(s) for Acme; inserted ${row.id}`);

  // Soft-deleting the customer honors the onDelete you declared: setnull clears the
  // reference instead of blocking the delete.
  await customers.delete(acme.id);
  grok.shell.info(`after delete: ${(await items.query().where('customer_id', '=', acme.id).top(1)).length} referrers`);
}
