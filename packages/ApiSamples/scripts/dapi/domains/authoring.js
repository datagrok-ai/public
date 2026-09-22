//api: DG.DomainsDataSource.draft, DG.DomainsDataSource.createSchema, DG.DomainSchemaClient.validate, DG.DomainSchemaClient.delete
// Authoring an EXTERNAL binding from a script: draft a manifest over a database you may
// query, dry-run it, create the schema, share it, delete it. Needs the CreateDomainSchema
// privilege and GetSchema + Query on the connection (here the Northwind demo container,
// core/docs/features/ems/external-bindings/fixtures/northwind). Nothing below writes the
// warehouse: the schema is registry-only and its delete leaves the connection untouched.

if (!(await grok.dapi.domains.schemas.list()).some((s) => s.name === 'northwind'))
  return grok.shell.info('Publish the Northwind binding fixture first (it ships the PostgresNorthwind connection)');

const connection = 'NorthwindBinding:PostgresNorthwind';
const name = `nw_${Math.random().toString(36).slice(2, 8)}`;

// 1. Draft: remote names kept, the primary key as the business key, one-column foreign keys
//    onto drafted tables as refs; `inventory` says what stayed out and why.
const draft = await grok.dapi.domains.draft({connection, schema: 'public', tables: ['orders', 'order_details']});
grok.shell.info(`drafted ${Object.keys(draft.manifest.tables).join(', ')}; ` +
  `${draft.inventory.relations.filter((r) => r.status === 'ref').length} ref(s), ` +
  `${draft.inventory.tables.filter((t) => !t.bindable).length} table(s) left out`);

// 2. Edit the draft as any manifest, then dry-run: every check, nothing registered.
draft.manifest.tables.orders.friendlyName = 'Sales orders';
const plan = await grok.dapi.domains.createSchema(name, {manifest: draft.manifest, dryRun: true});
grok.shell.info(`dry run: ${plan.status}`);

// 3. Create: validated live against the warehouse, then registered as ext_<name>.
const created = await grok.dapi.domains.createSchema(name, {friendlyName: 'Northwind sales', manifest: draft.manifest});
grok.shell.info(`created ${created.pgSchema}, binding ${created.binding.status}`);
const schema = grok.dapi.domains.schema(name);
try {
  // 4. Share: row access is granted per table (schema grants gate schema operations only).
  const allUsers = await grok.dapi.groups.filter('friendlyName = "All users"').first();
  await grok.dapi.domains.table(`${name}.orders`).grant(allUsers.id, 'View');
  grok.shell.info(`${await grok.dapi.domains.table(`${name}.orders`).count()} orders through the binding`);

  // 5. Re-validate later (a warehouse can drift): the recorded verdict, issues by manifest path.
  const validation = await schema.validate();
  grok.shell.info(`validated ${validation.validatedOn}: ${validation.status}`);
} finally {
  // 6. Delete: registry only — the connection and the warehouse stay.
  await schema.delete();
}
