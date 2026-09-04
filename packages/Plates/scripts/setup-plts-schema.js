//name: SetupPltsSchema
//description: One-time plts schema setup: grants on the securing tables to All users, then seeds the dictionaries. Idempotent; run by an admin after the first deploy.
//language: javascript

const table = (name) => grok.dapi.domains.table('plts.' + name);

// Only the securing tables need grants; master-mode tables delegate security to their template or plate.
const allUsers = await grok.dapi.groups.find(DG.Group.defaultGroupsIds['All users']);
for (const name of ['plate_type', 'property', 'template', 'plate'])
  for (const permission of ['View', 'Edit', 'Delete'])
    await table(name).grant(allUsers.id, permission);
grok.shell.info('plts: granted View/Edit/Delete to All users');

// insert() with a business key reports duplicates instead of failing, so seeding doubles as get-or-create.
const seed = async (name, rows) => {
  const reports = await table(name).insert(rows);
  const created = reports.filter((r) => r.created).length;
  grok.shell.info(`${name}: ${created} created, ${rows.length - created} already there`);
};

await seed('plate_type', [
  {name: 'Generic 96 wells', rows: 8, cols: 12},
  {name: 'Generic 384 wells', rows: 16, cols: 24},
  {name: 'Generic 1536 wells', rows: 32, cols: 48},
]);

await seed('property', [
  {name: 'Volume', type: 'double', scope: 'well', units: 'uL'},
  {name: 'Concentration', type: 'double', scope: 'well', units: 'uM'},
  {name: 'Sample', type: 'string', scope: 'well'},
  {name: 'Well Role', type: 'string', scope: 'well',
    choices: JSON.stringify(['Empty', 'Sample', 'DMSO', 'Low Control', 'High Control'])},
]);
