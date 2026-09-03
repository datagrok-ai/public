//name: SetupPltsSchema
//description: One-time setup of the plts domain schema through the server Domain API: grants View/Edit/Delete on the securing tables to All users, then seeds the dictionaries (plate types, default well-level properties) via business-key inserts. Idempotent — grants re-apply cleanly and every seed merges by business key. Run by an admin after the first deploy; replaces the legacy Plates:Plts connection grant and the seed rows of the retired 0000_init.sql.
//language: javascript

const table = (name) => grok.dapi.domains.table('plts.' + name);

// 1. Grants. Only the four securing tables need them: the junction, value,
// run, and result tables are master-mode and delegate security to their
// template or plate.
const allUsers = await grok.dapi.groups.find(DG.Group.defaultGroupsIds['All users']);
for (const name of ['plate_types', 'properties', 'templates', 'plates'])
  for (const permission of ['View', 'Edit', 'Delete'])
    await table(name).grant(allUsers.id, permission);
grok.shell.info('plts: granted View/Edit/Delete to All users');

// 2. Dictionaries. insert() on a table with a business key reports
// {status: 'duplicate'} instead of failing, so this doubles as a get-or-create.
const seed = async (name, rows) => {
  const reports = await table(name).insert(rows);
  const created = reports.filter((r) => r.created).length;
  grok.shell.info(`${name}: ${created} created, ${rows.length - created} already there`);
};

await seed('plate_types', [
  {name: 'Generic 96 wells', rows: 8, cols: 12},
  {name: 'Generic 384 wells', rows: 16, cols: 24},
  {name: 'Generic 1536 wells', rows: 32, cols: 48},
]);

await seed('properties', [
  {name: 'Volume', type: 'double', scope: 'well', units: 'uL'},
  {name: 'Concentration', type: 'double', scope: 'well', units: 'uM'},
  {name: 'Sample', type: 'string', scope: 'well'},
  {name: 'Well Role', type: 'string', scope: 'well',
    choices: JSON.stringify(['Empty', 'Sample', 'DMSO', 'Low Control', 'High Control'])},
]);
