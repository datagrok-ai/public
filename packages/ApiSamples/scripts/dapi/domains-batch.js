// Bulk upload into a domain table: a CSV string, a DG.DataFrame, row objects, or
// Uint8Array bytes ({format: 'd42' | 'parquet'}); upsert merges by the table's business key.

const items = grok.dapi.domains.table('apitests.item');
const stamp = `S-${Date.now()}`;

// Insert two rows from CSV
let report = await items.batch(`sku,name,quantity\n${stamp}-1,Bolt,10\n${stamp}-2,Nut,20`);
grok.shell.info(`inserted: ${report.inserted}`);

// Upsert: matched business keys are updated, new ones inserted
report = await items.batch([
  {sku: `${stamp}-1`, quantity: 11},
  {sku: `${stamp}-3`, name: 'Washer', quantity: 30},
], {mode: 'upsert'});
grok.shell.info(`inserted: ${report.inserted}, updated: ${report.updated}`);

// Cleanup: one filtered bulk delete (bounded; loop while hasMore for larger sets)
for (let guard = 0; guard < 100; guard++)
  if (!(await items.deleteWhere(`sku starts "${stamp}"`)).hasMore)
    break;
