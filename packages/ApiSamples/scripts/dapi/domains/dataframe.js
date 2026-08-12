// Query a domain table straight into a typed DataFrame (d42 wire format, 10M row cap
// vs 10k for query()). Columns carry db property tags, choices, and semantic types,
// so grids and property panels work out of the box.

const items = grok.dapi.domains.table('apitests.item');
const stamp = `DF-${Date.now()}`;
await items.batch([
  {sku: `${stamp}-1`, name: 'Bolt', quantity: 10},
  {sku: `${stamp}-2`, name: 'Nut', quantity: 7}]);

const df = await items.queryDf({filter: `sku starts "${stamp}"`, sort: '!created_on', limit: 100});
grok.shell.addTableView(df);

// Enrichment-style bulk fetch: rows for a set of ids as a typed DataFrame
// ('id' plus the requested fields; missing/invisible ids are simply absent rows)
const ids = [...Array(df.rowCount).keys()].map((i) => df.col('id').get(i));
const fields = await items.fetchFields(ids, ['sku', 'name']);
grok.shell.info(`fetchFields: ${fields.rowCount} rows, ${fields.columns.names().join(', ')}`);

// Cleanup: one bounded filtered bulk delete
for (let guard = 0; guard < 100; guard++)
  if (!(await items.deleteWhere(`sku starts "${stamp}"`)).hasMore)
    break;
