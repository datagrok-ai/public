// Server-side grouped aggregation over the domain-table rows visible to the caller.

const items = grok.dapi.domains.table('apitests.item');
const stamp = `AG-${Date.now()}`;
await items.batch([
  {sku: `${stamp}-1`, name: 'Bolt', quantity: 10},
  {sku: `${stamp}-2`, name: 'Bolt', quantity: 5},
  {sku: `${stamp}-3`, name: 'Nut', quantity: 7}]);

// JSON rows named by group column / measure alias
const totals = await items.aggregate({
  groupBy: ['name'],
  measures: [{fn: 'count'}, {fn: 'sum', column: 'quantity', as: 'total'}],
  filter: `sku starts "${stamp}"`,
  sort: '!total',
});
grok.shell.info(totals.map((t) => `${t.name}: ${t.count} rows, ${t.total} pcs`).join('\n'));

// The same aggregation as a typed DataFrame — feed grids and viewers directly
const df = await items.aggregateDf({
  groupBy: ['name'],
  measures: [{fn: 'sum', column: 'quantity', as: 'total'}],
  filter: `sku starts "${stamp}"`,
});
grok.shell.info(`aggregateDf: ${df.rowCount} rows, columns ${df.columns.names().join(', ')}`);

// Cleanup: one bounded filtered bulk delete
for (let guard = 0; guard < 100; guard++)
  if (!(await items.deleteWhere(`sku starts "${stamp}"`)).hasMore)
    break;
