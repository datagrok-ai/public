// Server-side grouped aggregation over the domain-table rows visible to the caller.

const items = grok.dapi.domains.table('apitests.item');
const stamp = `AG-${Date.now()}`;
await items.batch([
  {sku: `${stamp}-1`, name: 'Bolt', quantity: 10},
  {sku: `${stamp}-2`, name: 'Bolt', quantity: 5},
  {sku: `${stamp}-3`, name: 'Nut', quantity: 7}]);

const totals = await items.aggregate({
  groupBy: ['name'],
  measures: [{fn: 'count'}, {fn: 'sum', column: 'quantity', as: 'total'}],
  filter: `sku starts "${stamp}"`,
  sort: '!total',
});
grok.shell.info(totals.map((t) => `${t.name}: ${t.count} rows, ${t.total} pcs`).join('\n'));

// Cleanup
for (const row of await items.query({filter: `sku starts "${stamp}"`}))
  await items.delete(row.id);
