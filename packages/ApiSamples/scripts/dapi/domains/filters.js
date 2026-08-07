// Batched facets (filter-panel counts) and shareable saved filter presets on a domain table.

if (!(await grok.dapi.domains.schemas.list()).some((s) => s.name === 'apitests'))
  return grok.shell.info('Deploy the ApiTests package first (it declares the apitests domain schema)');

const items = grok.dapi.domains.table('apitests.item');
const stamp = `FL-${Date.now()}`;
await items.batch([
  {sku: `${stamp}-1`, name: 'Bolt', quantity: 10},
  {sku: `${stamp}-2`, name: 'Bolt', quantity: 5},
  {sku: `${stamp}-3`, name: 'Nut', quantity: 7}]);

// One round trip: category counts + a histogram + the row count, all computed
// under the filter (conditions on a facet's own column are stripped per facet).
const res = await items.facets({
  filter: [{property: 'sku', operator: 'like', value: `${stamp}%`}],
  facets: [
    {id: 'names', kind: 'categories', column: 'name'},
    {id: 'qty', kind: 'histogram', column: 'quantity', bins: 4, min: 0, max: 10},
    {id: 'n', kind: 'count'}],
});
const names = res.facets['names'].categories.map((c) => `${c.value}: ${c.filtered}/${c.total}`);
grok.shell.info(`${res.facets['n'].count} rows; ${names.join(', ')}; buckets ${res.facets['qty'].buckets}`);

// Saved presets are shareable entities carrying the panel's state maps.
const preset = await items.filters.save(`Bolts only ${stamp}`,
  {'Name': {type: 'categorical', column: 'name', selected: ['Bolt']}});
grok.shell.info((await items.filters.list()).map((f) => f.friendlyName).join(', '));
await items.filters.delete(preset.id);

// Cleanup: one bounded filtered bulk delete
for (let guard = 0; guard < 100; guard++)
  if (!(await items.deleteWhere(`sku starts "${stamp}"`)).hasMore)
    break;
