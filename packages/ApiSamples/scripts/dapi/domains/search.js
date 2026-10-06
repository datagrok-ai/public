// Substring search over a domain table's SEARCHABLE columns — the ones declared
// `"searchable": true` in schema.json, else the table's name column. Case-insensitive,
// ANDed with `filter`; a table with neither rejects with a filter error.

const items = grok.dapi.domains.table('apitests.item');

const info = await grok.dapi.domains.registry.tableInfo('apitests.item');
console.log(`searchable columns: ${info.searchableColumns.join(', ')}`);

// The query spec, the builder, and the count all take the same `search`.
const rows = await items.query({search: 'alp', limit: 20});
const same = await items.query().search('alp').top(20);
const total = await items.count(undefined, {search: 'alp'});
grok.shell.info(`${rows.length} of ${total} rows match "alp" (builder: ${same.length})`);

// The d42 path narrows the same way.
const df = await items.queryDf({search: 'alp', limit: 100});
grok.shell.addTableView(df);
