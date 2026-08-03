// Query a domain table straight into a typed DataFrame (d42 wire format, 1M row cap
// vs 10k for query()). Columns carry db property tags, choices, and semantic types,
// so grids and property panels work out of the box.

const df = await grok.dapi.domains.table('apitests.item').queryDf({sort: '!created_on', limit: 100});
grok.shell.addTableView(df);

// Enrichment-style bulk fetch: rows for a set of ids as a typed DataFrame
// ('id' plus the requested fields; missing/invisible ids are simply absent rows)
const ids = [...Array(Math.min(3, df.rowCount)).keys()].map((i) => df.col('id').get(i));
const fields = await grok.dapi.domains.table('apitests.item').fetchFields(ids, ['sku', 'name']);
grok.shell.info(`fetchFields: ${fields.rowCount} rows, ${fields.columns.names().join(', ')}`);
