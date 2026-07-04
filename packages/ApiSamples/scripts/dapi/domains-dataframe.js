// Query a domain table straight into a typed DataFrame (d42 wire format, 1M row cap
// vs 10k for query()). Columns carry db property tags, choices, and semantic types,
// so grids and property panels work out of the box.

const df = await grok.dapi.domains.table('apitests.item').queryDf({sort: '!created_on', limit: 100});
grok.shell.addTableView(df);
