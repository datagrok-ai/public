// Generated low-level REST client (grok.dapi2) for domain tables; init it once per session.

grok.dapi2Init(grok.dapi.root, grok.dapi.token);

const items = grok.dapi.domains.table('apitests.item');
const sku = `D2-${Date.now()}`;
const [ins] = await items.insert({sku: sku, name: 'Gen'});

const rows = await grok.dapi2.domains.queryRows('apitests', 'item', {filter: `sku = "${sku}"`});
grok.shell.info(`${rows.length} row(s), name: ${rows[0].name}`);

// Cleanup
await items.delete(ins.id);
