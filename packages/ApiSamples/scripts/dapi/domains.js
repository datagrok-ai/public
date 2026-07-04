// Row CRUD over entity-mapped domain schemas that plugins declare
// in databases/<schema>/schema.json — here the ApiTests package's 'apitests' schema.

// List registered domain schemas and their tables
const schemas = await grok.dapi.domains.schemas.list();
grok.shell.info(schemas.map((s) => `${s.name} v${s.version}: ${s.tables.map((t) => t.name).join(', ')}`).join('\n'));

const items = grok.dapi.domains.table('apitests.item');

// Insert a row (columns + dynamic jsonb keys in one object); business-key duplicates are reported, not duplicated
const [ins] = await items.insert({sku: `SKU-${Date.now()}`, name: 'Widget', quantity: 5, note: 'demo'});

// Query with a smart filter; update with optimistic concurrency; read the audit trail
const rows = await items.query({filter: 'name = "Widget"', sort: '!created_on', limit: 10});
await items.update(ins.id, {quantity: 6}, {version: rows[0].version});
const audit = await items.audit(ins.id);
grok.shell.info(`ops: ${audit.map((a) => a.op).join(', ')}`);

// Soft-delete
await items.delete(ins.id);
