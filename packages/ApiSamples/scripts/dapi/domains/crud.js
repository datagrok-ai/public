// Row CRUD over entity-mapped domain schemas that plugins declare
// in databases/<schema>/schema.json — here the ApiTests package's 'apitests' schema.
// In a package, `grok api` generates fully-typed clients (<schema>Db with Row/Insert/Column types).

// List registered domain schemas and their tables
const schemas = await grok.dapi.domains.schemas.list();
grok.shell.info(schemas.map((s) => `${s.name} v${s.version}: ${s.tables.map((t) => t.name).join(', ')}`).join('\n'));

const items = grok.dapi.domains.table('apitests.item');
const sku = `SKU-${Date.now()}`;

// Insert (columns + dynamic jsonb keys in one object); business-key duplicates are reported, not duplicated
const [ins] = await items.insert({sku: sku, name: 'Widget', quantity: 5, note: 'demo'});

// Reads: the fluent builder (awaitable) and the one-shot helpers
const rows = await items.query().where('name', '=', 'Widget').orderBy('created_on', true).top(10);
const byKey = await items.getByKey({sku: sku});        // business-key lookup; null on ambiguity
const n = await items.count({property: 'name', operator: '=', value: 'Widget'});
const has = await items.exists(`sku = "${sku}"`); // app-generated value — bind USER input via conditions instead
grok.shell.info(`${rows.length} rows, count ${n}, exists ${has}, byKey ${byKey?.id === ins.id}`);

// Update with optimistic concurrency; read the row and table audit trails
await items.update(ins.id, {quantity: 6}, {version: 1});
grok.shell.info(`ops: ${(await items.audit(ins.id)).map((a) => a.op).join(', ')}`);
grok.shell.info(`table-wide: ${(await items.auditLog({limit: 5})).length} recent events`);

// Promote the row to a full platform entity (sharing, favorites, comments)
await items.promote(ins.id);

// Watch the table (or one row) for change notifications
await items.watch();
grok.shell.info(`watching: ${await items.isWatching()}`);
await items.unwatch();

// Soft-delete
await items.delete(ins.id);
