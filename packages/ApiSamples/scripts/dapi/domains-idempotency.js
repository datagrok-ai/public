// Safe retries and optimistic concurrency on domain-table writes.
// 'plates.plate' (PlatesFixture package) declares "idempotency": true.

const plates = grok.dapi.domains.table('plates.plate');

// Retrying an insert with the same idempotencyKey (a UUID) never duplicates the row
const values = {barcode: `P-${Date.now()}`, row_count: 8, col_count: 12, idempotencyKey: crypto.randomUUID()};
const [first] = await plates.insert(values);
const [retry] = await plates.insert(values);
grok.shell.info(`created: ${first.created}; retry: ${retry.status}, same id: ${retry.id === first.id}`);

// Optimistic concurrency: pass the version you last read; a stale version fails the update
const plate = await plates.get(first.id);
await plates.update(first.id, {row_count: 16}, {version: plate.version});
try {
  await plates.update(first.id, {row_count: 24}, {version: plate.version}); // stale
} catch (e) {
  grok.shell.info(e.message); // Version conflict: expected 1, current 2
}

await plates.delete(first.id);
