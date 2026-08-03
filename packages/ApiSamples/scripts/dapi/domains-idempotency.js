// Safe retries and optimistic concurrency on domain-table writes.
// 'plates.plate' (PlatesFixture package) declares "idempotency": true.

const plates = grok.dapi.domains.table('plates.plate');

// Retrying an insert with the same idempotencyKey (a UUID) never duplicates the row
const values = {barcode: `P-${Date.now()}`, row_count: 8, col_count: 12, idempotencyKey: crypto.randomUUID()};
const [first] = await plates.insert(values);
const [retry] = await plates.insert(values);
grok.shell.info(`created: ${first.created}; retry: ${retry.status}, same id: ${retry.id === first.id}`);

// Optimistic concurrency: pass the version you last read; a stale version rejects
// with a TYPED DomainVersionConflictError (discriminate by class, not message text)
const plate = await plates.get(first.id);
await plates.update(first.id, {row_count: 16}, {version: plate.version});
try {
  await plates.update(first.id, {row_count: 24}, {version: plate.version}); // stale
} catch (e) {
  if (e instanceof DG.DomainVersionConflictError)
    grok.shell.info(`conflict: expected ${e.expectedVersion}, current ${e.currentVersion}`);
}
// ...or let the client retry with a fresh read:
await plates.updateWithRetry(first.id, (fresh) => ({row_count: fresh.row_count + 8}));

await plates.delete(first.id);
