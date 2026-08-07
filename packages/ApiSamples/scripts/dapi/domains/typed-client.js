// The typed domain surface: fluent builder queries, bound condition helpers,
// upsert, save/updateWithRetry, and typed error classes (never match error messages).
// In a package, `grok api` generates fully-typed clients over the same methods.

const items = grok.dapi.domains.table('apitests.item');
const stamp = `TC-${Date.now()}`;

// upsert: insert-or-merge ONE row by the table's business key
const up1 = await items.upsert({sku: `${stamp}-1`, name: "O'Brien's widget", quantity: 1});
const up2 = await items.upsert({sku: `${stamp}-1`, quantity: 2});
grok.shell.info(`${up1.status} then ${up2.status}`); // inserted then updated

// Bound conditions express ANY value — apostrophes are impossible in the
// smart-filter STRING grammar, but condition values are bound server-side:
const found = await items.query().where('name', '=', "O'Brien's widget").top(5);
const firstMatch = await items.query().where('sku', '=', `${stamp}-1`).first(); // row or null
// ...and compose into trees with DG.cond / DG.and / DG.or:
const trees = await items.query({filter: DG.or(
  DG.cond('name', '=', "O'Brien's widget"), DG.cond('quantity', '>', 100))});
grok.shell.info(`builder: ${found.length}, helpers: ${trees.length}`);

// save(): insert-or-update by row identity; updateWithRetry(): read-modify-write
// that retries on a version conflict instead of failing
const saved = await items.save({sku: `${stamp}-2`, name: 'Saved', quantity: 1});
await items.save({...saved, quantity: 3});           // version-checked update
await items.updateWithRetry(saved.id, (fresh) => ({quantity: (fresh.quantity ?? 0) + 1}));

// Typed errors: discriminate by class/code, never by message text
try {
  await items.update(saved.id, {quantity: 0}, {version: 1}); // stale version
} catch (e) {
  if (e instanceof DG.DomainVersionConflictError)
    grok.shell.info(`conflict: expected ${e.expectedVersion}, current ${e.currentVersion}`);
}

// Cleanup: one bounded filtered bulk delete (never a per-row loop)
for (let guard = 0; guard < 100; guard++)
  if (!(await items.deleteWhere(`sku starts "${stamp}"`)).hasMore)
    break;
