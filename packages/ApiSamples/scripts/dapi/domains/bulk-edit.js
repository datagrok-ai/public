// Bulk edit is one call: `updateWhere(filter, values)` patches every row the filter
// matches and the caller may edit, oldest first, in ONE transaction, capped at 1000
// rows (`limit` lowers it, `hasMore` says to loop). It rides the same per-row engine
// as update — writability, immutability, validation, audit — so a value the engine
// refuses refuses the WHOLE call and nothing is written.

const items = grok.dapi.domains.table('apitests.item');
const prefix = `BULK-${Date.now()}`;
const rows = await items.insert([0, 1, 2].map((i) =>
  ({sku: `${prefix}-${i}`, name: 'Bulk', quantity: 1, origin: 'seed'})));

// Over a selection — what a grid's "edit the selected rows" posts: the ids, bound
// server-side like any other filter value.
const picked = rows.slice(0, 2).map((r) => `"${r.id}"`).join(', ');
const selected = await items.updateWhere(`id in (${picked})`, {quantity: 7});
grok.shell.info(`updated ${selected.updated} selected rows, more matching: ${selected.hasMore}`);

// Over everything a filter matches, in pages: narrow the filter so a page of updated
// rows drops out of it, and loop while hasMore.
const pending = () => [{property: 'sku', operator: 'like', value: `${prefix}%`}, 'and',
  {property: 'quantity', operator: '=', value: 1}];
let updated = 0;
for (let report = {hasMore: true}; report.hasMore;) {
  report = await items.updateWhere(pending(), {quantity: 9}, {limit: 2});
  updated += report.updated;
}
grok.shell.info(`drained the rest in pages: ${updated} rows`);

// A refused value refuses the whole call: an immutable column that already has a
// different value, an unknown or service ('~') column, a relation name.
try {
  await items.updateWhere(`sku like "${prefix}%"`, {origin: 'rewritten'});
} catch (e) {
  grok.shell.info(`nothing written: ${e.rows?.[0]?.errors?.[0]?.message ?? e.message}`);
}

await items.deleteWhere(`sku like "${prefix}%"`);
