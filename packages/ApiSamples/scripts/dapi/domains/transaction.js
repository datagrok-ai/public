// Atomic multi-table transaction over a domain schema: later ops reference
// earlier ops' row ids via '$<ref>'; any failure rolls everything back.
// With generated clients, <schema>Db.transaction types each result by its op
// (insert → DomainInsertResult, update → DomainUpdateResult) — no positional casts.

const key = `TX-${Date.now()}`;
const [item, event, updated] = await grok.dapi.domains.transaction('apitests', [
  {op: 'insert', table: 'item', ref: 'i', values: {sku: key, name: 'Widget', quantity: 1}},
  {op: 'insert', table: 'item_event', values: {item_id: '$i', kind: 'created'}},
  {op: 'update', table: 'item', id: '$i', values: {quantity: 2}, expectedVersion: 1},
]);
grok.shell.info(`item ${item.id}: event ${event.id}, version ${updated.version}`);

// A failing op rejects with a TYPED error carrying the failing op's index:
try {
  await grok.dapi.domains.transaction('apitests', [
    {op: 'insert', table: 'item_event', values: {item_id: item.id}}, // missing required 'kind'
  ]);
} catch (e) {
  grok.shell.info(`op #${e.opIndex} failed: ${e.code}`); // DomainValidationError
}

await grok.dapi.domains.table('apitests.item').delete(item.id); // cascades item_event
