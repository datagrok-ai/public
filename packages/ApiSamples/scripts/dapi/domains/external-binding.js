//api: DG.DomainTableClient.access, DG.DomainTableClient.query, DG.DomainTableClient.get, DG.DomainUnsupportedError
// A domain table bound to an EXTERNAL database (the Northwind binding fixture,
// core/docs/features/ems/external-bindings/fixtures/northwind, published on a dev stand
// over the local demo container; the hosted Northwind has no primary keys). The registry, the
// grants, the query API and the u2 app are the same as for a platform-stored table; the
// differences arrive as declarations in
// access().support and access().fields — a client reads them, it never asks where the table lives.

if (!(await grok.dapi.domains.schemas.list()).some((s) => s.name === 'northwind'))
  return grok.shell.info('Publish the Northwind binding fixture first (it declares the northwind domain schema)');

const orders = grok.dapi.domains.table('northwind.order');
const access = await orders.access();
const s = access.support;

// Rows carry the REMOTE primary key as their id (percent-encoded, comma-joined for a composite key).
const [first] = await orders.query({limit: 1, sort: 'order_id'});
grok.shell.info(`First order: id "${first.id}" shipped to ${first.ship_country}`);
const lines = await grok.dapi.domains.table('northwind.order_detail').query({filter: `order_id = "${first.id}"`, limit: 3});
grok.shell.info(`${lines.length} lines; a line's id is its composite key, e.g. "${lines[0]?.id}"`);
const same = await orders.get(first.id);
grok.shell.info(same.id === first.id ? 'get(id) round-trips the key' : 'unexpected');

// What the table declares:
//   writes        false unless the binding says `"writable": true`
//   transaction   u2 saves through /transaction; false → the table is read-only in u2
//   concurrency   'expected' (old-value guard) instead of 'version' on a writable binding
//   captions      false: ref captions are resolved per row, not projected with the query
//   updateWhere   false: no set-based writes on a warehouse
//   filters       'basic': the offer leaves out under, regex, not-like, fuzzy, between, datetime !=, bool null tests
//   probe/deleted/audit/watch  false: no live poll, no trash, no history
grok.shell.info(`writes: ${s.writes}, transaction: ${s.transaction}, concurrency: ${s.concurrency}, ` +
  `captions: ${s.captions}, filters: ${s.filters}, probe: ${s.probe}`);

// An operation the storage cannot do answers DomainUnsupportedError (422) naming the op —
// no grant makes it succeed, so gate on `support` instead of catching.
try {
  await orders.query({limit: 1, deleted: 'only'});
} catch (e) {
  if (e instanceof DG.DomainUnsupportedError)
    grok.shell.info(`Refused by name: op "${e.op}" on ${e.body.table}`);
  else
    throw e;
}
