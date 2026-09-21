//api: DG.DomainTableClient.update, DG.DomainTableClient.insert, DG.DomainVersionConflictError
// Writing through a domain table bound to an EXTERNAL database (the extwlive scratch binding,
// core/docs/features/ems/external-bindings/fixtures/extwlive, published on a dev stand over the
// `test` demo database). The key is the row's id and the client supplies it on insert; an update
// is guarded by `expected` — the values of the changed columns as last read — instead of a row
// version (access().support.concurrency === 'expected'), and a mismatch names what moved.

if (!(await grok.dapi.domains.schemas.list()).some((s) => s.name === 'extwlive'))
  return grok.shell.info('Publish the extwlive scratch binding first (it declares the extwlive domain schema)');

const things = grok.dapi.domains.table('extwlive.thing');
const tid = 900000 + Math.floor(Math.random() * 100000);
const [{id}] = await things.insert({tid, note: `sample ${tid}`, n: 1});
grok.shell.info(`inserted: id "${id}" is the key the client supplied`);
try {
  // the guard matches what was last read, so the update lands
  await things.update(id, {n: 2}, {expected: {n: 1}});
  grok.shell.info(`n is now ${(await things.get(id)).n}`);

  // a stale guard (n is 2 by now) rejects with a 409 and writes nothing
  try {
    await things.update(id, {n: 3}, {expected: {n: 1}});
  } catch (e) {
    if (!(e instanceof DG.DomainVersionConflictError))
      throw e;
    grok.shell.info(`conflict: expected ${JSON.stringify(e.body.expected)}, current ${JSON.stringify(e.body.current)}`);
  }
} finally {
  await things.delete(id);
}
