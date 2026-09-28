import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';

// The per-table change token (`DomainTableClient.version`): what a live list polls instead of
// re-counting its rows. `seq` moves once per write TRANSACTION that touched the table's rows —
// not once per row, and not at all for a write that rolled back — so "did anything change?"
// costs one indexed read per open list rather than an aggregate. Real server, real writes.
// The fixture is a throwaway user-managed schema, not the shared `apitests.item`: every
// assertion here is an EQUALITY on `seq + 1`, which another category writing the same table
// would break. The category skips cleanly without the CreateDomainSchema privilege.
category('Dapi: domain live', () => {
  const name = `zzl${`${Date.now()}`.slice(-8)}`;
  const items = () => grok.dapi.domains.table(`${name}.item`);
  const sku = () => `LV-${Date.now()}-${Math.floor(Math.random() * 1e6)}`;
  const seq = async () => (await items().version()).seq;
  let skip: string | null = null;

  before(async () => {
    try {
      await grok.dapi.domains.createSchema(name, {friendlyName: 'Live probe'});
    } catch (e: any) {
      if (e instanceof DG.DomainError && (e.code === 'forbidden' || e.status === 403)) {
        skip = 'no CreateDomainSchema privilege';
        return;
      }
      throw e;
    }
    await grok.dapi.domains.schema(name).apply({tables: {
      item: {
        businessKey: ['sku'],
        columns: {
          sku: {type: 'string', required: true, unique: true},
          name: {type: 'string', isName: true, searchable: true},
          quantity: {type: 'int', min: 0},
        },
      },
    }});
  });

  after(async () => {
    if (skip == null)
      await grok.dapi.domains.schema(name).delete();
  });

  const skipped = (): boolean => {
    if (skip != null)
      console.log(`skipped: ${skip}`);
    return skip != null;
  };

  // NB: no braces in a test name — the runner's --test matcher drops them.
  test('shape: a seq counter and the instant it last moved', async () => {
    if (skipped())
      return;
    const v = await items().version();
    expect(typeof v.seq, 'number', `version() must answer a numeric seq: ${JSON.stringify(v)}`);
    expect(v.seq >= 0, true, `seq must be a counter: ${JSON.stringify(v)}`);
    expect(v.at === null || !isNaN(Date.parse(v.at)), true,
      `at must be an ISO instant or null before the first write: ${JSON.stringify(v)}`);
    // Reading the token is not a write.
    expect((await items().version()).seq, v.seq, 'version() must not move the token it reports');
  });

  test('one write is one bump, and a three-op transaction is also one', async () => {
    if (skipped())
      return;
    const start = await items().version();
    const [ins] = await items().insert({sku: sku(), name: 'Live probe'});
    try {
      const afterInsert = await items().version();
      expect(afterInsert.seq, start.seq + 1, 'one insert must bump the token by exactly one');
      expect(afterInsert.at != null, true, 'a write must stamp `at`');
      expect(Date.parse(afterInsert.at!) >= Date.parse(start.at ?? afterInsert.at!), true,
        `at must move forward: ${start.at} -> ${afterInsert.at}`);
      // Three ops, one transaction, one bump — the reason the token beats a row counter.
      const res = await grok.dapi.domains.transaction(name, [
        {op: 'update', table: 'item', id: ins.id, values: {quantity: 1}},
        {op: 'update', table: 'item', id: ins.id, values: {quantity: 2}},
        {op: 'update', table: 'item', id: ins.id, values: {quantity: 3}},
      ]);
      expect(res.length, 3);
      expect(await seq(), afterInsert.seq + 1, 'a 3-op transaction must bump the token once');
    } finally {
      await items().delete(ins.id);
    }
  });

  test('a rolled-back write does not bump', async () => {
    if (skipped())
      return;
    const key = sku();
    const [ins] = await items().insert({sku: key, name: 'Live guard'});
    try {
      const start = await seq();
      // allOrNothing aborts the whole batch, so nothing was written and nothing changed.
      const report = await items().batch([{sku: `${key}-ok`, quantity: 1}, {sku: `${key}-bad`, quantity: -5}]);
      expect(report.error, 'validation', `the guard payload must abort: ${JSON.stringify(report)}`);
      expect(await seq(), start, 'an aborted batch must not bump the token');
      let refused = false;
      try {
        await items().update(ins.id, {quantity: -1});
      } catch (_) {
        refused = true;
      }
      expect(refused, true, 'the fixture refuses a negative quantity');
      expect(await seq(), start, 'a refused update must not bump the token');
    } finally {
      await items().delete(ins.id);
    }
  });

  test('an all-duplicate insert batch does not bump', async () => {
    if (skipped())
      return;
    const key = sku();
    const [ins] = await items().insert({sku: key, name: 'Live dedup'});
    try {
      const start = await seq();
      // Every row collides with the business key that is already there, so the loader
      // skips all of them: nothing was written, and the token must not move.
      const report = await items().batch([{sku: key, name: 'Live dedup'}]);
      expect(report.inserted, 0, `the payload must dedup: ${JSON.stringify(report)}`);
      expect(report.skipped, 1, JSON.stringify(report));
      expect(await seq(), start, 'an all-duplicate batch must not bump the token');
    } finally {
      await items().delete(ins.id);
    }
  });

  test('a validate-only batch does not bump', async () => {
    if (skipped())
      return;
    const key = sku();
    try {
      const start = await seq();
      const preview = await items().batch(
        [{sku: `${key}-1`, name: 'Preview', quantity: 1}], {validateOnly: true});
      expect(preview.validateOnly, true, `a preview answers its own shape: ${JSON.stringify(preview)}`);
      expect(preview.willInsert, 1, JSON.stringify(preview));
      expect(await seq(), start, 'a validate-only batch must leave the token where it was');
      expect((await items().query({filter: `sku starts "${key}"`})).length, 0,
        'a validate-only batch must write no rows');
    } finally {
      // A server that ignored the flag committed the row — take it back either way.
      await items().deleteWhere(`sku starts "${key}"`);
    }
  });

  test('two readers see each other through the token', async () => {
    if (skipped())
      return;
    // What a live list actually does: hold a seq, poll, and re-read only when it moved.
    const held = await items().version();
    expect((await items().version()).seq, held.seq, 'an idle table must not move the token');
    const [ins] = await items().insert({sku: sku(), name: 'Live notify'});
    try {
      const polled = await items().version();
      expect(polled.seq > held.seq, true,
        `the writer's insert must be visible to the reader's poll: ${held.seq} -> ${polled.seq}`);
    } finally {
      await items().delete(ins.id);
    }
    expect((await items().version()).seq > held.seq, true, 'the delete moved it further');
  });
}, {owner: 'askalkin@datagrok.ai'});
