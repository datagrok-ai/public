import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';

// WO-2 (GROK-20599): parity methods on DomainTableClient — count/exists/first/
// getByKey/fetchFields/aggregateDf/upsert/auditLog/watch. Fixture: the package's
// own 'apitests' schema; the audit:false row-watch probe rides the Inventory
// schema (stock_movements declares "audit": false) and skips cleanly where it
// is not deployed — the domain-errors.ts / domain-handlers.ts pattern.
category('Dapi: domain parity', () => {
  const items = () => grok.dapi.domains.table('apitests.item');
  const marker = `PAR${Date.now()}${Math.floor(Math.random() * 1e4)}`;
  const sku = () => `${marker}-${Math.floor(Math.random() * 1e9)}`;
  const ids: string[] = [];

  async function thrown(action: () => Promise<any>): Promise<any> {
    try {
      await action();
      return null;
    } catch (e) {
      return e;
    }
  }

  before(async () => {
    // Two rows sharing a name (getByKey ambiguity), one distinct.
    const ins = await items().insert([
      {sku: sku(), name: `${marker} dup`, quantity: 10},
      {sku: sku(), name: `${marker} dup`, quantity: 20},
      {sku: sku(), name: `${marker} solo`, quantity: 30},
    ]);
    for (const r of ins)
      ids.push(r.id);
  });

  after(async () => {
    for (const id of ids)
      await items().delete(id).catch(() => {});
  });

  test('count and exists: string and tree filters agree', async () => {
    const byString = await items().count(`name = "${marker} dup"`);
    const byTree = await items().count([{property: 'name', operator: '=', value: `${marker} dup`}]);
    expect(byString, 2);
    expect(byTree, 2);
    expect(await items().exists([{property: 'name', operator: '=', value: `${marker} solo`}]), true);
    expect(await items().exists(`name = "${marker} nothing-here"`), false);
  });

  test('first: match honors sort, no match resolves null', async () => {
    const top = await items().first({filter: `name like "${marker}%"`, sort: '!quantity'});
    expect(top?.quantity, 30);
    expect(await items().first({filter: `name = "${marker} nothing-here"`}), null);
  });

  test('getByKey: composite equality hits, ambiguity resolves null', async () => {
    const hit = await items().getByKey({name: `${marker} dup`, quantity: 20});
    expect(hit?.id, ids[1]);
    expect(await items().getByKey({name: `${marker} dup`}), null, 'ambiguous key must resolve null');
  });

  test('upsert: same key twice → inserted then updated; typed validation failures', async () => {
    const s = sku();
    const first = await items().upsert({sku: s, name: 'up v1', quantity: 1});
    ids.push(first.id);
    expect(first.status, 'inserted');
    const second = await items().upsert({sku: s, name: 'up v2', quantity: 2});
    expect(second.status, 'updated');
    expect(second.id, first.id);
    const row = await items().get(first.id);
    expect(row.name, 'up v2');
    expect(row.quantity, 2);
    expect(row.version, 2, 'upsert update must bump the version');

    // Business-key column absent from the payload → report-less 400 validation.
    const noKey = await thrown(() => items().upsert({name: 'no sku'}));
    expect(noKey instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${noKey?.constructor?.name}: ${noKey?.message}`);

    // Invalid value → report-CARRYING failure (batch resolves it; upsert must still throw).
    const bad = await thrown(() => items().upsert({sku: sku(), name: 'bad', quantity: -5}));
    expect(bad instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${bad?.constructor?.name}: ${bad?.message}`);
    expect(bad.rows.length > 0, true, 'per-row errors must ride the envelope');
  });

  test('fetchFields: requested shape, absent ids, empty ids', async () => {
    const df = await items().fetchFields(
      [ids[0], ids[1], '00000000-0000-0000-0000-000000000000'], ['sku']);
    expect(df.rowCount, 2, 'missing ids must be absent rows');
    expect(df.columns.names().sort().join(','), 'id,sku', 'exactly id + requested fields');
    const empty = await items().fetchFields([]);
    expect(empty.rowCount, 0, 'empty ids must resolve to an empty DataFrame');
    // Current contract: an empty id list short-circuits BEFORE the projection,
    // so the frame is zero-column even with fields requested.
    const emptyWithFields = await items().fetchFields([], ['sku']);
    expect(emptyWithFields.rowCount, 0);
    expect(emptyWithFields.columns.length, 0, 'empty ids yield a zero-column frame');
  });

  test('aggregateDf: typed frame matches aggregate JSON', async () => {
    const spec = {
      groupBy: ['name'],
      measures: [{fn: 'count', as: 'n'}, {fn: 'sum', column: 'quantity', as: 'total'}],
      filter: `name = "${marker} dup"`,
    } as _DG.DomainAggregateSpec;
    const json = await items().aggregate(spec);
    const df = await items().aggregateDf(spec);
    expect(df.rowCount, 1);
    expect(df.col('n')!.type, DG.TYPE.INT);
    expect(df.col('n')!.get(0), json[0]['n']);
    expect(df.col('total')!.get(0), json[0]['total']);
    expect(json[0]['n'], 2);
    expect(json[0]['total'], 30);
  });

  test('watch round-trip: table and row scope', async () => {
    expect(await items().isWatching(), false);
    expect(await items().watch(), true);
    try {
      expect(await items().isWatching(), true);
    } finally {
      expect(await items().unwatch(), true);
    }
    expect(await items().isWatching(), false);

    expect(await items().isWatching(ids[0]), false);
    expect(await items().watch(ids[0]), true);
    try {
      expect(await items().isWatching(ids[0]), true);
    } finally {
      expect(await items().unwatch(ids[0]), true);
    }
    expect(await items().isWatching(ids[0]), false);
  });

  test('row watch on an audit:false table: DomainValidationError', async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === 'inventory')) {
      console.log('skipped: inventory schema (audit:false fixture) not deployed');
      return;
    }
    // The audit gate fires before row resolution, so a synthetic id suffices.
    const movements = grok.dapi.domains.table('inventory.stock_movements');
    const err = await thrown(() => movements.watch('00000000-0000-0000-0000-000000000000'));
    expect(err instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${err?.constructor?.name}: ${err?.message}`);
  });

  test('query builder: thenable chains match the spec form; terminals work', async () => {
    const p = `QB${Date.now()}${Math.floor(Math.random() * 1e4)}`;
    const ins = await items().insert(
      [1, 2, 3].map((i) => ({sku: `${p}-${i}`, name: `${p} b`, quantity: i * 10})));
    for (const r of ins)
      ids.push(r.id);
    const viaBuilder = await items().query()
      .where('name', '=', `${p} b`).orderBy('quantity', true).top(2);
    const viaSpec = await items().query({filter: [{property: 'name', operator: '=', value: `${p} b`}],
      sort: '!quantity', limit: 2});
    expect(viaBuilder.map((r: any) => r.id).join(','), viaSpec.map((r) => r.id).join(','),
      'builder and spec form must produce identical results');
    expect(viaBuilder[0].quantity, 30);
    expect(await items().query().where({name: `${p} b`}).count(), 3);
    expect(await items().query().where({name: `${p} b`, quantity: 20}).exists(), true);
    expect(await items().query().where('quantity', '>', 1000000000).first(), null);
    const df = await items().query().where({name: `${p} b`}).df();
    expect(df.rowCount, 3, '.df() must return a DataFrame over the same query');
    // skip + select combine: the middle row by quantity, projection narrowed.
    const [mid] = await items().query()
      .where({name: `${p} b`}).orderBy('quantity').select('quantity').skip(1).top(1);
    expect(mid.quantity, 20, 'skip(1) after orderBy must land on the middle row');
    expect((mid as any).name, undefined, 'unselected columns must be absent at runtime');
  });

  test('predicate helpers bind values the string grammar cannot express', async () => {
    const name = `O'Brien ${Date.now()}`;
    const [r] = await items().insert({sku: sku(), name: name});
    ids.push(r.id);
    // Apostrophes are inexpressible in the smart-filter string grammar — the
    // condition tree binds them server-side, so they just work.
    const viaWhere = await items().query().where('name', '=', name).top(2);
    expect(viaWhere.length, 1);
    expect(viaWhere[0].name, name);
    const viaHelpers = await items().query({filter: DG.or(DG.cond('name', '=', name),
      DG.cond('quantity', '>', 1000000000))});
    expect(viaHelpers.length, 1);
    // Mixing a raw string filter with conditions is a clear client-side error
    // — in BOTH orders (string→condition and condition→string).
    let err: any = null;
    try {
      items().query().where('name = "x"').where('quantity', '>', 1);
    } catch (e) {
      err = e;
    }
    expect(`${err}`.includes('cond()'), true, `${err}`);
    let err2: any = null;
    try {
      items().query().where('quantity', '>', 1).where('name = "x"');
    } catch (e) {
      err2 = e;
    }
    expect(`${err2}`.includes('cond()'), true, `${err2}`);
  });

  test('deleteWhere: filter forms, hasMore drain, empty filter rejects', async () => {
    const p = `DW${Date.now()}${Math.floor(Math.random() * 1e4)}`;
    const ins = await items().insert(
      [1, 2, 3, 4, 5].map((i) => ({sku: `${p}-${i}`, name: `${p} bulk`, quantity: i})));
    // The suite's after() drains leftovers by id if the probe fails mid-way.
    for (const r of ins)
      ids.push(r.id);
    // String form, limited: the two oldest go, more remain.
    const first = await items().deleteWhere(`name = "${p} bulk"`, {limit: 2});
    expect(first.deleted, 2);
    expect(first.hasMore, true);
    // Tree form drains the rest via the hasMore loop.
    let total = first.deleted;
    for (let guard = 0; guard < 5; guard++) {
      const r = await items().deleteWhere(
        [{property: 'name', operator: '=', value: `${p} bulk`}], {limit: 2});
      total += r.deleted;
      if (!r.hasMore)
        break;
    }
    expect(total, 5);
    expect(await items().count(`name = "${p} bulk"`), 0);

    const err = await thrown(() => items().deleteWhere([]));
    expect(err instanceof DG.DomainValidationError, true,
      `expected DomainValidationError, got ${err?.constructor?.name}: ${err?.message}`);
  });

  test('deleteWhere: restrict rejects the whole call with zero deletions', async () => {
    const schemas = await grok.dapi.domains.schemas.list();
    if (!schemas.some((s) => s.name === 'grit')) {
      console.log('skipped: grit schema (restrict fixture) not deployed');
      return;
    }
    const projects = grok.dapi.domains.table('grit.project');
    const issues = grok.dapi.domains.table('grit.issue');
    const key = `DW${Date.now() % 10000000}`;
    // The free project is OLDER than the blocked one, so the per-row loop
    // deletes it first — the restrict rollback must undo it too.
    await projects.insert({key: `${key}A`, name: 'DW free'});
    const [blocked] = await projects.insert({key: `${key}B`, name: 'DW blocked'});
    const [issue] = await issues.insert({project_id: blocked.id, number: 1, title: 'blocks'});
    try {
      // Assert the premise: the free project really is the older one, so the
      // per-row loop reaches (and must roll back) a completed delete.
      const ordered = await projects.query({filter: `key starts "${key}"`, sort: 'created_on'});
      expect(ordered[0].key, `${key}A`, 'premise: the free project must be older');
      const err = await thrown(() => projects.deleteWhere(`key starts "${key}"`));
      expect(err instanceof DG.DomainRestrictError, true,
        `expected DomainRestrictError, got ${err?.constructor?.name}: ${err?.message}`);
      expect(await projects.count(`key starts "${key}"`), 2, 'restrict must delete nothing');
    } finally {
      await issues.delete(issue.id);
      await projects.deleteWhere(`key starts "${key}"`);
    }
  });

  test('details-expand child datetimes are dayjs with the correct instant', async () => {
    // The details json_agg path used to emit NAIVE datetime strings (no Z);
    // the server now serializes the same treat-as-UTC Z-form as top-level
    // rows, and detailDatetimeColumns materializes them recursively.
    const withDetails = grok.dapi.domains.table('apitests.item', {
      datetimeColumns: ['created_on', 'updated_on'],
      detailDatetimeColumns: {'item_event': ['created_on', 'updated_on']},
    });
    const [item] = await withDetails.insert({sku: sku(), name: 'dt probe'});
    ids.push(item.id);
    const [ev] = await grok.dapi.domains.table('apitests.item_event')
      .insert({item_id: item.id, kind: 'dt', amount: 1});
    const [expanded] = await withDetails.query(
      {filter: [{property: 'id', operator: '=', value: item.id}], expand: ['details:item_event']});
    const child = (expanded as any).item_event[0];
    expect(typeof child.created_on, 'object', 'child created_on must be dayjs, not a string');
    // The INSTANT matches the same row fetched top-level — guards a tz shift
    // from naive-form handling.
    const events = grok.dapi.domains.table('apitests.item_event',
      {datetimeColumns: ['created_on', 'updated_on']});
    const top = await events.get(ev.id);
    expect(child.created_on.valueOf(), top.created_on.valueOf(),
      'expanded child instant must equal the top-level read');
  });

  test('save: insert-or-update by identity; stale saves reject typed', async () => {
    const saved = await items().save({sku: sku(), name: 'save probe', quantity: 1});
    ids.push(saved.id);
    expect(saved.version, 1, 'insert path must resolve version 1');
    const again = await items().save({...saved, quantity: 2});
    expect(again.version, 2, 'update path must carry the fresh version back');
    expect((await items().get(saved.id)).quantity, 2);
    // Re-saving the OLD object carries its stale version — typed conflict, not a lie.
    const stale = await thrown(() => items().save({...saved, quantity: 3}));
    expect(stale instanceof DG.DomainVersionConflictError, true,
      `expected DomainVersionConflictError, got ${stale?.constructor?.name}: ${stale?.message}`);
  });

  test('save: a business-key duplicate resolves the existing id with NO version', async () => {
    const s = sku();
    const first = await items().save({sku: s, name: 'dup v', quantity: 1});
    ids.push(first.id);
    await items().update(first.id, {quantity: 2}, {version: 1});   // real version is now 2
    const dup = await items().save({sku: s, name: 'dup v2'});      // id-less, key exists
    expect(dup.id, first.id, 'duplicate save must resolve the existing id');
    expect((dup as any).version == null, true,
      'no fabricated version on a duplicate (the server reports none)');
    // The follow-up save is an unversioned update — no phantom-v1 conflict.
    const after = await items().save({...dup, quantity: 5});
    expect(after.version, 3);
    expect((await items().get(first.id)).quantity, 5);
  });

  test('optimistic retry helpers converge under a forced conflict', async () => {
    const [ins] = await items().insert({sku: sku(), name: 'retry probe', quantity: 0});
    ids.push(ins.id);
    // retryOnVersionConflict: the first attempt loses to an out-of-band bump
    // fired between its fresh read and its write; the second converges.
    let attempts = 0;
    const res = await DG.retryOnVersionConflict(async () => {
      attempts++;
      const fresh = await items().get(ins.id);
      if (attempts === 1)
        await items().update(ins.id, {quantity: 1}, {version: fresh.version});
      return await items().update(ins.id, {quantity: 2}, {version: fresh.version});
    });
    expect(attempts, 2, 'the first attempt must conflict, the second converge');
    expect(res.version, 3);

    // updateWithRetry: mutate sees the fresh row; null skips; unknown id rejects typed.
    const upd = await items().updateWithRetry(ins.id, (fresh) => ({quantity: (fresh.quantity ?? 0) + 5}));
    expect(upd?.version, 4);
    expect((await items().get(ins.id)).quantity, 7);
    expect(await items().updateWithRetry(ins.id, () => null), null, 'null mutate must skip the write');
    const nf = await thrown(() =>
      items().updateWithRetry('00000000-0000-0000-0000-000000000000', () => ({quantity: 1})));
    expect(nf instanceof DG.DomainNotFoundError, true,
      `expected DomainNotFoundError, got ${nf?.constructor?.name}: ${nf?.message}`);
  });

  test('auditLog: newest first, typed entries', async () => {
    const [ins] = await items().insert({sku: sku(), name: 'audit probe', quantity: 1});
    ids.push(ins.id);
    await items().update(ins.id, {quantity: 2}, {version: 1});
    const log = await items().auditLog({limit: 200});
    expect(typeof log[0].id, 'number', 'audit id is a wire number');
    for (let i = 1; i < log.length; i++)
      expect(log[i - 1].id >= log[i].id, true, 'audit entries must be newest first');
    const mine = log.filter((a) => a.row_id === ins.id);
    expect(mine.length >= 2, true, `expected insert+update entries, got ${mine.length}`);
    expect(mine[0].op, 'update', 'the update must precede the insert (newest first)');
  });
});
