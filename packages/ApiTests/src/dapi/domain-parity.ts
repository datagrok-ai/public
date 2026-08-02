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

  test('deleteWhere: filter forms, hasMore drain, empty filter rejects', async () => {
    const p = `DW${Date.now()}${Math.floor(Math.random() * 1e4)}`;
    await items().insert([1, 2, 3, 4, 5].map((i) => ({sku: `${p}-${i}`, name: `${p} bulk`, quantity: i})));
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
      const err = await thrown(() => projects.deleteWhere(`key starts "${key}"`));
      expect(err instanceof DG.DomainRestrictError, true,
        `expected DomainRestrictError, got ${err?.constructor?.name}: ${err?.message}`);
      expect(await projects.count(`key starts "${key}"`), 2, 'restrict must delete nothing');
    } finally {
      await issues.delete(issue.id);
      await projects.deleteWhere(`key starts "${key}"`);
    }
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
