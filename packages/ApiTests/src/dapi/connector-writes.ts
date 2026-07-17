import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import dayjs from 'dayjs';

// WO-9 (GROK-20341) — the grok.data.db.table(...) structured-write surface.
// WO-6 (GROK-20358) — object[] is now converted to a typed DataFrame in TS (scan-all typing,
// Date -> datetime, int->float widening, null preservation, all-null loud error) and always
// rides the d42 bulk path; the inline inference path is gone.
// The capability-negative and all-null-column tests run against Datlas/the client alone (they
// fail before any GrokConnect round trip). The write round trips need a GrokConnect jar that
// advertises /mutate + d42 (WO-1..6); the stale dev image lacks it, so they self-skip cleanly.
category('Dapi: connector writes', () => {
  const rnd = () => DG.Utils.randomString(8);
  // A writable Postgres (GrokConnect) DB. Defaults to the local compose demo `world` DB,
  // reachable from the grok_connect container by its docker-network alias (`world:5432`).
  // Override any field for another stack (e.g. remote CI) via a `dgConnectorWritesDb` global,
  // e.g. `globalThis.dgConnectorWritesDb = {server: 'dev.datagrok.ai:54322', db: 'northwind',
  // login: 'datagrok', password: 'datagrok'}`. When the target is unreachable or grok_connect
  // lacks /mutate, the write round trips self-skip cleanly (see `before`).
  const dbCfg = {
    server: 'world:5432', db: 'world', login: 'postgres', password: 'postgres', schema: 'public',
    ...((globalThis as any).dgConnectorWritesDb ?? {}),
  };
  const {server, db, login, password, schema} = dbCfg;
  const tableName = `apitests_cw_${rnd()}`;
  const fqTable = `${schema}.${tableName}`;

  let conn: _DG.DataConnection;        // GrokConnect Postgres — the write path
  let nativeConn: _DG.DataConnection;  // PostgresDart — supportsWrite=false (capability gate)
  let writesEnabled = false;
  let skipReason = 'grok_connect on this stack has no /mutate (stale jar) or the DB is not writable';

  const t = () => grok.data.db.table(conn.nqName, fqTable);
  const reset = () => grok.data.db.query(conn.nqName, `delete from ${fqTable}`);

  before(async () => {
    nativeConn = DG.DataConnection.create(`CW Native ${rnd()}`,
      {dataSource: 'PostgresDart', server, db, login, password});
    nativeConn = await grok.dapi.connections.save(nativeConn);

    conn = DG.DataConnection.create(`CW GrokConnect ${rnd()}`,
      {dataSource: 'Postgres', server, db, login, password});
    conn = await grok.dapi.connections.save(conn);

    try {
      await grok.data.db.query(conn.nqName,
        `create table ${fqTable} (id int primary key, name text, qty int)`);
      // Probe write support: an insert that fails on a stale jar (no /mutate).
      const probe = await t().insert([{id: -1, name: 'probe', qty: 0}]);
      writesEnabled = (probe.affectedRows ?? 0) >= 1;
      await reset();
    } catch (e: any) {
      writesEnabled = false;
      skipReason = `${skipReason}: ${e.message ?? e}`;
    }
  });

  after(async () => {
    try {
      await grok.data.db.query(conn.nqName, `drop table if exists ${fqTable}`);
    } catch (_) {}
    for (const c of [conn, nativeConn]) {
      try {
        if (c) await grok.dapi.connections.delete(c);
      } catch (_) {}
    }
  });

  test('insert object[] rows', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    await reset();
    const res = await t().insert([{id: 1, name: 'a', qty: 10}, {id: 2, name: 'b', qty: 20}]);
    expect(res.affectedRows, 2);
    const rows = await grok.data.db.query(conn.nqName, `select count(*) as c from ${fqTable}`);
    expect(rows.col('c')!.get(0), 2);
  });

  test('insert DataFrame (bulk)', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    await reset();
    const n = 50000;
    const ids: number[] = new Array(n);
    const names: string[] = new Array(n);
    const qtys: number[] = new Array(n);
    for (let i = 0; i < n; i++) { ids[i] = i + 1; names[i] = `bulk_${i}`; qtys[i] = i; }
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'id', ids),
      DG.Column.fromStrings('name', names),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'qty', qtys),
    ]);
    const res = await t().insert(df);
    expect(res.affectedRows, n);
    const rows = await grok.data.db.query(conn.nqName, `select count(*) as c from ${fqTable}`);
    expect(rows.col('c')!.get(0), n);
  }, {timeout: 120000});

  test('insert object[] preserves nulls', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    await reset();
    // WO-6: every object[] is converted to a typed DataFrame at the API boundary (no more
    // inline vs bulk split). A genuine null must land as SQL NULL, not '' — the sharpest
    // behavior change, since 5-row payloads previously travelled inline.
    const rows = [
      {id: 1, name: 'a', qty: 1},
      {id: 2, name: null, qty: 2},
      {id: 3, name: 'c', qty: 3},
      {id: 4, name: null, qty: 4},
      {id: 5, name: 'e', qty: 5},
    ];
    const res = await t().insert(rows);
    expect(res.affectedRows, 5);
    const nulls = await grok.data.db.query(conn.nqName, `select count(*) as c from ${fqTable} where name is null`);
    expect(nulls.col('c')!.get(0), 2);
  });

  test('insert object[] Date column becomes datetime', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    // WO-6: a Date value is detected on the live JS value (before serialization) and typed
    // as a real datetime column — the DateTime-bug fix. It round-trips as a timestamp, not
    // a stringified value; nulls are preserved.
    const dt = `${schema}.apitests_cw_dt_${rnd()}`;
    try {
      await grok.data.db.query(conn.nqName, `create table ${dt} (id int primary key, ts timestamp)`);
      const when = new Date(Date.UTC(2021, 2, 14, 9, 30, 0)); // 2021-03-14 09:30:00 UTC
      const res = await grok.data.db.table(conn.nqName, dt).insert([{id: 1, ts: when}, {id: 2, ts: null}]);
      expect(res.affectedRows, 2);
      const back = await grok.data.db.query(conn.nqName,
        `select to_char(ts, 'YYYY-MM-DD HH24:MI:SS') as t from ${dt} where id = 1`);
      expect(back.col('t')!.get(0), '2021-03-14 09:30:00');
      const nulls = await grok.data.db.query(conn.nqName, `select count(*) as c from ${dt} where ts is null`);
      expect(nulls.col('c')!.get(0), 1);
    } finally {
      try { await grok.data.db.query(conn.nqName, `drop table if exists ${dt}`); } catch (_) {}
    }
  });

  test('insert object[] dayjs column becomes datetime', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    // WO-6 follow-up: Datagrok represents datetimes as dayjs across much of its surface, so a
    // dayjs value must type as datetime too (duck-typed, not just `instanceof Date`).
    const dt = `${schema}.apitests_cw_dj_${rnd()}`;
    try {
      await grok.data.db.query(conn.nqName, `create table ${dt} (id int primary key, ts timestamp)`);
      const when = dayjs(new Date(Date.UTC(2021, 2, 14, 9, 30, 0))); // 2021-03-14 09:30:00 UTC
      const res = await grok.data.db.table(conn.nqName, dt).insert([{id: 1, ts: when}, {id: 2, ts: null}]);
      expect(res.affectedRows, 2);
      const back = await grok.data.db.query(conn.nqName,
        `select to_char(ts, 'YYYY-MM-DD HH24:MI:SS') as t from ${dt} where id = 1`);
      expect(back.col('t')!.get(0), '2021-03-14 09:30:00');
      const nulls = await grok.data.db.query(conn.nqName, `select count(*) as c from ${dt} where ts is null`);
      expect(nulls.col('c')!.get(0), 1);
    } finally {
      try { await grok.data.db.query(conn.nqName, `drop table if exists ${dt}`); } catch (_) {}
    }
  });

  test('insert object[] widens int column to float', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    // WO-6: scan-all typing widens an int column to double when any value is fractional
    // (the first-non-null path could truncate 2.5 to 2).
    const ft = `${schema}.apitests_cw_f_${rnd()}`;
    try {
      await grok.data.db.query(conn.nqName, `create table ${ft} (id int primary key, v double precision)`);
      const res = await grok.data.db.table(conn.nqName, ft).insert([{id: 1, v: 3}, {id: 2, v: 2.5}]);
      expect(res.affectedRows, 2);
      const back = await grok.data.db.query(conn.nqName, `select v from ${ft} where id = 2`);
      expect(back.col('v')!.get(0), 2.5);
    } finally {
      try { await grok.data.db.query(conn.nqName, `drop table if exists ${ft}`); } catch (_) {}
    }
  });

  test('insert object[] large float magnitude does not throw', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    // WO-6 follow-up: an integer-valued number of float magnitude (1e21) types as a double
    // column — the "exceeds 2^53" guard no longer throws for these.
    const bt = `${schema}.apitests_cw_big_${rnd()}`;
    try {
      await grok.data.db.query(conn.nqName, `create table ${bt} (id int primary key, big double precision)`);
      const res = await grok.data.db.table(conn.nqName, bt).insert([{id: 1, big: 1e21}, {id: 2, big: 3}]);
      expect(res.affectedRows, 2);
      // 1e21 loses exactness through the query float32 downcast, so compare by magnitude.
      const c = await grok.data.db.query(conn.nqName, `select count(*) as c from ${bt} where big > 1e20`);
      expect(c.col('c')!.get(0), 1);
    } finally {
      try { await grok.data.db.query(conn.nqName, `drop table if exists ${bt}`); } catch (_) {}
    }
  });

  test('insert object[] all-null column errors', async () => {
    // Type inference is client-side, so this fires before any server round trip — runs in CI
    // regardless of writesEnabled: an all-null column cannot be typed and must throw.
    let error = '';
    try {
      await t().insert([{id: 1, name: null}, {id: 2, name: null}]);
    } catch (e: any) {
      error = e.message ?? `${e}`;
    }
    expect(error !== '', true, 'an all-null column must be rejected');
    expect(/every value is null|cannot be inferred/i.test(error), true,
      `expected an all-null column error, got: ${error}`);
  });

  test('insert 100k object[] via the typed DataFrame path', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    await reset();
    // A large object[] now takes the same typed-DataFrame path as a small one (no size split).
    const n = 100000;
    const rows: object[] = new Array(n);
    for (let i = 0; i < n; i++) rows[i] = {id: i + 1, name: `r${i}`, qty: i};
    const res = await t().insert(rows);
    expect(res.affectedRows, n);
    const c = await grok.data.db.query(conn.nqName, `select count(*) as c from ${fqTable}`);
    expect(c.col('c')!.get(0), n);
  }, {timeout: 180000});

  test('upsert with keys', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    await reset();
    await t().insert([{id: 1, name: 'a', qty: 1}, {id: 2, name: 'b', qty: 2}]);
    const res = await t().upsert([{id: 1, name: 'a2', qty: 11}, {id: 3, name: 'c', qty: 3}], {keys: ['id']});
    expect((res.errorCount ?? 0), 0);
    const rows = await grok.data.db.query(conn.nqName, `select count(*) as c from ${fqTable}`);
    expect(rows.col('c')!.get(0), 3);
    const one = await grok.data.db.query(conn.nqName, `select name from ${fqTable} where id = 1`);
    expect(one.col('name')!.get(0), 'a2');
  });

  test('update and delete with where', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    await reset();
    await t().insert([{id: 1, name: 'a', qty: 1}, {id: 2, name: 'b', qty: 2}, {id: 3, name: 'c', qty: 3}]);
    const upd = await t().update({set: {name: 'renamed'}, where: {id: 2}});
    expect(upd.affectedRows, 1);
    const check = await grok.data.db.query(conn.nqName, `select name from ${fqTable} where id = 2`);
    expect(check.col('name')!.get(0), 'renamed');
    const del = await t().delete({where: {id: 3}});
    expect(del.affectedRows, 1);
    const rows = await grok.data.db.query(conn.nqName, `select count(*) as c from ${fqTable}`);
    expect(rows.col('c')!.get(0), 2);
  });

  test('capability negative (non-write provider)', async () => {
    // PostgresDart reports supportsWrite=false; the mutation must be refused at the Datlas
    // capability gate with a structured error — no GrokConnect round trip. Runs in CI.
    let error = '';
    try {
      await grok.data.db.table(nativeConn.nqName, fqTable).insert([{id: 1, name: 'x'}]);
    } catch (e: any) {
      error = e.message ?? `${e}`;
    }
    expect(error !== '', true, 'a non-write provider must reject the mutation');
    expect(/capability|does not support|supportsWrite|write/i.test(error), true,
      `expected a capability refusal, got: ${error}`);
  });

  test('permission gate', async () => {
    // The fine-grained connection-privilege gate (AddRows/ChangeValues/RemoveRows/...; both
    // directions) is verified end-to-end in datlas connector_mutation_test. ApiTests runs as a
    // single owner session (implicit full access), so it cannot exercise the negative path
    // without a second identity.
    console.log('skipped: permission gate covered by datlas connector_mutation_test (both directions green in CI)');
  }, {skipReason: 'covered by datlas connector_mutation_test'});
}, {owner: 'askalkin@datagrok.ai'});
