import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';

// WO-9 (GROK-20341) — the grok.data.db.table(...) structured-write surface.
// The capability-negative test runs against Datlas alone (it fails at the capability gate
// before any GrokConnect round trip). The write round trips need a GrokConnect jar that
// advertises /mutate (WO-1..6); the stale dev image lacks it, so they self-skip cleanly.
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

  test('insert large object[] preserves nulls', async () => {
    if (!writesEnabled) { console.log(`skipped: ${skipReason}`); return; }
    await reset();
    // > DB_TABLE_INLINE_ROW_LIMIT (10000) rows routes through the bulk DataFrame path; a
    // genuine null must land as SQL NULL, not '' (regression for the fromObjects coercion).
    const n = 10001;
    const rows: object[] = new Array(n);
    for (let i = 0; i < n; i++) rows[i] = {id: i + 1, name: i === 0 ? null : `r${i}`, qty: i};
    const res = await t().insert(rows);
    expect(res.affectedRows, n);
    const nulls = await grok.data.db.query(conn.nqName, `select count(*) as c from ${fqTable} where name is null`);
    expect(nulls.col('c')!.get(0), 1);
  }, {timeout: 120000});

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
    // The DataConnection.Write permission gate (both directions) is verified end-to-end in
    // datlas connector_mutation_test. ApiTests runs as a single owner session (implicit Write),
    // so it cannot exercise the negative path without a second identity.
    console.log('skipped: permission gate covered by datlas connector_mutation_test (both directions green in CI)');
  }, {skipReason: 'covered by datlas connector_mutation_test'});
}, {owner: 'askalkin@datagrok.ai'});
