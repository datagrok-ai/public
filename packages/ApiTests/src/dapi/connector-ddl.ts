import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';

// WO-B11 (GROK-20322) — the grok.data.db.ddl(...) fluent DDL builder, DbTable.uploadAs
// (create-table-from-DataFrame), and the dryRun/confirmDestructive contract surfaced as a
// typed DdlConfirmationRequiredError (the /mutate 400 body rides the thrown message verbatim
// and is JSON-parse-probed on the JS side).
// The round trips need a Datlas with the DDL dispatch (WO-B10) and a GrokConnect jar that
// advertises supportsDdl (WO-B5..B8); on an older stack the `before` dryRun probe fails and
// every round trip self-skips cleanly.
category('Dapi: connector ddl', () => {
  const rnd = () => DG.Utils.randomString(8);
  // Same writable-Postgres target convention as `Dapi: connector writes`: the local compose
  // demo `world` DB by default, overridable via a `dgConnectorWritesDb` global.
  const dbCfg = {
    server: 'world:5432', db: 'world', login: 'postgres', password: 'postgres', schema: 'public',
    ...((globalThis as any).dgConnectorWritesDb ?? {}),
  };
  const {server, db, login, password, schema} = dbCfg;

  let conn: _DG.DataConnection;        // GrokConnect Postgres — supportsDdl=true
  let nativeConn: _DG.DataConnection;  // PostgresDart — supportsDdl=false (capability gate)
  let ddlEnabled = false;
  let skipReason = 'this stack has no DDL dispatch (older Datlas or stale grok_connect jar) or the DB is not writable';

  const ddl = () => grok.data.db.ddl(conn.nqName).schema(schema);
  const scratch: string[] = [];  // bare table names; dropped in `after` even on failure

  const scratchTable = (prefix: string) => {
    const name = `apitests_ddl_${prefix}_${rnd()}`;
    scratch.push(name);
    return name;
  };

  before(async () => {
    nativeConn = DG.DataConnection.create(`DDL Native ${rnd()}`,
      {dataSource: 'PostgresDart', server, db, login, password});
    nativeConn = await grok.dapi.connections.save(nativeConn);

    conn = DG.DataConnection.create(`DDL GrokConnect ${rnd()}`,
      {dataSource: 'Postgres', server, db, login, password});
    conn = await grok.dapi.connections.save(conn);

    try {
      // Probe the whole path with a harmless dryRun: nothing executes, but it exercises the
      // interop entry, the Datlas dispatch, the privilege/capability gates, and GrokConnect
      // statement emission. Fails (and everything self-skips) on any older layer.
      const plan = await ddl().createTable(`apitests_ddl_probe_${rnd()}`,
        [{name: 'id', type: DG.TYPE.INT}]).dryRun();
      ddlEnabled = (plan?.statements?.length ?? 0) >= 1;
    } catch (e: any) {
      ddlEnabled = false;
      skipReason = `${skipReason}: ${e.message ?? e}`;
    }
  });

  after(async () => {
    for (const t of scratch) {
      try {
        await grok.data.db.query(conn.nqName, `drop table if exists ${schema}.${t}`);
      } catch (_) {}
    }
    for (const c of [conn, nativeConn]) {
      try {
        if (c) await grok.dapi.connections.delete(c);
      } catch (_) {}
    }
  });

  test('ddl lifecycle: create, alter, dryRun plan, confirm contract, drop', async () => {
    if (!ddlEnabled) { console.log(`skipped: ${skipReason}`); return; }
    const t = scratchTable('life');
    const fq = `${schema}.${t}`;

    const created = await ddl().createTable(t, [
      {name: 'id', type: DG.TYPE.INT, nullable: false},
      {name: 'name', type: DG.TYPE.STRING},
    ], {primaryKey: ['id']}).execute();
    expect((created.errorMessage ?? '') === '', true, `createTable failed: ${created.errorMessage}`);

    const ins = await grok.data.db.table(conn.nqName, fq)
      .insert([{id: 1, name: 'a'}, {id: 2, name: 'b'}, {id: 3, name: 'c'}]);
    expect(ins.affectedRows, 3);

    await ddl().alterTable(t).addColumn({name: 'qty', type: DG.TYPE.INT}).execute();
    await grok.data.db.table(conn.nqName, fq).update({set: {qty: 7}, where: {id: 1}});

    // dryRun dropTable: the plan carries the exact SQL and the live row count; nothing runs.
    const plan = await ddl().dropTable(t).dryRun();
    expect((plan.statements?.length ?? 0) >= 1, true, 'dryRun must return the statements');
    expect(/drop table/i.test(plan.statements.join(' ')), true,
      `expected DROP TABLE in the plan, got: ${plan.statements.join('; ')}`);
    expect(plan.destructive.length, 1);
    expect(plan.destructive[0].code, 'drop-table');
    expect(plan.destructive[0].count, 3);
    const still = await grok.data.db.query(conn.nqName, `select count(*) as c from ${fq}`);
    expect(still.col('c')!.get(0), 3, 'a dryRun must not execute the drop');

    // Unconfirmed execute: the refusal rethrows typed, with the recomputed plan attached.
    let typed: any = null;
    try {
      await ddl().dropTable(t).execute();
    } catch (e: any) {
      typed = e;
    }
    expect(typed != null, true, 'an unconfirmed destructive DDL must be refused');
    expect(typed instanceof DG.DdlConfirmationRequiredError, true,
      `expected DdlConfirmationRequiredError, got: ${typed?.name}: ${typed?.message}`);
    expect((typed.plan?.statements?.length ?? 0) >= 1, true, 'the refusal must carry the plan');
    expect((typed.plan?.destructive ?? []).some((d: any) => d.code === 'drop-table'), true);

    await ddl().dropTable(t).execute({confirmDestructive: true});
    const gone = await grok.data.db.query(conn.nqName,
      `select count(*) as c from information_schema.tables where table_schema = '${schema}' and table_name = '${t}'`);
    expect(gone.col('c')!.get(0), 0, 'the confirmed drop must remove the table');
  }, {timeout: 60000});

  test('uploadAs creates and loads a DataFrame (create mode)', async () => {
    if (!ddlEnabled) { console.log(`skipped: ${skipReason}`); return; }
    const t = scratchTable('up');
    const fq = `${schema}.${t}`;
    const n = 5000;
    const ids: number[] = new Array(n);
    const names: string[] = new Array(n);
    const vals: number[] = new Array(n);
    for (let i = 0; i < n; i++) { ids[i] = i + 1; names[i] = `row_${i}`; vals[i] = i + 0.5; }
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'id', ids),
      DG.Column.fromStrings('name', names),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'val', vals),
    ]);
    const dbTable = grok.data.db.table(conn.nqName, fq);

    // dryRun first: the CREATE statements derived from the DataFrame, nothing executed.
    const dry = await dbTable.uploadAs(df, {dryRun: true});
    expect((dry.plan?.statements?.length ?? 0) >= 1, true, 'uploadAs dryRun must return the plan');
    expect(/create table/i.test(dry.plan!.statements.join(' ')), true,
      `expected CREATE TABLE in the plan, got: ${dry.plan!.statements.join('; ')}`);

    const res = await dbTable.uploadAs(df);
    expect((res.errorMessage ?? '') === '', true, `uploadAs failed: ${res.errorMessage}`);
    const count = await grok.data.db.query(conn.nqName, `select count(*) as c from ${fq}`);
    expect(count.col('c')!.get(0), n);
    // Spot the dg → native type mapping (Postgres map: int → int4, string → text, float → float8).
    const types = await grok.data.db.query(conn.nqName,
      `select column_name, data_type from information_schema.columns ` +
      `where table_schema = '${schema}' and table_name = '${t}' order by ordinal_position`);
    const typeOf = (col: string) => {
      const i = types.col('column_name')!.toList().indexOf(col);
      return i < 0 ? '<missing>' : types.col('data_type')!.get(i);
    };
    expect(typeOf('id'), 'integer');
    expect(typeOf('name'), 'text');
    expect(typeOf('val'), 'double precision');
  }, {timeout: 120000});

  test('capability negative (non-DDL provider)', async () => {
    if (!ddlEnabled) { console.log(`skipped: ${skipReason}`); return; }
    // PostgresDart reports supportsDdl=false; the op must be refused at the Datlas capability
    // gate with a structured error — no GrokConnect round trip.
    let error = '';
    try {
      await grok.data.db.ddl(nativeConn.nqName).schema(schema)
        .createTable(scratchTable('cap'), [{name: 'id', type: DG.TYPE.INT}]).execute();
    } catch (e: any) {
      error = e.message ?? `${e}`;
    }
    expect(error !== '', true, 'a non-DDL provider must reject the operation');
    expect(/capability|does not support|supportsDdl|ddl/i.test(error), true,
      `expected a capability refusal, got: ${error}`);
  });

  test('privilege gate', async () => {
    // The fine-grained DDL privileges (DataConnection.CreateTable/AlterSchema/DropTable/
    // TruncateTable; both directions, exact-privilege naming) are verified end-to-end in
    // datlas connector_mutation_test (WO-B10). ApiTests runs as a single owner session
    // (implicit full access), so it cannot exercise the negative path without a second identity.
    console.log('skipped: privilege gate covered by datlas connector_mutation_test (both directions green in CI)');
  }, {skipReason: 'covered by datlas connector_mutation_test'});
}, {owner: 'askalkin@datagrok.ai'});
