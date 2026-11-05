import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {Enrichment, deleteEnrichment, getEnrichConfigsNames, getFkRefs, readEnrichConfig, saveEnrichment}
  from '../enrichment';

/* DB Explorer's column enrichment over the platform's own database (System:Datagrok): events.session_id and
   func_calls.session_id refer to users_sessions.id, events.event_type_id to event_types.id. An enrichment is
   applied through PowerPack:runEnrichment, which the Enrich pane's link calls and a reopened project runs again
   from the table's history. Every enrichment a test saves is deleted when it ends, with the directories saving
   made. */
const ENRICHMENTS = 'System:AppData/PowerPack/enrichments';

let conn: DG.DataConnection;

async function queryTable(sql: string): Promise<DG.DataFrame> {
  return await conn.query('enrichment-test', sql).executeTable();
}

/** An enrichment of public.<table>.<column> that left-joins public.<joined> on the key pair and selects its fields. */
function enrichment(name: string, db: string | null, table: string, column: string,
  joined: string, keys: [string, string], fields: string[]): Enrichment {
  return {name, keyDb: db, keySchema: 'public', keyTable: table, keyColumn: column,
    fields: [`public.${table}.${column}`, ...fields.map((f) => `public.${joined}.${f}`)],
    joins: [{leftTableName: `public.${table}`, rightTableName: `public.${joined}`, joinType: DG.JOIN_TYPE.LEFT,
      leftTableKeys: [keys[0]], rightTableKeys: [keys[1]]}]};
}

function names(e: Enrichment): Promise<string[]> {
  return getEnrichConfigsNames(conn.nqName, e.keySchema, e.keyTable, e.keyColumn, e.keyDb);
}

function config(e: Enrichment): Promise<Enrichment | null> {
  return readEnrichConfig(conn.nqName, e.keySchema, e.keyTable, e.keyColumn, e.name, e.keyDb);
}

/** The query that gives the column of public.<table> for each id put in place of KEYS. */
const lookup = (table: string, column: string) =>
  `select id as key, ${column} as value from ${table} where id in (KEYS)`;

async function apply(e: Enrichment, df: DG.DataFrame, localColumn?: string): Promise<void> {
  await DG.Func.find({package: 'PowerPack', name: 'runEnrichment'})[0].prepare({conn, schema: e.keySchema,
    table: e.keyTable, column: e.keyColumn, name: e.name, df, db: e.keyDb, localColumn}).call();
}

/** Deletes the enrichments that are still there, then the directories down to their columns that are left empty. */
async function cleanUp(...enrichments: Enrichment[]): Promise<void> {
  const dirs = new Set<string>();
  for (const e of enrichments) {
    if (await config(e))
      await deleteEnrichment(conn.nqName, e.keySchema, e.keyTable, e.keyColumn, e.name, e.keyDb);
    let dir = `${ENRICHMENTS}/${conn.nqName.replace(':', '_')}`;
    dirs.add(dir);
    for (const part of [e.keyDb, e.keySchema, e.keyTable, e.keyColumn].filter((p) => p))
      dirs.add(dir = `${dir}/${part}`);
  }
  for (const dir of [...dirs].sort((a, b) => b.length - a.length)) {
    if (await grok.dapi.files.exists(dir) && (await grok.dapi.files.list(dir)).length === 0)
      await grok.dapi.files.delete(dir);
  }
}

const plain = (col: DG.Column, i: number): string =>
  col.isNone(i) ? '' : col.type === DG.TYPE.DATE_TIME ? `${col.get(i).valueOf()}` : `${col.get(i)}`;

/** Every row of the table holds, in the column, the value the query gives for the row's key. The query selects
 * "key" and "value" for the keys put in place of KEYS. */
async function expectJoined(df: DG.DataFrame, column: string, key: string, sql: string): Promise<void> {
  const target = df.col(column);
  expect(target != null, true, `the table has no "${column}" column: ${df.columns.names().join(', ')}`);
  const keys = df.col(key)!;
  const values = [...new Set(keys.toList().filter((v) => v != null).map((v) => `'${v}'`))];
  const source = await queryTable(sql.replace('KEYS', values.join(', ')));
  const byKey = new Map<string, string>();
  for (let i = 0; i < source.rowCount; i++)
    byKey.set(`${source.get('key', i)}`, plain(source.col('value')!, i));
  for (let i = 0; i < df.rowCount; i++) {
    if (!keys.isNone(i))
      expect(plain(target!, i), byKey.get(`${keys.get(i)}`) ?? '', `"${column}" of row ${i} (${key} ${keys.get(i)})`);
  }
}

category('Enrichment', () => {
  let events: () => Promise<DG.DataFrame>;
  let db: string | null;

  before(async () => {
    conn = await grok.dapi.connections.filter('namespace = "System:" and shortName = "Datagrok"').first();
    events = () => queryTable('select * from events where session_id is not null order by event_time desc limit 50');
    db = (await events()).col('session_id')!.tags.get(DG.Tags.Db) || null;
  });

  /** An enrichment of events.session_id that selects these columns of the event's session. */
  const ofSessions = (fields: string[]) => enrichment(`test-sessions-${Date.now()}`, db, 'events', 'session_id',
    'users_sessions', ['session_id', 'id'], fields);

  test('Applied, an enrichment joins the columns it selects, row by row', async () => {
    const e = ofSessions(['started', 'user_id']);
    try {
      await saveEnrichment(e, conn.nqName);
      expect((await names(e)).includes(e.name), true, `${e.name} is not listed for events.session_id`);
      const df = await events();
      const columns = df.columns.names();
      await apply(e, df);
      expectArray(df.columns.names().filter((n) => !columns.includes(n)).sort(), ['started', 'user_id']);
      expect(df.rowCount, 50, 'rows after the enrichment');
      await expectJoined(df, 'started', 'session_id', lookup('users_sessions', 'started'));
      await expectJoined(df, 'user_id', 'session_id', lookup('users_sessions', 'user_id'));
    }
    finally {
      await cleanUp(e);
    }
  });

  test('Edited, an enrichment adds the column it selects now when applied again', async () => {
    const e = ofSessions(['started']);
    try {
      await saveEnrichment(e, conn.nqName);
      const df = await events();
      await apply(e, df);
      expect(df.col('type'), null, 'the type column before the edit');
      const edited = (await config(e))!;
      edited.fields = [...edited.fields, 'public.users_sessions.type'];
      await saveEnrichment(edited, conn.nqName);
      expectArray((await config(e))!.fields,
        ['public.events.session_id', 'public.users_sessions.started', 'public.users_sessions.type']);
      await apply(e, df);
      await expectJoined(df, 'type', 'session_id', lookup('users_sessions', 'type'));
    }
    finally {
      await cleanUp(e);
    }
  });

  test('Deleted, an enrichment is gone, the columns it added stay, and applied again it is reported', async () => {
    const e = ofSessions(['is_short']);
    const balloons: string[] = [];
    const sub = grok.events.onEvent('d4-balloon-shown').subscribe((a: any) =>
      balloons.push(`${a?.args?.type}: ${a?.args?.message}`));
    try {
      await saveEnrichment(e, conn.nqName);
      const df = await events();
      await apply(e, df);
      expect(df.col('is_short') != null, true, 'the is_short column after the enrichment');
      expect(await deleteEnrichment(conn.nqName, 'public', 'events', 'session_id', e.name, db), true, 'deleted');
      expect((await names(e)).includes(e.name), false, `${e.name} is still listed for events.session_id`);
      expect(await config(e), null, 'the deleted configuration');
      expect(df.col('is_short') != null, true, 'the is_short column after the deletion');
      const columns = df.columns.names();
      await apply(e, df);
      expectArray(df.columns.names(), columns);
      expect(balloons.some((b) => b.startsWith('error') && b.includes('Could not find enrichment')), true,
        `the balloons shown: ${balloons.join('; ')}`);
    }
    finally {
      sub.unsubscribe();
      DG.Balloon.closeAll();
      await cleanUp(e);
    }
  });

  test('A column the table already has comes in as result.<name>', async () => {
    const e = enrichment(`test-types-${Date.now()}`, db, 'events', 'event_type_id', 'event_types',
      ['event_type_id', 'id'], ['friendly_name', 'source']);
    try {
      await saveEnrichment(e, conn.nqName);
      const df = await events();
      const friendlyNames = df.col('friendly_name')!.toList();
      await apply(e, df);
      expect(df.col('result.friendly_name') != null, true, `columns: ${df.columns.names().join(', ')}`);
      expectArray(df.col('friendly_name')!.toList(), friendlyNames);
      await expectJoined(df, 'result.friendly_name', 'event_type_id', lookup('event_types', 'friendly_name'));
      await expectJoined(df, 'source', 'event_type_id', lookup('event_types', 'source'));
    }
    finally {
      await cleanUp(e);
    }
  });

  test('A foreign key is offered the enrichments of the column it refers to, not those of its siblings', async () => {
    const sessions = enrichment(`test-users-${Date.now()}`, db, 'users_sessions', 'id', 'users', ['user_id', 'id'],
      ['login']);
    const sibling = ofSessions(['started']);
    try {
      await saveEnrichment(sessions, conn.nqName);
      await saveEnrichment(sibling, conn.nqName);
      const refs = await getFkRefs(conn, 'public', 'func_calls', 'session_id', db);
      expectArray(refs.map((r) => `${r.schema}.${r.table}.${r.column}`), ['public.users_sessions.id']);
      const offered = await getEnrichConfigsNames(conn.nqName, refs[0].schema, refs[0].table, refs[0].column, db);
      expect(offered.includes(sessions.name), true, `offered via users_sessions: ${offered.join(', ')}`);
      expect(offered.includes(sibling.name), false, `offered via users_sessions: ${offered.join(', ')}`);
      const own = await getEnrichConfigsNames(conn.nqName, 'public', 'func_calls', 'session_id', db);
      expect(own.includes(sibling.name), false, `the enrichments of func_calls.session_id: ${own.join(', ')}`);
      const calls = await queryTable('select * from func_calls where session_id is not null limit 20');
      await apply(sessions, calls, 'session_id');
      await expectJoined(calls, 'login', 'session_id', 'select s.id as key, u.login as value ' +
        'from users_sessions s join users u on u.id = s.user_id where s.id in (KEYS)');
    }
    finally {
      await cleanUp(sessions, sibling);
    }
  });
});
