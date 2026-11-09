/**
 * NCA EMS: platform — the S0 probes of nca-ems-storage-migration, kept as permanent platform-
 * regression tests: each pins a platform behaviour the repository relies on, so a core change that
 * breaks one fails here first, not in a scientist's Save.
 *
 * Each probe asserts the behaviour the repository relies on and returns what it measured (read it
 * with `grok test --csv`). They write throwaway rows labelled `zz-probe-…` into the real `nca`
 * schema: rows are soft-deleted at the end (a package schema has no row purge — decision 13), and
 * every snapshot file written here is deleted.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {ncaDb} from '../generated/db';
import {loadSampleByName} from '../views/nca-launch';
import {computeSnapshotHash, newId, snapshotDiff} from '../adapter/nca-repository';

const TIMEOUT = 180000;
const PROBE = 'zz-probe';
const FILE_ROOT = 'System:DomainFiles/nca/analysis';

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Bag = Record<string, any>;

/** One `=` condition. Loosely typed on purpose: the probes also filter through expanded refs
 *  (`study_id.kind`), which the generated column unions do not list. */
// eslint-disable-next-line @typescript-eslint/no-explicit-any
const eq = (property: string, value: string | number | boolean | null | string[]): any =>
  ({property, operator: '=' as const, value});

/** A study + one analysis + one result row, for the probes that need a tree. */
async function seed(label: string, kind: 'study' | 'standalone' = 'study'): Promise<{sid: string; aid: string}> {
  const [s] = await ncaDb.studies.insert({label: `${PROBE} ${label}`, kind});
  const [a] = await ncaDb.analyses.insert({
    study_id: s.id, label: `${PROBE} ${label} A`, source_kind: 'stored', recipe: {probe: true}, recipe_version: 1,
  });
  await ncaDb.analysisResults.insert({analysis_id: a.id, study_id: s.id, ord: 0});
  return {sid: s.id, aid: a.id};
}

/** The step-0 `string_list` values (ems-modification): empty, one item, many items, an empty
 *  string item, and items holding commas, quotes, a backslash, braces and non-ASCII — the shapes a
 *  list-literal encoding would mangle. */
const LISTS: string[][] = [
  [], ['one'], ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'], [''], ['x', '', 'y'],
  ['a,b', 'say "hi"', 'it\'s', 'back\\slash', '{braces}', 'Ü-12 日本 λz', 'NULL', ' padded '],
];

async function softDelete(studyIds: string[]): Promise<void> {
  for (const id of studyIds)
    try {await ncaDb.studies.delete(id);} catch {/* already deleted */}
}


/** The PR-9 / PR-18 edge-case table. */
function edgeTable(): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Subject', ['007', 'Ü-12', '日本', '007', 'x y']),
    DG.Column.fromInt32Array('ConcInt', Int32Array.from([1, 2, 0, -5, 7])),
    DG.Column.fromFloat64Array('Conc', Float64Array.from([1e-308, 1e308, 5e-324, -0, 0.1 + 0.2])),
    DG.Column.fromStrings('BLQ', ['<LLOQ', 'BLQ', '', '0.5', 'BLQ']),
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'Taken',
      [null, '2026-01-01T08:00:00Z', '2026-01-01T08:30:00Z', null, '2026-01-02T00:00:00Z']),
    // No null here: a BigInt null next to a zero is the known d42 gap (probe PR-18 below).
    DG.Column.fromBigInt64Array('Big',
      BigInt64Array.from([BigInt(7), BigInt('9007199254740993'), BigInt(1), BigInt(-2), BigInt(0)])),
  ]);
  df.col('Conc')!.set(3, null); // a missing float next to -0
  df.col('ConcInt')!.set(2, null);
  df.col('Conc')!.setTag('units', 'ng/mL');
  df.name = 'edge';
  return df;
}

category('NCA EMS: platform', () => {
  test('AC-0.1: the nca schema is deployed; a business-key row upserts and reads by key', async () => {
    const {sid, aid} = await seed('ac01');
    try {
      const first = await ncaDb.pipelineFlags.upsert({analysis_id: aid, profile_key: 'P1', verdict: 'Review'});
      const again = await ncaDb.pipelineFlags.upsert({analysis_id: aid, profile_key: 'P1', verdict: 'Accepted'});
      expect(again.status, 'updated', 'a second upsert must update, not insert');
      expect(again.id, first.id, 'the upsert merges into the same row');
      const row = await ncaDb.pipelineFlags.getByKey({analysis_id: aid, profile_key: 'P1'});
      expect(row?.verdict, 'Accepted', 'getByKey must return the upserted row');
      return `upsert → ${again.status}; getByKey ok`;
    } finally {await softDelete([sid]);}
  }, {timeout: TIMEOUT});

  // ems-modification step 0 (research G9.4 / core#746: `string_list` single-row writes were broken
  // on core 1.28). Every list column of the schema relies on this.
  test('step 0: string_list round-trips on single-row insert, update and in a transaction', async () => {
    const {sid, aid} = await seed('lists');
    const bad: string[] = [];
    try {
      const show = (v: unknown): string => JSON.stringify(v);
      // An absent list reads back as null or [] — the mappers read both as [] (recorded here).
      const [blank] = await ncaDb.pipelineFlags.insert({analysis_id: aid, profile_key: 'blank', verdict: 'Review'});
      const blankBack = (await ncaDb.pipelineFlags.get(blank.id) as Bag).reasons;
      for (let i = 0; i < LISTS.length; i++) {
        const v = LISTS[i];
        // (a) single-row insert
        const [ins] = await ncaDb.pipelineFlags.insert({analysis_id: aid, profile_key: `i${i}`, verdict: 'Review',
          reasons: v});
        const insBack = (await ncaDb.pipelineFlags.get(ins.id) as Bag).reasons ?? [];
        if (show(insBack) !== show(v)) bad.push(`insert ${show(v)} → ${show(insBack)}`);
        // (b) single-row update, from another value
        await ncaDb.pipelineFlags.update(ins.id, {info: ['seed']} as never);
        await ncaDb.pipelineFlags.update(ins.id, {info: v} as never);
        const updBack = (await ncaDb.pipelineFlags.get(ins.id) as Bag).info ?? [];
        if (show(updBack) !== show(v)) bad.push(`update ${show(v)} → ${show(updBack)}`);
        // (c) a transaction insert (the repository's write path)
        const [tx] = await grok.dapi.domains.transaction('nca', [{op: 'insert', table: 'pipeline_flag',
          values: {analysis_id: aid, profile_key: `t${i}`, verdict: 'Review', pipeline_only: v}}]) as Bag[];
        const txBack = (await ncaDb.pipelineFlags.get(tx.id) as Bag).pipeline_only ?? [];
        if (show(txBack) !== show(v)) bad.push(`transaction ${show(v)} → ${show(txBack)}`);
      }
      expect(bad.length, 0, `string_list does not round-trip: ${bad.join(' | ')}`);
      return `${LISTS.length} lists × insert / update / transaction ok; an unset list reads back as ${show(blankBack)}`;
    } finally {await softDelete([sid]);}
  }, {timeout: TIMEOUT});

  test('step 0 / AC-8.1: immutable columns (string_list included) refuse a change; the supersede update lands',
    async () => {
      const {sid, aid} = await seed('immutable');
      try {
        const [d] = await ncaDb.reviewDecisions.insert({analysis_id: aid, study_id: sid, unit_key: 'S1',
          action: 'exclude', reason_codes: ['A', 'B'], reason_texts: ['a, with comma', 'b'],
          decided_by: grok.shell.user.id, decided_by_label: 'u',
          decided_at: '2026-01-01T00:00:00.000Z', fingerprint: 'f'});
        const refused = async (values: Bag): Promise<boolean> => {
          try {await ncaDb.reviewDecisions.update(d.id, values as never); return false;} catch {return true;}
        };
        expect(await refused({reason_codes: ['A', 'C']}), true, 'a changed immutable string_list must be refused');
        expect(await refused({fingerprint: 'g'}), true, 'a changed immutable string must be refused');
        expect(await refused({decided_at: '2026-02-01T00:00:00.000Z'}), true, 'a changed immutable datetime refused');
        await ncaDb.reviewDecisions.update(d.id,
          {superseded_by: 'reset', superseded_at: '2026-01-02T00:00:00.000Z'} as never);
        // The same value again (what a full-row write would send) must not count as a change.
        await ncaDb.reviewDecisions.update(d.id, {reason_codes: ['A', 'B'], fingerprint: 'f'} as never);
        const back = await ncaDb.reviewDecisions.get(d.id) as Bag;
        expect(JSON.stringify(back.reason_codes), JSON.stringify(['A', 'B']), 'the list is unchanged');
        expect(back.superseded_by, 'reset', 'the supersede landed');
        return 'immutable list / string / datetime refused; supersede + same-value rewrite accepted';
      } finally {await softDelete([sid]);}
    }, {timeout: TIMEOUT});

  // Measured 2026-11-02 (core 1.28 bleeding-edge): a `user` column takes a uuid only, but does NOT
  // check that the user exists — a well-formed id of no user lands. The id comes from
  // grok.shell.user, so the gap does not bite the app; it is pinned here so a core change shows.
  test('IR-12: set_by / decided_by are user columns — a user id lands, a name is refused, no existence check',
    async () => {
      const {sid, aid} = await seed('user-ref');
      try {
        const me = grok.shell.user.id;
        const lz = (setBy: string, ordinal: number): Promise<unknown> => ncaDb.lzOverrides.insert({analysis_id: aid,
          profile_key: 'S1', time_key: '1', ordinal, conc_key: '2', action: 'force-out', set_by: setBy,
          set_by_label: 'x', set_at: '2026-01-01T00:00:00.000Z'});
        const [ok] = await lz(me, 0) as Bag[];
        expect((await ncaDb.lzOverrides.get(ok.id) as Bag).set_by, me, 'the user id reads back');
        const refused = async (f: () => Promise<unknown>): Promise<boolean> => {
          try {await f(); return false;} catch {return true;}
        };
        expect(await refused(() => lz('Analyst One', 1)), true, 'a display name is refused (a user column)');
        expect(await refused(() => ncaDb.reviewDecisions.insert({analysis_id: aid, study_id: sid, unit_key: 'S1',
          action: 'accept', reason_codes: [], reason_texts: [], decided_by: 'Analyst One', decided_by_label: 'x',
          decided_at: '2026-01-01T00:00:00.000Z', fingerprint: 'f'})), true, 'the same for a decision');
        const ghost = !(await refused(() => lz('00000000-0000-0000-0000-000000000000', 2)));
        return `a user id lands; a name is refused on both tables; an id of no user ${ghost ? 'lands (no FK)' :
          'is refused (FK now enforced)'}`;
      } finally {await softDelete([sid]);}
    }, {timeout: TIMEOUT});

  test('AC-0.2 / PR-3 / PR-9 / DG-B: a d42 snapshot in DomainFiles round-trips the frame and its hash',
    async () => {
      const frames: DG.DataFrame[] = [edgeTable()];
      for (const f of ['theoph.csv', 'indometh.csv', '04_mouse_sparse_destructive.csv',
        '05_multicompound_batch_screen-with-smiles.csv', 'demo/propranolol_multiroute.csv'])
        frames.push(await loadSampleByName(f));
      const lines: string[] = [];
      for (const df of frames) {
        const path = `${FILE_ROOT}/${newId()}/source.d42`;
        try {
          const t0 = performance.now();
          const bytes = df.toByteArray();
          await grok.dapi.files.write(path, bytes);
          const t1 = performance.now();
          const back = DG.DataFrame.fromByteArray(await grok.dapi.files.readAsBytes(path));
          const t2 = performance.now();
          const d = snapshotDiff(df, back);
          expect(d, null, `${df.name}: the reloaded frame differs — ${d}`);
          expect(computeSnapshotHash(back), computeSnapshotHash(df), `${df.name}: the snapshot hash changed on reload`);
          lines.push(`${df.name}: ${df.rowCount}×${df.columns.length}, d42 ${bytes.length} B vs csv ` +
          `${df.toCsv().length} B, write ${Math.round(t1 - t0)} ms, read ${Math.round(t2 - t1)} ms`);
        } finally {
          try {await grok.dapi.files.delete(path.slice(0, path.lastIndexOf('/')));} catch {/* gone */}
        }
      }
      return lines.join(' | ');
    }, {timeout: TIMEOUT});

  test('PR-18: BigInt through d42 — zeros survive; a nulled zero turns the other zeros missing (known gap)',
    async () => {
      const cells = (d: DG.DataFrame): string => [0, 1, 2, 3, 4].map((i) =>
      d.col('Big')!.isNone(i) ? 'null' : d.col('Big')!.getString(i)).join(' ');
      const make = (): DG.DataFrame => DG.DataFrame.fromColumns([DG.Column.fromBigInt64Array('Big',
        BigInt64Array.from([BigInt(0), BigInt('9007199254740993'), BigInt(1), BigInt(-2), BigInt(0)]))]);
      const plain = make();
      expect(cells(DG.DataFrame.fromByteArray(plain.toByteArray())), cells(plain), 'BigInt zeros must survive d42');
      // The shape the first AC-0.2 run failed on: a column holding 0 at row 0, then nulled there.
      const gappy = make();
    gappy.col('Big')!.set(0, null);
    const back = cells(DG.DataFrame.fromByteArray(gappy.toByteArray()));
    // Known platform defect: the zero at row 4 reads back as missing. The snapshot writer refuses
    // such a frame (verify-on-write), so no data is lost silently. When this assertion starts
    // failing, the platform has fixed it.
    expect(back !== cells(gappy), true,
      'the d42 BigInt null/zero defect is fixed — close TODO EXT-DEP-D42-BIGINT-NULL');
    return `before: ${cells(gappy)} | after d42: ${back}`;
    }, {timeout: TIMEOUT});

  test('AC-0.3: one query returns an analysis with its study, its pipeline run and access', async () => {
    const {sid, aid} = await seed('ac03');
    try {
      const [run] = await ncaDb.pipelineRuns.insert({study_id: sid, run_id: 'r1', run_user: 'probe',
        mapping: '{"b":1,"a":2}'});
      await ncaDb.analyses.update(aid, {pipeline_run: run.id, origin_kind: 'pipeline'} as never);
      const rows = await ncaDb.analyses.query({
        filter: eq('study_id', sid), expand: ['study_id', 'pipeline_run'], withAccess: true, limit: 10,
      }) as Bag[];
      expect(rows.length, 1);
      const r = rows[0];
      expect(r['study_id.label'], `${PROBE} ac03`, 'the study columns must come with the expand');
      expect(r['pipeline_run.run_id'], 'r1', 'the run columns must come with the second expand');
      expect(r['pipeline_run.mapping'], '{"b":1,"a":2}', 'a string column keeps JSON text verbatim');
      expect(r['~can_edit'], true, 'the author must be able to edit');
      // Child → parent filter through the expanded ref (the Ungrouped list reads standalone analyses so).
      const byKind = await ncaDb.analyses.query({
        filter: [eq('study_id.kind', 'study'), eq('id', aid)], expand: ['study_id'], limit: 10});
      expect(byKind.length, 1, 'a filter on study_id.kind must work');
      return `study + run expand ok, ~can_edit ${r['~can_edit']}, ~can_share ${r['~can_share']}, ` +
        'filter on study_id.kind ok';
    } finally {await softDelete([sid]);}
  }, {timeout: TIMEOUT});

  // ems-modification: the reads the repository issues IN PARALLEL rely on these filters.
  test('child rows filter through a ref path and by an id list; a deleted analysis\'s rows are not listed',
    async () => {
      const a = await seed('paths');
      const b = await seed('paths-b');
      try {
        await ncaDb.lzOverrides.insert({analysis_id: a.aid, profile_key: 'S1', time_key: '1', ordinal: 0,
          action: 'force-out'});
        // pipelineBase: the overrides of every analysis of a study, without first reading the analyses.
        const viaPath = await ncaDb.lzOverrides.query({filter: eq('analysis_id.study_id', a.sid), limit: 10});
        expect(viaPath.length, 1, 'a filter through analysis_id.study_id must work');
        // loadCurves: rows by an id list (= ANY).
        const [r2] = await ncaDb.analysisResults.insert({analysis_id: b.aid, study_id: b.sid, ord: 1,
          curve: {series: []}});
        const byIds = await ncaDb.analysisResults.query({filter: eq('id', [r2.id]), columns: ['curve'],
          limit: 10}) as Bag[];
        expect(byIds.length === 1 && Array.isArray(byIds[0].curve?.series), true, 'an id-list filter must work');
        // A deleted analysis: its result rows (no soft delete, not stamped by the cascade) are not listed.
        await ncaDb.analyses.delete(b.aid);
        const left = await ncaDb.analysisResults.query({filter: eq('study_id', b.sid), limit: 10});
        expect(left.length, 0, 'the result rows of a deleted analysis must not be listed by study');
        return 'ref-path filter ok; id list ok; a deleted analysis\'s result rows are not listed';
      } finally {await softDelete([a.sid, b.sid]);}
    }, {timeout: TIMEOUT});

  // The results write is diff-gated on the values read back: a float column that narrowed a double
  // would make every reopen look like a change (and bump the analysis version).
  test('float columns keep a double exactly; a datetime reads back as the same instant; "" reads back null',
    async () => {
      const {sid, aid} = await seed('float');
      try {
        const values = [0.1 + 0.2, 5e-324, 1e308, -0, 123456.789012345678, 1 / 3];
        const bad: string[] = [];
        for (const v of values) {
          const [r] = await ncaDb.analysisResults.insert({analysis_id: aid, study_id: sid, ord: 1, cmax_mean: v});
          const back = (await ncaDb.analysisResults.get(r.id) as Bag).cmax_mean;
          if (back !== v && !(v === 0 && back === 0)) bad.push(`${v} → ${back}`);
        }
        expect(bad.length, 0, `float columns do not keep a double: ${bad.join(', ')}`);
        const at = '2026-03-04T05:06:07.089Z';
        await ncaDb.analyses.update(aid, {confirmed_at: at} as never);
        const got = (await ncaDb.analyses.get(aid) as Bag).confirmed_at;
        const iso = typeof got?.toISOString === 'function' ? got.toISOString() : String(got);
        expect(iso, at, 'a datetime must read back as the same instant');
        // A scalar empty string is stored as null (the diff gate treats them as one value — ems-nca-rows
        // `norm`). When this starts failing, '' survives and `norm` can keep the two apart.
        await ncaDb.analyses.update(aid, {tier_flag_summary: ''} as never);
        const empty = (await ncaDb.analyses.get(aid) as Bag).tier_flag_summary;
        expect(empty === null || empty === undefined, true, `'' now reads back as ${JSON.stringify(empty)}`);
        return `${values.length} doubles exact; datetime reads back as ${typeof got} ${iso}; ` +
          `'' → ${JSON.stringify(empty)}`;
      } finally {await softDelete([sid]);}
    }, {timeout: TIMEOUT});

  test('a merged insert (business-key duplicate) in a transaction: its $ref names the existing row', async () => {
    const {sid, aid} = await seed('dupref');
    try {
      const [run] = await ncaDb.pipelineRuns.insert({study_id: sid, run_id: 'r1'});
      const res = await grok.dapi.domains.transaction('nca', [
        {op: 'insert', table: 'pipeline_run', ref: 'r', values: {study_id: sid, run_id: 'r1', run_user: 'again'}},
        {op: 'update', table: 'analysis', id: aid, values: {pipeline_run: '$r'}},
      ]) as Bag[];
      const back = await ncaDb.analyses.get(aid) as Bag;
      expect(back.pipeline_run, run.id, 'the $ref of a merged insert must resolve to the existing row');
      return `merged insert status ${res[0].status}; $ref → the existing run`;
    } finally {await softDelete([sid]);}
  }, {timeout: TIMEOUT});

  test('AC-0.4: update + delete of one row in one transaction; re-parent then delete the container', async () => {
    const a = await seed('ac04a', 'standalone');
    const b = await seed('ac04b');
    try {
      await grok.dapi.domains.transaction('nca', [
        {op: 'update', table: 'analysis', id: a.aid, values: {study_id: b.sid}},
        {op: 'update', table: 'study', id: a.sid, values: {moved_to: b.sid}},
        {op: 'delete', table: 'study', id: a.sid},
      ]);
      const moved = await ncaDb.analyses.get(a.aid);
      expect(moved !== null && moved.study_id === b.sid, true, 'the re-parented analysis must stay live in the target');
      const gone = (await ncaDb.studies.query({filter: eq('id', a.sid), deleted: 'only', limit: 1})) as Bag[];
      expect(gone.length === 1 && gone[0].moved_to === b.sid, true,
        'the deleted container must keep the values written in the same transaction');
      return 'analysis live in target; container deleted with moved_to';
    } finally {await softDelete([a.sid, b.sid]);}
  }, {timeout: TIMEOUT});

  // The restore probes (deleted children, restore order) were dropped with `delete_batch`: delete is
  // final since BD-16, so the app never restores.
  test('AC-0.5 / PR-14: paging, a light select without the json column, large json', async () => {
    const {sid, aid} = await seed('ac05');
    const out: string[] = [];
    try {
      // (c) Paging (explicit limit / offset) and a light select that leaves the curve out.
      for (let i = 1; i <= 3; i++)
        await ncaDb.analysisResults.insert({analysis_id: aid, study_id: sid, ord: i, curve: {series: []}});
      const page = await ncaDb.analysisResults.query({filter: eq('study_id', sid), sort: 'ord', limit: 2, offset: 2});
      expect(page.length, 2, 'offset paging must return the rest');
      const light = await ncaDb.analysisResults.query({filter: eq('study_id', sid), columns: ['ord', 'has_curve'],
        limit: 10}) as Bag[];
      expect(light.every((r) => r.curve === undefined), true, 'a light select must not carry the json column');
      out.push(`paging ok; light select keys: ${Object.keys(light[0]).sort().join(',')}`);

      // (d) A 5 MB json value.
      const big = {series: [{group: 'x'.repeat(5 * 1024 * 1024), time: [], geoMean: []}]};
      const t0 = performance.now();
      const [r] = await ncaDb.analysisResults.insert({analysis_id: aid, study_id: sid, ord: 9, curve: big});
      const back = await ncaDb.analysisResults.get(r.id) as Bag;
      expect((back?.curve?.series?.[0]?.group as string)?.length, 5 * 1024 * 1024, 'a 5 MB json must round-trip');
      out.push(`5 MB json round-trip ${Math.round(performance.now() - t0)} ms`);
      return out.join(' | ');
    } finally {await softDelete([sid]);}
  }, {timeout: TIMEOUT});

  test('S0: version() moves once per write; a stale expectedVersion conflicts', async () => {
    const {sid} = await seed('ver');
    try {
      const v0 = await ncaDb.studies.version();
      const row = await ncaDb.studies.get(sid);
      const u = await ncaDb.studies.update(sid, {updated_at: '2026-01-01T00:00:00Z'} as never, {version: row.version});
      const v1 = await ncaDb.studies.version();
      let conflict = false;
      try {
        await ncaDb.studies.update(sid, {updated_at: '2026-01-02T00:00:00Z'} as never, {version: row.version});
      } catch (e) {
        conflict = e instanceof DG.DomainVersionConflictError;
      }
      expect(conflict, true, 'a stale version must raise DomainVersionConflictError');
      return `table seq ${v0.seq} → ${v1.seq}; row version ${row.version} → ${u.version}; stale write → conflict`;
    } finally {await softDelete([sid]);}
  }, {timeout: TIMEOUT});
});
