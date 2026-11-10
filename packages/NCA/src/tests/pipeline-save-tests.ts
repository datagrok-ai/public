/**
 * `ncaSaveToStudy` — the PURE decision layer (no platform store): the "untouched" predicate, the
 * study plan, the pipeline flags, the source frame and the study-grid / results-grid wordings.
 * The study's identity is its EMS business key (owner, study key) and the plan's write is one
 * versioned transaction — both covered by the storage contract suite.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  makeAnalysisEntry, recipeOf, computeRecipeFingerprint, computeRecipeHash, RECIPE_HASH_VERSION,
  buildStudyGrid, FLAGS_COL, formatPipelineFlags,
  formatStudyAnalysesSummary, pipelineFlagCount, PER_ANIMAL_SPARSE_NOTE, sparseUnitKey,
} from '../adapter';
import type {AnalysisEntry, AnalysisRecipe, NcaStudy, PipelineOrigin, SummaryCache} from '../adapter';
import {
  pipelineUntouched, planStudyUpdate, pipelineFlagsOf, recipeFor, sourceFrame,
} from '../pipeline/pipeline-save';
import type {PlannedAnalysis} from '../pipeline/pipeline-save';
import {parseRules} from '../pipeline/pipeline-compute';
import {parseMapping} from '../pipeline/pipeline-mapping';
import {writePipelineColumns, PIPELINE_COL, LZ_PIPELINE_INFO_COL, pipelineTipAt} from '../views/results-view';

const M = parseMapping(JSON.stringify({
  subject: 'Subject', time: 'Time_hr', concentration: 'Conc_ng_mL', dose: 'Dose_mg_kg', route: 'Route',
  lloq: 'LLOQ', groupBy: ['Study', 'Compound'],
  units: {concentration: 'ng/mL', time: 'h', dose: 'mg/kg', bodyWeight: 'kg'},
}));
const SIG = {kind: 'content-hash' as const, rowCount: 1, mtime: null, hash: 'h', computedAt: 't'};
const CACHE: SummaryCache = {
  analysisTier: {nCompounds: 1, nSubjects: 1, route: 'PO', dose: '', status: 'ok',
    flagSummary: 'No flags'},
  compoundTier: [],
};

/** An untouched pipeline entry of `partition`, written by a run at `runAt`. */
function pipelineEntry(id: string, partition: string, runAt: string, sparse = false): AnalysisEntry {
  const recipe = recipeFor(M, '', sparse);
  const fp = computeRecipeFingerprint(recipe);
  const draft = makeAnalysisEntry({
    id, label: partition, source: {kind: 'platformPath', ref: `S/${id}.csv`}, sourceSignature: SIG, recipe,
    summaryCache: CACHE,
    workflow: {confirmed: {recipe, fingerprint: fp, at: runAt}, lzOverrides: [], reviewDecisions: [],
      stageCache: {stage: 'complete', pending: 0, notApplicable: 0}},
    origin: null, pipelineFlags: [],
  });
  const origin: PipelineOrigin = {
    kind: 'pipeline', runId: `run-${id}`, runAt, user: 'flow', engineVersion: 'x', mapping: {}, rules: '',
    partition, recipeHash: computeRecipeHash(recipeOf(draft)), recipeHashVersion: RECIPE_HASH_VERSION,
  };
  return {...draft, origin};
}

const study = (analyses: AnalysisEntry[]): NcaStudy => ({
  id: 'sid', label: 'HAL-1', kind: 'study', createdAt: 't', updatedAt: 't0',
  access: {canEdit: true, canDelete: true, canShare: true}, analyses, pipelineKey: {owner: 'o', studyKey: 'HAL-1'},
});

const planned = (partition: string, id: string, base = partition, sparse = false): PlannedAnalysis =>
  ({base, partition, sparse, entry: {...pipelineEntry(id, partition, '2026-10-15T10:00:00Z', sparse),
    source: {kind: 'platformPath', ref: `S/new-${id}.csv`}}});

const CTX = {studyLabel: 'HAL-1', pipelineKey: {owner: 'o', studyKey: 'HAL-1'}, now: '2026-10-15T10:00:00Z'};

category('Pipeline save: untouched (pure)', () => {
  test('AC-4.2: mutating ANY recipe field makes the entry touched (compile-enforced over AnalysisRecipe)', async () => {
    const e = pipelineEntry('a1', 'DGX-1001', '2026-10-14T10:00:00Z');
    expect(pipelineUntouched(e, 'DGX-1001'), true, 'as written');
    // A Record over EVERY key of AnalysisRecipe: a field added to the recipe fails to compile here
    // until its mutation is listed — so the predicate's coverage cannot silently fall behind.
    const mutate: Record<keyof AnalysisRecipe, (r: AnalysisRecipe) => AnalysisRecipe> = {
      columns: (r) => ({...r, columns: {...r.columns, subject: 'Animal'}}),
      modes: (r) => ({...r, modes: {...r.modes, nonmem: true}}),
      ncaRules: (r) => ({...r, ncaRules: {...r.ncaRules, extrapWarnPct: r.ncaRules.extrapWarnPct + 1}}),
      units: (r) => ({...r, units: {...r.units, time: 'min'}}),
      groupByCols: (r) => ({...r, groupByCols: ['Sex']}),
      passthroughCols: (r) => ({...r, passthroughCols: [{name: 'Smiles', isStructure: true}]}),
      acceptanceRSquared: (r) => ({...r, acceptanceRSquared: 0.5}),
    };
    for (const [field, fn] of Object.entries(mutate)) {
      const edited: AnalysisEntry = {...e, ...fn(recipeOf(e))};
      expect(pipelineUntouched(edited, 'DGX-1001'), false, `a changed '${field}' is the scientist's work`);
    }
  });

  test('AC-4.3 / 6d: overrides, any ledger, origin and version make the entry touched', async () => {
    const e = pipelineEntry('a1', 'DGX-1001', '2026-10-14T10:00:00Z');
    const override = {profileKey: 'R1', timeKey: '1', ordinal: 0, concKey: '1', action: 'force-out' as const,
      user: 'u', at: 't'};
    expect(pipelineUntouched({...e, lzOverrides: [override]}, 'DGX-1001'), false, 'λz override');
    expect(pipelineUntouched({...e, reviewDecisions: [{unitKey: 'R1', current: null, history: []}]}, 'DGX-1001'),
      false, 'a ledger reset to To review is still the scientist\'s record');
    expect(pipelineUntouched({...e, origin: null}, 'DGX-1001'), false, 'made in the app');
    expect(pipelineUntouched({...e, origin: {kind: 'standalone', label: 'x', containerId: 'c'}}, 'DGX-1001'),
      false, 'moved from Ungrouped');
    expect(pipelineUntouched(e, 'DGX-1002'), false, 'another partition');
    const o = e.origin as PipelineOrigin;
    expect(pipelineUntouched({...e, origin: {...o, recipeHashVersion: RECIPE_HASH_VERSION + 1}}, 'DGX-1001'), false,
      'a hash of another version cannot be compared');
    expect(pipelineUntouched({...e, confirmed: {...e.confirmed!, fingerprint: 'other'}}, 'DGX-1001'), false,
      'a saved, unconfirmed Setup edit');
    expect(pipelineUntouched({...e, confirmed: {...e.confirmed!, at: '2026-10-15T00:00:00Z'}}, 'DGX-1001'), true,
      'open + re-confirm (only the time changes) is not work');
  });
});

category('Pipeline save: study plan (pure)', () => {
  test('AC-4.1: an untouched entry is replaced in place, keeping its id and the STORED label', async () => {
    const old = {...pipelineEntry('a1', 'DGX-1001', '2026-10-14T10:00:00Z'), label: 'Renamed by Jane'};
    const plan = planStudyUpdate(study([old]), [planned('DGX-1001', 'n1')], CTX);
    expect(plan.study.analyses.length, 1);
    expect(plan.study.analyses[0].id, 'a1', 'same id');
    expect(plan.study.analyses[0].label, 'Renamed by Jane', 'a rename is kept');
    expect(plan.study.analyses[0].source.ref, 'S/new-n1.csv', 'the new source');
    expect(plan.outcomes[0].status, 'updated');
    expect(plan.replacedRefs.join(','), 'S/a1.csv', 'the old source is up for removal');
  });

  test('AC-3.5: a touched entry gets the run saved beside it', async () => {
    const old = {...pipelineEntry('a1', 'DGX-1001', '2026-10-14T10:00:00Z'),
      reviewDecisions: [{unitKey: 'R1', current: null, history: []}]};
    const plan = planStudyUpdate(study([old]), [planned('DGX-1001', 'n1')], CTX);
    expect(plan.study.analyses.length, 2);
    expect(plan.study.analyses[0], old, 'the scientist\'s entry is untouched');
    expect(plan.study.analyses[1].label, 'DGX-1001-2');
    expect(plan.outcomes[0].status, `saved beside 'DGX-1001' (has review work)`);
    expect(plan.replacedRefs.length, 0, 'nothing removed');
  });

  test('AC-3.4c / 6d-1: other entries of the stable read are kept; a hand-made same label is ignored', async () => {
    const handMade = {...pipelineEntry('h1', 'DGX-1001', '2026-10-14T10:00:00Z'), origin: null};
    const other = {...pipelineEntry('o1', 'DGX-2002', '2026-10-14T10:00:00Z'), label: 'renamed concurrently'};
    const plan = planStudyUpdate(study([handMade, other]), [planned('DGX-1001', 'n1')], CTX);
    expect(plan.study.analyses.map((a) => a.id).join(','), 'h1,o1,n1');
    expect(plan.study.analyses[1].label, 'renamed concurrently', 'a concurrent rename survives');
    expect(plan.study.analyses[2].label, 'DGX-1001-2', 'beside the hand-made one');
    expect(plan.outcomes[0].status, 'created', 'no pipeline candidate → created');
  });

  test('AC-3.7: a new study carries the pipeline key (its business key) and the label', async () => {
    const plan = planStudyUpdate(undefined, [planned('DGX-1001', 'n1')], CTX);
    expect(plan.study.label, 'HAL-1');
    expect(plan.study.pipelineKey?.owner, 'o');
    expect(plan.study.pipelineKey?.studyKey, 'HAL-1');
    expect(plan.outcomes[0].status, 'created');
  });

  test('AC-4.1: two untouched candidates → the latest run is replaced', async () => {
    const older = pipelineEntry('a1', 'DGX-1001', '2026-10-13T10:00:00Z');
    const newer = {...pipelineEntry('a2', 'DGX-1001', '2026-10-14T10:00:00Z'), label: 'DGX-1001-2'};
    const plan = planStudyUpdate(study([older, newer]), [planned('DGX-1001', 'n1')], CTX);
    expect(plan.outcomes[0].analysisId, 'a2');
  });

  test('AC-4.0 / F6: route partitions never match each other; a design flip is superseded, left in place', async () => {
    const po = pipelineEntry('p1', 'A · PO', '2026-10-14T10:00:00Z', true);
    const plan1 = planStudyUpdate(study([po]), [planned('A · PO', 'n1', 'A', true), planned('A · IV', 'n2', 'A', true)],
      CTX);
    expect(plan1.outcomes.map((o) => o.status).join(','), 'updated,created');
    expect(plan1.superseded.length, 0, 'a new route is not a design flip');

    const dense = pipelineEntry('d1', 'A', '2026-10-14T10:00:00Z');
    const plan2 = planStudyUpdate(study([dense]), [planned('A · PO', 'n3', 'A', true)], CTX);
    expect(plan2.study.analyses.length, 2, 'the old dense entry is left in place');
    expect(plan2.superseded.map((s) => `${s.label}→${s.by}`).join(','), 'A→A · PO');
    const added = plan2.study.analyses[1];
    expect(added.origin?.kind === 'pipeline' ? added.origin.supersedes : '', 'd1', 'the new entry records it');

    // A dense partition whose label merely starts with `A · ` (another groupBy) is not a design flip.
    const day1 = pipelineEntry('x1', 'A · Day1', '2026-10-14T10:00:00Z');
    const plan3 = planStudyUpdate(study([day1]), [planned('A', 'n4')], CTX);
    expect(plan3.superseded.length, 0, 'same design, prefix only');
  });
});

/** A minimal `ncaCompute` results frame for `pipelineFlagsOf`. */
function resultsFrame(rows: Array<{subject: string; verdict: string; reasons: string[]; dose: number;
  nPoints?: number}>): DG.DataFrame {
  const f64 = (name: string, v: number[]): DG.Column =>
    DG.Column.fromFloat64Array(name, Float64Array.from(v.map((x) => (Number.isNaN(x) ? DG.FLOAT_NULL : x))));
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('SUBJECT', rows.map((r) => r.subject)),
    DG.Column.fromStrings('VERDICT', rows.map((r) => r.verdict)),
    DG.Column.fromStrings('REASONS_JSON', rows.map((r) => JSON.stringify(r.reasons))),
    DG.Column.fromStrings('INFO_FLAGS_JSON', rows.map(() => '["routine note"]')),
    f64('DOSE', rows.map((r) => r.dose)),
    f64('SPAN_RATIO', rows.map(() => 3)), f64('LAMBDAZ_ADJ_R2', rows.map(() => 0.99)),
    f64('LAMBDAZ_NPOINTS', rows.map((r) => r.nPoints ?? 5)), f64('PCT_EXTRAP', rows.map(() => 5)),
  ]);
}

category('Pipeline save: flags, source frame, wordings (pure)', () => {
  test('AC-8.2: flags parse the JSON (a " | " survives); pipelineOnly drops the app\'s own reasons', async () => {
    const acc = parseRules('').acceptance;
    const nPts = `N points (λz) ≤ ${acc.lambdaMinPoints}`;
    const res = resultsFrame([
      {subject: 'R1', verdict: 'Review', reasons: [nPts, 'CRO comment: "spat | regurgitated"'], dose: 1,
        nPoints: acc.lambdaMinPoints},
      {subject: 'R2', verdict: 'Accepted', reasons: [], dose: 1},
    ]);
    const src = DG.DataFrame.fromColumns([DG.Column.fromFloat64Array('Dose_mg_kg', Float64Array.from([1]))]);
    const flags = pipelineFlagsOf(res, [0, 1], src, M, acc);
    expect(flags[0].profileKey, 'R1');
    expect(flags[0].reasons.length, 2);
    expect(flags[0].pipelineOnly.join(';'), 'CRO comment: "spat | regurgitated"',
      'the acceptance reason is the app\'s');
    expect(flags[0].info.join(';'), 'routine note');
    expect(flags[1].pipelineOnly.length, 0);
    expect(pipelineFlagCount(flags), 1);
  });

  test('AC-8.2: composite keys follow the app (raw dose): fractional, empty, 2 spellings; a miss throws', async () => {
    const acc = parseRules('').acceptance;
    const comp = (dose: number): DG.DataFrame =>
      resultsFrame([{subject: '(composite)', verdict: 'Review', reasons: ['x'], dose}]);
    const src = DG.DataFrame.fromColumns([DG.Column.fromStrings('Dose_mg_kg', ['0.3', '0.3', '10', '10.0', ''])]);
    expect(pipelineFlagsOf(comp(0.3), [0], src, M, acc).map((f) => f.profileKey).join(','),
      sparseUnitKey('DOSE=0.3'), 'fractional dose');
    // Two spellings of one dose exist only in a genuinely STRING dose column (a non-numeric cell
    // keeps DG from inferring the column numeric — guardrails, numeric-looking string columns).
    const spelled = DG.DataFrame.fromColumns([DG.Column.fromStrings('Dose_mg_kg', ['10', '10.0', 'not dosed'])]);
    expect(pipelineFlagsOf(comp(10), [0], spelled, M, acc).map((f) => f.profileKey).join(','),
      `${sparseUnitKey('DOSE=10')},${sparseUnitKey('DOSE=10.0')}`, 'two spellings → one flag per app key');
    expect(pipelineFlagsOf(comp(NaN), [0], src, M, acc)[0].profileKey, sparseUnitKey('DOSE='), 'missing dose');
    let threw = false;
    try {
      pipelineFlagsOf(comp(7), [0], src, M, acc);
    } catch {
      threw = true;
    }
    expect(threw, true, 'no orphan flag');
  });

  test('AC-7.1: an all-empty MAPPED column stays in the CSV, an unmapped one is dropped', async () => {
    const t = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Subject', ['S1', 'S1']), DG.Column.fromFloat64Array('Time_hr', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('Conc_ng_mL', Float64Array.from([1, 2])),
      DG.Column.fromFloat64Array('LLOQ', Float64Array.from([DG.FLOAT_NULL, DG.FLOAT_NULL])),
      DG.Column.fromStrings('Extra', ['', '']),
    ]);
    const src = sourceFrame(t, [0, 1], M);
    expect(src.col('LLOQ') !== null, true, 'mapped LLOQ kept');
    expect(src.col('Extra') === null, true, 'unmapped empty dropped');
  });

  test('AC-6.2: the sparse sampling type follows the mapping design', async () => {
    expect(recipeFor(M, '', true).modes.sparse.samplingType, '(auto)');
    expect(recipeFor({...M, design: 'sparse-batch'}, '', true).modes.sparse.samplingType, 'batch');
    expect(recipeFor({...M, design: 'sparse-destructive'}, '', true).modes.sparse.samplingType, 'destructive');
    expect(recipeFor(M, '', false).passthroughCols.length, 0, 'no pipeline passthrough (Q1)');
  });

  test('F9: the study-grid Warnings lines, in order, and the studies-list overlay', async () => {
    expect(formatPipelineFlags(0), '');
    expect(formatPipelineFlags(1), '1 pipeline flag');
    expect(formatPipelineFlags(3), '3 pipeline flags');
    expect(formatStudyAnalysesSummary([{stage: 'complete', pending: 0, notApplicable: 0}], 2),
      '1 analysis · 2 pipeline flags');
    const old = {...pipelineEntry('d1', 'A', 't'), label: 'A'};
    const sparse: AnalysisEntry = {...pipelineEntry('s1', 'A · PO', 't', true),
      stageCache: {stage: 'review', pending: 2, notApplicable: 0},
      pipelineFlags: [{profileKey: 'sparse:DOSE=1', verdict: 'Review', reasons: ['x'], info: [], pipelineOnly: ['x']}]};
    const s2 = {...sparse, origin: {...(sparse.origin as PipelineOrigin), supersedes: 'd1'}};
    const grid = buildStudyGrid(study([old, s2]));
    const flags = grid.col(FLAGS_COL)!;
    const lines = (i: number): string[] => flags.getString(i).split('\n');
    expect(lines(1).slice(0, 4).join(' / '),
      `2 not reviewed / 1 pipeline flag / ${PER_ANIMAL_SPARSE_NOTE} / supersedes 'A'`);
    expect(lines(0)[0], `superseded by 'A · PO'`);
  });

  test('F10: the Pipeline column is written only when it carries a reason; tooltip lists every reason', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('SUBJECT', ['R1', 'R2'])]);
    const pv = {runAt: '2026-10-15T09:12:00Z', user: 'flow', byKey: new Map([['R1', ['a reason', 'another']]])};
    expect(writePipelineColumns(df, pv, (i) => df.col('SUBJECT')!.getString(i)), true);
    expect(df.col(PIPELINE_COL)!.getString(0), 'a reason (+1)');
    expect(df.col(PIPELINE_COL)!.getString(1), '');
    expect(df.col(LZ_PIPELINE_INFO_COL) !== null, true);
    expect(pipelineTipAt(df, 0)?.extra?.length, 2, 'one bullet per reason');
    expect(pipelineTipAt(df, 1), null);
    writePipelineColumns(df, null, (i) => df.col('SUBJECT')!.getString(i));
    expect(df.col(PIPELINE_COL) === null && df.col(LZ_PIPELINE_INFO_COL) === null, true, 'removed with no reasons');
  });
});
