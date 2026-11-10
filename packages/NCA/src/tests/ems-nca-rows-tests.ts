/**
 * NCA EMS: row mapping — the pure row ↔ domain functions of the 0.2.0 layout (`ems-nca-rows.ts`),
 * without a server. Every round trip goes through `JSON.parse(JSON.stringify(...))`, the transport
 * the rows really take (NaN becomes null there), and compares with a NaN-aware serialisation — a
 * plain JSON comparison would read NaN and null as equal and hide exactly what these pin.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  resultRowsOf, summaryCacheOf, tierCols, groupKeysOf, lzRowsOf, lzOverridesOf, normLz, decisionRowsOf, ledgersOf,
  normLedgers, decisionSteps, storedDecisionOf, flagRowsOf, flagsOf, confirmedCols, confirmedOf,
  recipeSearchCols, originCols, originOf, pipelineRunCols, signatureCols, signatureOf, stageCols, stageOf,
} from '../adapter/ems-nca-rows';
import type {Bag} from '../adapter/ems-nca-rows';
import {computeRecipeFingerprint} from '../adapter/fingerprint';
import {isCurvePending, pendingCurveIds, withCurves} from '../adapter/nca-repository';
import type {LoadedStudy} from '../adapter/nca-repository';
import type {
  AnalysisRecipe, CompoundTierRow, LzOverride, MeanCv, PipelineOrigin, ReviewLedger, SummaryCache,
} from '../adapter/study-schema';

/** Canonical JSON (sorted keys, `undefined` dropped) with NaN / ±Infinity as tokens — a NaN that
 *  became null must NOT compare equal. */
const show = (x: unknown): string => {
  const walk = (v: unknown): unknown => {
    if (Array.isArray(v)) return v.map(walk);
    if (v !== null && typeof v === 'object') {
      return Object.fromEntries(Object.keys(v as Bag).sort().filter((k) => (v as Bag)[k] !== undefined)
        .map((k) => [k, walk((v as Bag)[k])]));
    }
    return typeof v === 'number' && !Number.isFinite(v) ? `#${v}` : v;
  };
  return JSON.stringify(walk(x));
};
/** The wire: what a row looks like after the JSON transport. */
const wire = <T>(x: T): T => JSON.parse(JSON.stringify(x)) as T;

const mc = (mean: number, cvPct: number, n: number, nExcluded?: number): MeanCv =>
  (nExcluded !== undefined ? {mean, cvPct, n, nExcluded} : {mean, cvPct, n});
const NA = mc(NaN, NaN, 0);

function tierRow(over: Partial<CompoundTierRow>): CompoundTierRow {
  return {
    compoundId: 'C1', isWholeAnalysis: false, groupValues: [{key: 'Compound', value: 'C1'}, {key: 'Sex', value: 'F'}],
    keyInventory: [{key: 'Dose', values: [{value: '1', count: 3}]}], passthrough: [], n: 3, route: 'PO',
    doseNorm: 'per_kg', cmax: mc(10.5, 12.3, 3), tmax: mc(1, NaN, 3), aucInf: mc(55.25, 8.1, 3, 1),
    thalf: mc(2.2, 5, 2, 1), cl: mc(0.123456789012345, 7, 3), vz: mc(1.5, 9, 3), mrt: mc(3.3, 4, 3),
    vss: NA, fPct: mc(45.6, 10, 3), pctExtrap: mc(12, 30, 3),
    curve: {series: [{group: 'C1 · F', time: [0, 1, 2], geoMean: [0, 5.5, 2.25], geoLower: [0, 4, 1],
      geoUpper: [0, 7, 3]}]},
    warnings: [{code: 'W1', message: 'first, with a comma'}, {code: 'W2', message: 'second "quoted"'}],
    ...over,
  };
}

const CACHE: SummaryCache = {
  analysisTier: {nCompounds: 2, nSubjects: 9, route: 'PO', dose: '1 mg/kg', status: 'warn', flagSummary: '0 failed, 2'},
  compoundTier: [
    tierRow({}),
    // %F computed, not estimable (no IV reference); every λz contributor excluded (FIELD-02); no curve.
    tierRow({compoundId: 'C2', groupValues: [{key: 'Compound', value: 'C2'}, {key: 'Sex', value: ''}],
      fPct: NA, aucInf: mc(NaN, NaN, 0, 3), curve: null,
      passthrough: [{key: 'SMILES', values: [{value: 'CCO', count: 3}], isStructure: true}]}),
    // IV: %F does not apply; several structures; an empty curve.
    tierRow({compoundId: 'C3', groupValues: [{key: 'Compound', value: 'C3'}, {key: 'Sex', value: 'M'}],
      route: 'IV-bolus', fPct: null, curve: {series: []},
      passthrough: [{key: 'SMILES', values: [{value: 'CCO', count: 1}, {value: 'CCN', count: 2}], isStructure: true},
        {key: 'Name', values: [{value: 'x', count: 3}], isStructure: false}]}),
  ],
};

const STAMPS = {engineVersion: '0.12.0', resultsBasis: 'per-kg-1'};

/** The analysis-row columns a summary is read back from. */
const analysisOf = (cache: SummaryCache | null): Bag => wire({id: 'a1', ...tierCols(cache)});

const RECIPE = {
  columns: {subject: 'SUBJ', time: 'TIME', concentration: 'CONC'},
  modes: {nonmem: false, confirmedRoute: null, sparse: {on: true, nominalTimeCol: 'NT', animalIdCol: 'A',
    samplingType: 'destructive'}},
  ncaRules: {aucMethod: 'linear-up-log-down', blq: {preFirstMeasurable: 'set-zero', embedded: 'missing',
    afterLast: 'set-half-lloq', consecutiveAfterLast: 'exclude'}},
  units: {}, groupByCols: ['Compound'], passthroughCols: [], acceptanceRSquared: 0.9,
} as unknown as AnalysisRecipe;

const LZ: LzOverride[] = [
  {profileKey: 'S2', timeKey: '8', ordinal: 0, concKey: '1.25', action: 'force-out', user: 'Ann',
    userId: 'u-ann', at: '2026-01-02T03:04:05Z'},
  {profileKey: 'S1', timeKey: '0.5', ordinal: 1, concKey: '', action: 'force-in', user: '',
    at: '2026-01-02T03:04:06.007Z'},
];

const rec = (action: 'accept' | 'acknowledge' | 'exclude', at: string, extra: Bag = {}) => ({action,
  reasons: [{code: 'R1', text: 'one, two'}, {code: 'OTHER', text: 'free "text" — λz'}], user: 'Bob', at,
  fingerprint: `fp-${at}`, ...extra});

const LEDGERS: ReviewLedger[] = [
  {unitKey: 'sparse:DOSE=1', current: null, history: [rec('exclude', '2026-01-01T00:00:00.000Z',
    {supersededBy: 'reset', supersededAt: '2026-01-01T01:00:00.000Z'})]},
  {unitKey: 'S1', current: rec('accept', '2026-01-03T00:00:00.000Z', {userId: 'u-bob'}), history: [
    rec('exclude', '2026-01-01T00:00:00.000Z', {supersededBy: 'auto-reset', supersededAt: '2026-01-02T00:00:00.000Z'}),
    rec('acknowledge', '2026-01-02T00:00:00.000Z',
      {supersededBy: 'replaced', supersededAt: '2026-01-03T00:00:00.000Z'}),
  ]},
];

category('NCA EMS: row mapping', () => {
  test('SummaryCache ↔ rows: every field round-trips (NaN, nExcluded, %F null / NE, curves)', async () => {
    const rows = wire(resultRowsOf('a1', 's1', CACHE, STAMPS));
    expect(rows.length, 3, 'one row per tier row');
    const back = summaryCacheOf(analysisOf(CACHE), [...rows].reverse()); // any order: sorted by ord
    expect(show(back), show(CACHE), 'the summary comes back exactly');
    expect(Number.isNaN(back!.compoundTier[1].aucInf.mean) && back!.compoundTier[1].aucInf.n === 0, true,
      'an all-excluded mean reads back as NaN with n = 0');
    expect(back!.compoundTier[1].fPct !== null && Number.isNaN(back!.compoundTier[1].fPct.mean), true,
      'the %F "not estimable" sentinel survives');
    expect(back!.compoundTier[2].fPct === null, true, '%F that does not apply stays null');
    expect(rows[1].f_pct_ne, true, 'f_pct_ne marks the NE sentinel');
    expect(rows.every((r) => r.engine_version === STAMPS.engineVersion && r.results_basis === STAMPS.resultsBasis),
      true, 'the rows carry the analysis stamps (decision 15)');
  });

  test('a light read + the lazily loaded curves give the same summary; curve null ≠ pending', async () => {
    const rows: Bag[] = wire(resultRowsOf('a1', 's1', CACHE, STAMPS)).map((r, i) => ({...r, id: `r${i}`}));
    const light = rows.map(({curve: _c, ...rest}) => rest);
    const curves = new Map(rows.filter((r) => r.has_curve).map((r) => [r.id as string, r.curve]));
    expect(show(summaryCacheOf(analysisOf(CACHE), light, curves)), show(CACHE), 'curves filled by id');
    const none = summaryCacheOf(analysisOf(CACHE), light)!;
    expect(none.compoundTier.every((r) => r.curve === null), true, 'not loaded yet → null (the caller tracks pending)');
    expect(show(rows.map((r) => r.has_curve)), show([true, false, true]), 'has_curve: a curve, none, an empty one');
  });

  test('withCurves: an id the load did not return reads as no curve, and nothing stays pending', async () => {
    const light: SummaryCache = {...CACHE, compoundTier: CACHE.compoundTier.map((r) => ({...r, curve: null}))};
    const study = {analyses: [{id: 'a1', summaryCache: light}],
      pendingCurves: new Map([['a1', ['r0', null, 'r2']]])} as unknown as LoadedStudy;
    expect(show(pendingCurveIds(study)), show(['r0', 'r2']));
    // r2 was deleted between the study read and the curve read: the load returns r0 only.
    const after = withCurves(study, new Map([['r0', CACHE.compoundTier[0].curve]]));
    const tier = after.analyses[0].summaryCache!.compoundTier;
    expect(show(tier[0].curve), show(CACHE.compoundTier[0].curve), 'a returned curve is filled in');
    expect(tier[2].curve === null && !isCurvePending(after, 'a1', 2), true, 'a missing one is no curve, not pending');
    expect(pendingCurveIds(after).length, 0, 'nothing left to load: the overview does not load again');
  });

  test('has_cache: false → no summary; true with zero rows → an empty tier, not null', async () => {
    expect(summaryCacheOf(analysisOf(null), []) === null, true, 'never confirmed → null (FIELD-04)');
    const empty: SummaryCache = {...CACHE, compoundTier: []};
    const back = summaryCacheOf(analysisOf(empty), []);
    expect(back !== null && back.compoundTier.length === 0, true, 'computed, zero rows → an empty summary');
  });

  test('structure: the SMILES of a single-structure row; null for several or none', async () => {
    const rows = resultRowsOf('a1', 's1', CACHE, STAMPS);
    expect(show(rows.map((r) => r.structure)), show([null, 'CCO', null]));
  });

  test('group keys: one list per analysis; rows that disagree refuse to store', async () => {
    expect(show(groupKeysOf(CACHE)), show(['Compound', 'Sex']));
    const whole: SummaryCache = {...CACHE, compoundTier: [tierRow({isWholeAnalysis: true, groupValues: []})]};
    expect(show(groupKeysOf(whole)), show([]), 'a whole-analysis summary has no keys');
    const bad: SummaryCache = {...CACHE, compoundTier: [tierRow({}),
      tierRow({groupValues: [{key: 'Dose', value: '1'}]})]};
    let threw = false;
    try {groupKeysOf(bad);} catch {threw = true;}
    expect(threw, true, 'disagreeing keys must not be stored under one label list');
  });

  test('λz overrides ↔ rows: the set round-trips (keyed, times as ISO)', async () => {
    const back = lzOverridesOf(wire(lzRowsOf('a1', LZ)));
    expect(show(back), show(normLz(LZ)));
    expect(back[0].profileKey, 'S1', 'ordered by key');
    expect(back[1].at, '2026-01-02T03:04:05.000Z', 'a time reads back in ISO with milliseconds');
    const rows = lzRowsOf('a1', LZ);
    expect(rows[0].set_by === 'u-ann' && rows[0].set_by_label === 'Ann', true, 'the user id and the name as recorded');
    expect(rows[1].set_by === null && back[0].userId === undefined, true, 'no platform user → null, read back absent');
  });

  test('review ledgers ↔ rows: current + history round-trip, history oldest first', async () => {
    const back = ledgersOf(wire(decisionRowsOf('a1', 's1', LEDGERS)));
    expect(show(back), show(normLedgers(LEDGERS)));
    const s1 = back.find((l) => l.unitKey === 'S1')!;
    expect(s1.current?.action, 'accept');
    expect(s1.current?.userId === 'u-bob' && s1.current?.user === 'Bob', true, 'the user id beside the name');
    expect(s1.history.every((h) => h.userId === undefined), true, 'a record without an id stays without one');
    expect(show(s1.history.map((h) => h.supersededBy)), show(['auto-reset', 'replaced']));
    expect(back.find((l) => l.unitKey === 'sparse:DOSE=1')!.current === null, true, 'a reset unit has no current');
  });

  test('decisionSteps: inserts new records, updates only superseded_*, never deletes', async () => {
    const stored = wire(decisionRowsOf('a1', 's1', LEDGERS)).map((r, i) => storedDecisionOf({...r, id: `d${i}`}));
    const decisionOps = (target: ReviewLedger[]): Bag[] => decisionSteps(stored, 'a1', 's1', target).steps
      .map((s) => s.op).filter((op): op is Bag => op !== null);
    expect(decisionOps([...LEDGERS]).length, 0, 'the same ledgers → no ops');
    // S1's current is reset, and a new decision is made on S2.
    const s1 = LEDGERS[1];
    const next: ReviewLedger[] = [LEDGERS[0],
      {...s1, current: null, history: [...s1.history, {...s1.current!, supersededBy: 'reset',
        supersededAt: '2026-01-04T00:00:00.000Z'}]},
      {unitKey: 'S2', current: rec('exclude', '2026-01-05T00:00:00.000Z'), history: []}];
    const ops = decisionOps(next);
    expect(show(ops.map((o) => o.op).sort()), show(['insert', 'update']));
    const upd = ops.find((o) => o.op === 'update')!;
    expect(show(Object.keys(upd.values).sort()), show(['superseded_at', 'superseded_by']),
      'an update touches the supersede columns only (the core is immutable)');
    // A ledger that vanished from the target: nothing is deleted.
    const dropped = decisionOps([LEDGERS[1]]);
    expect(dropped.every((o) => o.op !== 'delete'), true, 'append-only: no delete');
    expect(decisionSteps(stored, 'a1', 's1', [LEDGERS[1]]).unmatched.length, 1, 'the dropped record stays tracked');
  });

  test('pipeline flags: rows round-trip; no rows → [] for a pipeline entry, undefined otherwise', async () => {
    const flags = [{profileKey: 'S1', verdict: 'Review' as const, reasons: ['a', 'b, c'], info: [],
      pipelineOnly: ['b, c']}];
    expect(show(flagsOf(wire(flagRowsOf('a1', flags)), 'pipeline')), show(flags));
    expect(show(flagsOf([], 'pipeline')), show([]), 'a pipeline entry with nothing flagged');
    expect(flagsOf([], null) === undefined, true, 'an app-made entry has no flags');
  });

  test('confirmed: the snapshot is stored only when its fingerprint differs from the recipe', async () => {
    const fp = computeRecipeFingerprint(RECIPE);
    const same = confirmedCols({recipe: RECIPE, fingerprint: fp, at: '2026-01-01T00:00:00Z'}, fp);
    expect(same.confirmed_recipe === null, true, 'Confirm + Save → no second copy of the recipe');
    expect(show(confirmedOf(wire(same), RECIPE)), show({recipe: RECIPE, fingerprint: fp,
      at: '2026-01-01T00:00:00.000Z'}), 'rebuilt from the recipe on read');
    const edited = {...RECIPE, units: {time: 'min'}} as unknown as AnalysisRecipe;
    const kept = confirmedCols({recipe: RECIPE, fingerprint: fp, at: '2026-01-01T00:00:00Z'},
      computeRecipeFingerprint(edited));
    expect(kept.confirmed_recipe !== null, true, 'an unconfirmed Setup edit was saved → the snapshot is kept');
    expect(show(confirmedOf(wire(kept), edited)?.recipe), show(RECIPE), 'Discard restores the confirmed recipe');
    expect(confirmedOf(wire(confirmedCols(null, fp)), RECIPE) === null, true, 'never confirmed');
  });

  test('recipe search columns follow the recipe; a dense recipe has no sampling type', async () => {
    expect(show(recipeSearchCols(RECIPE)), show({is_sparse: true, sparse_sampling_type: 'destructive',
      auc_method: 'linear-up-log-down', blq_pre_first: 'set-zero', blq_embedded: 'missing',
      blq_after_last: 'set-half-lloq', blq_consecutive: 'exclude'}));
    const dense = {...RECIPE, modes: {...RECIPE.modes, sparse: {...RECIPE.modes.sparse, on: false}}};
    expect(recipeSearchCols(dense).sparse_sampling_type === null && recipeSearchCols(dense).is_sparse === false, true);
  });

  test('origin: standalone and pipeline round-trip; mapping and rules come back byte-identical', async () => {
    expect(show(originOf(wire(originCols({kind: 'standalone', label: 'L', containerId: 'c'})))),
      show({kind: 'standalone', label: 'L', containerId: 'c'}));
    const o: PipelineOrigin = {kind: 'pipeline', runId: 'run-1', runAt: '2026-01-01T00:00:00.000Z', user: 'svc',
      engineVersion: '0.12.0', mapping: {subject: 'ID', time: 'T', nested: {b: 1, a: [2, 'x']}},
      rules: '{"z": 1,  "a": "keep  spacing"}', partition: 'P · PO', recipeHash: 'h', recipeHashVersion: 1,
      supersedes: 'old'};
    const run = wire(pipelineRunCols(o, 's1'));
    const expanded = Object.fromEntries(Object.entries(run).map(([k, v]) => [`pipeline_run.${k}`, v]));
    const row = {...wire(originCols(o)), ...expanded};
    const back = originOf(row) as PipelineOrigin;
    expect(show(back), show(o));
    expect(back.rules, o.rules, 'rules verbatim');
    expect(run.mapping, JSON.stringify(o.mapping), 'mapping stored as JSON text');
    expect(originOf(wire(originCols(null))) === null, true, 'made in the app');
  });

  test('signature and stage columns round-trip', async () => {
    const sig = {kind: 'content-hash' as const, rowCount: 7, mtime: null, hash: 'abc',
      computedAt: '2026-01-01T00:00:00.000Z'};
    expect(show(signatureOf(wire(signatureCols(sig)))), show(sig));
    const stage = {stage: 'review' as const, pending: 3, notApplicable: 1};
    expect(show(stageOf(wire(stageCols(stage)))), show(stage));
    expect(stageOf(wire(stageCols(null))) === null, true);
  });
});
