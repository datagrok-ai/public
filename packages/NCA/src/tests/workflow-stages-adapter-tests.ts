/**
 * nca-workflow-stages — S1 adapter foundations (pure / headless).
 *
 * Covers F-S1.1 (acceptance move: stable flag codes), F-S1.2 (`reportableResults`),
 * F-S1.3 (review-unit classifier on three REAL fixtures — rule 16: the counts are measured,
 * not mirrored), F-S1.4 (fingerprints), F-S1.5 (stage + badge), F-S1.6 (reason registry),
 * F-S1.7 (schema v15), and the F-S4.8 adapter half (sparse composite acceptance).
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {DEFAULT_RULES, DEFAULT_ACCEPTANCE_RSQUARED} from '../adapter/defaults';
import {
  evaluateAcceptance, ACCEPT_ACCEPTED, ACCEPT_FLAGGED, ACCEPT_NA, ACCEPT_REASONS_COL, reasonsCellValue,
} from '../adapter/acceptance';
import type {AcceptanceThresholds, AcceptanceReading} from '../adapter/acceptance';
import {
  reportableResults, LZ_DEPENDENT_PARAMS, SPARSE_LZ_DEPENDENT_COLS, isSparseLzDependentCol, excludedCompositeRows,
} from '../adapter/reportable';
import {
  classifyReviewUnits, pendingUnits, profileMeasurability, denseProfilesViable, sparseUnitKey,
} from '../adapter/review-units';
import type {ProfileUnitInput, CompositeUnitInput} from '../adapter/review-units';
import {computeRecipeFingerprint, computeDecisionFingerprint} from '../adapter/fingerprint';
import type {DecisionFingerprintInput} from '../adapter/fingerprint';
import {
  deriveStage, formatStageBadge, formatStudyAnalysesSummary, formatNotReviewed,
} from '../adapter/workflow-stage';
import {
  judgementReasons, criterionReasons, reasonsValid, formatReasons, ACTIONS_BY_CATEGORY,
} from '../adapter/review-reasons';
import type {AnalysisRecipe, ReviewLedger} from '../adapter/study-schema';
import {aggregateNcaResults} from '../adapter/results-aggregator';
import {numAt} from '../adapter/adapter-utils';
import {computeSparseNca} from '../adapter/sparse-nca';
import {
  normalizeEventDrivenObservations, SUBJECT_COL, TAD_COL, CONC_COL, DOSE_COL, OCCASION_COL,
  type EventContract,
} from '../adapter/event-driven';
import {poolCompositeByCell} from '../adapter/composite-pool';
import {DERIVED_ROUTE_CODE_COL} from '../adapter/route-normalize';
import type {ColumnContract} from '../adapter/column-contract';
import {
  loadFixture, runPipelineOrThrow, withDoseMg, theophContract, theophScalars,
} from './test-utils';

const TIMEOUT = 60000;

const T: AcceptanceThresholds = {
  lambdaAcceptanceRSquared: DEFAULT_ACCEPTANCE_RSQUARED,
  lambdaMinPoints: DEFAULT_RULES.lambdaZ.minPoints,
  extrapWarnPct: DEFAULT_RULES.extrapWarnPct,
};

/** Acceptance of one results row, read the way the grid reads it (numAt, isNone-first). */
function rowAcceptance(results: DG.DataFrame, i: number): AcceptanceReading {
  const v = (n: string): number => {const c = results.col(n); return c === null ? NaN : numAt(c, i);};
  return evaluateAcceptance({
    spanRatio: v('SPAN_RATIO'), adjR2: v('LAMBDAZ_ADJ_R2'),
    nPoints: v('LAMBDAZ_NPOINTS'), pctExtrap: v('PCT_EXTRAP'),
  }, T);
}

/** Classifier inputs for every results row of a dense run. */
function profileInputs(
  results: DG.DataFrame, meas: Map<string, {measurable: number; blqTimes: string[]}>,
): ProfileUnitInput[] {
  const subj = results.col('SUBJECT')!;
  const status = results.col('STATUS')!;
  const out: ProfileUnitInput[] = [];
  for (let i = 0; i < results.rowCount; i++) {
    const key = subj.getString(i);
    const m = meas.get(key);
    out.push({
      key, status: status.getString(i), acceptance: rowAcceptance(results, i),
      measurableCount: m?.measurable ?? 0, blqTimes: m?.blqTimes ?? [],
    });
  }
  return out;
}

/** The real UC-04 sparse path on the mouse destructive fixture (BLQ flag column already added). */
function mouseSparseNca(df: DG.DataFrame): ReturnType<typeof computeSparseNca> {
  const conc = new Float64Array(df.rowCount);
  const cc = df.col('Conc_ng_mL')!;
  for (let i = 0; i < df.rowCount; i++) conc[i] = (cc.isNone(i) ? 0 : Number(cc.get(i))) / 1000;
  df.columns.add(DG.Column.fromFloat64Array('Conc', conc));
  return computeSparseNca(df, {
    nominalTimeCol: 'NominalTime_hr', animalIdCol: 'MouseID', blqCol: 'BLQflag',
    doseCol: 'Dose_mg_kg', groupCol: 'Dose_mg_kg', bootstrapIterations: 50, acceptance: T,
  });
}

function baseRecipe(): AnalysisRecipe {
  return {
    columns: {
      subject: 'Subject', time: 'Time', concentration: 'conc',
      doseSource: 'column', doseCol: 'Dose', doseScalar: 0,
      routeSource: 'value', routeScalar: 'PO', routeCol: null,
      infusionSource: 'none', infusionDurationCol: null, infusionDurationScalar: 0,
      bodyWeightSource: 'column', bodyWeightCol: 'Wt', bodyWeightScalar: 0,
      blqFlagCol: null, lloqSource: 'none', lloqCol: null, lloqScalar: 0,
      sexCol: null, sexValueMap: {},
    },
    modes: {
      nonmem: false, confirmedRoute: null,
      sparse: {on: false, nominalTimeCol: null, animalIdCol: null, samplingType: '(auto)'},
    },
    ncaRules: DEFAULT_RULES,
    units: {time: 'h', concentration: 'mg/L', dose: 'mg/kg', bodyWeight: 'kg'},
    groupByCols: [],
    passthroughCols: [],
    acceptanceRSquared: 0.9,
  };
}

function fpInput(over: Partial<DecisionFingerprintInput> = {}): DecisionFingerprintInput {
  return {
    points: [['0', '0.74'], ['0.25', '2.84'], ['24.37', '3.28']],
    dose: '320', route: 'PO', pointsUsedTimes: ['9.05', '12.12', '24.37'],
    lambdaZ: 0.04847, adjR2: 0.99861, aucInf: 216.63, pctExtrap: 14.23,
    verdict: ACCEPT_FLAGGED, flagCodes: ['SPAN_LT'],
    warnings: [{code: 'LAMBDAZ_LOW_SPAN', severity: 'info'}],
    ...over,
  };
}

category('Workflow stages: acceptance core (F-S1.1)', () => {
  test('flagCodes are stable ids in criterion order; reasons unchanged short form', async () => {
    const r = evaluateAcceptance({spanRatio: 1.5, adjR2: 0.8, nPoints: 3, pctExtrap: 25}, T);
    expect(r.verdict, ACCEPT_FLAGGED);
    expect(r.flagCodes.join(','), 'SPAN_LT,ADJR2_LT,NPTS_LE,EXTRAP_GT');
    expect(r.reasons.join(' · '), 'span ratio < 2 · adj. R² < 0.9 · N points (λz) ≤ 3 · %AUCextrap > 20%');
  });
  test('all four null → Not estimated, no codes; clean fit → Accepted', async () => {
    const na = evaluateAcceptance({spanRatio: NaN, adjR2: NaN, nPoints: NaN, pctExtrap: NaN}, T);
    expect(na.verdict, ACCEPT_NA);
    expect(na.flagCodes.length, 0);
    const ok = evaluateAcceptance({spanRatio: 4, adjR2: 0.99, nPoints: 5, pctExtrap: 5}, T);
    expect(ok.verdict, ACCEPT_ACCEPTED);
  });
});

category('Workflow stages: reportable results (F-S1.2)', () => {
  test('excluding theoph subject 1: λz-dependent N −1, Cmax/AUClast N unchanged, source untouched', async () => {
    const theoph = await loadFixture('theoph.csv');
    const {results} = runPipelineOrThrow(
      withDoseMg(theoph), theophContract(theoph), theophScalars(), DEFAULT_RULES);
    const subj = results.col('SUBJECT')!.getString(0);
    const before = results.col('AUC_INF')!.get(0);
    const rep = reportableResults(results, new Set([subj]));
    expect(rep.excludedCount, 1);
    expect(rep.df !== results, true, 'a clone, not the source frame');
    expect(results.col('AUC_INF')!.get(0), before, 'source frame not mutated');
    expect(rep.df.col('AUC_INF')!.isNone(0), true, 'AUC_INF blanked on the clone');
    expect(rep.df.col('CMAX')!.isNone(0), false, 'CMAX kept');
    expect(rep.df.col('LAMBDAZ_ADJ_R2')!.isNone(0), false, 'diagnostics kept (AC-S1.2.2)');
    // Whole-analysis aggregate: add a constant group column on both frames.
    const addAll = (df: DG.DataFrame): DG.DataFrame => {
      const c = df.clone();
      c.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'ALL', new Array(c.rowCount).fill('all')));
      return c;
    };
    const aggRaw = aggregateNcaResults(addAll(results), 'ALL');
    const aggEx = aggregateNcaResults(addAll(rep.df), 'ALL');
    const nOf = (agg: DG.DataFrame, p: string): number => Number(agg.col(`N_${p}`)?.get(0) ?? NaN);
    for (const p of ['HALF_LIFE', 'AUC_INF', 'CL']) {
      if (aggHas(aggRaw, p))
        expect(nOf(aggEx, p), nOf(aggRaw, p) - 1, `N_${p} drops by one`);
    }
    expect(nOf(aggEx, 'CMAX'), nOf(aggRaw, 'CMAX'), 'N_CMAX unchanged');
    expect(LZ_DEPENDENT_PARAMS.includes('CMAX'), false);
    expect(rep.excludedCountBy(['ALL']).size, 0, 'missing key column → empty map');
  }, {timeout: TIMEOUT});
});

function aggHas(agg: DG.DataFrame, p: string): boolean {
  return agg.col(`N_${p}`) !== null;
}

category('Workflow stages: review-unit classifier (F-S1.3)', () => {
  test('AC-S1.3.3 theoph (dense): 12 units, each classified by its shipped verdict', async () => {
    const theoph = await loadFixture('theoph.csv');
    const {results} = runPipelineOrThrow(
      withDoseMg(theoph), theophContract(theoph), theophScalars(), DEFAULT_RULES);
    const meas = profileMeasurability(theoph, 'Subject', 'Time', 'conc');
    const units = classifyReviewUnits({
      profiles: profileInputs(results, meas), viable: denseProfilesViable(theoph, 'Subject'),
      composites: [], ledgers: new Map(), fingerprints: new Map(),
    });
    expect(units.length, 12, 'one unit per subject');
    for (let i = 0; i < results.rowCount; i++) {
      const u = units[i];
      const v = rowAcceptance(results, i).verdict;
      const want = v === ACCEPT_ACCEPTED ? 'not-needed' : 'review';
      expect(u.need, want, `subject ${u.key} (${v})`);
    }
    const flagged = units.filter((u) => u.category === 'flagged').length;
    grok.shell.info(`theoph: ${flagged} flagged, ${pendingUnits(units).length} pending`);
  }, {timeout: TIMEOUT});

  test('AC-S1.3.1 mouse destructive: 48 animals N/A, 2 composites, 0 no-exposure; 10 BLQ-only disclosed', async () => {
    const df = await loadFixture('04_mouse_sparse_destructive.csv');
    const blq = df.col('BLQ')!;
    const blqBool = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'BLQflag',
      Array.from({length: df.rowCount}, (_, i) => blq.getString(i).trim().toLowerCase() === 'true'));
    df.columns.add(blqBool);
    expect(denseProfilesViable(df, 'MouseID'), false, 'one sample per animal → non-viable');
    const meas = profileMeasurability(df, 'MouseID', 'NominalTime_hr', 'Conc_ng_mL', 'BLQflag');
    const na: AcceptanceReading = evaluateAcceptance({spanRatio: NaN, adjR2: NaN, nPoints: NaN, pctExtrap: NaN}, T);
    const profiles: ProfileUnitInput[] = [...meas.entries()].map(([key, m]) => ({
      key, status: 'failed', acceptance: na, measurableCount: m.measurable, blqTimes: m.blqTimes,
    }));
    // Composites from the real sparse path.
    const sparse = mouseSparseNca(df);
    const composites: CompositeUnitInput[] = sparse.groups.map((g) => ({
      groupKey: g.groupKey, acceptance: g.acceptance,
      measurableTimepoints: g.composite.timepoints.filter((t) => t.mean > 0).length,
    }));
    const units = classifyReviewUnits({
      profiles, viable: false, composites, ledgers: new Map(), fingerprints: new Map(),
    });
    const animals = units.filter((u) => u.kind === 'profile');
    expect(animals.length, 48, '48 animals');
    expect(animals.every((u) => u.need === 'not-applicable'), true, 'all per-animal rows N/A (SF-3)');
    expect(units.filter((u) => u.kind === 'composite').length, 2, 'two dose-group composites');
    expect(units.filter((u) => u.need === 'review' && u.category === 'no-exposure').length, 0,
      'no counted no-exposure unit from late BLQ bleeds (C-1)');
    const blqOnly = animals.filter((u) => u.blqOnly).map((u) => u.key).sort();
    expect(blqOnly.join(','), 'M017,M020,M022,M023,M024,M043,M044,M045,M046,M048',
      'the 10 BLQ animals are disclosed (AC-S1.3.1a)');
    expect(animals.find((u) => u.key === 'M017')!.blqTimes.length > 0, true, 'BLQ time carried');
    for (const g of sparse.groups) {
      const want = evaluateAcceptance({spanRatio: g.spanRatio, adjR2: g.lambdaZAdjR2,
        nPoints: g.lambdaZNPoints, pctExtrap: g.pctExtrap}, T).verdict;
      expect(g.acceptance.verdict, want, `AC-S4.8.1 composite ${g.groupKey} verdict`);
      expect(units.find((u) => u.key === sparseUnitKey(g.groupKey)) !== undefined, true);
    }
    const frag = sparse.parameterTable.col('FRAGILITY')!;
    for (let i = 0; i < frag.length; i++)
      expect(frag.getString(i).includes('R²'), false, 'no R²<0.90 token in FRAGILITY (U2)');
    expect(sparse.parameterTable.col('ACCEPTANCE') !== null, true, 'ACCEPTANCE column');
    // acceptance-reason-pills AC-7: the composite carries the same `Violated` pills column.
    const pills = sparse.parameterTable.col(ACCEPT_REASONS_COL)!;
    expect(pills !== null, true, 'Violated column on the composite table');
    expect(pills.getTag('cell.renderer'), 'Tags', 'drawn by the platform Tags renderer');
    sparse.groups.forEach((g, i) =>
      expect(pills.getString(i), reasonsCellValue(g.acceptance), `composite ${g.groupKey} pills = its reasons`));
  }, {timeout: TIMEOUT});

  test('AC-S1.3.2 rat TK28d NONMEM composite: control cells are no-exposure, dosed cells by verdict', async () => {
    const df = await loadFixture('06_rat_tk28d_nonmem.csv');
    const ec: EventContract = {
      id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', rate: 'RATE', mdv: 'MDV',
      addl: 'ADDL', ii: 'II', route: 'ROUTE', bodyWeight: 'BW_kg', covariates: ['Sex', 'DOSE_MGKG'],
      day: 'DAY',
    };
    const norm = normalizeEventDrivenObservations(df, ec);
    expect(norm.errors.length, 0);
    expect(denseProfilesViable(norm.observations, SUBJECT_COL), false, 'destructive → pooled');
    const {pooled, memberHashByCell} = poolCompositeByCell(norm.observations,
      {doseGroupCol: 'DOSE_MGKG', sexCol: 'Sex'});
    const contract: ColumnContract = {observations: {
      subjectId: {name: SUBJECT_COL, required: true, rawType: 'string'},
      time: {name: TAD_COL, required: true, rawType: 'float64', units: 'h'},
      concentration: {name: CONC_COL, required: true, rawType: 'float64', units: 'ng/mL'},
      route: {name: DERIVED_ROUTE_CODE_COL, required: false, rawType: 'int32'},
    }};
    const {results} = runPipelineOrThrow(pooled, contract, {dose: 0, route: nca.ROUTE_OTHER,
      doseUnits: 'mg', concentrationUnits: 'ng/mL', timeUnits: 'h', bodyWeight: null,
      bodyWeightUnit: 'kg', infusionDuration: null, lloq: 0}, DEFAULT_RULES);
    const meas = profileMeasurability(pooled, SUBJECT_COL, TAD_COL, CONC_COL);
    const units = classifyReviewUnits({
      profiles: profileInputs(results, meas), viable: denseProfilesViable(pooled, SUBJECT_COL),
      composites: [], ledgers: new Map(), fingerprints: new Map(),
    });
    const control = units.filter((u) => u.key.startsWith('D0_'));
    expect(control.length > 0, true, 'control cells present');
    expect(control.every((u) => u.category === 'no-exposure' && u.need === 'review'), true,
      'every control cell is a counted no-exposure unit');
    expect(control.every((u) => u.allowedActions.join() === 'acknowledge'), true, 'Acknowledge only');
    const dosed = units.filter((u) => !u.key.startsWith('D0_'));
    expect(dosed.some((u) => u.category === 'no-exposure'), false, 'no dosed cell is no-exposure');
    expect(memberHashByCell.size, results.rowCount, 'one member hash per cell');
    grok.shell.info(`rat tk28d: ${control.length} control cells, ${pendingUnits(units).length} pending`);
  }, {timeout: TIMEOUT});

  test('a ledger decision is valid only while its fingerprint matches (D13)', async () => {
    const acc = evaluateAcceptance({spanRatio: 1.5, adjR2: 0.99, nPoints: 5, pctExtrap: 5}, T);
    const ledger: ReviewLedger = {unitKey: 'S1', history: [], current: {
      action: 'accept', reasons: [{code: 'ACC_VISUAL_OK', text: 'x'}], user: 'u', at: 't', fingerprint: 'abc'}};
    const run = (fp: string) => classifyReviewUnits({
      profiles: [{key: 'S1', status: 'ok', acceptance: acc, measurableCount: 6, blqTimes: []}],
      viable: true, composites: [], ledgers: new Map([['S1', ledger]]), fingerprints: new Map([['S1', fp]]),
    })[0];
    expect(run('abc').decision !== null, true, 'match → decided');
    expect(pendingUnits([run('abc')]).length, 0);
    const stale = run('zzz');
    expect(stale.decision, null, 'mismatch → not decided');
    expect(stale.staleDecision !== null, true, 'reported as stale');
  });
});

category('Workflow stages: fingerprints (F-S1.4)', () => {
  test('AC-S1.4.1 copy-edit invariance: warning order / reason text never move the hash', async () => {
    const a = computeDecisionFingerprint(fpInput());
    const b = computeDecisionFingerprint(fpInput({warnings: [{code: 'LAMBDAZ_LOW_SPAN', severity: 'info'}]}));
    expect(a, b);
    const two = [{code: 'A', severity: 'info'}, {code: 'B', severity: 'warning'}];
    expect(computeDecisionFingerprint(fpInput({warnings: two})),
      computeDecisionFingerprint(fpInput({warnings: [...two].reverse()})), 'warning order irrelevant');
  });
  test('nca-cl-bw-normalised-units: the per-kg reporting notes never move a decision fingerprint', async () => {
    // A saved λz decision must survive the per-kg basis change: these notes describe CL / Vz / Vss
    // only, and the fit the analyst judged is unchanged.
    const before = computeDecisionFingerprint(fpInput({warnings: []}));
    for (const code of ['PER_KG_NOT_FORMED', 'BW_VARIES_WITHIN_PROFILE', 'SPARSE_PER_KG_WEIGHT_BASED',
      'COMPOSITE_BW_PARTIAL', 'COMPOSITE_NOT_WEIGHT_BASED']) {
      expect(computeDecisionFingerprint(fpInput({warnings: [{code, severity: 'info'}]})), before,
        `${code} is not hashed`);
    }
    // A warning that is NOT a reporting note still moves it.
    expect(computeDecisionFingerprint(fpInput({warnings: [{code: 'CL_IMPLAUSIBLE_HEURISTIC_WARNING',
      severity: 'warning'}]})) === before, false, 'a unit-error heuristic still resets the decision');
  });
  test('fit changes move the hash; sub-6-digit noise does not', async () => {
    const a = computeDecisionFingerprint(fpInput());
    expect(a === computeDecisionFingerprint(fpInput({adjR2: 0.95})), false, 'adj-R² change');
    expect(a === computeDecisionFingerprint(fpInput({pointsUsedTimes: ['12.12', '24.37']})), false, 'points used');
    expect(a === computeDecisionFingerprint(fpInput({verdict: ACCEPT_ACCEPTED, flagCodes: []})), false, 'verdict');
    expect(a === computeDecisionFingerprint(fpInput({points: [['0', '0.75']]})), false, 'source data');
    expect(a, computeDecisionFingerprint(fpInput({lambdaZ: 0.0484700001})), '6 significant digits');
  });
  test('AC-S1.4.2 recipe fingerprint: no-op toggle equal; groupBy excluded; rule change moves', async () => {
    const r = baseRecipe();
    const fp = computeRecipeFingerprint(r);
    const toggled = {...r, modes: {...r.modes, nonmem: true}};
    const back = {...toggled, modes: {...toggled.modes, nonmem: false}};
    expect(computeRecipeFingerprint(back), fp, 'toggle-and-untoggle');
    expect(computeRecipeFingerprint({...r, groupByCols: ['Sex']}), fp, 'Report input not hashed');
    expect(computeRecipeFingerprint({...r, acceptanceRSquared: undefined}), fp, 'absent R_acc = default');
    expect(computeRecipeFingerprint({...r, acceptanceRSquared: 0.8}) === fp, false, 'R_acc is Setup');
    const rules = {...r.ncaRules, lambdaZ: {...r.ncaRules.lambdaZ, minPoints: 4}};
    expect(computeRecipeFingerprint({...r, ncaRules: rules}) === fp, false, 'rules are Setup');
  });
  test('AC-S1.4.3 composite membership: reassigning an animal changes the cell hash', async () => {
    const mk = (sexOfA2: string): DG.DataFrame => DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, SUBJECT_COL, ['A1', 'A2', 'A3']),
      DG.Column.fromFloat64Array(TAD_COL, Float64Array.from([1, 1, 2])),
      DG.Column.fromFloat64Array(CONC_COL, Float64Array.from([5, 6, 3])),
      DG.Column.fromFloat64Array(DOSE_COL, Float64Array.from([1, 1, 1])),
      DG.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, Int32Array.from([2, 2, 2])),
      DG.Column.fromFloat64Array(OCCASION_COL, Float64Array.from([1, 1, 1])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', sexOfA2, 'M']),
    ]);
    const h1 = poolCompositeByCell(mk('M'), {doseGroupCol: null, sexCol: 'Sex'}).memberHashByCell;
    const h2 = poolCompositeByCell(mk('F'), {doseGroupCol: null, sexCol: 'Sex'}).memberHashByCell;
    const cell = [...h1.keys()].find((k) => k.includes('_M_'))!;
    expect(h2.has(cell), true, 'the M cell exists in both');
    expect(h1.get(cell) === h2.get(cell), false, 'membership change flips the hash');
  });
});

category('Workflow stages: stage + badge (F-S1.5)', () => {
  const u = (need: 'review' | 'not-needed' | 'not-applicable', decided = false) => ({need, decided});
  test('AC-S1.5.1 stage table', async () => {
    expect(deriveStage({confirmedFingerprint: null, currentFingerprint: 'a', hasResults: false,
      units: [u('review')]}).stage, 'setup');
    expect(deriveStage({confirmedFingerprint: 'a', currentFingerprint: 'a', hasResults: true,
      units: [u('review'), u('not-applicable'), u('not-needed')]}).stage, 'review');
    const done = deriveStage({confirmedFingerprint: 'a', currentFingerprint: 'a', hasResults: true,
      units: [u('review', true), u('not-needed'), u('not-applicable')]});
    expect(done.stage, 'complete');
    expect(done.notApplicable, 1);
    const edited = deriveStage({confirmedFingerprint: 'a', currentFingerprint: 'b', hasResults: true,
      units: [u('review', true)]});
    expect(edited.stage, 'setup', 'un-confirmed Setup edit → setup (P6)');
    expect(edited.setupChanged, true);
  });
  test('AC-S1.5.2 a stale entry never renders a bare Complete', async () => {
    expect(formatStageBadge({stage: 'complete', pending: 0, notApplicable: 0}, false), 'Complete');
    expect(formatStageBadge({stage: 'complete', pending: 0, notApplicable: 0}, true),
      'Complete · stale (refreshes on open)');
    expect(formatStageBadge({stage: 'review', pending: 4, notApplicable: 0}, false), 'Review λz · 4 left');
    expect(formatStageBadge(null, false), 'Setup');
    expect(formatNotReviewed({stage: 'review', pending: 3, notApplicable: 0}), '3 not reviewed');
    expect(formatNotReviewed({stage: 'complete', pending: 0, notApplicable: 0}), '');
    expect(formatStudyAnalysesSummary([
      {stage: 'complete', pending: 0, notApplicable: 0},
      {stage: 'review', pending: 2, notApplicable: 0},
      {stage: 'complete', pending: 0, notApplicable: 0},
    ]), '3 analyses · 1 in review');
  });
});

category('Workflow stages: reasons (F-S1.6)', () => {
  test('AC-S1.6.0 reason sets keyed by (action, category); all-BLQ only under no-exposure', async () => {
    const noExp = judgementReasons('acknowledge', 'no-exposure');
    const notEst = judgementReasons('acknowledge', 'not-estimated');
    expect(noExp[0].code, 'NOEXP_TRUE_FINDING', 'true-finding reason listed first (R2/F8)');
    const codes = new Set(noExp.map((r) => r.code));
    expect(notEst.some((r) => codes.has(r.code)), false, 'disjoint');
    expect(notEst.some((r) => /LLOQ throughout/.test(r.text)), false);
    expect(judgementReasons('accept', 'not-estimated').length, 0, 'accept not allowed on no fit');
    expect(ACTIONS_BY_CATEGORY['flagged'].join(), 'accept,exclude');
  });
  test('AC-S1.6.1/2 display join + ≥1 non-blank reason required', async () => {
    const crit = criterionReasons(['ADJR2_LT', 'SPAN_LT'], T);
    expect(formatReasons(crit), 'adj. R² < 0.9; span ratio < 2');
    expect(reasonsValid([]), false);
    expect(reasonsValid([{code: 'OTHER', text: '  '}]), false);
    expect(reasonsValid([{code: 'OTHER', text: 'odd profile'}]), true);
  });
});

category('Workflow stages: composite exclusion (review SF-1, option b)', () => {
  test('Exclude λz on a composite marks exactly that row; the λz-dependent set matches the real table', async () => {
    const df = await loadFixture('04_mouse_sparse_destructive.csv');
    const blq = df.col('BLQ')!;
    df.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'BLQflag',
      Array.from({length: df.rowCount}, (_, i) => blq.getString(i).trim().toLowerCase() === 'true')));
    const sparse = mouseSparseNca(df);
    const keys = sparse.groups.map((g) => g.groupKey);
    expect(keys.length, 2, 'two dose-group composites');
    // Drift guard: every column the exclusion withholds exists on the table buildParameterTable makes.
    for (const c of SPARSE_LZ_DEPENDENT_COLS)
      expect(sparse.parameterTable.col(c) !== null, true, `sparse table carries ${c}`);
    // Observed composite data + fit diagnostics are never withheld (AC-S1.2.2 analogue).
    for (const c of ['AUC_LAST', 'CMAX', 'TMAX', 'LAMBDAZ_ADJ_R2', 'ACCEPTANCE'])
      expect(isSparseLzDependentCol(c), false, `${c} is not λz-dependent`);
    expect(excludedCompositeRows(keys, new Set()).join(), 'false,false', 'nothing excluded');
    expect(excludedCompositeRows(keys, new Set([sparseUnitKey(keys[1])])).join(), 'false,true',
      'the excluded composite row, and only it');
    // A dense subject key that happens to equal the group label never excludes the composite.
    expect(excludedCompositeRows(keys, new Set([keys[0]])).join(), 'false,false', 'keys are namespaced');
  }, {timeout: TIMEOUT});
});
