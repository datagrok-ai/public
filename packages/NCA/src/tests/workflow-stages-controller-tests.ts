/**
 * nca-workflow-stages S2 / S5 — controller-level persistence and exclusion behaviour.
 *
 * - AC-S2.1.2 (D22): a force-out on a NONMEM destructive-composite CELL survives (i) a full
 *   Compute that REBUILDS the pooled frame and (ii) a persist → fresh-controller "reopen".
 *   Measured on the real `06_rat_tk28d_nonmem.csv` through the production pooling path.
 * - AC-S2.2.1 / AC-S2.2.2 (D13): decisions persist with history; a fit change auto-resets;
 *   a λz point edit replaces; every event is audited.
 * - AC-S5.1.3: an Exclude λz on a Day-1 TK cell suppresses the λz-based Formula B / RAC_THEO for
 *   that cell exactly like a not-estimated Day-1 λz, and leaves Formula A unchanged.
 */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {NcaController} from '../controller/nca-controller';
import {DEFAULT_RULES, DEFAULT_ACCEPTANCE_RSQUARED} from '../adapter/defaults';
import {AuditLog} from '../adapter/audit-log';
import {reportableResults} from '../adapter/reportable';
import {appendGroupingColumns} from '../adapter/grouping-columns';
import {computeAccumulationRatios} from '../adapter/toxicokinetics';
import {
  normalizeEventDrivenObservations, SUBJECT_COL, TAD_COL, CONC_COL, type EventContract,
} from '../adapter/event-driven';
import {poolCompositeByCell, STUDY_DAY_COL} from '../adapter/composite-pool';
import {applySexMap} from '../adapter/sex-mapping';
import {DERIVED_ROUTE_CODE_COL} from '../adapter/route-normalize';
import type {ColumnContract} from '../adapter/column-contract';
import {decisionRowsOf, ledgersOf} from '../adapter/ems-nca-rows';
import type {ProfileScalars} from '../adapter/profile-inputs';
import {loadFixture, makeSyntheticObs, singleSubjectContract, singleSubjectScalars} from './test-utils';

const TIMEOUT = 60000;
const T = {
  lambdaAcceptanceRSquared: DEFAULT_ACCEPTANCE_RSQUARED,
  lambdaMinPoints: DEFAULT_RULES.lambdaZ.minPoints,
  extrapWarnPct: DEFAULT_RULES.extrapWarnPct,
};

const EC06: EventContract = {
  id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', rate: 'RATE', mdv: 'MDV',
  addl: 'ADDL', ii: 'II', route: 'ROUTE', bodyWeight: 'BW_kg', covariates: ['Sex', 'DOSE_MGKG'], day: 'DAY',
};
const POOLED_CONTRACT: ColumnContract = {observations: {
  subjectId: {name: SUBJECT_COL, required: true, rawType: 'string'},
  time: {name: TAD_COL, required: true, rawType: 'float64', units: 'h'},
  concentration: {name: CONC_COL, required: true, rawType: 'float64', units: 'ng/mL'},
  route: {name: DERIVED_ROUTE_CODE_COL, required: false, rawType: 'int32'},
}};
const POOLED_SCALARS: ProfileScalars = {
  dose: 0, route: nca.ROUTE_OTHER, doseUnits: 'mg', concentrationUnits: 'ng/mL', timeUnits: 'h',
  bodyWeight: null, bodyWeightUnit: 'kg', infusionDuration: null, lloq: 0,
};

async function pooledTk(): Promise<DG.DataFrame> {
  const df = await loadFixture('06_rat_tk28d_nonmem.csv');
  const norm = normalizeEventDrivenObservations(df, EC06);
  return poolCompositeByCell(norm.observations, {doseGroupCol: 'DOSE_MGKG', sexCol: 'Sex'}).pooled;
}

function rowsOf(df: DG.DataFrame, subject: string): number[] {
  const c = df.col(SUBJECT_COL)!;
  const out: number[] = [];
  for (let i = 0; i < df.rowCount; i++) if (c.getString(i) === subject) out.push(i);
  return out;
}

category('Workflow stages: controller persistence (S2)', () => {
  test('AC-S2.1.2 D22: a NONMEM composite-cell override survives a rebuilt frame and a reopen', async () => {
    const cell = 'D10_F_Day28';
    const pooled1 = await pooledTk();
    const ctrl = new NcaController(pooled1, POOLED_CONTRACT, POOLED_SCALARS, DEFAULT_RULES);
    ctrl.computeAll();
    const rows1 = rowsOf(pooled1, cell);
    expect(rows1.length, 6, 'the Day-28 cell has 6 composite points');
    // Force the FIRST point out (a real analyst action on a composite cell).
    ctrl.overrideLambdaZPoints(cell, [{rowIdx: rows1[0], action: 'force-out'}]);
    expect(ctrl.overriddenSubjects().has(cell), true, 'overridden');
    const saved = JSON.parse(JSON.stringify(ctrl.getLzOverrides()));

    // (i) a full Compute REBUILDS the pooled frame (a new instance) — the regression C-3.
    const pooled2 = await pooledTk();
    ctrl.setInputs(pooled2, POOLED_CONTRACT, POOLED_SCALARS, DEFAULT_RULES);
    ctrl.computeAll();
    expect(ctrl.consumeReconciled(), 0, 'nothing dropped');
    expect(ctrl.currentLzAction(rowsOf(pooled2, cell)[0]), 'force-out', 'still forced out after Compute');
    expect(ctrl.overriddenSubjects().has(cell), true, '✎ source survives a Compute');

    // (ii) persist → a fresh controller on a fresh frame (Save → reopen).
    const pooled3 = await pooledTk();
    const reopened = new NcaController(pooled3, POOLED_CONTRACT, POOLED_SCALARS, DEFAULT_RULES);
    reopened.setLzOverrides(saved);
    reopened.computeAll();
    expect(reopened.currentLzAction(rowsOf(pooled3, cell)[0]), 'force-out', 'still forced out after reopen');
    expect(reopened.fitRowsBySubject().get(cell)!.has(rowsOf(pooled3, cell)[0]), false,
      'the facet marker source agrees: the point is not in the fit');
  }, {timeout: TIMEOUT});

  test('a λz override updates the outcome\'s per-profile result the session summary is built from', async () => {
    // The view keeps `outcome.perProfile` and builds the saved summary (warnings, status) from
    // it; a stale entry saved the pre-override warnings until the next open recomputed them.
    const df = makeSyntheticObs({
      subjects: ['S1', 'S1', 'S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S2', 'S2'],
      times: Float64Array.from([0, 1, 2, 4, 8, 0, 1, 2, 4, 8]),
      concs: Float64Array.from([10, 5, 2.5, 0.625, 0.039, 20, 10, 5, 1.25, 0.078]),
    });
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES,
      {log: new AuditLog()});
    const session = ctrl.computeAll();
    if (session.status !== 'completed') throw new Error('compute failed');
    const i = session.perProfile.findIndex((p) => p.index.subjectKey === 'S1');
    const before = session.perProfile[i].result;
    // Force out two terminal points: the fit falls below the minimum and λz is not estimated.
    ctrl.overrideLambdaZPoints('S1', [{rowIdx: 3, action: 'force-out'}, {rowIdx: 4, action: 'force-out'}]);
    const live = session.perProfile[i].result;
    expect(live !== before, true, 'the held per-profile entry is the recomputed one');
    const fresh = ctrl.computeAll(); // what a reopen computes (overrides are sticky)
    if (fresh.status !== 'completed') throw new Error('recompute failed');
    const codes = (r: typeof live): string => r.provenance.warnings.map((w) => w.code).sort().join(',');
    expect(codes(live), codes(fresh.perProfile[i].result), 'the session warnings = a full recompute\'s');
    expect(live.status, fresh.perProfile[i].result.status);
    expect(codes(live) !== codes(before), true, 'the override changed the warnings (the case is not vacuous)');
  }, {timeout: TIMEOUT});

  test('AC-S2.2.1 / AC-S2.2.2 decisions: record, keep, auto-reset, replace — history + audit', async () => {
    const df = makeSyntheticObs({
      subjects: ['S1', 'S1', 'S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S2', 'S2'],
      times: Float64Array.from([0, 1, 2, 4, 8, 0, 1, 2, 4, 8]),
      concs: Float64Array.from([10, 5, 2.5, 0.625, 0.039, 20, 10, 5, 1.25, 0.078]),
    });
    const log = new AuditLog();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES, {log});
    ctrl.computeAll();
    const snap = ctrl.reviewSnapshot(T, {viable: true});
    // Decide on S1 with the fingerprint of what is on screen (whatever its verdict).
    const fp = ctrl.reviewSnapshot(T, {viable: true}).units.find((u) => u.key === 'S1');
    expect(fp !== undefined, true, 'S1 is a unit');
    const allFps = new Map<string, string>([['S1', 'fp-S1']]);
    ctrl.decide('S1', 'exclude', [{code: 'EXC_ERRATIC', text: 'Non-monotonic'}], 'fp-S1');
    expect(ctrl.excludedKeys().has('S1'), true, 'exclude decision feeds the exclusion set');
    expect(ctrl.reconcileDecisions(allFps), 0, 'same fingerprint → kept');
    expect(ctrl.reconcileDecisions(new Map([['S1', 'changed']])), 1, 'fit changed → auto-reset');
    expect(ctrl.decisionFor('S1'), null, 'back to To review');
    expect(ctrl.excludedKeys().has('S1'), false, 'a reset decision excludes nothing');
    expect(ctrl.ledgerFor('S1')!.history.length, 1, 'superseded decision kept in history');
    expect(ctrl.ledgerFor('S1')!.history[0].supersededBy, 'auto-reset');

    // A point edit REPLACES a current decision (D13).
    ctrl.decide('S1', 'accept', [{code: 'ACC_VISUAL_OK', text: 'ok'}], 'fp2');
    ctrl.overrideLambdaZPoints('S1', [{rowIdx: 0, action: 'force-out'}]);
    expect(ctrl.decisionFor('S1'), null, 'edit after Accept as is → replaced');
    expect(ctrl.ledgerFor('S1')!.history.map((h) => h.supersededBy).join(), 'auto-reset,replaced');
    // IR-12: both records carry the analyst's user id beside the display name.
    const me = grok.shell.user.id;
    expect(ctrl.ledgerFor('S1')!.history.every((h) => h.userId === me), true, 'each decision carries the user id');
    expect(ctrl.getLzOverrides().every((o) => o.userId === me && o.user !== ''), true,
      'the override carries the user id and the name');

    // Manual reset.
    ctrl.decide('S2', 'accept', [{code: 'ACC_VISUAL_OK', text: 'ok'}], 'fpS2');
    ctrl.resetDecision('S2');
    expect(ctrl.ledgerFor('S2')!.history[0].supersededBy, 'reset');

    const actions = log.entries().map((e) => e.action);
    for (const a of ['review-exclude', 'review-auto-reset', 'review-accept', 'review-reset'])
      expect(actions.includes(a), true, `audited: ${a}`);

    // AC-S2.2.1: the ledgers round-trip through storage (append-only `review_decision` rows, one per
    // record — ems-modification) with history intact.
    const stored = ledgersOf(decisionRowsOf('a', 's', ctrl.getReviewLedgers()));
    const back = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    back.setReviewLedgers(stored);
    expect(back.ledgerFor('S1')!.history.length, 2, 'superseded decisions readable after reopen');
    expect(snap.units.length, 2, 'two review units');
  });
});

category('Workflow stages: Confirm dry run (F-S3.5 / AD-7)', () => {
  // The pane's dry run = a THROW-AWAY controller + AuditLog seeded with the live overrides and
  // ledgers, run on the edited Setup, reconciled against the same snapshot the real Confirm uses.
  // These pin the two contracts at the controller level (the pane code is the thin wrapper).
  function setup(): {df: DG.DataFrame; ctrl: NcaController; fpS1: string; log: AuditLog} {
    const df = makeSyntheticObs({
      subjects: ['S1', 'S1', 'S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S2', 'S2'],
      times: Float64Array.from([0, 1, 2, 4, 8, 0, 1, 2, 4, 8]),
      concs: Float64Array.from([10, 5, 2.5, 0.625, 0.039, 20, 10, 5, 1.25, 0.078]),
    });
    const log = new AuditLog();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES, {log});
    ctrl.computeAll();
    const snap = ctrl.reviewSnapshot(T, {viable: true});
    // Record decisions on every unit on screen with its REAL current fingerprint.
    const all = ctrl.reviewSnapshot(T, {viable: true});
    for (const u of all.units) {
      const fp = (u.need === 'review' ? snap.fingerprints.get(u.key) : undefined) ?? 'fp';
      ctrl.decide(u.key, 'accept', [{code: 'ACC_VISUAL_OK', text: 'ok'}], fp);
    }
    return {df, ctrl, fpS1: snap.fingerprints.get('S1') ?? '', log};
  }
  const dryRun = (live: NcaController, df: DG.DataFrame, rules: typeof DEFAULT_RULES,
    log: AuditLog): number => {
    const tmp = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), rules, {log});
    tmp.setLzOverrides(live.getLzOverrides());
    tmp.setReviewLedgers(live.getReviewLedgers());
    tmp.computeAll();
    return tmp.reconcileDecisions(tmp.reviewSnapshot(T, {viable: true}).fingerprints);
  };

  test('AC-S3.5.1 the dry run predicts exactly the resets the real Confirm makes', async () => {
    const {df, ctrl} = setup();
    // An edit that changes the fits: AUC method → linear (λz/AUC∞ move → fingerprints move).
    const edited = {...DEFAULT_RULES, lambdaZ: {...DEFAULT_RULES.lambdaZ, minPoints: 4}};
    const predicted = dryRun(ctrl, df, edited, new AuditLog());
    ctrl.setInputs(df, singleSubjectContract(df), singleSubjectScalars(), edited);
    ctrl.computeAll();
    const actual = ctrl.reconcileDecisions(ctrl.reviewSnapshot(T, {viable: true}).fingerprints);
    expect(predicted, actual, `predicted ${predicted} = actual ${actual}`);
  });

  test('AC-S3.5.2 a dry run followed by Cancel leaves the session audit log and the live state alone', async () => {
    const {df, ctrl, log} = setup();
    const before = log.entries().length;
    const ledgersBefore = JSON.stringify(ctrl.getReviewLedgers());
    const throwAway = new AuditLog();
    dryRun(ctrl, df, {...DEFAULT_RULES, aucMethod: 'linear'}, throwAway);
    expect(throwAway.entries().length > 0, true, 'the dry run DID compute (into the throw-away log)');
    expect(log.entries().length, before, 'no phantom compute entries in the live audit log');
    expect(JSON.stringify(ctrl.getReviewLedgers()), ledgersBefore, 'live decisions untouched');
  });
});

category('Workflow stages: exclusion consumers (S5)', () => {
  // AC-S5.1.3. Measured 2026-09-26 (rule 16): on `06_rat_tk28d_nonmem.csv` EVERY Day-1 composite
  // cell has λz = NA — the PKNCA per-cell reference agrees (toxicokinetics-tests AC-6: all 8
  // cells fail the gate) — so Formula B is unreachable on the bundled TK fixtures and cannot be
  // "suppressed" there. The λz-path contract is therefore pinned on a synthetic pass-cell (the
  // toxicokinetics-tests AC-7 pattern, here in results-frame column names so the production
  // `reportableResults` applies unchanged); the real fixture pins that exclusion never moves
  // Formula A.
  const TK_OPTS = {tau: 24, doseCol: 'DOSE_MGKG', studyDayCol: STUDY_DAY_COL, stratifyCols: ['Sex']};

  test('AC-S5.1.3 (synthetic pass-cell) excluded Day-1 λz → Formula B / RAC_THEO suppressed like NA', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('SUBJECT', ['D10_F_Day1', 'D10_F_Day28']),
      DG.Column.fromStrings('Sex', ['F', 'F']),
      DG.Column.fromFloat64Array('DOSE_MGKG', Float64Array.from([10, 10])),
      DG.Column.fromInt32Array(STUDY_DAY_COL, Int32Array.from([1, 28])),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from([2, 8])),
      DG.Column.fromFloat64Array('AUC_LAST', Float64Array.from([40, 200])),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([50, DG.FLOAT_NULL])),
      DG.Column.fromFloat64Array('LAMBDA_Z', Float64Array.from([0.05, DG.FLOAT_NULL])),
      DG.Column.fromFloat64Array('LAMBDAZ_ADJ_R2', Float64Array.from([0.95, DG.FLOAT_NULL])),
      DG.Column.fromFloat64Array('PCT_EXTRAP', Float64Array.from([10, DG.FLOAT_NULL])),
      DG.Column.fromFloat64Array('TLAST', Float64Array.from([24, 24])),
    ]);
    const optsB = {...TK_OPTS, formula: 'B' as const};
    const base = computeAccumulationRatios(df, optsB).table;
    expect(Number.isFinite(base.col('RAC_THEO')!.get(0) as number), true, 'gate passes before exclusion');
    expect(String(base.col('RAC_FORMULA')!.get(0)).includes('fallback'), false, 'Formula B used');

    const ex = computeAccumulationRatios(reportableResults(df, new Set(['D10_F_Day1'])).df, optsB).table;
    expect(ex.col('RAC_THEO')!.isNone(0), true, 'RAC_THEO suppressed once the Day-1 λz is excluded');
    expect(String(ex.col('RAC_FORMULA')!.get(0)).includes('fallback'), true, 'Formula B falls back to A');
    // The same state as a genuinely not-estimated Day-1 λz (AD-6: excluded ≡ not estimated).
    const naDf = df.clone();
    for (const c of ['LAMBDA_Z', 'LAMBDAZ_ADJ_R2', 'PCT_EXTRAP', 'AUC_INF']) naDf.col(c)!.set(0, DG.FLOAT_NULL);
    const na = computeAccumulationRatios(naDf, optsB).table;
    expect(String(ex.col('RAC_FORMULA')!.get(0)), String(na.col('RAC_FORMULA')!.get(0)), 'excluded ≡ not estimated');
    expect(ex.col('RAC_AUC')!.get(0), na.col('RAC_AUC')!.get(0), 'same Rac value as not estimated');
  });

  test('AC-S5.1.3 (real fixture) excluding a Day-1 cell leaves Formula A unchanged', async () => {
    const pooled = await pooledTk();
    const ctrl = new NcaController(pooled, POOLED_CONTRACT, POOLED_SCALARS, DEFAULT_RULES);
    const outcome = ctrl.computeAll();
    if (outcome.status !== 'completed') throw new Error('compute failed');
    const results = outcome.results;
    appendGroupingColumns(results, pooled, outcome.perProfile.map((p) => p.index), POOLED_CONTRACT);
    const a = computeAccumulationRatios(results, TK_OPTS).table;
    const b = computeAccumulationRatios(reportableResults(results, new Set(['D10_F_Day1'])).df, TK_OPTS).table;
    expect(a.rowCount, b.rowCount, 'same cells');
    for (let i = 0; i < a.rowCount; i++) {
      expect(b.col('RAC_AUC')!.get(i), a.col('RAC_AUC')!.get(i), `row ${i} RAC_AUC unchanged`);
      expect(b.col('RAC_CMAX')!.get(i), a.col('RAC_CMAX')!.get(i), `row ${i} RAC_CMAX unchanged`);
    }
  }, {timeout: TIMEOUT});
});

category('Workflow stages: override identity on derived frames (TEST-GAP-WS-01, AC-S2.1.2)', () => {
  /** Serial multi-occasion NONMEM (DAY 1 / 28 per animal) → AD-8 occasion-scoped pseudo-subjects. */
  const EC08: EventContract = {
    id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', mdv: 'MDV', day: 'DAY',
  };
  /** Serial multi-route NONMEM carrying a Sex covariate — the sex-map path's input. */
  const EC03: EventContract = {
    id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', rate: 'RATE', mdv: 'MDV',
    route: 'ROUTE', covariates: ['Sex'],
  };
  async function frame(file: string, ec: EventContract): Promise<DG.DataFrame> {
    return normalizeEventDrivenObservations(await loadFixture(file), ec).observations;
  }
  /** The first (pseudo-)subject with ≥ 4 measurable points, and its LAST measurable row. */
  function pickSubject(df: DG.DataFrame): {subject: string; row: number} {
    const conc = df.col(CONC_COL)!;
    const subj = df.col(SUBJECT_COL)!;
    const seen = new Map<string, number[]>();
    for (let i = 0; i < df.rowCount; i++) {
      if (conc.isNone(i) || !(Number(conc.get(i)) > 0)) continue;
      const k = subj.getString(i);
      const arr = seen.get(k) ?? [];
      arr.push(i);
      seen.set(k, arr);
    }
    for (const [k, rows] of seen) if (rows.length >= 4) return {subject: k, row: rows[rows.length - 1]};
    throw new Error('no subject with ≥ 4 measurable points');
  }
  const rowOfSame = (df: DG.DataFrame, subject: string, like: DG.DataFrame, likeRow: number): number => {
    const t = Number(like.col(TAD_COL)!.get(likeRow));
    const subj = df.col(SUBJECT_COL)!;
    const time = df.col(TAD_COL)!;
    for (let i = 0; i < df.rowCount; i++)
      if (subj.getString(i) === subject && Number(time.get(i)) === t) return i;
    return -1;
  };

  test('an occasion pseudo-subject override survives a rebuilt frame and a reopen', async () => {
    const raw = await loadFixture('08_perday_reset_nonmem.csv');
    const rawIds = new Set(Array.from({length: raw.rowCount}, (_, i) => raw.col('ID')!.getString(i)));
    const f1 = await frame('08_perday_reset_nonmem.csv', EC08);
    const {subject, row} = pickSubject(f1);
    expect(rawIds.has(subject), false, `an occasion-scoped pseudo-key, not a raw ID (${subject})`);
    const ctrl = new NcaController(f1, POOLED_CONTRACT, POOLED_SCALARS, DEFAULT_RULES);
    ctrl.computeAll();
    ctrl.overrideLambdaZPoints(subject, [{rowIdx: row, action: 'force-out'}]);
    const saved = JSON.parse(JSON.stringify(ctrl.getLzOverrides()));
    const f2 = await frame('08_perday_reset_nonmem.csv', EC08); // a fresh instance, as a re-Compute
    ctrl.setInputs(f2, POOLED_CONTRACT, POOLED_SCALARS, DEFAULT_RULES);
    ctrl.computeAll();
    expect(ctrl.consumeReconciled(), 0, 'nothing dropped');
    expect(ctrl.currentLzAction(rowOfSame(f2, subject, f1, row)), 'force-out', 'survives the rebuilt frame');
    const f3 = await frame('08_perday_reset_nonmem.csv', EC08);
    const reopened = new NcaController(f3, POOLED_CONTRACT, POOLED_SCALARS, DEFAULT_RULES);
    reopened.setLzOverrides(saved);
    reopened.computeAll();
    expect(reopened.currentLzAction(rowOfSame(f3, subject, f1, row)), 'force-out', 'survives a reopen');
  }, {timeout: TIMEOUT});

  test('an override survives a sex-remapped derived frame (the Sex column replaced, keys intact)', async () => {
    const f1 = await frame('03_rat_multiroute_nonmem.csv', EC03);
    expect(f1.col('Sex') !== null, true, 'the Sex covariate is carried');
    const {subject, row} = pickSubject(f1);
    const ctrl = new NcaController(f1, POOLED_CONTRACT, POOLED_SCALARS, DEFAULT_RULES);
    ctrl.computeAll();
    ctrl.overrideLambdaZPoints(subject, [{rowIdx: row, action: 'force-out'}]);
    const f2 = (await frame('03_rat_multiroute_nonmem.csv', EC03)).clone();
    applySexMap(f2, 'Sex', new Map([['M', 'F'], ['F', 'M']])); // the sex-map path's derived frame
    ctrl.setInputs(f2, POOLED_CONTRACT, POOLED_SCALARS, DEFAULT_RULES);
    ctrl.computeAll();
    expect(ctrl.consumeReconciled(), 0, 'nothing dropped');
    expect(ctrl.currentLzAction(rowOfSame(f2, subject, f1, row)), 'force-out', 'survives the remap');
  }, {timeout: TIMEOUT});
});
