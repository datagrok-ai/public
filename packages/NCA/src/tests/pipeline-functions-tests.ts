/**
 * Pipeline functions (ncaDetectMapping / ncaCompute / ncaSummarize).
 *
 * Synthetic profiles are noise-free one-compartment curves, so their expectations
 * are about routing (verdicts, keys, blocking), not about parameter precision.
 * Precision is asserted against the PKNCA references in the
 * `Pipeline functions: reference` category, through `ncaCompute` itself.
 */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';

import {computePipeline, VERDICT_ACCEPTED, VERDICT_BLOCKED, VERDICT_REVIEW, COMPOSITE_SUBJECT}
  from '../pipeline/pipeline-compute';
import {detectPipelineMapping} from '../pipeline/pipeline-mapping';
import {summarizePipeline} from '../pipeline/pipeline-summary';
import {saveToStudies, recipeFor} from '../pipeline/pipeline-save';
import {parseMapping} from '../pipeline/pipeline-mapping';
import {nca} from '@datagrok-libraries/sci-comp';
import {
  newId, recipeOf, SourceReplacedError, getNcaRepository, setNcaRepository, computeSnapshotHash, RESULTS_BASIS,
  DEFAULT_RULES, pendingCurveIds, withCurves,
} from '../adapter';
import type {LzOverride, NcaStudy, AnalysisRecipe, ColumnContract} from '../adapter';
import {runNcaPipeline} from '../controller/nca-controller';
import {ncaDb} from '../generated/db';
import {openAnalysis} from '../views/nca-analysis-open';
import {updateAnalysisInPlace} from '../views/nca-add-to-study';
import {currentRecipeFingerprint, liveWorkflowState, RESULTS_GRID_TAG, isAnalysisDirty}
  from '../views/nca-mapping-pane';
import type {StudyBinding} from '../views/nca-mapping-pane';
import {MemoryNcaRepository} from './memory-nca-repository';
import {signatureCols} from '../adapter/ems-nca-rows';
import {loadFixture, loadReference, indexFixtureProfiles, type FixtureProfile} from './test-utils';

const TIMEOUT = 60000;
const FIXTURE_PATH_PREFIX = 'System:AppData/Nca';

const T_IV = [0.083, 0.25, 0.5, 1, 2, 4, 6, 8, 24];
const T_PO = [0, 0.25, 0.5, 1, 2, 4, 6, 8, 24];
const T_SPARSE = [0.25, 0.5, 1, 2, 4, 6, 8, 24];

interface Arm {
  study: string; compound: string; route: 'IV' | 'PO' | 'IV-infusion'; dose: number; subjects: string[];
  cl?: number; v?: number; ka?: number; f?: number; period?: string;
  /** Sex of every subject in the arm (default `M`). */
  sex?: string;
  /** Per-subject multiplier on the concentrations (exposure outlier). */
  scale?: Record<string, number>;
  /** Per-subject {time: conc} overrides (spikes). */
  override?: Record<string, Record<number, number>>;
  comment?: Record<string, string>;
  times?: number[];
}

/** One-compartment concentration in ng/mL; CL in L/h/kg, V in L/kg, dose in mg/kg. */
function conc(a: Arm, t: number): number {
  const cl = a.cl ?? 0.72; const v = a.v ?? 1.4; const k = cl / v;
  if (a.route !== 'PO') return 1000 * a.dose / v * Math.exp(-k * t);
  const ka = a.ka ?? 6; const f = a.f ?? 0.55;
  return 1000 * f * a.dose * ka / (v * (ka - k)) * (Math.exp(-k * t) - Math.exp(-ka * t));
}

function table(arms: Arm[], opts: {withPeriod?: boolean; compoundCol?: string} = {}): DG.DataFrame {
  const rows: Array<Record<string, unknown>> = [];
  for (const a of arms) {
    for (const s of a.subjects) {
      for (const t of a.times ?? (a.route === 'PO' ? T_PO : T_IV)) {
        let c = conc(a, t) * (a.scale?.[s] ?? 1);
        if (a.override?.[s]?.[t] !== undefined) c = a.override[s][t];
        const blq = c < 1;
        rows.push({Study: a.study, [opts.compoundCol ?? 'Compound']: a.compound, Subject: s, Sex: a.sex ?? 'M',
          ...(opts.withPeriod ? {Period: a.period ?? ''} : {}),
          Route: a.route, Dose_mg_kg: a.dose, BW_kg: 0.25, Time_hr: t, Conc_ng_mL: blq ? 0 : c, BLQ: blq,
          LLOQ: 1, Comment: a.comment?.[s] ?? ''});
      }
    }
  }
  const names = Object.keys(rows[0]);
  return DG.DataFrame.fromColumns(names.map((n) => {
    const vals = rows.map((r) => r[n]);
    if (typeof vals[0] === 'boolean') return DG.Column.fromList(DG.COLUMN_TYPE.BOOL, n, vals);
    if (typeof vals[0] === 'number') return DG.Column.fromFloat64Array(n, Float64Array.from(vals as number[]));
    return DG.Column.fromStrings(n, vals as string[]);
  }));
}

const ivPo = (study: string, compound: string, extra: Partial<Arm> = {}, cl?: number): Arm[] => [
  {study, compound, route: 'IV', dose: 1, subjects: ['R1001', 'R1002', 'R1003'], cl, ...extra},
  {study, compound, route: 'PO', dose: 5, subjects: ['R2001', 'R2002', 'R2003'], cl, ...extra},
];

function rowsWhere(df: DG.DataFrame, col: string, value: string): number[] {
  const c = df.col(col)!;
  return Array.from({length: df.rowCount}, (_, i) => i).filter((i) => c.get(i) === value);
}
const s = (df: DG.DataFrame, col: string, i: number): string => String(df.col(col)!.get(i) ?? '');
/** A numeric result cell, NaN for a platform null (`get` returns the raw FLOAT_NULL sentinel). */
const num = (df: DG.DataFrame, col: string, i: number): number =>
  df.col(col)!.isNone(i) ? NaN : df.col(col)!.get(i) as number;
const reasonsOf = (df: DG.DataFrame): string =>
  Array.from({length: df.rowCount}, (_, i) => `${s(df, 'SUBJECT', i)}: ${s(df, 'REASONS', i)}`).join(' | ');

const UNITS = {concentration: 'ng/mL', time: 'h', dose: 'mg/kg', bodyWeight: 'kg'};
/** A manual mapping of the synthetic table, every role and unit declared; `extra` overrides. */
const manual = (extra: Record<string, unknown> = {}): string => JSON.stringify({
  subject: 'Subject', time: 'Time_hr', concentration: 'Conc_ng_mL', dose: 'Dose_mg_kg', route: 'Route',
  bodyWeight: 'BW_kg', blq: 'BLQ', lloq: 'LLOQ', comment: 'Comment', sex: 'Sex', groupBy: ['Study', 'Compound'],
  units: UNITS, ...extra});

/** The error message `fn` throws ('' when it does not throw). */
function thrown(fn: () => unknown): string {
  try {
    fn();
  } catch (e) {
    return e instanceof Error ? e.message : String(e);
  }
  return '';
}

/** One-sample-per-animal arms (destructive design): 3 animals per nominal time. */
const sparseArms = (study: string, route: Arm['route'], dose: number, prefix: string,
  extra: Partial<Arm> = {}): Arm[] => T_SPARSE.map((t, i) => ({
  study, compound: 'DGX-1008', route, dose, times: [t],
  subjects: [0, 1, 2].map((k) => `${prefix}${String(3 * i + k + 1).padStart(3, '0')}`),
  cl: 0.3, v: 0.3, ka: 0.6, f: 0.12, ...extra,
}));

/** A copy of `df` without `col`. */
function without(df: DG.DataFrame, col: string): DG.DataFrame {
  const c = df.clone();
  c.columns.remove(col);
  return c;
}

/** The mouse destructive fixture with its detected mapping (plus the text BLQ flag). */
async function mouseFixture(): Promise<{df: DG.DataFrame; mapping: Record<string, unknown>}> {
  const df = await loadFixture('04_mouse_sparse_destructive.csv');
  const det = detectPipelineMapping(df);
  expect(det.confidence, 'high', det.issues.join('; '));
  return {df, mapping: {...det.mapping!, blq: 'BLQ'} as unknown as Record<string, unknown>};
}

const compositeRows = (df: DG.DataFrame): number[] => rowsWhere(df, 'SUBJECT', COMPOSITE_SUBJECT);

category('Pipeline functions', () => {
  test('detect: well-formed table → high confidence with study and compound keys', async () => {
    const det = detectPipelineMapping(table(ivPo('HAL-1', 'DGX-1001')));
    expect(det.confidence, 'high', det.issues.join('; '));
    expect(det.mapping!.subject, 'Subject');
    expect(det.mapping!.concentration, 'Conc_ng_mL');
    expect(det.mapping!.units.concentration, 'ng/mL');
    expect(det.mapping!.units.time, 'h');
    expect(det.mapping!.groupBy.join(','), 'Study,Compound');
    expect(det.mapping!.comment, 'Comment');
  });

  test('detect: reassay column → low confidence naming it', async () => {
    const df = table(ivPo('MPS-1', 'DGX-1002'));
    df.columns.add(DG.Column.fromFloat64Array('Conc_ng_mL_reassay', new Float64Array(df.rowCount).fill(5)));
    const det = detectPipelineMapping(df);
    expect(det.confidence, 'low');
    expect(det.issues.some((i) => i.includes('Conc_ng_mL_reassay')), true, det.issues.join('; '));
  });

  test('compute: clean IV/PO study → 6 accepted profiles with engine version', async () => {
    const res = computePipeline(table(ivPo('HAL-1', 'DGX-1001')), '', '');
    expect(res.rowCount, 6);
    for (let i = 0; i < res.rowCount; i++) {
      expect(s(res, 'VERDICT', i), VERDICT_ACCEPTED, `${s(res, 'SUBJECT', i)}: ${s(res, 'REASONS', i)}`);
      for (const p of ['CMAX', 'AUC_LAST', 'AUC_INF', 'HALF_LIFE'])
        expect(Number.isFinite(res.col(p)!.get(i)) && !res.col(p)!.isNone(i), true, `${p} row ${i}`);
      expect(s(res, 'ENGINE_VERSION', i) !== '', true);
    }
  });

  test('compute: canonical unit tags and per-kg CL in mL/min/kg', async () => {
    const res = computePipeline(table(ivPo('HAL-1', 'DGX-1001')), '', '');
    expect(res.col('CMAX')!.getTag('units'), 'mg/L');
    expect(res.col('CL_ML_MIN_KG')!.getTag('units'), 'mL/min/kg');
    // Simulated CL 0.72 L/h/kg = 12 mL/min/kg; Cmax of a 1 mg/kg IV bolus with V 1.4 L/kg ≈ 0.68 mg/L.
    const iv = rowsWhere(res, 'ROUTE', 'IV-bolus');
    for (const i of iv) {
      const cl = res.col('CL_ML_MIN_KG')!.get(i);
      expect(Math.abs(cl - 12) < 1.5, true, `CL ${cl}`);
      expect(res.col('CMAX')!.get(i) < 1, true, `CMAX ${res.col('CMAX')!.get(i)} should be mg/L`);
    }
  });

  test('compute: ambiguous concentration blocks only the study that has it', async () => {
    const df = table([...ivPo('MPS-1', 'DGX-1002'), ...ivPo('HAL-1', 'DGX-1001')]);
    const re = new Float64Array(df.rowCount).fill(DG.FLOAT_NULL);
    rowsWhere(df, 'Study', 'MPS-1').slice(0, 3).forEach((r) => re[r] = 100);
    df.columns.add(DG.Column.fromFloat64Array('Conc_ng_mL_reassay', re));
    const res = computePipeline(df, '', '');
    const mps = rowsWhere(res, 'Study', 'MPS-1');
    expect(mps.length, 1);
    expect(s(res, 'VERDICT', mps[0]), VERDICT_BLOCKED);
    expect(s(res, 'REASONS', mps[0]).includes('Conc_ng_mL_reassay'), true);
    expect(rowsWhere(res, 'Study', 'HAL-1').every((i) => s(res, 'VERDICT', i) === VERDICT_ACCEPTED), true);
  });

  test('compute: sampling too short for t½ → %AUCextrap review, AUClast kept', async () => {
    const res = computePipeline(table(ivPo('HAL-2', 'DGX-1004', {v: 6.0}, 0.18)), '', '');
    const flagged = rowsWhere(res, 'VERDICT', VERDICT_REVIEW);
    expect(flagged.length > 0, true);
    expect(flagged.some((i) => s(res, 'REASONS', i).includes('%AUCextrap > 20%')), true,
      flagged.map((i) => s(res, 'REASONS', i)).join(' | '));
    expect(flagged.every((i) => !res.col('AUC_LAST')!.isNone(i)), true);
  });

  test('compute: exposure outlier and CRO comment → review, still computed', async () => {
    const arms = ivPo('HAL-3', 'DGX-1003');
    arms[1].scale = {R2003: 0.04};
    arms[1].comment = {R2003: 'Partial regurgitation observed at dosing'};
    const res = computePipeline(table(arms), '', '');
    const r = rowsWhere(res, 'SUBJECT', 'R2003')[0];
    expect(s(res, 'VERDICT', r), VERDICT_REVIEW);
    expect(s(res, 'REASONS', r).includes('exposure outlier'), true, s(res, 'REASONS', r));
    expect(s(res, 'REASONS', r).includes('Partial regurgitation'), true);
    expect(res.col('AUC_LAST')!.isNone(r), false);
    for (const other of ['R2001', 'R2002'])
      expect(s(res, 'VERDICT', rowsWhere(res, 'SUBJECT', other)[0]), VERDICT_ACCEPTED);
  });

  test('compute: late spike → terminal rise review', async () => {
    const arms = ivPo('MPS-2', 'DGX-1006');
    arms[1].override = {R2002: {24: 410}};
    const res = computePipeline(table(arms), '', '');
    const r = rowsWhere(res, 'SUBJECT', 'R2002')[0];
    expect(s(res, 'VERDICT', r), VERDICT_REVIEW);
    expect(s(res, 'REASONS', r).includes('terminal rise'), true, s(res, 'REASONS', r));
  });

  test('compute: cassette with compound left out of groupBy → split, never pooled', async () => {
    const arms: Arm[] = ['DGX-1001', 'DGX-1002', 'DGX-1003'].map((c, i) => ({
      study: 'HAL-4', compound: c, route: 'PO', dose: 2, subjects: ['R6001', 'R6002', 'R6003'], cl: 0.4 + 0.3 * i,
    }));
    const df = table(arms, {compoundCol: 'Test Material'});
    const mapping = JSON.stringify({subject: 'Subject', time: 'Time_hr', concentration: 'Conc_ng_mL',
      dose: 'Dose_mg_kg', route: 'Route', bodyWeight: 'BW_kg', blq: 'BLQ', lloq: 'LLOQ', groupBy: ['Study'],
      units: {concentration: 'ng/mL', time: 'h', dose: 'mg/kg', bodyWeight: 'kg'}});
    const res = computePipeline(df, mapping, '');
    expect(res.rowCount, 9);
    expect(res.col('Test Material') !== null, true);
    for (let i = 0; i < res.rowCount; i++) {
      expect(s(res, 'REASONS', i).startsWith(`split by 'Test Material'`), true, s(res, 'REASONS', i));
      expect(res.col('N_SAMPLES')!.get(i) <= 9, true);
    }
  });

  test('compute: crossover with period left out of groupBy → split by period', async () => {
    const arms: Arm[] = [
      {study: 'KL-1', compound: 'DGX-1005', route: 'IV', dose: 0.5, subjects: ['D01', 'D02'], period: 'PERIOD 1'},
      {study: 'KL-1', compound: 'DGX-1005', route: 'PO', dose: 2, subjects: ['D01', 'D02'], period: 'PERIOD 2'},
    ];
    const df = table(arms, {withPeriod: true});
    const mapping = JSON.stringify({subject: 'Subject', time: 'Time_hr', concentration: 'Conc_ng_mL',
      dose: 'Dose_mg_kg', route: 'Route', bodyWeight: 'BW_kg', blq: 'BLQ', groupBy: ['Study', 'Compound'],
      units: {concentration: 'ng/mL', time: 'h', dose: 'mg/kg', bodyWeight: 'kg'}});
    const res = computePipeline(df, mapping, '');
    expect(res.rowCount, 4);
    expect(rowsWhere(res, 'ROUTE', 'PO').length, 2);
    const reasons = Array.from({length: 4}, (_, i) => s(res, 'REASONS', i));
    expect(reasons.every((r) => r.startsWith(`split by 'Period'`)), true);
    // With Period detected into groupBy, the same table needs no review for keying.
    const auto = computePipeline(df, '', '');
    expect(auto.rowCount, 4);
    expect(Array.from({length: 4}, (_, i) => s(auto, 'REASONS', i)).some((r) => r.includes('split by')), false);
  });

  test('compute: unexplained duplicate samples → blocked', async () => {
    const df = table([{study: 'X-1', compound: 'C', route: 'IV', dose: 1, subjects: ['S1', 'S1']}]);
    const res = computePipeline(df, '', '');
    expect(res.rowCount, 1);
    expect(s(res, 'VERDICT', 0), VERDICT_BLOCKED);
    expect(s(res, 'REASONS', 0).includes('duplicate samples'), true);
  });

  test('compute: sparse destructive design → composite row per arm', async () => {
    const times = [0.25, 0.5, 1, 2, 4, 6, 8, 24];
    const arms: Arm[] = times.map((t, i) => ({
      study: 'KL-2', compound: 'DGX-1008', route: 'PO' as const, dose: 10, times: [t],
      subjects: [0, 1, 2].map((k) => `M${String(3 * i + k + 1).padStart(3, '0')}`), cl: 0.3, v: 0.3, ka: 0.6, f: 0.12,
    }));
    const res = computePipeline(table(arms), '', '');
    expect(res.rowCount, 1, Array.from({length: res.rowCount}, (_, i) => s(res, 'REASONS', i)).join(' | '));
    expect(s(res, 'SUBJECT', 0), COMPOSITE_SUBJECT);
    expect(Number.isFinite(res.col('AUC_LAST')!.get(0)) && !res.col('AUC_LAST')!.isNone(0), true);
    // Same canonical units as the dense path: Cmax in mg/L (true ≈ 0.4 mg/L for this curve).
    const cmax = res.col('CMAX')!.get(0);
    expect(cmax > 0.05 && cmax < 5, true, `composite CMAX ${cmax} should be mg/L`);
  });

  test('compute: rules JSON overrides the %AUCextrap threshold', async () => {
    const df = table(ivPo('HAL-2', 'DGX-1004', {v: 6.0}, 0.18));
    const strict = computePipeline(df, '', '');
    const lax = computePipeline(df, '', '{"extrapWarnPct": 95}');
    const n = (r: DG.DataFrame): number => Array.from({length: r.rowCount}, (_, i) => s(r, 'REASONS', i))
      .filter((x) => x.includes('%AUCextrap')).length;
    expect(n(strict) > 0, true);
    expect(n(lax), 0);
  });

  test('compute: acceptanceIgnore moves the listed criteria to INFO_FLAGS', async () => {
    // A 4-point measurable terminal phase over 2 h: span ratio < 2 under the app default.
    const arms = ivPo('HAL-6', 'DGX-1002', {v: 1.8}, 0.48);
    const df = table(arms);
    const app = computePipeline(df, '', '');
    const spanFlagged = rowsWhere(app, 'VERDICT', VERDICT_REVIEW)
      .filter((i) => s(app, 'REASONS', i).includes('span ratio'));
    const policy = computePipeline(df, '', '{"acceptanceIgnore": ["SPAN_LT", "NPTS_LE"]}');
    for (let i = 0; i < policy.rowCount; i++)
      expect(/span ratio|N points/.test(s(policy, 'REASONS', i)), false);
    for (const i of spanFlagged) expect(s(policy, 'INFO_FLAGS', i).includes('span ratio'), true);
  });

  test('summarize: arm rows with N, CV%, N_REVIEW and F% from the IV arm', async () => {
    const arms = [...ivPo('HAL-1', 'DGX-1001'), ...ivPo('HAL-5', 'DGX-1002', {f: 0.62}, 0.48)];
    const res = computePipeline(table(arms), '', '');
    const sum = summarizePipeline(res);
    expect(sum.rowCount, 4);
    for (const c of ['N', 'MEAN_CMAX', 'CV_PCT_AUC_INF', 'MEAN_HALF_LIFE', 'MEAN_CL', 'MEAN_VZ', 'N_REVIEW', 'F_PCT'])
      expect(sum.col(c) !== null, true, `missing ${c}`);
    const po = rowsWhere(sum, 'ROUTE', 'PO');
    expect(po.length, 2);
    const f = (study: string): number => sum.col('F_PCT')!.get(po.find((i) => s(sum, 'Study', i) === study)!);
    expect(Math.abs(f('HAL-1') - 55) < 3, true, `F HAL-1 = ${f('HAL-1')}`);
    expect(Math.abs(f('HAL-5') - 62) < 3, true, `F HAL-5 = ${f('HAL-5')}`);
  });

  test('registered: functions are Flow-listed and callable by name', async () => {
    for (const name of ['ncaDetectMapping', 'ncaCompute', 'ncaSummarize']) {
      const f = DG.Func.find({package: 'Nca', name})[0];
      expect(f !== undefined, true, `${name} not registered`);
      expect(String(f.options['includeInFlow']), 'true');
    }
    const out: DG.DataFrame = await grok.functions.call('Nca:ncaCompute',
      {table: table(ivPo('HAL-1', 'DGX-1001')), mapping: '', rules: ''});
    expect(out.rowCount, 6);
  });

  // ---------- Phase 1: strict mapping and configuration ----------

  test('AC-1.1: a mapping without the concentration unit blocks every study, nothing computed', async () => {
    const df = table([...ivPo('HAL-1', 'DGX-1001'), ...ivPo('HAL-5', 'DGX-1002')]);
    const res = computePipeline(df, manual({units: {time: 'h', dose: 'mg/kg', bodyWeight: 'kg'}}), '');
    expect(res.rowCount, 2, reasonsOf(res));
    for (let i = 0; i < res.rowCount; i++) {
      expect(s(res, 'VERDICT', i), VERDICT_BLOCKED);
      expect(s(res, 'REASONS', i).includes('concentration unit'), true, s(res, 'REASONS', i));
      expect(res.col('AUC_LAST')!.isNone(i), true);
    }
  });

  test('AC-1.1: every missing item is listed in one reason', async () => {
    const res = computePipeline(without(table(ivPo('HAL-1', 'DGX-1001')), 'Route'),
      manual({route: undefined, units: {dose: 'mg/kg', bodyWeight: 'kg'}}), '');
    const reason = s(res, 'REASONS', 0);
    expect(reason, 'mapping: missing concentration unit, time unit, route (no route column and no defaultRoute)');
  });

  test('AC-1.2: no route column and no defaultRoute → blocked, reason names route', async () => {
    const res = computePipeline(without(table(ivPo('HAL-1', 'DGX-1001')), 'Route'), manual({route: undefined}), '');
    expect(res.rowCount, 1);
    expect(s(res, 'VERDICT', 0), VERDICT_BLOCKED);
    expect(s(res, 'REASONS', 0).includes('route'), true, s(res, 'REASONS', 0));
  });

  test('AC-1.3: defaultRoute IV-bolus on a table without a route column → computed as IV bolus', async () => {
    const df = without(table([ivPo('HAL-1', 'DGX-1001')[0]]), 'Route');
    const res = computePipeline(df, manual({route: undefined, defaultRoute: 'IV-bolus'}), '');
    expect(res.rowCount, 3, reasonsOf(res));
    for (let i = 0; i < res.rowCount; i++) {
      expect(s(res, 'ROUTE', i), 'IV-bolus');
      expect(s(res, 'CL_BASIS', i), 'CL');
      expect(s(res, 'VERDICT', i) !== VERDICT_BLOCKED, true, s(res, 'REASONS', i));
    }
  });

  test('AC-1.4: auto-detect on a table without a route column → low confidence naming the route', async () => {
    const det = detectPipelineMapping(without(table(ivPo('HAL-1', 'DGX-1001')), 'Route'));
    expect(det.confidence, 'low');
    expect(det.issues.some((i) => i.includes('route')), true, det.issues.join('; '));
    expect(det.mapping!.units.concentration, 'ng/mL');
  });

  test('detect: an unread body-weight unit is an issue and is not defaulted', async () => {
    const df = table(ivPo('HAL-1', 'DGX-1001'));
    df.col('BW_kg')!.name = 'BW';
    const det = detectPipelineMapping(df);
    expect(det.confidence, 'low');
    expect(det.issues.some((i) => i.includes('body-weight unit')), true, det.issues.join('; '));
    expect(det.mapping!.units.bodyWeight === undefined, true);
  });

  test('AC-1.5: unknown or malformed configuration throws with a readable message', async () => {
    const df = table(ivPo('HAL-1', 'DGX-1001'));
    const typo = thrown(() => computePipeline(df, '', '{"extrapWarnPc": 30}'));
    expect(typo.includes('extrapWarnPc'), true, typo);
    expect(thrown(() => computePipeline(df, '', '{')).startsWith('rules: invalid JSON'), true);
    expect(thrown(() => computePipeline(df, '', '{"lambdaZ": {"minPts": 4}}')).includes('lambdaZ.minPts'), true);
    expect(thrown(() => computePipeline(df, '', '{"extrapWarnPct": "30"}')).includes('must be a number'), true);
    expect(thrown(() => computePipeline(df, manual({conc: 'Conc_ng_mL'}), '')).includes(`unknown key 'conc'`), true);
    expect(thrown(() => computePipeline(df, manual({units: {...UNITS, concentration: 'ng/ml'}}), ''))
      .includes(`unknown concentration unit 'ng/ml'`), true);
    expect(thrown(() => computePipeline(df, '{"subject": ', '')).startsWith('mapping: invalid JSON'), true);
    // The ncaDetectMapping wrapper is still accepted.
    const wrapped = JSON.stringify(detectPipelineMapping(df));
    expect(computePipeline(df, wrapped, '').rowCount, 6);
  });

  // ---------- Phase 2: sex in the results ----------

  test('AC-2.1: rat multi-route auto-detects Sex and carries SEX M/F on every profile', async () => {
    const df = await loadFixture('03_rat_multiroute.csv');
    const det = detectPipelineMapping(df);
    expect(det.mapping!.sex, 'Sex');
    const res = computePipeline(df, '', '');
    expect(res.col('SEX') !== null, true);
    const values = new Set(Array.from({length: res.rowCount}, (_, i) => s(res, 'SEX', i)));
    expect([...values].sort().join(','), 'F,M');
  }, {timeout: TIMEOUT});

  test('AC-2.2: a Gender column gives the same SEX values', async () => {
    const sex = computePipeline(await loadFixture('03_rat_multiroute.csv'), '', '');
    const gender = computePipeline(await loadFixture('03_rat_multiroute_gender.csv'), '', '');
    expect(gender.rowCount, sex.rowCount);
    for (let i = 0; i < sex.rowCount; i++) {
      expect(s(gender, 'SUBJECT', i), s(sex, 'SUBJECT', i));
      expect(s(gender, 'SEX', i), s(sex, 'SEX', i));
    }
  }, {timeout: TIMEOUT});

  test('AC-2.3: an unrecognised sex value → no SEX column, NCA rows unaffected', async () => {
    const arms = ivPo('HAL-1', 'DGX-1001');
    const base = computePipeline(table(arms), '', '');
    arms[1].sex = 'X';
    const res = computePipeline(table(arms), '', '');
    expect(res.col('SEX'), null);
    expect(res.getTag('nca.sexReason'), 'sex values not recognised: X');
    expect(res.rowCount, base.rowCount);
    for (let i = 0; i < res.rowCount; i++) {
      expect(num(res, 'AUC_LAST', i), num(base, 'AUC_LAST', i));
      expect(s(res, 'VERDICT', i), s(base, 'VERDICT', i));
    }
  });

  test('SEX: a subject whose sex varies across its rows gets an empty SEX', async () => {
    const df = table(ivPo('HAL-1', 'DGX-1001'));
    df.col('Sex')!.set(rowsWhere(df, 'Subject', 'R1001')[0], 'F');
    const res = computePipeline(df, '', '');
    expect(s(res, 'SEX', rowsWhere(res, 'SUBJECT', 'R1001')[0]), '');
    expect(s(res, 'SEX', rowsWhere(res, 'SUBJECT', 'R1002')[0]), 'M');
  });

  // ---------- Phase 3: F% and the arm summary ----------

  // AC-3.1 (F% values vs PKNCA) is in the `Pipeline functions: reference` category.
  test('AC-3.4: rat multi-route summary has one arm per Sex × Route × Dose, F% on every oral arm', async () => {
    const sum = summarizePipeline(computePipeline(await loadFixture('03_rat_multiroute.csv'), '', ''));
    // IV 3 doses, PO 3 doses, SC 2 doses, × 2 sexes.
    expect(sum.rowCount, 16);
    for (let i = 0; i < sum.rowCount; i++) {
      const tag = `${s(sum, 'ROUTE', i)}/${s(sum, 'SEX', i)}/${num(sum, 'DOSE', i)}`;
      expect(s(sum, 'F_PCT_REASON', i), '', tag);
      if (s(sum, 'ROUTE', i) === 'IV-bolus') {
        expect(sum.col('F_PCT')!.isNone(i), true, tag);
        expect(sum.col('F_PCT_IV_N')!.isNone(i), true, tag);
      } else {
        expect(Number.isFinite(num(sum, 'F_PCT', i)), true, tag);
        expect(sum.col('F_PCT_IV_N')!.get(i) > 0, true, tag);
      }
    }
  }, {timeout: TIMEOUT});

  test('AC-3.2: without a sex column there is no F%, and the reason says so', async () => {
    const sum = summarizePipeline(computePipeline(without(await loadFixture('03_rat_multiroute.csv'), 'Sex'), '', ''));
    expect(sum.col('SEX'), null);
    for (let i = 0; i < sum.rowCount; i++) {
      expect(sum.col('F_PCT')!.isNone(i), true);
      if (s(sum, 'ROUTE', i) !== 'IV-bolus') expect(s(sum, 'F_PCT_REASON', i), 'no sex column');
    }
  }, {timeout: TIMEOUT});

  test('AC-3.3: a PO arm with no IV arm in its study → no IV reference', async () => {
    const sum = summarizePipeline(computePipeline(table([ivPo('HAL-1', 'DGX-1001')[1]]), '', ''));
    expect(sum.rowCount, 1);
    expect(sum.col('F_PCT')!.isNone(0), true);
    expect(s(sum, 'F_PCT_REASON', 0), 'no IV reference');
  });

  test('AC-3.5: a crossover (IV period 1, PO period 2, same subjects) still gets F%', async () => {
    const arms: Arm[] = [
      {study: 'KL-1', compound: 'DGX-1005', route: 'IV', dose: 0.5, subjects: ['D01', 'D02', 'D03'],
        period: 'PERIOD 1'},
      {study: 'KL-1', compound: 'DGX-1005', route: 'PO', dose: 2, subjects: ['D01', 'D02', 'D03'],
        period: 'PERIOD 2', f: 0.4},
    ];
    // The period key is recognised by the same normalised name rule as detection (`Visit Day` too).
    for (const periodName of ['Period', 'Visit Day']) {
      const df = table(arms, {withPeriod: true});
      df.col('Period')!.name = periodName;
      const sum = summarizePipeline(computePipeline(df, '', ''));
      const po = rowsWhere(sum, 'ROUTE', 'PO');
      expect(po.length, 1, periodName);
      expect(Math.abs(num(sum, 'F_PCT', po[0]) - 40) < 3, true,
        `${periodName}: F ${num(sum, 'F_PCT', po[0])} ${s(sum, 'F_PCT_REASON', po[0])}`);
      expect(sum.col('F_PCT_IV_N')!.get(po[0]), 3);
    }
  });

  test('F%: a sparse composite pooling sexes gets no F%, with the reason', async () => {
    const arms = [...sparseArms('KL-3', 'IV', 1, 'I'), ...sparseArms('KL-3', 'PO', 5, 'P')];
    arms[0].sex = 'F';
    arms[arms.length - 1].sex = 'F';
    const res = computePipeline(table(arms), '', '');
    expect(compositeRows(res).map((i) => s(res, 'SEX', i)).join(','), 'mixed,mixed', reasonsOf(res));
    const sum = summarizePipeline(res);
    const po = rowsWhere(sum, 'ROUTE', 'PO');
    expect(sum.col('F_PCT')!.isNone(po[0]), true);
    expect(s(sum, 'F_PCT_REASON', po[0]), 'composite pools sexes');
  }, {timeout: TIMEOUT});

  test('F%: single-sex composites match on sex like dense arms', async () => {
    const res = computePipeline(table([...sparseArms('KL-3', 'IV', 1, 'I'), ...sparseArms('KL-3', 'PO', 5, 'P')]),
      '', '');
    const sum = summarizePipeline(res);
    const po = rowsWhere(sum, 'ROUTE', 'PO');
    expect(s(sum, 'SEX', po[0]), 'M');
    expect(Number.isFinite(num(sum, 'F_PCT', po[0])), true, s(sum, 'F_PCT_REASON', po[0]));
  }, {timeout: TIMEOUT});

  // ---------- Phase 4: sparse path ----------

  // AC-4.1 (mouse composite vs PKNCA) is in the `Pipeline functions: reference` category.
  // The mouse fixture is PO: Vss is route-gated (intravascular only), so the
  // body-weight-dependent volume compared here is Vz; Vss must be empty in both runs.
  test('AC-4.2: body weight in g gives the same CL / Vz as in kg', async () => {
    const {df, mapping} = await mouseFixture();
    const g = computePipeline(df, JSON.stringify(mapping), '');
    const bw = df.col('BodyWeight_g')!;
    df.columns.add(DG.Column.fromFloat64Array('BodyWeight_kg',
      Float64Array.from({length: df.rowCount}, (_, i) => bw.get(i) / 1000)));
    const kgMap = {...mapping, bodyWeight: 'BodyWeight_kg',
      units: {...(mapping.units as object), bodyWeight: 'kg'}};
    const kg = computePipeline(without(df, 'BodyWeight_g'), JSON.stringify(kgMap), '');
    const rg = compositeRows(g); const rk = compositeRows(kg);
    expect(rg.length, 2);
    for (let j = 0; j < rg.length; j++) {
      for (const p of ['CL', 'VZ', 'CL_ML_MIN_KG', 'VZ_L_KG', 'CL_ABS', 'VZ_ABS', 'BW_KG']) {
        const a = num(g, p, rg[j]); const b = num(kg, p, rk[j]);
        expect(Number.isFinite(a) && a > 0, true, `${p} ${a}`);
        expectFloat(a, b, Math.abs(b) * 1e-9, `${p} g vs kg`);
      }
      expect(g.col('VSS')!.isNone(rg[j]) && kg.col('VSS')!.isNone(rk[j]), true, 'VSS empty for PO');
    }
  }, {timeout: TIMEOUT});

  test('AC-4.3: nominal time in minutes gives the same t½ and AUC as in hours', async () => {
    const {df, mapping} = await mouseFixture();
    const h = computePipeline(df, JSON.stringify(mapping), '');
    for (const c of [mapping.time as string, mapping.nominalTime as string]) {
      const col = df.col(c)!;
      const v = Float64Array.from({length: df.rowCount}, (_, i) => col.get(i) * 60);
      df.columns.remove(c);
      df.columns.add(DG.Column.fromFloat64Array(c, v));
    }
    const minMap = {...mapping, units: {...(mapping.units as object), time: 'min'}};
    const min = computePipeline(df, JSON.stringify(minMap), '');
    const rh = compositeRows(h); const rm = compositeRows(min);
    expect(rm.length, rh.length);
    for (let j = 0; j < rh.length; j++) {
      for (const p of ['HALF_LIFE', 'AUC_LAST', 'AUC_INF', 'TMAX']) {
        const a = num(min, p, rm[j]); const b = num(h, p, rh[j]);
        expect(Number.isFinite(b), true, `${p} hours run`);
        expectFloat(a, b, Math.abs(b) * 1e-9, `${p} min vs h`);
      }
    }
  }, {timeout: TIMEOUT});

  test('AC-4.4: sparse IV and PO arms → one composite per route, also at the same dose', async () => {
    for (const poDose of [5, 1]) {
      const arms = [...sparseArms('KL-4', 'IV', 1, 'I'), ...sparseArms('KL-4', 'PO', poDose, 'P')];
      const res = computePipeline(table(arms), '', '');
      const rows = compositeRows(res);
      expect(rows.length, 2, reasonsOf(res));
      expect(rows.map((i) => s(res, 'ROUTE', i)).sort().join(','), 'IV-bolus,PO');
      for (const i of rows) expect(s(res, 'CL_BASIS', i), s(res, 'ROUTE', i) === 'PO' ? 'CL/F' : 'CL');
    }
  }, {timeout: TIMEOUT});

  test('AC-4.5: a sparse IV-infusion arm → one blocked row naming infusion', async () => {
    const res = computePipeline(table(sparseArms('KL-5', 'IV-infusion', 1, 'I')), '', '');
    expect(res.rowCount, 1, reasonsOf(res));
    expect(s(res, 'VERDICT', 0), VERDICT_BLOCKED);
    expect(s(res, 'REASONS', 0).includes('infusion'), true, s(res, 'REASONS', 0));
  });

  test('AC-4.6: batch design → composite; common 9-point schedule → dense; mixed → blocked', async () => {
    // Batch: six animals, two samples each, schedules differ between animals.
    const schedule: Record<string, number[]> = {A: [1, 2], B: [1, 2], C: [1, 4], D: [2, 4], E: [2, 4], F: [1, 4]};
    const batch = computePipeline(table(Object.entries(schedule).map(([a, times]) => ({
      study: 'KL-6', compound: 'DGX-1008', route: 'PO' as const, dose: 10, subjects: [a], times}))), '', '');
    expect(compositeRows(batch).length, 1, reasonsOf(batch));

    const serial = computePipeline(table(ivPo('KL-7', 'DGX-1001')), '', '');
    expect(compositeRows(serial).length, 0);
    expect(serial.rowCount, 6);

    const mixed = computePipeline(table([...sparseArms('KL-8', 'PO', 10, 'M'),
      {study: 'KL-8', compound: 'DGX-1008', route: 'PO', dose: 10, subjects: ['S01', 'S02', 'S03']}]), '', '');
    expect(mixed.rowCount, 1, reasonsOf(mixed));
    expect(s(mixed, 'VERDICT', 0), VERDICT_BLOCKED);
    expect(s(mixed, 'REASONS', 0).includes(`declare 'design'`), true, s(mixed, 'REASONS', 0));
    // Declaring the design resolves it.
    const declared = computePipeline(table([...sparseArms('KL-8', 'PO', 10, 'M'),
      {study: 'KL-8', compound: 'DGX-1008', route: 'PO', dose: 10, subjects: ['S01', 'S02', 'S03']}]),
    manual({design: 'serial'}), '');
    expect(rowsWhere(declared, 'VERDICT', VERDICT_BLOCKED).length, 0, reasonsOf(declared));
  }, {timeout: TIMEOUT});

  test('sparse: fewer than 3 single-sample animals → blocked, not a 1–2 point composite', async () => {
    const res = computePipeline(table([{study: 'KL-9', compound: 'DGX-1008', route: 'PO', dose: 10,
      subjects: ['M001', 'M002'], times: [2]}]), '', '');
    expect(res.rowCount, 1, reasonsOf(res));
    expect(s(res, 'VERDICT', 0), VERDICT_BLOCKED);
    expect(s(res, 'REASONS', 0).includes('fewer than 3 animals'), true, s(res, 'REASONS', 0));
  });

  test('sparse: a declared topology that contradicts the data is a review reason', async () => {
    const schedule: Record<string, number[]> = {A: [1, 2], B: [1, 2], C: [1, 4], D: [2, 4], E: [2, 4], F: [1, 4]};
    const df = table(Object.entries(schedule).map(([a, times]) => ({
      study: 'KL-6', compound: 'DGX-1008', route: 'PO' as const, dose: 10, subjects: [a], times})));
    const res = computePipeline(df, manual({design: 'sparse-destructive'}), '');
    const r = compositeRows(res);
    expect(r.length, 1, reasonsOf(res));
    expect(s(res, 'REASONS', r[0]).includes('sampling design'), true, s(res, 'REASONS', r[0]));
  }, {timeout: TIMEOUT});

  // ---------- Phase 5: heuristic thresholds and CRO comments ----------

  test('AC-5.1: exposureOutlierFold raises the bar of the screening heuristic', async () => {
    const arms = ivPo('HAL-3', 'DGX-1003');
    arms[1].scale = {R2003: 0.04};
    const df = table(arms);
    const r = (res: DG.DataFrame): number => rowsWhere(res, 'SUBJECT', 'R2003')[0];
    const dflt = computePipeline(df, '', '');
    expect(s(dflt, 'REASONS', r(dflt)).includes('screening heuristic: exposure outlier'), true,
      s(dflt, 'REASONS', r(dflt)));
    const lax = computePipeline(df, '', '{"exposureOutlierFold": 50}');
    expect(s(lax, 'REASONS', r(lax)).includes('exposure outlier'), false);
  });

  test('AC-5.2: a fold of 1 or less is a configuration error', async () => {
    const df = table(ivPo('HAL-1', 'DGX-1001'));
    expect(thrown(() => computePipeline(df, '', '{"terminalRiseFold": 1}')).includes('terminalRiseFold'), true);
    expect(thrown(() => computePipeline(df, '', '{"exposureOutlierFold": "4"}')).includes('exposureOutlierFold'), true);
  });

  test('terminalRiseFold: a laxer threshold drops the terminal-rise flag', async () => {
    const arms = ivPo('MPS-2', 'DGX-1006');
    arms[1].override = {R2002: {24: 410}};
    const df = table(arms);
    const r = (res: DG.DataFrame): number => rowsWhere(res, 'SUBJECT', 'R2002')[0];
    const dflt = computePipeline(df, '', '');
    expect(s(dflt, 'REASONS', r(dflt)).includes('screening heuristic: terminal rise'), true);
    const lax = computePipeline(df, '', '{"terminalRiseFold": 1000}');
    expect(s(lax, 'REASONS', r(lax)).includes('terminal rise'), false);
  });

  const commented = (): DG.DataFrame => {
    const arms = ivPo('HAL-7', 'DGX-1001');
    arms[1].comment = {R2001: 'Partial regurgitation observed at dosing', R2002: 'Terminal bleed'};
    return table(arms);
  };

  test('AC-5.3: a dosing-incident comment → REASONS, verdict Review (default pattern)', async () => {
    const res = computePipeline(commented(), '', '');
    const r = rowsWhere(res, 'SUBJECT', 'R2001')[0];
    expect(s(res, 'VERDICT', r), VERDICT_REVIEW);
    expect(s(res, 'REASONS', r), 'CRO comment: "Partial regurgitation observed at dosing"');
  });

  test('AC-5.4: a routine comment → INFO_FLAGS only, verdict unaffected', async () => {
    const res = computePipeline(commented(), '', '');
    const r = rowsWhere(res, 'SUBJECT', 'R2002')[0];
    expect(s(res, 'VERDICT', r), VERDICT_ACCEPTED, s(res, 'REASONS', r));
    expect(s(res, 'REASONS', r), '');
    // INFO_FLAGS may also carry engine warnings (nca-pipeline-followups phase 3).
    expect(s(res, 'INFO_FLAGS', r).includes('CRO comment: "Terminal bleed"'), true, s(res, 'INFO_FLAGS', r));
  });

  test('AC-5.5: reviewCommentPattern replaces the default pattern', async () => {
    const res = computePipeline(commented(), '', '{"reviewCommentPattern": "bleed"}');
    const bleed = rowsWhere(res, 'SUBJECT', 'R2002')[0];
    const regurg = rowsWhere(res, 'SUBJECT', 'R2001')[0];
    expect(s(res, 'REASONS', bleed).includes('Terminal bleed'), true);
    expect(s(res, 'VERDICT', bleed), VERDICT_REVIEW);
    expect(s(res, 'INFO_FLAGS', regurg).includes('Partial regurgitation'), true);
    expect(s(res, 'VERDICT', regurg), VERDICT_ACCEPTED);
  });

  test('AC-5.6: an invalid reviewCommentPattern is a configuration error', async () => {
    const msg = thrown(() => computePipeline(commented(), '', '{"reviewCommentPattern": "("}'));
    expect(msg.startsWith('rules: invalid reviewCommentPattern'), true, msg);
  });

  // ---------- nca-pipeline-followups ----------

  test('followups AC-1.1: a serial study with one dropout animal stays serial', async () => {
    const arms: Arm[] = [
      {study: 'FU-1', compound: 'DGX-1001', route: 'PO', dose: 5, subjects: ['R2001', 'R2002', 'R2003', 'R2004',
        'R2005']},
      {study: 'FU-1', compound: 'DGX-1001', route: 'PO', dose: 5, subjects: ['R2006'], times: [0.25, 0.5]},
    ];
    const res = computePipeline(table(arms), '', '');
    expect(res.rowCount, 6, reasonsOf(res));
    expect(rowsWhere(res, 'VERDICT', VERDICT_BLOCKED).length, 0, reasonsOf(res));
    const short = rowsWhere(res, 'SUBJECT', 'R2006')[0];
    expect(/λz not estimable|NCA failed/.test(s(res, 'REASONS', short)), true, s(res, 'REASONS', short));
  });

  const withSingles = (single: Record<string, number[]>): DG.DataFrame => table([
    {study: 'FU-2', compound: 'DGX-1001', route: 'PO', dose: 5, subjects: ['R2001', 'R2002', 'R2003']},
    ...Object.entries(single).map(([a, times]) =>
      ({study: 'FU-2', compound: 'DGX-1001', route: 'PO' as const, dose: 5, subjects: [a], times})),
  ]);

  test('followups AC-1.2: complete animals in the minority (3 of 7) → mixed design, blocked', async () => {
    const res = computePipeline(withSingles({S1: [1], S2: [2], S3: [4], S4: [8]}), '', '');
    expect(res.rowCount, 1, reasonsOf(res));
    expect(s(res, 'VERDICT', 0), VERDICT_BLOCKED);
    expect(s(res, 'REASONS', 0).includes(`declare 'design'`), true, s(res, 'REASONS', 0));
  });

  test('followups AC-1.3: exactly half complete → blocked (strict majority required)', async () => {
    const res = computePipeline(withSingles({S1: [1, 2], S2: [2, 4], S3: [4, 8]}), '', '');
    expect(res.rowCount, 1, reasonsOf(res));
    expect(s(res, 'VERDICT', 0), VERDICT_BLOCKED);
  });

  const sexArm = (withSex: boolean): DG.DataFrame => {
    const df = table([
      {study: 'FU-3', compound: 'DGX-1001', route: 'PO', dose: 5, subjects: ['R2001', 'R2002', 'R2003'], sex: 'M'},
      {study: 'FU-3', compound: 'DGX-1001', route: 'PO', dose: 5, subjects: ['R2004'], sex: 'F', scale: {R2004: 5}},
    ]);
    return withSex ? df : without(df, 'Sex');
  };

  test('followups AC-2.1: a sex difference in exposure is not an exposure outlier', async () => {
    const res = computePipeline(sexArm(true), '', '');
    for (let i = 0; i < res.rowCount; i++) {
      expect(s(res, 'REASONS', i).includes('exposure outlier'), false,
        `${s(res, 'SUBJECT', i)}: ${s(res, 'REASONS', i)}`);
    }
  });

  test('followups AC-2.3: without a sex column the same arm still flags the 5x profile', async () => {
    const res = computePipeline(sexArm(false), '', '');
    expect(res.col('SEX'), null);
    const r = rowsWhere(res, 'SUBJECT', 'R2004')[0];
    expect(s(res, 'REASONS', r).includes('exposure outlier'), true, s(res, 'REASONS', r));
  });

  test('followups AC-3.1: per-kg dose without body weight → PER_KG_NO_BW in INFO_FLAGS only', async () => {
    const res = computePipeline(without(table(ivPo('FU-4', 'DGX-1001')), 'BW_kg'), '', '');
    for (let i = 0; i < res.rowCount; i++) {
      expect(s(res, 'INFO_FLAGS', i).includes('PER_KG_NO_BW: '), true, s(res, 'INFO_FLAGS', i));
      expect(s(res, 'REASONS', i).includes('PER_KG_NO_BW'), false);
      // The clean simulated study is Accepted; the per-kg note does not change that.
      expect(s(res, 'VERDICT', i), VERDICT_ACCEPTED, s(res, 'REASONS', i));
    }
  });

  test('followups AC-3.2: ng/mL data mapped as mg/L → CONC_MAGNITUDE heuristic in REASONS, Review', async () => {
    const df = table(ivPo('FU-5', 'DGX-1001', {dose: 10}));
    const res = computePipeline(df, manual({units: {...UNITS, concentration: 'mg/L'}}), '');
    const iv = rowsWhere(res, 'ROUTE', 'IV-bolus');
    expect(iv.length, 3, reasonsOf(res));
    for (const i of iv) {
      expect(s(res, 'REASONS', i).includes('CONC_MAGNITUDE_HEURISTIC_WARNING: '), true, s(res, 'REASONS', i));
      expect(s(res, 'VERDICT', i), VERDICT_REVIEW);
    }
  });

  test('followups AC-3.3: body weight 250 declared as kg → BW_UNIT heuristic in REASONS, Review', async () => {
    const df = table(ivPo('FU-6', 'DGX-1001'));
    df.columns.remove('BW_kg');
    df.columns.add(DG.Column.fromFloat64Array('BW_kg', new Float64Array(df.rowCount).fill(250)));
    const res = computePipeline(df, '', '');
    for (let i = 0; i < res.rowCount; i++) {
      expect(s(res, 'REASONS', i).includes('BW_UNIT_HEURISTIC_WARNING: '), true, s(res, 'REASONS', i));
      expect(s(res, 'VERDICT', i), VERDICT_REVIEW);
    }
  });

  test('followups AC-3.4: bare IV route → BARE_IV_DEFAULTED_TO_BOLUS in INFO_FLAGS, verdict unchanged', async () => {
    const res = computePipeline(table(ivPo('FU-7', 'DGX-1001')), '', '');
    for (const i of rowsWhere(res, 'ROUTE', 'IV-bolus')) {
      expect(s(res, 'INFO_FLAGS', i).includes('BARE_IV_DEFAULTED_TO_BOLUS: '), true, s(res, 'INFO_FLAGS', i));
      expect(s(res, 'VERDICT', i), VERDICT_ACCEPTED, s(res, 'REASONS', i));
    }
  });

  test('followups AC-3.5: composite warnings in INFO_FLAGS, bootstrap notes left out', async () => {
    const {df, mapping} = await mouseFixture();
    const res = computePipeline(df, JSON.stringify(mapping), '');
    expect(compositeRows(res).length, 2, reasonsOf(res));
    for (const i of compositeRows(res)) {
      expect(s(res, 'INFO_FLAGS', i).includes('SPARSE_COMPOSITE_BW: '), true, s(res, 'INFO_FLAGS', i));
      expect(s(res, 'INFO_FLAGS', i).includes('SPARSE_BOOTSTRAP_'), false);
    }
  }, {timeout: TIMEOUT});

  test('followups AC-3.6: warnings that restate an acceptance criterion are not repeated', async () => {
    const res = computePipeline(table(ivPo('HAL-2', 'DGX-1004', {v: 6.0}, 0.18)), '', '');
    let extrap = 0;
    for (let i = 0; i < res.rowCount; i++) {
      const all = `${s(res, 'REASONS', i)} ${s(res, 'INFO_FLAGS', i)}`;
      expect(/AUC_EXTRAP_HIGH:|LAMBDAZ_FEW_POINTS:/.test(all), false, all);
      const times = s(res, 'REASONS', i).split('%AUCextrap > 20%').length - 1;
      expect(times <= 1, true, `stated ${times}x: ${s(res, 'REASONS', i)}`);
      if (times === 1) extrap++;
    }
    expect(extrap > 0, true, reasonsOf(res));
  });

  test('followups AC-3.8: implied CL far below the plausible band → CL_IMPLAUSIBLE heuristic in REASONS', async () => {
    // A low dose keeps the mis-scaled concentrations under the 500 mg/L magnitude ceiling.
    const df = table(ivPo('FU-8', 'DGX-1001', {dose: 0.1}));
    const res = computePipeline(df, manual({units: {...UNITS, concentration: 'mg/L'}}), '');
    const iv = rowsWhere(res, 'ROUTE', 'IV-bolus');
    for (const i of iv) {
      expect(s(res, 'REASONS', i).includes('CL_IMPLAUSIBLE_HEURISTIC_WARNING: '), true, s(res, 'REASONS', i));
      expect(s(res, 'REASONS', i).includes('CONC_MAGNITUDE_HEURISTIC_WARNING'), false, s(res, 'REASONS', i));
      expect(s(res, 'VERDICT', i), VERDICT_REVIEW);
    }
  });

  test('followups AC-4.1: an infusion duration in the wrong time unit → Review', async () => {
    const df = table([{study: 'FU-9', compound: 'DGX-1001', route: 'IV-infusion', dose: 1,
      subjects: ['R1001', 'R1002', 'R1003']}]);
    // 100 h of infusion over a 24 h window: minutes entered as hours, or the wrong column.
    df.columns.add(DG.Column.fromFloat64Array('InfDur', new Float64Array(df.rowCount).fill(100)));
    const res = computePipeline(df, manual({infusionDuration: 'InfDur'}), '');
    expect(res.rowCount, 3, reasonsOf(res));
    for (let i = 0; i < res.rowCount; i++) {
      expect(s(res, 'REASONS', i).includes('INFUSION_DURATION_HEURISTIC_WARNING: '), true, s(res, 'REASONS', i));
      expect(s(res, 'VERDICT', i), VERDICT_REVIEW);
    }
  });

  test('followups AC-4.2: REASONS / INFO_FLAGS text uses " | ", the _JSON columns hold the arrays', async () => {
    const arms = ivPo('FU-10', 'DGX-1001');
    arms[0].comment = {R1001: 'Terminal bleed'};
    const res = computePipeline(table(arms), '', '');
    for (let i = 0; i < res.rowCount; i++) {
      for (const [text, json] of [['REASONS', 'REASONS_JSON'], ['INFO_FLAGS', 'INFO_FLAGS_JSON']]) {
        const arr = JSON.parse(s(res, json, i)) as string[];
        expect(Array.isArray(arr), true, `${json} row ${i}`);
        expect(arr.join(' | '), s(res, text, i), `${text} row ${i}`);
      }
    }
    // R1001 carries a comment and the bare-IV note: two entries, one per array element.
    const r = rowsWhere(res, 'SUBJECT', 'R1001')[0];
    const info = JSON.parse(s(res, 'INFO_FLAGS_JSON', r)) as string[];
    expect(info.length >= 2, true, s(res, 'INFO_FLAGS', r));
    expect(info.some((x) => x.startsWith('BARE_IV_DEFAULTED_TO_BOLUS: ')), true, info.join(' || '));
  });

  test('followups AC-4.3: sparse — ng/mL declared as mg/L → unit-error heuristics, Review', async () => {
    const {df, mapping} = await mouseFixture();
    const res = computePipeline(df, JSON.stringify({...mapping, units: {...(mapping.units as object),
      concentration: 'mg/L'}}), '');
    const rows = compositeRows(res);
    expect(rows.length, 2, reasonsOf(res));
    for (const i of rows) {
      expect(s(res, 'REASONS', i).includes('CL_IMPLAUSIBLE_HEURISTIC_WARNING: '), true, s(res, 'REASONS', i));
      expect(s(res, 'VERDICT', i), VERDICT_REVIEW);
    }
    expect(rows.some((i) => s(res, 'REASONS', i).includes('CONC_MAGNITUDE_HEURISTIC_WARNING: ')), true,
      rows.map((i) => s(res, 'REASONS', i)).join(' || '));
  }, {timeout: TIMEOUT});

  test('followups AC-4.4: sparse — body weight in g declared as kg → BW_UNIT heuristic on the composites', async () => {
    const {df, mapping} = await mouseFixture();
    const res = computePipeline(df, JSON.stringify({...mapping, units: {...(mapping.units as object),
      bodyWeight: 'kg'}}), '');
    const rows = compositeRows(res);
    expect(rows.length, 2, reasonsOf(res));
    for (const i of rows) {
      expect(s(res, 'REASONS', i).includes('BW_UNIT_HEURISTIC_WARNING: '), true, s(res, 'REASONS', i));
      expect(s(res, 'VERDICT', i), VERDICT_REVIEW);
    }
  }, {timeout: TIMEOUT});

  test('followups AC-4.5: sparse — a correctly mapped composite raises no unit-error heuristic', async () => {
    const {df, mapping} = await mouseFixture();
    const res = computePipeline(df, JSON.stringify(mapping), '');
    const rows = compositeRows(res);
    expect(rows.length, 2);
    for (const i of rows)
      expect(/HEURISTIC_WARNING/.test(s(res, 'REASONS', i)), false, s(res, 'REASONS', i));
  }, {timeout: TIMEOUT});

  test('CRO comments on a composite: routine notes are information', async () => {
    const {df, mapping} = await mouseFixture();
    const res = computePipeline(df, JSON.stringify(mapping), '');
    for (const i of compositeRows(res)) {
      expect(s(res, 'REASONS', i).includes('CRO comment'), false, s(res, 'REASONS', i));
      expect(s(res, 'INFO_FLAGS', i).includes('CRO comment: "Terminal bleed"'), true, s(res, 'INFO_FLAGS', i));
    }
  }, {timeout: TIMEOUT});
});

/** The mouse fixture's two composites against `mouse_sparse-reference.json`. */
async function assertMouseComposite(res: DG.DataFrame): Promise<void> {
  const golden = JSON.parse(await grok.dapi.files.readAsText(`${FIXTURE_PATH_PREFIX}/mouse_sparse-reference.json`)) as
    {destructive: Array<{dose_mg_kg: number; sparse_auclast: number;
      composite: {nominal_time: number[]; mean: number[]}}>};
  const rows = compositeRows(res);
  expect(rows.length, golden.destructive.length, reasonsOf(res));
  for (const g of golden.destructive) {
    const r = rows.find((i) => num(res, 'DOSE', i) === g.dose_mg_kg);
    expect(r !== undefined, true, `composite for ${g.dose_mg_kg} mg/kg`);
    expect(s(res, 'ROUTE', r!), 'PO');
    // The PKNCA mouse reference pools both sexes per dose group, as the composite does.
    expect(s(res, 'SEX', r!), 'mixed');
    expectFloat(num(res, 'AUC_LAST', r!), g.sparse_auclast, Math.abs(g.sparse_auclast) * 0.001,
      `${g.dose_mg_kg} mg/kg: AUClast (0.1%)`);
    // The oracle carries the composite mean profile; Cmax / Tmax are its maximum.
    let iMax = 0;
    g.composite.mean.forEach((v, i) => {if (v > g.composite.mean[iMax]) iMax = i;});
    expectFloat(num(res, 'CMAX', r!), g.composite.mean[iMax], g.composite.mean[iMax] * 0.001,
      `${g.dose_mg_kg} mg/kg: Cmax (0.1%)`);
    expect(num(res, 'TMAX', r!), g.composite.nominal_time[iMax], `${g.dose_mg_kg} mg/kg: Tmax`);
  }
}

/** One dense result row against its PKNCA reference profile (tolerances of reference-suite-tests.ts).
 *  The references' CL / Vz / Vss are absolute (mg doses) → asserted on the companions; `bwKg` (the
 *  subject's weight as recorded in the fixture) adds the per-kg IDENTITY (nca-cl-bw-normalised-units). */
function assertParity(res: DG.DataFrame, i: number, ref: FixtureProfile, tag: string, bwKg?: number): void {
  const p = ref.parameters;
  expect(s(res, 'VERDICT', i) !== VERDICT_BLOCKED, true, `${tag}: ${s(res, 'REASONS', i)}`);
  expect(s(res, 'STATUS', i), 'ok', `${tag}: STATUS`);
  if (p.cmax !== null) expect(num(res, 'CMAX', i), p.cmax, `${tag}: Cmax`);
  if (p.tmax !== null) expect(num(res, 'TMAX', i), p.tmax, `${tag}: Tmax`);
  const rel = (col: string, want: number | null | undefined, tol: number): void => {
    if (want == null) return;
    expectFloat(num(res, col, i), want, Math.abs(want) * tol, `${tag}: ${col} out of ${tol * 100}%`);
  };
  rel('AUC_LAST', p.auclast, 0.001);
  rel('AUC_INF', p.aucinf, 0.001);
  rel('LAMBDA_Z', p.lambda_z, 0.005);
  rel('HALF_LIFE', p.half_life, 0.01);
  rel('CL_ABS', p.cl, 0.01);
  rel('VZ_ABS', p.vz, 0.01);
  if (bwKg !== undefined && p.cl != null) rel('CL', p.cl / bwKg * 1000 / 60, 0.01);
  if (bwKg !== undefined && p.vz != null) rel('VZ', p.vz / bwKg, 0.01);
  rel('MRT', p.mrt, 0.01);
  if (p.vss === null) expect(res.col('VSS')!.isNone(i), true, `${tag}: VSS empty for extravascular`);
  else rel('VSS_ABS', p.vss, 0.01);
  if (p.pct_aucextrap !== null)
    expectFloat(num(res, 'PCT_EXTRAP', i), p.pct_aucextrap, 0.5, `${tag}: %AUCextrap out of 0.5 pp`);
}

/** Every reference profile has exactly one matching result row, and it passes parity. */
async function assertFixture(res: DG.DataFrame, refName: string, n: number,
  bwBySubject?: ReadonlyMap<string, number>): Promise<void> {
  const byKey = indexFixtureProfiles(await loadReference(refName));
  expect(res.rowCount, n, reasonsOf(res));
  for (let i = 0; i < res.rowCount; i++) {
    const ref = byKey.get(s(res, 'SUBJECT', i));
    expect(ref !== undefined, true, `${refName}: subject '${s(res, 'SUBJECT', i)}' not in the reference`);
    assertParity(res, i, ref!, `${refName} ${s(res, 'SUBJECT', i)}`, bwBySubject?.get(s(res, 'SUBJECT', i)));
  }
}

/** A constant numeric column (a scalar the fixture's reference run takes, as a mapped column). */
function addConst(df: DG.DataFrame, name: string, v: number): DG.DataFrame {
  df.columns.add(DG.Column.fromFloat64Array(name, new Float64Array(df.rowCount).fill(v)));
  return df;
}

category('Pipeline functions: reference', () => {
  test('AC-6.1: theoph — 12 subjects within PKNCA tolerances through ncaCompute', async () => {
    const df = await loadFixture('theoph.csv');
    const wt = new Map<string, number>();
    for (let r = 0; r < df.rowCount; r++) wt.set(String(df.col('Subject')!.get(r)), df.col('Wt')!.get(r) as number);
    const res = computePipeline(df, JSON.stringify({
      subject: 'Subject', time: 'Time', concentration: 'conc', dose: 'Dose', bodyWeight: 'Wt', groupBy: [],
      defaultRoute: 'PO', units: {concentration: 'mg/L', time: 'h', dose: 'mg/kg', bodyWeight: 'kg'}}), '');
    await assertFixture(res, 'theoph-reference.json', 12, wt);
  }, {timeout: TIMEOUT});

  test('AC-6.1: indometh — 6 IV-bolus subjects within PKNCA tolerances', async () => {
    // The reference run's scalar 25 mg dose, carried as a mapped column.
    const res = computePipeline(addConst(await loadFixture('indometh.csv'), 'Dose', 25), JSON.stringify({
      subject: 'Subject', time: 'time', concentration: 'conc', dose: 'Dose', groupBy: [],
      defaultRoute: 'IV-bolus', units: {concentration: 'mg/L', time: 'h', dose: 'mg'}}), '');
    await assertFixture(res, 'indometh-reference.json', 6);
  }, {timeout: TIMEOUT});

  test('AC-6.1: IV infusion — 100 mg over 1 h, MRT / Vss with the −T_inf/2 correction', async () => {
    const df = addConst(addConst(await loadFixture('04_iv_infusion.csv'), 'Dose', 100), 'InfDur', 1);
    const res = computePipeline(df, JSON.stringify({
      subject: 'Subject', time: 'time', concentration: 'conc', dose: 'Dose', infusionDuration: 'InfDur',
      groupBy: [], defaultRoute: 'IV-infusion', units: {concentration: 'mg/L', time: 'h', dose: 'mg'}}), '');
    expect(s(res, 'ROUTE', 0), 'IV-infusion');
    await assertFixture(res, '04_iv_infusion-reference.json', 1);
  }, {timeout: TIMEOUT});

  test('AC-6.1 / AC-3.1: rat multi-route F% per Sex × route × dose (auto mapping) within 1% of PKNCA', async () => {
    const res = computePipeline(await loadFixture('03_rat_multiroute.csv'), '', '');
    expect(rowsWhere(res, 'VERDICT', VERDICT_BLOCKED).length, 0, reasonsOf(res));
    const sum = summarizePipeline(res);
    const ref = DG.DataFrame.fromCsv(
      await grok.dapi.files.readAsText(`${FIXTURE_PATH_PREFIX}/rat_multiroute_bioavailability.csv`));
    expect(ref.rowCount, 10);
    for (let j = 0; j < ref.rowCount; j++) {
      const row = Array.from({length: sum.rowCount}, (_, i) => i).find((i) =>
        s(sum, 'ROUTE', i) === s(ref, 'Route', j) && s(sum, 'SEX', i) === s(ref, 'Sex', j) &&
        num(sum, 'DOSE', i) === num(ref, 'Dose_mg_kg', j));
      expect(row !== undefined, true, `arm ${s(ref, 'Route', j)}/${s(ref, 'Sex', j)}/${num(ref, 'Dose_mg_kg', j)}`);
      const want = num(ref, 'F_pct', j);
      expectFloat(num(sum, 'F_PCT', row!), want, Math.abs(want) * 0.01,
        `F% ${s(ref, 'Route', j)}/${s(ref, 'Sex', j)}/${num(ref, 'Dose_mg_kg', j)}`);
    }
  }, {timeout: TIMEOUT});

  test('AC-6.1 / AC-4.1: mouse destructive composite (auto mapping) vs PKNCA', async () => {
    const {df, mapping} = await mouseFixture();
    await assertMouseComposite(computePipeline(df, JSON.stringify(mapping), ''));
  }, {timeout: TIMEOUT});
});

/** The runner's login — the default `owner` must be the runner or one of its groups (PR-11). */
async function myLogin(): Promise<string> {
  const u = await grok.dapi.users.current();
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  return ((u as any).login as string | undefined) ?? u.name;
}

/** Turn a recipe into an AnalysisRecipe-typed value (the pipeline recipe already is one). */
const asRecipe = (r: AnalysisRecipe): AnalysisRecipe => r;

/**
 * nca-cl-bw-normalised-units F5 (AD-6): the dense pipeline has no normalisation of its own — it reads
 * the app's results frame. E4 (BW in grams divided raw: 1000× low) and E5 (no CL unit on a mixed
 * absolute / per-kg frame) are asserted on columns that already existed before the fix
 * (`CL_ML_MIN_KG`, the `CL` tag), so each test is the failing reproduction AND the regression guard.
 */
category('Pipeline functions: per-kg output contract', () => {
  /** 03_rat_iv through ncaCompute, BW in `bwUnit`; optionally one subject's BW blanked. */
  async function ratIv(bwUnit: 'g' | 'kg', blank?: string): Promise<DG.DataFrame> {
    const df = await loadFixture('03_rat_iv.csv');
    const bw = df.col('BodyWeight_g')!;
    const subj = df.col('SubjectID')!;
    df.columns.add(DG.Column.fromFloat64Array('BW', Float64Array.from({length: df.rowCount}, (_, r) =>
      String(subj.get(r)) === blank ? DG.FLOAT_NULL : (bw.get(r) as number) / (bwUnit === 'kg' ? 1000 : 1))));
    return computePipeline(df, JSON.stringify({
      subject: 'SubjectID', time: 'Time_hr', concentration: 'Conc_ng_mL', dose: 'Dose_mg_kg', bodyWeight: 'BW',
      groupBy: [], defaultRoute: 'IV-bolus',
      units: {concentration: 'ng/mL', time: 'h', dose: 'mg/kg', bodyWeight: bwUnit}}), '');
  }
  const rowOf = (res: DG.DataFrame, subj: string): number => rowsWhere(res, 'SUBJECT', subj)[0];

  test('AC-5.2 (E4): BW mapped in g gives the same per-kg CL as in kg — R001 5.1548 mL/min/kg', async () => {
    const g = await ratIv('g');
    const kg = await ratIv('kg');
    for (const p of ['CL_ML_MIN_KG', 'VZ_L_KG']) {
      const a = num(g, p, rowOf(g, 'R001')); const b = num(kg, p, rowOf(kg, 'R001'));
      expectFloat(a, b, Math.abs(b) * 1e-9, `${p}: g ≡ kg (E4 divided grams raw → 1000× low)`);
    }
    expectFloat(num(g, 'CL_ML_MIN_KG', rowOf(g, 'R001')), 5.1548, 5e-4, 'R001 per-kg CL (REGRESSION_PIN)');
  }, {timeout: TIMEOUT});

  test('AC-5.3 (E5): a mixed absolute / per-kg frame tags CL in mL/min/kg on every row', async () => {
    const res = await ratIv('g', 'R001'); // R001 has no BW → the per-kg engine path; the others absolute
    const norms = new Set(rowsWhere(res, 'VERDICT', VERDICT_ACCEPTED).concat(rowsWhere(res, 'VERDICT', VERDICT_REVIEW))
      .map((i) => s(res, 'DOSE_NORM', i)));
    expect(norms.size, 2, `both engine bases present: ${[...norms].join(', ')}`);
    expect(res.col('CL')!.getTag('units'), 'mL/min/kg', 'CL tag on the mixed frame (E5: was unset)');
    expect(res.col('VZ')!.getTag('units'), 'L/kg', 'VZ tag on the mixed frame');
    expectFloat(num(res, 'CL', rowOf(res, 'R001')), 5.1548, 5e-4, 'per-kg R001 on the per-kg engine path');
  }, {timeout: TIMEOUT});

  test('AC-5.1: dense rows equal the app results frame (CL / VZ / VSS / companions / BW_KG) to 1e-12', async () => {
    const df = await loadFixture('03_rat_iv.csv');
    const pipe = await ratIv('g');
    const app = runNcaPipeline(df, {observations: {
      subjectId: {name: 'SubjectID', required: true, rawType: 'string'},
      time: {name: 'Time_hr', required: true, rawType: 'float64', units: 'h'},
      concentration: {name: 'Conc_ng_mL', required: true, rawType: 'float64', units: 'ng/mL'},
      dose: {name: 'Dose_mg_kg', required: false, rawType: 'float64', units: 'mg/kg'},
      bodyWeight: {name: 'BodyWeight_g', required: false, rawType: 'float64'},
    }} as ColumnContract, {dose: 0, route: nca.ROUTE_IV_BOLUS, doseUnits: 'mg/kg', concentrationUnits: 'ng/mL',
      timeUnits: 'h', bodyWeight: null, bodyWeightUnit: 'g', infusionDuration: null, lloq: 0}, DEFAULT_RULES);
    expect(app.status, 'completed', 'app run');
    if (app.status !== 'completed') return;
    for (let i = 0; i < app.results.rowCount; i++) {
      const subj = String(app.results.col('SUBJECT')!.get(i));
      const j = rowOf(pipe, subj);
      for (const p of ['CL', 'VZ', 'VSS', 'CL_ABS', 'VZ_ABS', 'VSS_ABS', 'BW_KG']) {
        const a = app.results.col(p)!; const b = pipe.col(p)!;
        if (a.isNone(i)) {expect(b.isNone(j), true, `${subj} ${p}: both missing`); continue;}
        expectFloat(num(pipe, p, j), a.get(i) as number, Math.abs(a.get(i) as number) * 1e-12, `${subj} ${p}`);
      }
    }
  }, {timeout: TIMEOUT});

  test('F-2 / USER-CONFIRM-1: RESULTS_BASIS column + tag; explicit-unit aliases are deprecated', async () => {
    const res = await ratIv('g');
    expect(res.getTag('nca.resultsBasis'), RESULTS_BASIS, 'frame tag');
    for (let i = 0; i < res.rowCount; i++) expect(s(res, 'RESULTS_BASIS', i), RESULTS_BASIS, `row ${i} basis`);
    expect(res.col('CL_ML_MIN_KG')!.getTag('nca.deprecated'), 'use CL', 'alias marked deprecated');
    expect(res.col('VZ_L_KG')!.getTag('nca.deprecated'), 'use VZ', 'alias marked deprecated');
    for (let i = 0; i < res.rowCount; i++) {
      if (res.col('CL')!.isNone(i)) continue;
      expect(num(res, 'CL_ML_MIN_KG', i), num(res, 'CL', i), `alias = CL on row ${i}`);
    }
    // AC-5.5: the summary covers the per-kg family and the companions, and carries the basis.
    const sum = summarizePipeline(res);
    for (const p of ['MEAN_CL', 'MEAN_VZ', 'MEAN_VSS', 'MEAN_CL_ABS', 'MEAN_VZ_ABS', 'MEAN_VSS_ABS'])
      expect(sum.col(p) !== null, true, `summary has ${p}`);
    expect(sum.getTag('nca.resultsBasis'), RESULTS_BASIS, 'summary carries the basis');
  }, {timeout: TIMEOUT});
});

/**
 * The save, live (nca-pipeline-save-rework §1c; nca-ems-storage-migration S5) — on the EMS store,
 * whose reads are fresh: a write is visible on the next read (S0 AC-0.2).
 */
category('Pipeline functions: save to NCA Studio', () => {
  const repo = () => getNcaRepository();
  const stable = async (id: string): Promise<NcaStudy | undefined> => repo().loadStudy(id).catch(() => undefined);
  const idOf = (saved: DG.DataFrame): string => {
    const c = saved.col('StudyId')!;
    for (let i = 0; i < saved.rowCount; i++) if (c.get(i) !== '') return c.get(i);
    return '';
  };
  const cleanup = async (...ids: string[]): Promise<void> => {
    for (const id of ids) if (id !== '') await repo().deleteStudy(id).catch(() => {});
  };
  const statusOf = (saved: DG.DataFrame, analysis: string): string =>
    s(saved, 'Status', rowsWhere(saved, 'Analysis', analysis)[0]);

  test('AC-3.3 / AC-5 / Q1: one analysis per compound, new snapshots, provenance, stage, flags, clean data',
    async () => {
      const study = `TEST-SAVE-${Date.now()}`;
      const arms = ivPo(study, 'DGX-1003');
      arms[1].scale = {R2003: 0.04};
      const cassette: Arm[] = ['DGX-1001', 'DGX-1002'].map((c, i) => ({
        study, compound: c, route: 'PO' as const, dose: 2, subjects: ['R6001', 'R6002', 'R6003'], cl: 0.4 + 0.3 * i,
      }));
      const df = table([...arms, ...cassette]);
      const res = computePipeline(df, '', '');
      let id = '';
      try {
        const saved = await saveToStudies(df, res);
        id = idOf(saved);
        expect(saved.rowCount, 3, saved.toCsv());
        const p = await stable(id);
        expect(p !== undefined, true, 'study not found');
        expect(p!.label, study);
        expect(p!.analyses.map((a) => a.label).sort().join(','), 'DGX-1001,DGX-1002,DGX-1003');
        const a = p!.analyses.find((x) => x.label === 'DGX-1003')!;
        const runId = s(saved, 'RunId', 0);
        expect(a.origin?.kind === 'pipeline' && a.origin.runId === runId, true, 'run provenance');
        expect(a.source.kind, 'stored', 'a stored snapshot');
        expect(a.passthroughCols.length, 0, 'no pipeline passthrough columns (Q1)');
        expect(a.summaryCache !== null, true, 'summary cache');
        expect(a.stageCache !== null && a.stageCache !== undefined, true, 'stage (review 5)');
        expect(a.confirmed?.fingerprint !== undefined && a.confirmed.fingerprint !== '', true);
        // The overview carries the stored count; the flags themselves come with an open (ems-modification).
        expect((a.pipelineFlagCount ?? 0) > 0, true, 'the overview counts the pipeline flags (FIELD-11)');
        const flags = (await repo().loadAnalysis(a.id)).entry.pipelineFlags ?? [];
        expect(flags.some((f) => f.pipelineOnly.some((r) => r.includes('exposure outlier'))), true,
          'the outlier reason is the pipeline\'s own');
        // The overview loads the curves after its first read (ems-modification step 5).
        const loaded = withCurves(p!, await repo().loadCurves(pendingCurveIds(p!)));
        const curves = loaded.analyses.find((x) => x.id === a.id)!;
        expect((curves.summaryCache?.compoundTier ?? []).some((r) => r.curve !== null), true, 'curve (AC-5.4)');
        const data = (await repo().readSource(a.source)).toCsv();
        expect(data.includes('PIPELINE_') || data.includes('exposure outlier'), false,
          'the snapshot is the delivered data');
        expect(s(saved, 'Url', rowsWhere(saved, 'Analysis', 'DGX-1003')[0]), `/apps/Nca/${study}/DGX-1003`);
        expect(s(saved, 'Verdict', rowsWhere(saved, 'Analysis', 'DGX-1003')[0]), VERDICT_REVIEW);
        const me = (await grok.dapi.users.current()).id;
        expect((await repo().pipelineBase(me, study)).study?.id, id, 'found by its business key (owner, study key)');
      } finally {
        await cleanup(id);
      }
    }, {timeout: TIMEOUT});

  test('AC-4.1 / AC-3.6: a re-run updates an untouched analysis in place and removes the replaced snapshot',
    async () => {
      const study = `TEST-SAVE-${Date.now()}`;
      const df = table(ivPo(study, 'DGX-1001'));
      const res = computePipeline(df, '', '');
      let id = '';
      try {
        id = idOf(await saveToStudies(df, res));
        const first = (await stable(id))!.analyses[0];
        const again = await saveToStudies(df, res);
        expect(statusOf(again, 'DGX-1001'), 'updated');
        const p = (await stable(id))!;
        expect(p.analyses.length, 1);
        expect(p.analyses[0].id, first.id, 'same analysis');
        expect(p.analyses[0].source.ref !== first.source.ref, true, 'a new snapshot');
        expect(await grok.dapi.files.exists(first.source.ref), false, 'the replaced snapshot is removed');
      } finally {
        await cleanup(id);
      }
    }, {timeout: TIMEOUT});

  test('7d / AC-3.5: a re-run never overwrites review work — saved beside, the scientist\'s entry and file intact',
    async () => {
      const study = `TEST-SAVE-${Date.now()}`;
      const df = table(ivPo(study, 'DGX-1001'));
      const res = computePipeline(df, '', '');
      let id = '';
      try {
        id = idOf(await saveToStudies(df, res));
        const aid = (await stable(id))!.analyses[0].id;
        const {entry: a, df: stored} = await repo().loadAnalysis(aid);
        const lz: LzOverride = {profileKey: 'R1001', timeKey: '0.5', ordinal: 0, concKey: '1',
          action: 'force-out', user: 'test', at: new Date().toISOString()};
        // The scientist's λz edit, saved through the app's own Save.
        await repo().saveAnalysis(aid, a.source.ref, {df: stored, recipe: asRecipe(recipeOf(a)),
          summaryCache: a.summaryCache, workflow: {confirmed: a.confirmed ?? null, lzOverrides: [lz],
            reviewDecisions: a.reviewDecisions ?? [], stageCache: a.stageCache ?? null}});
        expect(((await repo().loadAnalysis(aid)).entry.lzOverrides ?? []).length, 1, 'the scientist edit reads back');
        const third = await saveToStudies(df, res);
        expect(statusOf(third, 'DGX-1001-2'), `saved beside 'DGX-1001' (has review work)`);
        const after = (await stable(id))!;
        expect(after.analyses.length, 2, after.analyses.map((x) => x.label).join(','));
        const kept = (await repo().loadAnalysis(aid)).entry;
        expect((kept.lzOverrides ?? []).length, 1, 'the override survives');
        expect(kept.source.ref, a.source.ref, 'its snapshot untouched');
        expect(computeSnapshotHash(await repo().readSource(kept.source)), computeSnapshotHash(stored),
          'its data untouched');
      } finally {
        await cleanup(id);
      }
    }, {timeout: TIMEOUT});

  test('AC-3.1b / PR-11: a same-label study elsewhere → our own study with a unique label; a foreign owner → refused',
    async () => {
      const study = `TEST-SAVE-${Date.now()}`;
      const df = table(ivPo(study, 'DGX-1001'));
      const res = computePipeline(df, '', '');
      let appId = '';
      let id = '';
      try {
        appId = (await repo().createStudy(study)).id;
        const saved = await saveToStudies(df, res);
        id = idOf(saved);
        expect(id !== appId, true, 'never the app-made study');
        expect(s(saved, 'Status', 0), `created study '${study}-2' (a study named '${study}' already exists)`);
        expect(s(saved, 'Url', 0), `/apps/Nca/${encodeURIComponent(`${study}-2`)}/DGX-1001`);
        expect(statusOf(await saveToStudies(df, res), 'DGX-1001'), 'updated', 'found again by its key');
        // PR-11: an owner that is neither the runner nor one of its groups → refused before any read.
        let msg = '';
        try {
          await saveToStudies(df, res, '', `not-me-${newId().slice(0, 6)}`);
        } catch (e) {
          msg = e instanceof Error ? e.message : String(e);
        }
        expect(msg.includes('is not you or one of your groups'), true, msg);
      } finally {
        await cleanup(id, appId);
      }
    }, {timeout: TIMEOUT});

  test('Q2: owner is part of the study identity; share grants a NEW study View; an unknown group throws',
    async () => {
      const study = `TEST-SAVE-${Date.now()}`;
      const df = table(ivPo(study, 'DGX-1001'));
      const res = computePipeline(df, '', '');
      const ids: string[] = [];
      try {
        // Two owners the runner may use: its id (the default) and its login.
        const byId = await saveToStudies(df, res);
        ids.push(idOf(byId));
        const byLogin = await saveToStudies(df, res, '', await myLogin());
        ids.push(idOf(byLogin));
        expect(idOf(byLogin) !== idOf(byId), true, 'another owner → another study');
        const all = (await grok.dapi.groups.list())
          .find((g) => g.friendlyName === 'All users' || g.name === 'All users');
        expect(all !== undefined, true, 'the stand has an "All users" group');
        const study2 = `TEST-SAVE-${Date.now()}-b`;
        const df2 = table(ivPo(study2, 'DGX-1001'));
        const shared = await saveToStudies(df2, computePipeline(df2, '', ''),
          JSON.stringify({view: ['All users'], edit: ['All users']}));
        ids.push(idOf(shared));
        // Decision 10 / SD-16: an Edit list is granted View, and the run says so.
        expect(s(shared, 'Status', 0).includes('share: edit not available yet (View granted)'), true,
          s(shared, 'Status', 0));
        const [entity] = await grok.dapi.getEntities([idOf(shared)]);
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        const perms = await grok.dapi.permissions.get(entity) as any;
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
        expect((perms.view as any[]).some((g) => g.id === all!.id), true, 'All users can view the new study');
        let msg = '';
        try {
          await saveToStudies(df, res, JSON.stringify({view: ['no-such-group-xyz']}));
        } catch (e) {
          msg = e instanceof Error ? e.message : String(e);
        }
        expect(msg, `share: unknown group 'no-such-group-xyz'`);
      } finally {
        await cleanup(...ids);
      }
    }, {timeout: TIMEOUT});

  test('F6: an IV + PO sparse partition → one analysis per route, each snapshot holds its route', async () => {
    const study = `TEST-SAVE-${Date.now()}`;
    const df = table([...sparseArms(study, 'IV', 1, 'A'), ...sparseArms(study, 'PO', 5, 'B')]);
    const res = computePipeline(df, '', '');
    let id = '';
    try {
      const saved = await saveToStudies(df, res);
      id = idOf(saved);
      const p = (await stable(id))!;
      const labels = p.analyses.map((a) => a.label).sort();
      expect(labels.length, 2, saved.toCsv());
      expect(labels.every((l) => l.startsWith('DGX-1008 · ')), true, labels.join(','));
      for (const a of p.analyses) {
        expect(a.modes.sparse.on && a.modes.sparse.samplingType === '(auto)', true, a.label);
        const data = await repo().readSource(a.source);
        expect(new Set(data.col('Route')!.toList()).size, 1, `${a.label}: one route`);
      }
    } finally {
      await cleanup(id);
    }
  }, {timeout: TIMEOUT});

  test('AC-3.7a: two overlapping runs on one study in one session both land', async () => {
    const study = `TEST-SAVE-${Date.now()}`;
    const df1 = table(ivPo(study, 'DGX-1001'));
    const df2 = table(ivPo(study, 'DGX-1002'));
    let id = '';
    try {
      const [a, b] = await Promise.all([
        saveToStudies(df1, computePipeline(df1, '', '')), saveToStudies(df2, computePipeline(df2, '', ''))]);
      id = idOf(a);
      expect(idOf(b), id);
      expect(s(a, 'Status', 0).startsWith('conflict') || s(b, 'Status', 0).startsWith('conflict'), false);
      expect((await stable(id))!.analyses.map((x) => x.label).sort().join(','), 'DGX-1001,DGX-1002');
    } finally {
      await cleanup(id);
    }
  }, {timeout: TIMEOUT});

  test('AC-5.6 / AC-5.7: the app opens a saved analysis on the saved stage (dense and sparse)', async () => {
    const study = `TEST-SAVE-${Date.now()}`;
    const dense = table(ivPo(study, 'DGX-1001'));
    const {df: mouse, mapping} = await mouseFixture();
    const ids: string[] = [];
    try {
      ids.push(idOf(await saveToStudies(dense, computePipeline(dense, '', ''))));
      ids.push(idOf(await saveToStudies(mouse, computePipeline(mouse, JSON.stringify(mapping), ''))));
      for (const id of ids) {
        const p = (await stable(id))!;
        for (const entry of p.analyses) {
          const view = await openAnalysis(entry.id) as DG.TableView;
          try {
            await DG.delay(800); // the replay settles, the deferred compute starts
            expect(currentRecipeFingerprint(view), entry.confirmed!.fingerprint,
              `${entry.label}: the replayed Setup fingerprints as saved (AC-5.7)`);
            const bag = view.temp as Record<string, unknown>;
            const t0 = Date.now();
            while (bag[RESULTS_GRID_TAG] == null && Date.now() - t0 < 30000) await DG.delay(200);
            const live = liveWorkflowState(view);
            expect(live.stage.stage, entry.stageCache!.stage, `${entry.label}: stage`);
            expect(live.stage.pending, entry.stageCache!.pending, `${entry.label}: pending`);
            const composites = new Set(live.units.filter((u) => u.kind === 'composite').map((u) => u.key));
            // The flags come with an open, not with the study (ems-modification).
            for (const f of (await repo().loadAnalysis(entry.id)).entry.pipelineFlags ?? []) {
              if (f.profileKey.startsWith('sparse:'))
                expect(composites.has(f.profileKey), true, `${entry.label}: flag key ${f.profileKey} is an app unit`);
            }
          } finally {
            view.close();
          }
        }
      }
    } finally {
      await cleanup(...ids);
    }
  }, {timeout: 240000});

  test('blocked partitions are reported, not saved', async () => {
    const study = `TEST-SAVE-${Date.now()}`;
    const df = table([{study, compound: 'C', route: 'IV', dose: 1, subjects: ['S1', 'S1']}]);
    const res = computePipeline(df, '', '');
    const saved = await saveToStudies(df, res);
    expect(saved.rowCount, 1);
    expect(s(saved, 'Verdict', 0), VERDICT_BLOCKED);
    expect(s(saved, 'Status', 0).startsWith('not saved'), true);
    const me = (await grok.dapi.users.current()).id;
    expect((await repo().pipelineBase(me, study)).study === undefined, true, 'no study was created');
  });
});

/** The app's source-replaced guard (nca-pipeline-save-rework F11, review 3a). */
category('Pipeline functions: source-replaced guard (app)', () => {
  const repo = () => getNcaRepository();
  const workflow = {confirmed: null, lzOverrides: [], reviewDecisions: [], stageCache: null};
  const CACHE = {analysisTier: {nCompounds: 1, nSubjects: 6, route: 'IV', dose: '',
    status: 'ok' as const, flagSummary: 'No flags'}, compoundTier: []};

  /** A one-analysis study / standalone container whose analysis another writer has re-pointed to a
   *  NEW snapshot (a pipeline re-run) after this view opened it. */
  async function replacedStudy(container: 'study' | 'standalone'): Promise<{
    id: string; binding: StudyBinding; df: DG.DataFrame; storedRef: string;
  }> {
    const df = table(ivPo('S', 'DGX-1001'));
    const recipe = recipeFor(parseMapping(manual()), '', false);
    const input = {label: 'DGX-1001', df, recipe, summaryCache: null, workflow, origin: null};
    let id: string;
    let aid: string;
    let openedRef: string;
    if (container === 'study') {
      id = (await repo().createStudy(`GUARD-${newId().slice(0, 6)}`)).id;
      const a = await repo().appendAnalysis(id, input);
      aid = a.analysisId;
      openedRef = a.source.ref;
    } else {
      const a = await repo().createStandalone({...input, label: `GUARD-${newId().slice(0, 6)}`});
      id = a.studyId;
      aid = a.analysisId;
      openedRef = a.source.ref;
    }
    await repo().loadAnalysis(aid); // the view opened it
    // Another session (a pipeline run) writes new data and points the analysis at it.
    const replaced = await repo().writeSnapshot(df, new Date().toISOString());
    await ncaDb.analyses.update(aid, {source_file: `file://${replaced.source.ref}`,
      ...signatureCols(replaced.signature)} as never);
    // …and removes the snapshot it replaced, as a run does (else the file outlives every row).
    await repo().deleteSnapshot({kind: 'stored', ref: openedRef});
    const binding: StudyBinding = {container, studyId: id, analysisId: aid, label: 'DGX-1001',
      source: {kind: 'stored', ref: openedRef}, sourceHash: computeSnapshotHash(df), passthroughCols: []};
    return {id, binding, df, storedRef: replaced.source.ref};
  }

  test('AC-11.1 / AC-11.3: the repository refuses a Save onto a replaced source; the writeback skips', async () => {
    const {id, binding, df} = await replacedStudy('study');
    try {
      let replaced = false;
      try {
        await repo().saveAnalysis(binding.analysisId, binding.source.ref, {df,
          recipe: recipeFor(parseMapping(manual()), '', false), summaryCache: null, workflow});
      } catch (e) {
        replaced = e instanceof SourceReplacedError;
      }
      expect(replaced, true, 'SourceReplacedError');
      await repo().writeSummary(binding.analysisId, binding.source.ref, {summaryCache: CACHE, workflow,
        frameHash: computeSnapshotHash(df)});
      const after = (await getNcaRepository().loadAnalysis(binding.analysisId)).entry;
      expect(after.summaryCache, null, 'the stale writeback did not land');
    } finally {
      await repo().deleteStudy(id).catch(() => {});
    }
  }, {timeout: TIMEOUT});

  test('AC-11.2: a SAVE on a replaced source saves beside, on a new snapshot, as the scientist\'s', async () => {
    const {id, binding, df, storedRef} = await replacedStudy('study');
    try {
      const recipe = recipeFor(parseMapping(manual()), '', false);
      const res = await updateAnalysisInPlace({df, recipe, computed: null, workflow, binding});
      expect(res.ok && res.rebound !== undefined, true, 'rebound');
      const p = await repo().loadStudy(id);
      expect(p.analyses.length, 2);
      const beside = p.analyses.find((a) => a.id === res.rebound!.analysisId)!;
      expect(beside.label, 'DGX-1001-2');
      expect(beside.origin, null, 'the scientist\'s — no pipeline run will update it');
      expect(beside.source.ref !== binding.source.ref && beside.source.ref !== storedRef, true, 'a new snapshot');
      expect(await grok.dapi.files.exists(beside.source.ref), true);
      expect(p.analyses.find((a) => a.id === binding.analysisId)!.source.ref, storedRef, 'the replaced entry kept');

      // AC-11.2a: the view is rebound — a second SAVE updates the new entry (no third one), and a
      // compute writeback through the new pointer lands on the new entry, never on the old one.
      const res2 = await updateAnalysisInPlace({df, recipe, computed: null, workflow, binding: res.rebound!});
      expect(res2.ok && res2.rebound === undefined, true, 'the second save is an in-place update');
      expect((await repo().loadStudy(id)).analyses.length, 2, 'still two entries');
      await repo().writeSummary(res.rebound!.analysisId, res.rebound!.source.ref, {summaryCache: CACHE, workflow,
        frameHash: computeSnapshotHash(df)});
      const p2 = await repo().loadStudy(id);
      expect(p2.analyses.find((a) => a.id === res.rebound!.analysisId)!.summaryCache !== null, true,
        'the writeback reached the new entry');
      expect(p2.analyses.find((a) => a.id === binding.analysisId)!.summaryCache, null, 'the old entry untouched');
    } finally {
      await repo().deleteStudy(id).catch(() => {});
    }
  }, {timeout: TIMEOUT});

  test('AC-11.4: a standalone analysis is never split — the save is refused', async () => {
    const {id, binding, df} = await replacedStudy('standalone');
    try {
      const res = await updateAnalysisInPlace({df, recipe: recipeFor(parseMapping(manual()), '', false),
        computed: null, workflow, binding});
      expect(res.ok, false);
      expect((await repo().listStudies()).find((s) => s.id === id)?.analysisCount, 1, 'still one analysis');
    } finally {
      await repo().deleteStudy(id).catch(() => {});
    }
  }, {timeout: TIMEOUT});
});

// ---------- The save protocol on an in-memory store (AC-3.4d) ----------

/**
 * The save protocol against another writer, staged deterministically on the in-memory repository
 * (decision 7): its hook runs right before each plan's write — between the base read and the
 * transaction, where a scientist's save or another run would land. (The pre-EMS cases "a store that
 * never reads stable" and "a save another writer overwrites" are gone with the store that had them:
 * reads are fresh and the plan's write is one versioned transaction.)
 */
category('Pipeline functions: save protocol (in-memory store)', () => {
  /** A first, clean run: the study exists with one pipeline analysis and one snapshot. */
  const seeded = async (): Promise<{mem: MemoryNcaRepository; df: DG.DataFrame; res: DG.DataFrame;
    owner: string; studyId: string}> => {
    const study = `TEST-MEM-${Date.now()}`;
    const df = table(ivPo(study, 'DGX-1001'));
    const res = computePipeline(df, '', '');
    const mem = new MemoryNcaRepository();
    const owner = await myLogin();
    const first = await saveToStudies(df, res, '', owner, {repo: mem});
    expect(s(first, 'Status', 0), 'created');
    return {mem, df, res, owner, studyId: s(first, 'StudyId', 0)};
  };

  test('AC-3.4d: a writer between the base read and the write → re-planned, the writer\'s change kept', async () => {
    const {mem, df, res, owner, studyId} = await seeded();
    const firstWrite = mem.pipelineWrites;
    mem.beforePipelineWrite = (n) => {
      // Another writer adds an analysis to the study right before this run's first write.
      if (n === firstWrite + 1) mem.addForeignAnalysis(studyId, 'Hand-made');
    };
    const out = await saveToStudies(df, res, '', owner, {repo: mem});
    expect(s(out, 'Status', 0), 'updated', out.toCsv());
    const p = (await mem.pipelineBase(owner, s(out, 'Study', 0))).study!;
    expect(p.analyses.some((a) => a.label === 'Hand-made'), true, 'the other writer\'s analysis survives');
    expect(p.analyses.length, 2);
  }, {timeout: TIMEOUT});

  /** Save a pipeline run (optionally with a `.`-marker annotation column), open it in the app, and
   *  check the open is clean — not dirty — and its summary write-back lands. */
  const opensClean = async (marked: boolean): Promise<void> => {
    const df = table(ivPo(`TEST-MEM-${Date.now()}`, 'DGX-1001'));
    if (marked) {
      df.columns.add(DG.Column.fromStrings('Batch',
        Array.from({length: df.rowCount}, (_, i) => (i % 3 === 0 ? '.' : String(100 + i)))));
    }
    const mem = new MemoryNcaRepository();
    const owner = await myLogin();
    const out = await saveToStudies(df, computePipeline(df, '', ''), '', owner, {repo: mem});
    expect(s(out, 'Status', 0), 'created', out.toCsv());
    const entry = (await mem.pipelineBase(owner, s(out, 'Study', 0))).study!.analyses[0];
    const stored = await mem.readSource(entry.source);
    if (marked)
      expect(stored.col('Batch')!.type, DG.COLUMN_TYPE.STRING, 'the pipeline stores the frame un-normalized');
    const prev = setNcaRepository(mem);
    try {
      const before = mem.analyses.get(entry.id)!.version;
      const view = await openAnalysis(entry.id) as DG.TableView;
      try {
        const bag = view.temp as Record<string, unknown>;
        const t0 = Date.now();
        while (bag[RESULTS_GRID_TAG] == null && Date.now() - t0 < 30000) await DG.delay(200);
        await DG.delay(500); // the write-back is async
        if (marked) expect(view.dataFrame.col('Batch')!.type, DG.COLUMN_TYPE.FLOAT, 'the open normalized the markers');
        expect(computeSnapshotHash(view.dataFrame) === entry.sourceSignature.hash, !marked,
          'the derived ~RouteCode the compute added does not move the hash');
        expect(isAnalysisDirty(view), false, 'the open is not an edit');
        expect(mem.analyses.get(entry.id)!.version > before, true, 'the summary write-back landed');
      } finally {
        view.close();
      }
    } finally {
      setNcaRepository(prev);
    }
  };

  // Every open of a pipeline-saved analysis was marked dirty (the recipe replay fired onChanged
  // before the replay guard was up), so its write-back never ran.
  test('a pipeline-saved analysis opens clean and its write-back lands', async () => opensClean(false),
    {timeout: TIMEOUT});

  // The pipeline stores its frame un-normalized; the open retypes a `.`-marker column. That frame is
  // still the stored data.
  test('a pipeline snapshot with a `.`-marker column opens clean and its write-back lands', async () =>
    opensClean(true), {timeout: TIMEOUT});

  test('AC-3.4d: a base that moves on every attempt → nothing written, this run\'s snapshots deleted', async () => {
    const {mem, df, res, owner} = await seeded();
    const snapshotsBefore = mem.snapshots.size;
    const study = s(await saveToStudies(df, res, '', owner, {repo: mem}), 'Study', 0); // a clean re-run
    const base = (await mem.pipelineBase(owner, study)).study!;
    const before = JSON.stringify(base.analyses);
    const snapshots = mem.snapshots.size;
    // Every write attempt finds the study moved (another writer touches it each time).
    mem.beforePipelineWrite = () => {void mem.renameStudy(base.id, base.label);};
    const out = await saveToStudies(df, res, '', owner, {repo: mem});
    expect(s(out, 'Status', 0), `not saved: study '${base.label}' kept changing during the save — re-run`, out.toCsv());
    expect(JSON.stringify((await mem.pipelineBase(owner, study)).study!.analyses), before, 'nothing written');
    expect(mem.snapshots.size, snapshots, 'this run\'s snapshots were removed');
    expect(snapshotsBefore > 0, true);
  }, {timeout: TIMEOUT});
});
