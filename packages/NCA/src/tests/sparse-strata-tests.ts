/**
 * BUG-19 — one sparse composite compute for the app and the pipeline (`computeSparseStrata`).
 *
 * Layers:
 *  - **Canonicalization** (F1): every unit converted on synthetic frames with hard-coded values —
 *    time, conc, BW, per-kg vs absolute dose, LLOQ; molar and text-conc refusals; route strata.
 *  - **Real fixture** (rule 18): the mouse sparse fixture's nominal time DERIVED into min / s / d
 *    and run through the app's entry point (`computeSparseForRecipe`) equals the hours run, and the
 *    PKNCA golden AUClast / SE / df. λz / t½ / CL / MRT are identity-covered only
 *    (DATA-GAP-PIPE-SPARSE-LZ: no PKNCA sparse λz reference).
 *  - **App ↔ pipeline**: both paths agree, key and BLQ and missing-concentration rules included.
 *  - **Decisions**: a non-hour study's composite fingerprint changes (the disclosed reset).
 */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect, expectFloat, before} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {
  computeSparseStrata, computeSparseForRecipe, computeSparseNca, sparseCompositeUnits, sparseCompositeKey,
  assembleSparseResult,
} from '../adapter/sparse-nca';
import {
  normalizeEventDrivenObservations, sparseDoseSource, SUBJECT_COL, TAD_COL, CONC_COL, DOSE_COL,
  type EventContract,
} from '../adapter/event-driven';
import {CONC_TO_MG_L, TIME_TO_HOUR, BW_TO_KG, DOSE_TO_MG} from '../adapter/unit-conversion';
import type {
  SparseUnits, SparseStrataSpec, SparseComputeMapping, SparseGroupResult, SparseStrata, SparseStratumCell,
} from '../adapter/sparse-nca';
import {DEFAULT_RULES, DEFAULT_ACCEPTANCE_RSQUARED} from '../adapter/defaults';
import {ensureBoolBlq} from '../adapter/adapter-utils';
import {DERIVED_ROUTE_CODE_COL} from '../adapter/route-normalize';
import {sparseUnitKey} from '../adapter/review-units';
import {getNcaRepository, recipeOf} from '../adapter';
import {NcaController} from '../controller/nca-controller';
import {computePipeline, parseRules, COMPOSITE_SUBJECT} from '../pipeline/pipeline-compute';
import {detectPipelineMapping, parseMapping} from '../pipeline/pipeline-mapping';
import {pipelineFlagsOf, saveToStudies} from '../pipeline/pipeline-save';
import {openAnalysis} from '../views/nca-analysis-open';
import {liveWorkflowState, RESULTS_GRID_TAG} from '../views/nca-mapping-pane';
import {loadFixture, makeSyntheticObs, singleSubjectContract, singleSubjectScalars} from './test-utils';

const TIMEOUT = 120000;
const FIXTURE_PATH_PREFIX = 'System:AppData/Nca';
const ACC = {
  lambdaAcceptanceRSquared: DEFAULT_ACCEPTANCE_RSQUARED,
  lambdaMinPoints: DEFAULT_RULES.lambdaZ.minPoints,
  extrapWarnPct: DEFAULT_RULES.extrapWarnPct,
};
const H: SparseUnits = {time: 'h', concentration: 'mg/L', dose: 'mg/kg', bodyWeight: 'kg'};
/** The destructive TK NONMEM fixture's event contract (as in workflow-stages-controller-tests). */
const EC06: EventContract = {
  id: 'ID', time: 'TIME', dv: 'DV', amt: 'AMT', evid: 'EVID', cmt: 'CMT', rate: 'RATE', mdv: 'MDV',
  addl: 'ADDL', ii: 'II', route: 'ROUTE', bodyWeight: 'BW_kg', covariates: ['Sex', 'DOSE_MGKG'], day: 'DAY',
};

// ─── Synthetic destructive frames ─────────────────────────────────────────────

const T_SPARSE = [0.25, 0.5, 1, 2, 4, 6, 8, 24];

/** One-compartment concentration (mg/L) per mg/kg: extravascular unless `iv`. */
function pk(t: number, dose: number, iv: boolean): number {
  const cl = 0.3; const v = 0.3; const k = cl / v; const ka = 0.6; const f = 0.12;
  return iv ? dose / v * Math.exp(-k * t) : f * dose * ka / (v * (ka - k)) * (Math.exp(-k * t) - Math.exp(-ka * t));
}

interface Arm {
  dose: number | string; route?: number | null; prefix: string; sex?: string; comment?: string;
  /** The mg/kg the curve is generated with (default: the dose value) — an absolute-dose arm keeps
   *  the per-kg arm's concentrations. */
  curveDose?: number;
  /** Per-animal concentration overrides (null = a missing cell). */
  conc?: Record<string, number | null>;
  /** Per-animal BLQ flag overrides (bool, or a text token for a text column). */
  blq?: Record<string, boolean | string>;
  lloq?: number;
}

interface FrameOpts {
  /** Multiply the nominal time / conc / BW / dose by these factors (a unit change of the data). */
  timeF?: number; concF?: number; bwF?: number; doseF?: number;
  /** Text BLQ column ('Y' / 'N' unless an arm overrides it) instead of bool. */
  textBlq?: boolean;
  /** Dose as a STRING column. */
  textDose?: boolean;
}

/** A destructive frame: 3 animals per nominal time per arm, one sample each, BW 0.25 kg. */
function frame(arms: Arm[], o: FrameOpts = {}): DG.DataFrame {
  const animal: string[] = []; const time: number[] = []; const conc: number[] = [];
  const blqB: boolean[] = []; const blqT: string[] = []; const dose: Array<number | string> = [];
  const route: number[] = []; const bw: number[] = []; const lloq: number[] = []; const sex: string[] = [];
  const comment: string[] = [];
  for (const a of arms) {
    T_SPARSE.forEach((t, i) => {
      for (let k = 0; k < 3; k++) {
        const id = `${a.prefix}${String(3 * i + k + 1).padStart(3, '0')}`;
        // A text dose label ("High") has no number: its curve uses 1 mg/kg.
        const d = a.curveDose ?? (Number.isFinite(Number(a.dose)) ? Number(a.dose) : 1);
        let c: number | null = pk(t, d, a.route === nca.ROUTE_IV_BOLUS);
        if (a.conc !== undefined && id in a.conc) c = a.conc[id];
        const flag = a.blq?.[id];
        animal.push(id); time.push(t * (o.timeF ?? 1));
        conc.push(c === null ? DG.FLOAT_NULL : c * (o.concF ?? 1));
        blqB.push(flag === true); blqT.push(typeof flag === 'string' ? flag : (flag === true ? 'Y' : 'N'));
        dose.push(typeof a.dose === 'string' ? a.dose : a.dose * (o.doseF ?? 1));
        route.push(a.route === undefined ? nca.ROUTE_PO : (a.route === null ? DG.INT_NULL : a.route));
        bw.push(0.25 * (o.bwF ?? 1)); lloq.push(a.lloq ?? 0.001); sex.push(a.sex ?? 'M');
        comment.push(a.comment ?? '');
      }
    });
  }
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Animal', animal),
    DG.Column.fromFloat64Array('Time', Float64Array.from(time)),
    DG.Column.fromFloat64Array('Conc', Float64Array.from(conc)),
    // Explicit STRING columns: `fromStrings` infers a numeric type from numeric-looking text.
    o.textBlq ? DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'BLQ', blqT) :
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'BLQ', blqB),
    // A dose column as a CSV load types it: integers → INT (getString '10'), else FLOAT.
    o.textDose ? DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Dose', dose.map(String)) :
      (dose as number[]).every((v) => Number.isInteger(v)) ?
        DG.Column.fromInt32Array('Dose', Int32Array.from(dose as number[])) :
        DG.Column.fromFloat64Array('Dose', Float64Array.from(dose as number[])),
    DG.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, Int32Array.from(route)),
    DG.Column.fromStrings('Route', route.map((r) => (r === DG.INT_NULL ? '' : nca.ROUTE_IV_BOLUS === r ? 'IV-bolus' :
      r === nca.ROUTE_IV_INFUSION ? 'IV-infusion' : r === nca.ROUTE_SC ? 'SC' : 'PO'))),
    DG.Column.fromFloat64Array('BW', Float64Array.from(bw)),
    DG.Column.fromFloat64Array('LLOQ', Float64Array.from(lloq)),
    DG.Column.fromStrings('Sex', sex),
    DG.Column.fromStrings('Comment', comment),
  ]);
}

const spec = (units: Partial<SparseUnits> = {}, extra: Partial<SparseStrataSpec> = {}): SparseStrataSpec => ({
  concCol: 'Conc', nominalTimeCol: 'Time', animalIdCol: 'Animal', blqCol: 'BLQ', doseCol: 'Dose', bodyWeightCol: 'BW',
  lloq: {colName: null, scalar: 0}, route: {codesCol: DERIVED_ROUTE_CODE_COL, scalar: 'PO'},
  units: {...H, ...units}, acceptance: ACC, bootstrapIterations: 50, ...extra,
});

function cellsOf(s: SparseStrata): SparseStratumCell[] {
  if (!s.ok) throw new Error(`strata refused: ${s.reason}`);
  return s.cells;
}
const groupsOf = (s: SparseStrata): SparseGroupResult[] => cellsOf(s).map((c) => c.group);

/** The composite parameters a unit change must leave unchanged once canonicalized. */
const PARAMS = ['aucLast', 'aucSe', 'aucDf', 'cmax', 'tmax', 'aucInf', 'lambdaZ', 'halfLife', 'cl', 'vz', 'mrt',
  'vss', 'pctExtrap', 'lambdaZAdjR2', 'spanRatio', 'lambdaZNPoints'] as const;

function expectSame(a: number, b: number, tag: string, rel = 1e-9): void {
  if (Number.isNaN(b)) {
    expect(Number.isNaN(a), true, `${tag}: NaN in the reference, ${a} here`);
    return;
  }
  expectFloat(a, b, Math.abs(b) * rel + 1e-12, tag);
}

function expectSameGroup(a: SparseGroupResult, b: SparseGroupResult, tag: string): void {
  for (const p of PARAMS) expectSame(a[p], b[p], `${tag} ${p}`);
  for (const k of ['bootstrapAucExploratory', 'bootstrapCmax', 'bootstrapHalfLife', 'bootstrapCl'] as const) {
    expectSame(a[k].ciLow, b[k].ciLow, `${tag} ${k} low`);
    expectSame(a[k].ciHigh, b[k].ciHigh, `${tag} ${k} high`);
  }
}

// ─── The mouse fixture ────────────────────────────────────────────────────────

interface SparseGolden {
  destructive: Array<{dose_mg_kg: number; sparse_auclast: number; sparse_auc_se: number; sparse_auc_df: number}>;
}

/** The app's sparse mapping of the mouse fixture with `timeCol` declared in `time`. */
const mouseMapping = (time: nca.TimeUnit, timeCol = 'NominalTime_hr'): SparseComputeMapping => ({
  subjectName: 'MouseID', timeName: timeCol, concName: 'Conc_ng_mL', doseColName: 'Dose_mg_kg',
  bwColName: 'BodyWeight_g', units: {time, concentration: 'ng/mL', dose: 'mg/kg', bodyWeight: 'g'},
  route: {codesCol: null, scalar: 'PO'}, lloq: {colName: 'LLOQ_ng_mL', scalar: 0}, blqColName: 'BLQ',
  nominalTimeName: timeCol, animalIdName: 'MouseID', declaredTopology: null, acceptance: ACC,
  bootstrapIterations: 200,
});

/** A copy of the hours column `src` × `f` named `name` (a derived fixture: never hand-typed values). */
function derive(df: DG.DataFrame, src: string, name: string, f: number): DG.DataFrame {
  const c = df.col(src)!;
  df.columns.add(DG.Column.fromFloat64Array(name, Float64Array.from({length: df.rowCount}, (_, i) => c.get(i) * f)));
  return df;
}

function appGroups(df: DG.DataFrame, m: SparseComputeMapping): SparseGroupResult[] {
  const r = computeSparseForRecipe(df, m);
  if (r.result === null) throw new Error(`no sparse result: ${r.error}`);
  return r.result.groups;
}

const compositeRows = (res: DG.DataFrame): number[] =>
  Array.from({length: res.rowCount}, (_, i) => i).filter((i) => res.col('SUBJECT')!.get(i) === COMPOSITE_SUBJECT);
const num = (df: DG.DataFrame, col: string, i: number): number =>
  df.col(col)!.isNone(i) ? NaN : df.col(col)!.get(i) as number;
const str = (df: DG.DataFrame, col: string, i: number): string => String(df.col(col)!.get(i) ?? '');

/** The mouse fixture's pipeline mapping (auto-detected, plus the BLQ column). */
function mousePipelineMapping(df: DG.DataFrame): Record<string, unknown> {
  const det = detectPipelineMapping(df);
  expect(det.confidence, 'high', det.issues.join('; '));
  return {...det.mapping!, blq: 'BLQ'} as unknown as Record<string, unknown>;
}

/** The mouse fixture with its observed AND nominal time columns rescaled to minutes in place. */
function toMinutes(df: DG.DataFrame, m: Record<string, unknown>): Record<string, unknown> {
  for (const c of new Set([m.time as string, m.nominalTime as string])) {
    const col = df.col(c)!;
    const v = Float64Array.from({length: df.rowCount}, (_, i) => col.get(i) * 60);
    df.columns.remove(c);
    df.columns.add(DG.Column.fromFloat64Array(c, v));
  }
  return {...m, units: {...(m.units as object), time: 'min'}};
}

/** The pipeline mapping of a synthetic {@link frame} (every role a column, mg/kg / h / mg/L / kg). */
const synthMapping = (extra: Record<string, unknown> = {}): string => JSON.stringify({
  subject: 'Animal', time: 'Time', concentration: 'Conc', dose: 'Dose', route: 'Route', bodyWeight: 'BW',
  blq: 'BLQ', lloq: 'LLOQ', comment: 'Comment', sex: 'Sex', groupBy: [], design: 'sparse-destructive',
  units: {concentration: 'mg/L', time: 'h', dose: 'mg/kg', bodyWeight: 'kg'}, ...extra,
});

category('Adapter: sparse strata (BUG-19)', () => {
  let golden: SparseGolden;
  before(async () => {
    golden = JSON.parse(await grok.dapi.files.readAsText(`${FIXTURE_PATH_PREFIX}/mouse_sparse-reference.json`));
  });

  // ─── F1: canonicalization ───────────────────────────────────────────────────

  test('AC-1.1 / AC-1.3: min, ng/mL, g, µg/kg read as h, mg/L, kg, mg/kg (hard-coded values)', async () => {
    const arm: Arm = {dose: 10, prefix: 'A'};
    const h = groupsOf(computeSparseStrata(frame([arm]), spec()));
    const x = groupsOf(computeSparseStrata(frame([arm], {timeF: 60, concF: 1000, bwF: 1000, doseF: 1000}),
      spec({time: 'min', concentration: 'ng/mL', bodyWeight: 'g', dose: 'µg/kg'})));
    expect(h.length, 1); expect(x.length, 1);
    expectSameGroup(x[0], h[0], 'min/ng/g/µg vs h/mg/kg/mg');
    expectFloat(x[0].rawPoints[0].nominalTime, 0.25, 1e-12, 'time in h (15 min)');
    expectFloat(x[0].rawPoints[23].nominalTime, 24, 1e-12, 'time in h (1440 min)');
    expectFloat(x[0].rawPoints[0].conc, h[0].rawPoints[0].conc, 1e-12, 'conc in mg/L');
    expectFloat(x[0].meanBodyWeightKg!, 0.25, 1e-12, 'BW in kg');
    expectFloat(x[0].doseMgKg, 10, 1e-12, 'dose in mg/kg (canonical compute dose)');
    expectFloat(x[0].doseMg, 2.5, 1e-12, 'mg dose = 10 mg/kg × 0.25 kg');
    expect(x[0].doseLevel.key, '10000', 'the key keeps the raw µg/kg spelling');
    expect(x[0].groupKey, 'DOSE=10000');
    expect(cellsOf(computeSparseStrata(frame([arm]), spec())).length, 1);
  }, {timeout: TIMEOUT});

  test('AC-1.1: every unit of the four tables — data in that unit equals the canonical run', async () => {
    const arm: Arm = {dose: 10, prefix: 'A'};
    const base = groupsOf(computeSparseStrata(frame([arm]), spec()))[0];
    const check = (o: FrameOpts, u: Partial<SparseUnits>, tag: string): void => {
      const g = groupsOf(computeSparseStrata(frame([arm], o), spec(u)))[0];
      expectSameGroup(g, base, tag);
      expectFloat(g.rawPoints[0].nominalTime, 0.25, 1e-12, `${tag}: time in h`);
    };
    for (const [u, f] of Object.entries(TIME_TO_HOUR))
      check({timeF: 1 / f}, {time: u as nca.TimeUnit}, `time ${u}`);
    for (const [u, f] of Object.entries(CONC_TO_MG_L)) {
      if (!Number.isFinite(f)) continue; // molar — refused (AC-1.2)
      check({concF: 1 / f}, {concentration: u as nca.ConcentrationUnit}, `conc ${u}`);
    }
    for (const [u, f] of Object.entries(BW_TO_KG)) check({bwF: 1 / f}, {bodyWeight: u as 'kg' | 'g'}, `BW ${u}`);
    for (const u of ['mg/kg', 'µg/kg'] as const) check({doseF: 1 / DOSE_TO_MG[u]}, {dose: u}, `dose ${u}`);
    // Absolute units: 2.5 mg ≡ 10 mg/kg × 0.25 kg, each spelled in its own unit.
    for (const u of ['mg', 'µg', 'ng'] as const) {
      const g = groupsOf(computeSparseStrata(frame([{dose: 2.5 / DOSE_TO_MG[u], curveDose: 10, prefix: 'A'}]),
        spec({dose: u})))[0];
      expectFloat(g.doseMg, 2.5, 1e-9, `dose ${u} → mg`);
      expectSame(g.cl, base.cl, `dose ${u}: CL`);
    }
  }, {timeout: TIMEOUT});

  test('AD-8: ensureBoolBlq reads bool, 0/1 and the text tokens (the moved truthy set)', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'F',
      ['1', 'y', 'YES', 'true', 'BLQ', 'bloq', '<LLOQ', ' Y ', '0', 'N', 'false', '', 'ND'])]);
    const name = ensureBoolBlq(df, 'F');
    expect(name, '~BLQ');
    expect(Array.from({length: df.rowCount}, (_, i) => df.col(name)!.get(i) ? 1 : 0).join(''), '1111111100000');
    const bool = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'B', [true, false])]);
    expect(ensureBoolBlq(bool, 'B'), 'B', 'a bool column is read as is');
  });

  test('F4 (IV): an IV-bolus composite — Vss, MRT, CL identical in min and h (R2-N1)', async () => {
    // indometh is a serial study, so the IV-bolus composite is a synthetic destructive arm.
    const arm: Arm = {dose: 1, prefix: 'I', route: nca.ROUTE_IV_BOLUS};
    const h = groupsOf(computeSparseStrata(frame([arm]), spec()))[0];
    const min = groupsOf(computeSparseStrata(frame([arm], {timeF: 60}), spec({time: 'min'})))[0];
    expect(h.route, 'IV-bolus');
    expect(Number.isFinite(h.vss) && h.vss > 0, true, `Vss reported for IV-bolus (${h.vss})`);
    expectSameGroup(min, h, 'IV-bolus min vs h');
  }, {timeout: TIMEOUT});

  test('AC-1.3 / F3: an absolute mg or µg dose → CL on that dose, no BW factor; µg/kg → CL ÷ 1000', async () => {
    const perKg = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A'}]), spec()))[0];
    // 2.5 mg absolute ≡ 10 mg/kg × 0.25 kg — the same CL, without the group-mean BW convention.
    const mg = groupsOf(computeSparseStrata(frame([{dose: 2.5, curveDose: 10, prefix: 'A'}]), spec({dose: 'mg'})))[0];
    const ug = groupsOf(computeSparseStrata(frame([{dose: 2500, curveDose: 10, prefix: 'A'}]),
      spec({dose: 'µg'})))[0];
    expect(Number.isNaN(mg.doseMgKg), true, 'absolute: no per-kg compute dose');
    expectFloat(mg.doseMg, 2.5, 1e-12, 'mg dose'); expectFloat(ug.doseMg, 2.5, 1e-12, 'µg → mg');
    expectSame(mg.cl, perKg.cl, 'CL mg vs mg/kg × BW');
    expectSame(ug.cl, perKg.cl, 'CL µg vs mg/kg × BW');
    expect(mg.warnings.some((w) => w.code === 'SPARSE_COMPOSITE_BW'), false, 'no composite-BW note on an mg dose');
    const ugkg = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A'}]), spec({dose: 'µg/kg'})))[0];
    expectSame(ugkg.cl, perKg.cl / 1000, 'µg/kg → CL ÷ 1000');
    expectSame(ugkg.aucLast, perKg.aucLast, 'AUC unaffected by the dose unit');
    // Same on the pipeline path (S-3).
    const df = frame([{dose: 10, prefix: 'A'}]);
    const pMgKg = computePipeline(df, synthMapping(), '');
    const pUgKg = computePipeline(df, synthMapping({units: {concentration: 'mg/L', time: 'h', dose: 'µg/kg',
      bodyWeight: 'kg'}}), '');
    expectSame(num(pUgKg, 'CL', compositeRows(pUgKg)[0]), num(pMgKg, 'CL', compositeRows(pMgKg)[0]) / 1000,
      'pipeline µg/kg → CL ÷ 1000');
    const dfMg = frame([{dose: 2.5, curveDose: 10, prefix: 'A'}]);
    const pMg = computePipeline(dfMg, synthMapping({units: {concentration: 'mg/L', time: 'h', dose: 'mg',
      bodyWeight: 'kg'}}), '');
    expectSame(num(pMg, 'CL', compositeRows(pMg)[0]), num(pMgKg, 'CL', compositeRows(pMgKg)[0]),
      'pipeline mg vs mg/kg × BW');
  }, {timeout: TIMEOUT});

  test('AC-1.2 / AC-7.5: a molar unit and a text concentration column are refused', async () => {
    const molar = computeSparseStrata(frame([{dose: 10, prefix: 'A'}]), spec({concentration: 'nmol/L'}));
    expect(molar.ok ? '' : molar.reason, `sparse: molar concentration unit 'nmol/L' is not supported`);
    const df = frame([{dose: 10, prefix: 'A'}]);
    const c = df.col('Conc')!;
    df.columns.remove('Conc');
    df.columns.add(DG.Column.fromStrings('Conc', Array.from({length: df.rowCount},
      (_, i) => (i === 0 ? 'BLQ' : i === 1 ? 'ND' : String(c.get(i))))));
    const text = computeSparseStrata(df, spec());
    expect(text.ok ? '' : text.reason, `concentration column 'Conc' is not numeric: "BLQ", "ND"`);
    // The pipeline blocks the partition with that reason (R2-N4).
    const res = computePipeline(df, synthMapping(), '');
    expect(str(res, 'VERDICT', 0), 'Blocked');
    expect(str(res, 'REASONS', 0).includes(`is not numeric: "BLQ", "ND"`), true, str(res, 'REASONS', 0));
  }, {timeout: TIMEOUT});

  test('AC-1.4 / AD-9: LLOQ = the largest value × the conc factor; a varying column is noted, not fingerprinted',
    async () => {
      const df = frame([{dose: 10, prefix: 'A'}]);
      const constant = groupsOf(computeSparseStrata(df, spec({}, {lloq: {colName: 'LLOQ', scalar: 0}})));
      expect(constant[0].warnings.some((w) => w.code === 'SPARSE_LLOQ_VARIES'), false, 'one LLOQ value → no note');
      const lloq = df.col('LLOQ')!;
      lloq.set(0, 0.005); lloq.set(1, 0.002);
      const varying = groupsOf(computeSparseStrata(df, spec({concentration: 'µg/L'},
        {lloq: {colName: 'LLOQ', scalar: 0}})));
      const note = varying[0].warnings.find((w) => w.code === 'SPARSE_LLOQ_VARIES');
      expect(note?.message ?? '', 'The LLOQ column holds 3 distinct values; the composite uses the largest, ' +
        '0.000005 mg/L.');
      // The note changes no number, so it must not reset a decision (AD-9).
      const fp = (g: SparseGroupResult[]): string =>
        sparseCompositeUnits({groups: g, parameterTable: null!, compositeTable: null!, warnings: []})
          .map((u) => u.fingerprint).join();
      const noLloq = groupsOf(computeSparseStrata(df, spec({concentration: 'µg/L'})));
      expect(fp(varying), fp(noLloq), 'a varying LLOQ leaves the fingerprint byte-identical');
    }, {timeout: TIMEOUT});

  test('AC-1.5 / AD-3: routes — unprefixed for one route, prefixed for two; blank = the scalar; ' +
    'PO + SC two strata; infusion blocked per dose', async () => {
    const one = cellsOf(computeSparseStrata(frame([{dose: 1, prefix: 'A'}, {dose: 5, prefix: 'B'}]), spec()));
    expect(one.map((c) => c.group.groupKey).join(','), 'DOSE=1,DOSE=5', 'one route → unprefixed');
    const two = cellsOf(computeSparseStrata(frame([{dose: 1, prefix: 'I', route: nca.ROUTE_IV_BOLUS},
      {dose: 5, prefix: 'P'}]), spec()));
    expect(two.map((c) => c.group.groupKey).join(','), 'ROUTE=IV-bolus|DOSE=1,ROUTE=PO|DOSE=5');
    expect(two.map((c) => c.group.route).join(','), 'IV-bolus,PO');
    expect(two[0].rows.length + two[1].rows.length, 48, 'rows index the input frame, every row in one cell');
    // A blank route cell takes the scalar.
    const blank = cellsOf(computeSparseStrata(frame([{dose: 1, prefix: 'A', route: null}]),
      spec({}, {route: {codesCol: DERIVED_ROUTE_CODE_COL, scalar: 'IV-bolus'}})));
    expect(blank.map((c) => c.group.groupKey).join(','), 'DOSE=1', 'blank → the scalar route, one route');
    expect(blank[0].group.route, 'IV-bolus');
    const poSc = cellsOf(computeSparseStrata(frame([{dose: 1, prefix: 'P'}, {dose: 1, prefix: 'S',
      route: nca.ROUTE_SC}]), spec()));
    expect(poSc.map((c) => c.group.groupKey).join(','), 'ROUTE=PO|DOSE=1,ROUTE=SC|DOSE=1', 'PO and SC split');
    // IV-bolus + IV-infusion: the infusion arm is blocked per dose group; the bolus keys still carry
    // the route, so a later infusion fix cannot move them (R2-N5).
    const inf = computeSparseStrata(frame([{dose: 1, prefix: 'B', route: nca.ROUTE_IV_BOLUS},
      {dose: 1, prefix: 'I', route: nca.ROUTE_IV_INFUSION}, {dose: 3, prefix: 'J', route: nca.ROUTE_IV_INFUSION}]),
    spec());
    if (!inf.ok) throw new Error(inf.reason);
    expect(inf.cells.map((c) => c.group.groupKey).join(','), 'ROUTE=IV-bolus|DOSE=1');
    expect(inf.blocked.map((b) => `${b.route}:${b.doseKey}`).join(','), 'IV-infusion:1,IV-infusion:3');
    expect(inf.blocked[0].reason.includes('infusion'), true);
    // An IV arm of ONE row still forms its own stratum.
    const single = frame([{dose: 1, prefix: 'P'}]);
    single.col(DERIVED_ROUTE_CODE_COL)!.set(0, nca.ROUTE_IV_BOLUS);
    const s = cellsOf(computeSparseStrata(single, spec()));
    // Strata come in first-seen route order: row 0 is the IV one.
    expect(s.map((c) => `${c.group.groupKey}:${c.rows.length}`).join(','),
      'ROUTE=IV-bolus|DOSE=1:1,ROUTE=PO|DOSE=1:23');
  }, {timeout: TIMEOUT});

  test('AD-3: "10" and "10.0" are one composite keyed "10" with a note (both paths); 2 < 10 numerically',
    async () => {
      const df = frame([{dose: '10', prefix: 'A'}, {dose: '10.0', prefix: 'B'}, {dose: '2', prefix: 'C'}],
        {textDose: true});
      const g = groupsOf(computeSparseStrata(df, spec()));
      expect(g.map((x) => x.groupKey).join(','), 'DOSE=2,DOSE=10', 'numeric order, merged spellings');
      expect(g[1].warnings.some((w) => w.code === 'SPARSE_DOSE_SPELLING_MERGED'), true, 'merge note');
      const res = computePipeline(df, synthMapping(), '');
      expect(compositeRows(res).length, 2, 'the pipeline pools the same two composites');
      const r10 = compositeRows(res).find((i) => num(res, 'DOSE', i) === 10)!;
      expect(str(res, 'INFO_FLAGS', r10).includes('SPARSE_DOSE_SPELLING_MERGED: '), true,
        `pipeline INFO carries the merge note: ${str(res, 'INFO_FLAGS', r10)}`);
      // Absolute basis: "2" still sorts before "10".
      const abs = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A'}, {dose: 2, prefix: 'B'}]),
        spec({dose: 'mg'})));
      expect(abs.map((x) => x.doseLevel.key).join(','), '2,10');
      // A text label groups by its string and sorts last.
      const lbl = groupsOf(computeSparseStrata(frame([{dose: 'High', prefix: 'A'}, {dose: '5', prefix: 'B'}],
        {textDose: true}), spec()));
      expect(lbl.map((x) => x.groupKey).join(','), 'DOSE=5,DOSE=High');
    }, {timeout: TIMEOUT});

  // ─── AD-7 / AD-8: missing concentrations and BLQ flags ──────────────────────

  test('AD-7: an unflagged missing conc is dropped, a flagged blank kept; the note counts per time; ' +
    'n = 1 keeps a finite SE', async () => {
    // A001–A003 sit at 0.25 h, A022–A024 at 24 h. Drop two at 0.25 h, all three at 24 h.
    const missing = {A001: null, A002: null, A022: null, A023: null, A024: null};
    const g = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A', conc: missing}]), spec()))[0];
    const note = g.warnings.find((w) => w.code === 'SPARSE_CONC_MISSING_DROPPED');
    expect(note?.message ?? '', '5 sample(s) with no concentration and no BLQ flag dropped from the composite: ' +
      '2 at 0.25 h, 3 at 24 h. No sample left at 24 h: that time is not in the composite.');
    expect(g.composite.timepoints.some((t) => t.time === 24), false, 'the all-dropped time leaves the composite');
    expect(g.composite.timepoints.find((t) => t.time === 0.25)?.n, 1, '0.25 h left with n = 1');
    expect(g.varianceModeled, true, 'the singleton variance is modeled (AC-7.4)');
    expect(g.minTimepointN, 1, 'drives the n<4 flag');
    expect(Number.isFinite(g.aucSe) && g.aucSe > 0, true, `finite SE > 0 (${g.aucSe})`);
    const fragility = String(assembleSparseResult([g], 'mg/kg').parameterTable.col('FRAGILITY')!.get(0));
    expect(fragility.includes('variance-modeled') && fragility.includes('n<4'), true,
      `the existing FRAGILITY flags fire (${fragility})`);
    expect(g.rawPoints.length, 19, '24 rows − 5 dropped');
    expect(g.rows.length, 24, 'member rows keep every row (the pipeline attributes per cell on them)');
    // The group-mean BW (CL / Vz / Vss dose) is over the pooled rows only (review-gate decision):
    // the frame's BW column is 0.25 kg everywhere, so give a dropped animal another weight.
    const heavyDf = frame([{dose: 10, prefix: 'A', conc: {A001: null}}]);
    heavyDf.col('BW')!.set(0, 1.0); // A001 is row 0
    const heavy = groupsOf(computeSparseStrata(heavyDf, spec()))[0];
    expectFloat(heavy.meanBodyWeightKg!, 0.25, 1e-12, 'the dropped animal\'s 1.0 kg is not in the mean');
    // A BLQ-flagged blank is kept (read as 0).
    const flagged = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A', conc: {A022: null},
      blq: {A022: true}}]), spec()))[0];
    expect(flagged.rawPoints.length, 24, 'flagged blank kept');
    expect(flagged.warnings.some((w) => w.code === 'SPARSE_CONC_MISSING_DROPPED'), false);
    // AC-7.6: the last time keeps only a BLQ sample — Tlast moves, and the note names that time.
    const lastBlq = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A',
      conc: {A022: null, A023: null, A024: null}, blq: {A022: true}}]), spec()))[0];
    const msg = lastBlq.warnings.find((w) => w.code === 'SPARSE_CONC_MISSING_DROPPED')?.message ?? '';
    expect(msg.endsWith('No measurable sample left at 24 h.'), true, msg);
    // AD-9: a drop changes the composite, so it changes the fingerprint.
    const full = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A'}]), spec()));
    const fp = (x: SparseGroupResult[]): string =>
      sparseCompositeUnits({groups: x, parameterTable: null!, compositeTable: null!, warnings: []})[0]
        .fingerprint;
    expect(fp([g]) !== fp(full), true, 'a dropped conc changes the fingerprint');
    // Every sample of a dose group missing → no composite, reported as blocked.
    const all: Record<string, null> = {};
    for (let k = 1; k <= 24; k++) all[`A${String(k).padStart(3, '0')}`] = null;
    const empty = computeSparseStrata(frame([{dose: 10, prefix: 'A', conc: all}, {dose: 5, prefix: 'B'}]), spec());
    if (!empty.ok) throw new Error(empty.reason);
    expect(empty.cells.map((c) => c.group.groupKey).join(','), 'DOSE=5');
    expect(empty.blocked.map((b) => b.doseKey).join(','), '10');
  }, {timeout: TIMEOUT});

  test('AD-8: text BLQ flags read as the bool ones (both paths); a flagged reported value is zeroed and counted',
    async () => {
      const arm: Arm = {dose: 10, prefix: 'A', conc: {A022: null, A023: null}, blq: {A022: true, A023: true}};
      const bool = groupsOf(computeSparseStrata(frame([arm]), spec()))[0];
      const textArm: Arm = {...arm, blq: {A022: 'BLQ', A023: '<LLOQ'}};
      const text = groupsOf(computeSparseStrata(frame([textArm], {textBlq: true}), spec()))[0];
      expectSameGroup(text, bool, 'text vs bool flags');
      expect(text.rawPoints.length, 24, 'text-flagged blanks kept, not dropped');
      const pBool = computePipeline(frame([arm]), synthMapping(), '');
      const pText = computePipeline(frame([textArm], {textBlq: true}), synthMapping(), '');
      expectSame(num(pText, 'AUC_LAST', compositeRows(pText)[0]), num(pBool, 'AUC_LAST', compositeRows(pBool)[0]),
        'pipeline text vs bool');
      expectSame(num(pText, 'AUC_LAST', compositeRows(pText)[0]), bool.aucLast, 'pipeline = app');
      // A flagged row with a reported value (user decision (a)): read as 0, with a count.
      const valued = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A', conc: {A022: 0.5},
        blq: {A022: 'Y'}}], {textBlq: true}), spec()))[0];
      const zero = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A', conc: {A022: 0},
        blq: {A022: true}}]), spec()))[0];
      expectSame(valued.aucLast, zero.aucLast, 'the flagged 0.5 counts as 0');
      expect(valued.warnings.find((w) => w.code === 'SPARSE_BLQ_FLAG_WITH_VALUE')?.message.startsWith('1 '), true);
      expect(zero.warnings.some((w) => w.code === 'SPARSE_BLQ_FLAG_WITH_VALUE'), false, 'a flagged 0 is no "value"');
    }, {timeout: TIMEOUT});

  // ─── The mouse fixture (rule 18) ────────────────────────────────────────────

  test('AD-7: the mouse fixture has no empty concentration cell (the drop rule changes nothing on it)', async () => {
    const df = await loadFixture('04_mouse_sparse_destructive.csv');
    expect(df.col('Conc_ng_mL')!.stats.missingValueCount, 0);
  }, {timeout: TIMEOUT});

  test('P3: mouse nominal time derived into min / s / d → every parameter equals the hours run; AUClast / SE / df ' +
    'match PKNCA — λz / t½ / CL / MRT identity-covered only', async () => {
    const df = await loadFixture('04_mouse_sparse_destructive.csv');
    derive(df, 'NominalTime_hr', 'NominalTime_min', 60);
    derive(df, 'NominalTime_hr', 'NominalTime_s', 3600);
    derive(df, 'NominalTime_hr', 'NominalTime_d', 1 / 24);
    const h = appGroups(df, mouseMapping('h'));
    expect(h.map((g) => g.groupKey).join(','), 'DOSE=10,DOSE=50');
    for (const [unit, col] of [['min', 'NominalTime_min'], ['s', 'NominalTime_s'], ['d', 'NominalTime_d']] as const) {
      const x = appGroups(df, mouseMapping(unit, col));
      expect(x.length, h.length, unit);
      x.forEach((g, i) => expectSameGroup(g, h[i], `${unit} ${g.groupKey}`));
      expect(x.map((g) => g.groupKey).join(','), 'DOSE=10,DOSE=50', `${unit}: keys unchanged`);
    }
    const min = appGroups(df, mouseMapping('min', 'NominalTime_min'));
    for (const ref of golden.destructive) {
      const g = min.find((x) => x.doseLevel.value === ref.dose_mg_kg)!;
      expectFloat(g.aucLast, ref.sparse_auclast, Math.abs(ref.sparse_auclast) * 0.001, `${ref.dose_mg_kg}: AUClast`);
      expectFloat(g.aucSe, ref.sparse_auc_se, Math.abs(ref.sparse_auc_se) * 0.01, `${ref.dose_mg_kg}: SE`);
      expectFloat(g.aucDf, ref.sparse_auc_df, Math.abs(ref.sparse_auc_df) * 0.01, `${ref.dose_mg_kg}: df`);
    }
  }, {timeout: TIMEOUT});

  test('F2 / AC-2.1: hours, mg/kg, one route — the app path equals the pre-BUG-19 direct compute ' +
    '(numbers, keys, fingerprints)', async () => {
    const df = await loadFixture('04_mouse_sparse_destructive.csv');
    const app = computeSparseForRecipe(df, mouseMapping('h')).result!;
    // The pre-change app path: conc and BW scaled by hand, `computeSparseNca` with the raw dose column
    // as both dose and group, the fingerprint's dose = String(doseMgKg).
    const ref = df.clone();
    derive(ref, 'Conc_ng_mL', 'Conc_mgL', 1 / 1000);
    derive(ref, 'BodyWeight_g', 'BW_kg', 1 / 1000);
    // A bool BLQ column (a no-op when the CSV parsed it bool; the pre-change path read bool / 0-1 only).
    const blq = ensureBoolBlq(ref, 'BLQ');
    const old = computeSparseNca(ref, {concCol: 'Conc_mgL', nominalTimeCol: 'NominalTime_hr', animalIdCol: 'MouseID',
      blqCol: blq, doseCol: 'Dose_mg_kg', bodyWeightCol: 'BW_kg', groupCol: 'Dose_mg_kg', route: 'PO',
      samplingTypeCol: 'SamplingType', concUnit: 'mg/L', acceptance: ACC, bootstrapIterations: 200});
    expect(app.groups.map((g) => g.groupKey).join(','), old.groups.map((g) => g.groupKey).join(','));
    app.groups.forEach((g, i) => {
      expectSameGroup(g, old.groups[i], g.groupKey);
      expect(String(g.doseLevel.value), String(old.groups[i].doseMgKg), `${g.groupKey}: fingerprint dose byte`);
    });
    expect(sparseCompositeUnits(app).map((u) => u.fingerprint).join(),
      sparseCompositeUnits(old).map((u) => u.fingerprint).join(), 'fingerprints byte-identical');
    expect(app.parameterTable.col('DOSE_LEVEL')!.getTag('label'), 'Dose (mg/kg)', 'AC-2.3 header from the unit');
  }, {timeout: TIMEOUT});

  test('AC-7.1: the mouse 10 mg/kg group with one measurable 1 h sample emptied → AUClast 3.717087 ' +
    '(PKNCA, row removed) on both paths', async () => {
    const df = await loadFixture('04_mouse_sparse_destructive.csv');
    const id = df.col('MouseID')!;
    const row = Array.from({length: df.rowCount}, (_, i) => i).find((i) => id.get(i) === 'M007')!;
    expect(df.col('NominalTime_hr')!.get(row), 1, 'M007 is the first measurable 1 h sample of 10 mg/kg');
    df.col('Conc_ng_mL')!.set(row, DG.FLOAT_NULL);
    const g = appGroups(df, mouseMapping('h')).find((x) => x.doseLevel.value === 10)!;
    expectFloat(g.aucLast, 3.717087, 3.717087 * 0.001, 'app AUClast (0.1 %)');
    expect(g.warnings.some((w) => w.code === 'SPARSE_CONC_MISSING_DROPPED'), true, 'the drop is disclosed');
    const res = computePipeline(df, JSON.stringify(mousePipelineMapping(df)), '');
    const r = compositeRows(res).find((i) => num(res, 'DOSE', i) === 10)!;
    expectFloat(num(res, 'AUC_LAST', r), 3.717087, 3.717087 * 0.001, 'pipeline AUClast (0.1 %)');
    expect(str(res, 'INFO_FLAGS', r).includes('SPARSE_CONC_MISSING_DROPPED: '), true, 'pipeline INFO, not REASONS');
    expect(str(res, 'REASONS', r).includes('SPARSE_CONC_MISSING_DROPPED'), false);
  }, {timeout: TIMEOUT});

  test('S-2 parity: the mouse fixture in minutes — ncaCompute and the app agree on every composite parameter',
    async () => {
      const df = await loadFixture('04_mouse_sparse_destructive.csv');
      const m = toMinutes(df, mousePipelineMapping(df));
      const res = computePipeline(df, JSON.stringify(m), '');
      const app = appGroups(df, mouseMapping('min', m.nominalTime as string));
      const cols: Array<[string, keyof SparseGroupResult]> = [['AUC_LAST', 'aucLast'], ['AUC_INF', 'aucInf'],
        ['CMAX', 'cmax'], ['TMAX', 'tmax'], ['LAMBDA_Z', 'lambdaZ'], ['HALF_LIFE', 'halfLife'], ['MRT', 'mrt'],
        ['PCT_EXTRAP', 'pctExtrap'], ['CL', 'cl'], ['VZ', 'vz']];
      expect(compositeRows(res).length, app.length);
      for (const g of app) {
        const r = compositeRows(res).find((i) => num(res, 'DOSE', i) === g.doseLevel.value)!;
        for (const [col, key] of cols) expectSame(num(res, col, r), g[key] as number, `${g.groupKey} ${col}`, 1e-6);
      }
    }, {timeout: TIMEOUT});

  // ─── The pipeline's per-cell attribution and keys ───────────────────────────

  test('F5 / N-1: two dose groups with their own sexes and comments keep their per-cell SEX, reasons and DOSE',
    async () => {
      const df = frame([{dose: 1, prefix: 'A', sex: 'F', comment: 'Terminal bleed'},
        {dose: 5, prefix: 'B', sex: 'M', comment: 'regurgitated dose'}]);
      const res = computePipeline(df, synthMapping(), '');
      const rows = compositeRows(res);
      expect(rows.length, 2, str(res, 'REASONS', 0));
      const at = (d: number): number => rows.find((i) => num(res, 'DOSE', i) === d)!;
      expect(str(res, 'SEX', at(1)), 'F'); expect(str(res, 'SEX', at(5)), 'M');
      expect(str(res, 'INFO_FLAGS', at(1)).includes('CRO comment: "Terminal bleed"'), true);
      expect(str(res, 'REASONS', at(1)).includes('regurgitated'), false, 'the other cell\'s comment stays there');
      expect(str(res, 'REASONS', at(5)).includes('CRO comment: "regurgitated dose"'), true);
      expect(str(res, 'ROUTE', at(1)), 'PO');
    }, {timeout: TIMEOUT});

  test('N-2: pipeline flag keys of a one-route analysis equal the app\'s composite keys', async () => {
    const df = frame([{dose: 10, prefix: 'A'}, {dose: 50, prefix: 'B'}]);
    const res = computePipeline(df, synthMapping(), '');
    const flags = pipelineFlagsOf(res, compositeRows(res), df, parseMapping(synthMapping()), parseRules('').acceptance);
    const app = cellsOf(computeSparseStrata(df, spec())).map((c) => sparseUnitKey(c.group.groupKey));
    expect(flags.map((f) => f.profileKey).sort().join(','), [...app].sort().join(','));
    expect(app.join(','), `${sparseUnitKey(sparseCompositeKey(null, '10'))},${sparseUnitKey('DOSE=50')}`);
  }, {timeout: TIMEOUT});

  // ─── AD-6: the NONMEM dose basis ────────────────────────────────────────────

  test('AD-6: NONMEM — non-composite reads Dose_mg as absolute mg (no BW factor); composite keeps the per-kg ' +
    'group column, unchanged', async () => {
    const base = {doseColName: 'Dose', doseUnit: 'µg/kg' as const};
    expect(JSON.stringify(sparseDoseSource({...base, mode: 'standard', compositeDoseGroupCol: null})),
      '{"doseColName":"Dose","doseUnit":"µg/kg"}', 'standard: the mapped column, its declared unit');
    expect(JSON.stringify(sparseDoseSource({...base, mode: 'nonmem', compositeDoseGroupCol: 'DOSE_MGKG'})),
      `{"doseColName":"${DOSE_COL}","doseUnit":"mg"}`, 'non-composite: AMT as absolute mg');
    expect(JSON.stringify(sparseDoseSource({...base, mode: 'nonmem-composite', compositeDoseGroupCol: 'DOSE_MGKG'})),
      '{"doseColName":"DOSE_MGKG","doseUnit":"mg/kg"}', 'composite: the per-kg group column');
    expect(JSON.stringify(sparseDoseSource({...base, mode: 'nonmem-composite', compositeDoseGroupCol: null})),
      `{"doseColName":"${DOSE_COL}","doseUnit":"mg"}`, 'composite fallback: absolute');
    // Non-composite: a constant 2.5 mg AMT ≡ 10 mg/kg × 0.25 kg — CL on the dose of record. Before
    // BUG-19 the app passed Dose_mg as a per-kg dose: CL from 2.5 × 0.25 = 0.625 "mg".
    const df = frame([{dose: 2.5, curveDose: 10, prefix: 'A'}]);
    df.columns.add(DG.Column.fromFloat64Array(DOSE_COL, new Float64Array(df.rowCount).fill(2.5)));
    const nm = sparseDoseSource({...base, mode: 'nonmem', compositeDoseGroupCol: null});
    const abs = groupsOf(computeSparseStrata(df, spec({dose: nm.doseUnit}, {doseCol: nm.doseColName})))[0];
    const perKg = groupsOf(computeSparseStrata(frame([{dose: 10, prefix: 'A'}]), spec()))[0];
    expectSame(abs.cl, perKg.cl, 'CL on the absolute AMT');
    expectFloat(abs.doseMg, 2.5, 1e-12, 'dose of record');
    // Composite: the real destructive TK NONMEM fixture — the per-kg branch is unchanged (direct
    // pre-change compute: DOSE_MGKG as the per-kg dose and group, conc scaled by hand).
    const norm = normalizeEventDrivenObservations(await loadFixture('06_rat_tk28d_nonmem.csv'), EC06).observations;
    const comp = sparseDoseSource({...base, mode: 'nonmem-composite', compositeDoseGroupCol: 'DOSE_MGKG'});
    const app = computeSparseForRecipe(norm, {
      subjectName: SUBJECT_COL, timeName: TAD_COL, concName: CONC_COL, doseColName: comp.doseColName,
      bwColName: null, units: {time: 'h', concentration: 'ng/mL', dose: comp.doseUnit, bodyWeight: 'kg'},
      route: {codesCol: DERIVED_ROUTE_CODE_COL, scalar: 'PO'}, lloq: {colName: null, scalar: 0}, blqColName: null,
      nominalTimeName: TAD_COL, animalIdName: SUBJECT_COL, declaredTopology: null, acceptance: ACC,
      bootstrapIterations: 50,
    });
    expect(app.error, null);
    const ref = norm.clone();
    derive(ref, CONC_COL, 'Conc_mgL', 1 / 1000);
    const routes = new Set(app.result!.groups.map((g) => g.route));
    expect(routes.size, 1, `06 is one route (${[...routes].join(',')})`);
    const old = computeSparseNca(ref, {concCol: 'Conc_mgL', nominalTimeCol: TAD_COL, animalIdCol: SUBJECT_COL,
      doseCol: 'DOSE_MGKG', groupCol: 'DOSE_MGKG', route: [...routes][0], concUnit: 'mg/L', acceptance: ACC,
      bootstrapIterations: 50});
    expect(app.result!.groups.map((g) => g.groupKey).join(','), old.groups.map((g) => g.groupKey).join(','));
    app.result!.groups.forEach((g, i) => expectSameGroup(g, old.groups[i], `06 ${g.groupKey}`));
    // TAD in minutes (the Time unit picker reads it on the NONMEM path too) = TAD in hours.
    derive(norm, TAD_COL, 'TAD_min', 60);
    const appMin = computeSparseForRecipe(norm, {
      subjectName: SUBJECT_COL, timeName: 'TAD_min', concName: CONC_COL, doseColName: comp.doseColName,
      bwColName: null, units: {time: 'min', concentration: 'ng/mL', dose: comp.doseUnit, bodyWeight: 'kg'},
      route: {codesCol: DERIVED_ROUTE_CODE_COL, scalar: 'PO'}, lloq: {colName: null, scalar: 0}, blqColName: null,
      nominalTimeName: 'TAD_min', animalIdName: SUBJECT_COL, declaredTopology: null, acceptance: ACC,
      bootstrapIterations: 50,
    });
    expect(appMin.result!.groups.map((g) => g.groupKey).join(','), app.result!.groups.map((g) => g.groupKey).join(','),
      'TAD min: keys unchanged');
    appMin.result!.groups.forEach((g, i) => expectSameGroup(g, app.result!.groups[i], `06 TAD min ${g.groupKey}`));
  }, {timeout: TIMEOUT});

  // ─── Decisions ──────────────────────────────────────────────────────────────

  test('P4: a minute study\'s composite fingerprint changes with the fix, so its decision resets', async () => {
    const df = await loadFixture('04_mouse_sparse_destructive.csv');
    derive(df, 'NominalTime_hr', 'NominalTime_min', 60);
    const fpOf = (m: SparseComputeMapping): string =>
      sparseCompositeUnits(computeSparseForRecipe(df, m).result)[0].fingerprint;
    // Pre-fix behaviour: minutes read as hours (the time unit ignored).
    const buggy = fpOf(mouseMapping('h', 'NominalTime_min'));
    const fixed = fpOf(mouseMapping('min', 'NominalTime_min'));
    expect(fixed, fpOf(mouseMapping('h')), 'the fixed minute run fingerprints as the hours run');
    expect(buggy !== fixed, true, 'the pre-fix fingerprint differs');
    const obs = makeSyntheticObs({subjects: ['S1', 'S1', 'S1'], times: Float64Array.from([0, 1, 2]),
      concs: Float64Array.from([10, 5, 2.5])});
    const ctrl = new NcaController(obs, singleSubjectContract(obs), singleSubjectScalars(), DEFAULT_RULES);
    const key = sparseUnitKey('DOSE=10');
    ctrl.decide(key, 'exclude', [{code: 'EXC_ERRATIC', text: 'x'}], buggy);
    expect(ctrl.reconcileDecisions(new Map([[key, fixed]]), (k) => k.startsWith('sparse:')), 1, 'auto-reset');
    expect(ctrl.decisionFor(key), null);
  }, {timeout: TIMEOUT});

  test('P4 live: a saved minute study reopens with its stale composite decision reset', async () => {
    const df = await loadFixture('04_mouse_sparse_destructive.csv');
    const study = `TEST-BUG19-${Date.now()}`;
    df.col('StudyID')!.init(() => study);
    const m = toMinutes(df, mousePipelineMapping(df));
    const repo = getNcaRepository();
    let id = '';
    try {
      const saved = await saveToStudies(df, computePipeline(df, JSON.stringify(m), ''));
      id = String(saved.col('StudyId')!.get(0));
      const aid = String(saved.col('AnalysisId')!.get(0));
      const {entry, df: stored} = await repo.loadAnalysis(aid);
      const key = sparseUnitKey('DOSE=10');
      const stale = {action: 'exclude' as const, reasons: [{code: 'EXC_ERRATIC', text: 'x'}], user: 'test',
        at: new Date().toISOString(), fingerprint: 'pre-bug19-fingerprint'};
      // A decision recorded before BUG-19, saved through the app's own Save.
      await repo.saveAnalysis(aid, entry.source.ref, {df: stored, recipe: recipeOf(entry),
        summaryCache: entry.summaryCache, workflow: {confirmed: entry.confirmed ?? null,
          lzOverrides: entry.lzOverrides ?? [], reviewDecisions: [{unitKey: key, current: stale, history: []}],
          stageCache: entry.stageCache ?? null}});
      const view = await openAnalysis(aid) as DG.TableView;
      try {
        const bag = view.temp as Record<string, unknown>;
        const t0 = Date.now();
        while (bag[RESULTS_GRID_TAG] == null && Date.now() - t0 < 60000) await DG.delay(200);
        const unit = liveWorkflowState(view).units.find((u) => u.key === key);
        expect(unit !== undefined, true, `composite ${key} is a unit`);
        expect(unit!.decision, null, 'the stale decision no longer holds');
      } finally {
        view.close();
      }
    } finally {
      if (id !== '') await repo.deleteStudy(id);
    }
  }, {timeout: 240000});
});
