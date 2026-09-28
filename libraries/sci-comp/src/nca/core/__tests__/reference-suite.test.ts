/**
 * Full validation gate for `computeNca` against PKNCA fixtures
 * (Task 1.9.5 in `docs/nca_development_plan_v2.md`).
 *
 * Loads each Phase 0 dataset (theoph, indometh, rat_simple) and the
 * matching PKNCA fixture, runs the entire NCA pipeline through
 * `computeNca`, and asserts the eight Phase 1 parameters fall within the
 * tolerances defined in §9.2 of the plan:
 *
 *   Cmax, AUClast, AUCinf:  < 0.1% relative
 *   Tmax:                   exact match
 *   lambda_z:               < 0.5% relative
 *   t½, CL, Vz:             < 1% relative
 *   %AUCextrap:             < 0.5% absolute
 *
 * If this gate fails, the offending parameter / subject / dataset is
 * surfaced — investigate before promoting Task 1.9.
 */

import {readFileSync} from 'fs';
import {join} from 'path';
import {computeNca} from '../compute-nca';
import {augmentProfile} from '../augment';
import {sparseAuc} from '../sparse';
import type {SparseInput} from '../sparse';
import {ROUTE_IV_BOLUS, ROUTE_IV_INFUSION, ROUTE_PO} from '../types';
import type {ProfileInputs, NcaRules, RouteCode} from '../types';

const PKNCA_RULES: NcaRules = {
  aucMethod: 'linear-up-log-down',
  blq: {
    preFirstMeasurable: 'set-zero',
    embedded: 'set-zero',
    afterLast: 'set-zero',
    consecutiveAfterLast: 'set-zero',
  },
  lambdaZ: {
    mode: 'auto-best-fit',
    minPoints: 3,
    minRSquared: 0.85,
    excludeCmax: true,
    adjRSquaredFactor: 1e-4,
  },
  extrapWarnPct: 20,
  extrapErrorPct: 50,
  extrapWarnPctAumc: 20,
  compensatedSummation: false,
};

interface FixtureProfile {
  profile_key: {subject: string; route: string};
  parameters: {
    cmax: number | null;
    tmax: number | null;
    auclast: number | null;
    aucinf: number | null;
    pct_aucextrap: number | null;
    lambda_z: number | null;
    half_life: number | null;
    cl: number | null;
    vz: number | null;
    // FR-200 moment / lag params (null = not applicable for the route).
    aumclast?: number | null;
    aumcinf_obs?: number | null;
    mrt?: number | null;
    vss?: number | null; // present for IV; null for extravascular (gate)
    tlag?: number | null; // present for extravascular; null for IV (gate)
    pct_aumcextrap?: number | null;
  };
  /** PKNCA fit provenance, committed alongside the parameters. */
  provenance?: {
    lambda_z_n_points?: number | null;
    lambda_z_r_squared?: number | null;
    lambda_z_adj_r_squared?: number | null;
    lambda_z_time_first?: number | null;
    lambda_z_time_last?: number | null;
    clast_obs?: number | null;
    c0_extrapolated?: number | null;
    /** PKNCA's own `span.ratio` — the rule-18 oracle for `LambdaZResult.spanRatio`. */
    span_ratio?: number | null;
    /** PKNCA's own `c0` PPTESTCD on the RAW IV-bolus profile — the independent
     *  oracle for `provenance.c0.value` (closes the c0 circularity, F8). */
    c0_pknca?: number | null;
    /** Back-extrapolated share of AUCinf — hand formula on PKNCA's `c0` + first
     *  observation + `aucinf.obs` (Phoenix `AUC_%Back_Ext`, no PKNCA equivalent). */
    pct_auc_back_extrap?: number | null;
  };
}

interface FixtureFile {
  dataset: string;
  profiles: FixtureProfile[];
}

const DATASETS_DIR = join(__dirname, '..', '..', '__tests__', 'datasets');
const FIXTURES_DIR = join(__dirname, '..', '..', '__tests__', 'fixtures');

function parseCsv(text: string): {headers: string[]; rows: string[][]} {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0);
  const stripQuotes = (s: string) => s.replace(/^"|"$/g, '');
  const headers = lines[0].split(',').map(stripQuotes);
  const rows = lines.slice(1).map((l) => l.split(',').map(stripQuotes));
  return {headers, rows};
}

function loadFixture(name: string): FixtureFile {
  return JSON.parse(readFileSync(join(FIXTURES_DIR, name), 'utf-8'));
}

interface SubjectRows {
  subject: string;
  rows: Record<string, number>[];
}

function loadAndGroup(csvName: string): SubjectRows[] {
  const text = readFileSync(join(DATASETS_DIR, csvName), 'utf-8');
  const {headers, rows} = parseCsv(text);
  const subjCol = headers.indexOf('Subject');
  const grouped = new Map<string, Record<string, number>[]>();
  for (const r of rows) {
    const s = r[subjCol];
    if (!grouped.has(s)) grouped.set(s, []);
    const obj: Record<string, number> = {};
    for (let i = 0; i < headers.length; i++)
      if (i !== subjCol) obj[headers[i]] = parseFloat(r[i]);
    grouped.get(s)!.push(obj);
  }
  return Array.from(grouped.entries()).map(([subject, rows]) => ({subject, rows}));
}

function buildInputs(
  rows: Record<string, number>[], timeCol: string,
  dose: number, route: RouteCode, infusionDuration: number | null = null,
): ProfileInputs {
  return {
    time: Float64Array.from(rows.map((r) => r[timeCol])),
    conc: Float64Array.from(rows.map((r) => r.conc)),
    blqMask: new Uint8Array(rows.length),
    lloq: 0.001,
    dose,
    doseUnits: 'mg',
    concentrationUnits: 'mg/L',
    timeUnits: 'h',
    route,
    infusionDuration,
    bodyWeight: null,
  };
}

interface Tolerances {
  cmax: number; // relative
  auclast: number; // relative
  aucinf: number; // relative
  lambdaZ: number; // relative
  halfLife: number; // relative
  cl: number; // relative
  vz: number; // relative
  pctExtrap: number; // absolute (percentage points)
  aumcLast: number; // relative (area-class, AD-12)
  aumcInf: number; // relative
  mrt: number; // relative (derived-class)
  vss: number; // relative (derived volume)
  pctExtrapAumc: number; // absolute (percentage points)
}

// Tolerance classes by dimensional analogy to NFR-04/05 (AD-12), set only
// after MEASURING the actual core-vs-PKNCA deviation across all four fixtures
// in both summation modes (see REGEN.md → "Validation output"). The measured
// max deviations land ~9 orders of magnitude inside these gates — aumcLast
// 7e-12, aumcInf/mrt/vss ~2e-11 (rel), %AUMCextrap 1e-9 (pp), tlag exact —
// so the analogy values hold with enormous margin (the core matches PKNCA to
// floating-point, not merely to the gate).
const TOL: Tolerances = {
  cmax: 0.001,
  auclast: 0.001,
  aucinf: 0.001,
  lambdaZ: 0.005,
  halfLife: 0.01,
  cl: 0.01,
  vz: 0.01,
  pctExtrap: 0.5,
  aumcLast: 0.001,
  aumcInf: 0.001,
  mrt: 0.01,
  vss: 0.01,
  pctExtrapAumc: 0.5,
};

function relErr(got: number, expected: number): number {
  return Math.abs(got - expected) / Math.abs(expected);
}

const NO_SKIP: ReadonlySet<string> = new Set();

/**
 * Assert one profile against its fixture. `rules` defaults to the PKNCA
 * configuration every committed fixture was produced under; `skip` names
 * fixture parameter keys the caller asserts separately (the 06 BLQ-rules
 * suite, where a few cells are documented divergences rather than parity).
 */
function assertProfile(
  subject: string, dataset: string,
  inputs: ProfileInputs, fx: FixtureProfile,
  rules: NcaRules = PKNCA_RULES, skip: ReadonlySet<string> = NO_SKIP,
) {
  const r = computeNca(inputs, rules);
  const p = fx.parameters;
  const tag = `${dataset} subject ${subject}`;
  const on = (key: string): boolean => !skip.has(key);

  if (on('cmax') && p.cmax !== null)
    expect(relErr(r.values.cmax, p.cmax)).toBeLessThan(TOL.cmax);

  if (on('tmax') && p.tmax !== null) {
    // Tmax — exact match (per §9.2)
    expect(r.values.tmax).toBe(p.tmax);
  }
  if (on('auclast') && p.auclast !== null)
    expect(relErr(r.values.aucLast, p.auclast)).toBeLessThan(TOL.auclast);

  if (on('aucinf') && p.aucinf !== null)
    expect(relErr(r.values.aucInf, p.aucinf)).toBeLessThan(TOL.aucinf);

  if (on('lambda_z') && p.lambda_z !== null)
    expect(relErr(r.values.lambdaZ, p.lambda_z)).toBeLessThan(TOL.lambdaZ);

  if (on('half_life') && p.half_life !== null)
    expect(relErr(r.values.halfLife, p.half_life)).toBeLessThan(TOL.halfLife);

  if (on('cl') && p.cl !== null)
    expect(relErr(r.values.cl, p.cl)).toBeLessThan(TOL.cl);

  if (on('vz') && p.vz !== null)
    expect(relErr(r.values.vz, p.vz)).toBeLessThan(TOL.vz);

  if (on('pct_aucextrap') && p.pct_aucextrap !== null) {
    expect(Math.abs(r.values.pctExtrap - p.pct_aucextrap))
      .toBeLessThan(TOL.pctExtrap);
  }

  // ── FR-200 moment / lag parameters ──────────────────────────────────────
  if (on('aumclast') && p.aumclast != null)
    expect(relErr(r.values.aumcLast, p.aumclast)).toBeLessThan(TOL.aumcLast);

  if (on('aumcinf_obs') && p.aumcinf_obs != null)
    expect(relErr(r.values.aumcInf, p.aumcinf_obs)).toBeLessThan(TOL.aumcInf);

  if (on('mrt') && p.mrt != null)
    expect(relErr(r.values.mrt, p.mrt)).toBeLessThan(TOL.mrt);

  // Vss: gate-first. null (extravascular) ⇒ assert the route gate (NaN);
  // a number (IV) ⇒ assert parity. `== null` also covers a missing key.
  if (on('vss')) {
    if (p.vss == null)
      expect(Number.isNaN(r.values.vss)).toBe(true);
    else
      expect(relErr(r.values.vss, p.vss)).toBeLessThan(TOL.vss);
  }

  // Tlag: gate-first. null (IV) ⇒ assert the route gate (NaN); a number
  // (extravascular) ⇒ assert the exact observed time.
  if (on('tlag')) {
    if (p.tlag == null)
      expect(Number.isNaN(r.values.tlag)).toBe(true);
    else
      expect(r.values.tlag).toBe(p.tlag);
  }

  if (on('pct_aumcextrap') && p.pct_aumcextrap != null) {
    expect(Math.abs(r.values.pctExtrapAumc - p.pct_aumcextrap))
      .toBeLessThan(TOL.pctExtrapAumc);
  }

  // Quiet helper — used to avoid unused-arg lint
  void tag;
}

describe('reference suite (Task 1.9.5) — full computeNca pipeline vs PKNCA', () => {
  describe('01 Theoph (12 subjects, oral)', () => {
    const subjects = loadAndGroup('01_theoph.csv');
    const fixture = loadFixture('01_theoph.json');
    const fxBySubject = new Map(
      fixture.profiles.map((f) => [f.profile_key.subject, f]));

    it.each(subjects.map((s) => [s.subject, s] as const))(
      'subject %s — all parameters within tolerance',
      (subjectId, group) => {
        const fx = fxBySubject.get(subjectId);
        expect(fx).toBeDefined();
        const dose = group.rows[0].Dose * group.rows[0].Wt; // mg/kg × kg
        const inputs = buildInputs(group.rows, 'Time', dose, ROUTE_PO);
        assertProfile(subjectId, '01_theoph', inputs, fx!);
      });
  });

  describe('02 Indomethacin (6 subjects, IV bolus, dose = 25 mg)', () => {
    const subjects = loadAndGroup('02_indometh.csv');
    const fixture = loadFixture('02_indometh.json');
    const fxBySubject = new Map(
      fixture.profiles.map((f) => [f.profile_key.subject, f]));

    it.each(subjects.map((s) => [s.subject, s] as const))(
      'subject %s — all parameters within tolerance',
      (subjectId, group) => {
        const fx = fxBySubject.get(subjectId);
        expect(fx).toBeDefined();
        const inputs = buildInputs(group.rows, 'time', 25, ROUTE_IV_BOLUS);
        assertProfile(subjectId, '02_indometh', inputs, fx!);
      });
  });

  describe('03 Rat synthetic (8 subjects, oral, dose = 2.5 mg)', () => {
    const subjects = loadAndGroup('03_rat_simple.csv');
    const fixture = loadFixture('03_rat_simple.json');
    const fxBySubject = new Map(
      fixture.profiles.map((f) => [f.profile_key.subject, f]));

    it.each(subjects.map((s) => [s.subject, s] as const))(
      'subject %s — all parameters within tolerance',
      (subjectId, group) => {
        const fx = fxBySubject.get(subjectId);
        expect(fx).toBeDefined();
        const dose = group.rows[0].Dose;
        const inputs = buildInputs(group.rows, 'Time', dose, ROUTE_PO);
        assertProfile(subjectId, '03_rat_simple', inputs, fx!);
      });
  });

  // The only fixture that exercises the −T_inf/2 infusion correction (AC-D3).
  // Synthetic 1-compartment model (CL=2, V=10, k=0.2, dose=100, T_inf=1h) →
  // analytic AUCinf=50, MRT_iv=5, Vss=10; PKNCA references validate the wiring.
  describe('04 IV infusion (1 subject, dose = 100 mg, T_inf = 1 h)', () => {
    const subjects = loadAndGroup('04_iv_infusion.csv');
    const fixture = loadFixture('04_iv_infusion.json');
    const fxBySubject = new Map(
      fixture.profiles.map((f) => [f.profile_key.subject, f]));

    it.each(subjects.map((s) => [s.subject, s] as const))(
      'subject %s — all parameters within tolerance (incl. −T_inf/2 MRT/Vss)',
      (subjectId, group) => {
        const fx = fxBySubject.get(subjectId);
        expect(fx).toBeDefined();
        const inputs = buildInputs(
          group.rows, 'time', 100, ROUTE_IV_INFUSION, 1);
        assertProfile(subjectId, '04_iv_infusion', inputs, fx!);
      });
  });

  // UC-04 sparse / destructive-sampling AUC parity (AC-D1, rule 18). Consumed in
  // the gate-aware reference suite — like the 04 IV-infusion fixture — so a sparse
  // regression fails the same gate as a dense one. Destructive doses come from
  // PKNCA 0.12.1; the batch case from the hand-derived Holder oracle (df=NA in
  // PKNCA). Tolerances: AUC ≤ 0.1 % (NFR-04); SE, df ≤ 1 % (NFR-05).
  describe('05 mouse sparse (PKNCA pk.calc.sparse_auclast + batch Holder oracle)', () => {
    interface SparseFx {
      destructive: Array<{
        dose_mg_kg: number; sparse_auclast: number; sparse_auc_se: number; sparse_auc_df: number;
      }>;
      batch: {
        times: number[]; samples: Array<{time: number; animal: string[]; conc: number[]}>;
        sparse_auclast: number; sparse_auc_se: number; sparse_auc_df: number;
      };
    }
    const fx: SparseFx = JSON.parse(
      readFileSync(join(FIXTURES_DIR, '05_mouse_sparse.json'), 'utf8'));
    const csv = readFileSync(
      join(DATASETS_DIR, '04_mouse_sparse_destructive.csv'), 'utf8').trim().split(/\r?\n/);
    const rows = csv.slice(1).map((l) => l.split(','));
    const rel = (a: number, b: number): number => Math.abs(a - b) / Math.abs(b);

    it.each(fx.destructive.map((d) => [d.dose_mg_kg, d] as const))(
      'dose %s mg/kg — sparse AUClast/SE/df match PKNCA within NFR-04/05',
      (dose, golden) => {
        const sub = rows.filter((f) => Number(f[7]) === dose);
        const input: SparseInput = {
          nominalTime: Float64Array.from(sub.map((f) => Number(f[8]))),
          conc: Float64Array.from(sub.map((f) => Number(f[10]) / 1000)),
          blqMask: Uint8Array.from(sub.map((f) => (f[11] === 'True' ? 1 : 0))),
          lloq: 1 / 1000,
          animalId: Int32Array.from(sub.map((f) => Number(f[1].replace(/\D/g, '')))),
        };
        const res = sparseAuc(input, {blqRule: 'set-zero'});
        expect(res.topology).toBe('destructive');
        expect(rel(res.auc, golden.sparse_auclast)).toBeLessThan(1e-3);
        expect(rel(res.se, golden.sparse_auc_se)).toBeLessThan(1e-2);
        expect(rel(res.df, golden.sparse_auc_df)).toBeLessThan(1e-2);
      });

    it('batch — Holder A1+A3 + Nedelman-Jia df match the hand-derived oracle', () => {
      const b = fx.batch;
      const nominalTime: number[] = []; const conc: number[] = []; const animalId: number[] = [];
      const code = new Map<string, number>();
      for (const s of b.samples) {
        s.animal.forEach((a, k) => {
          if (!code.has(a)) code.set(a, code.size + 1);
          nominalTime.push(s.time); conc.push(s.conc[k]); animalId.push(code.get(a)!);
        });
      }
      const input: SparseInput = {
        nominalTime: Float64Array.from(nominalTime),
        conc: Float64Array.from(conc),
        blqMask: new Uint8Array(nominalTime.length),
        lloq: 0,
        animalId: Int32Array.from(animalId),
      };
      const res = sparseAuc(input, {blqRule: 'set-zero'});
      expect(res.topology).toBe('batch');
      expect(rel(res.auc, b.sparse_auclast)).toBeLessThan(1e-3);
      expect(rel(res.se, b.sparse_auc_se)).toBeLessThan(1e-2);
      expect(rel(res.df, b.sparse_auc_df)).toBeLessThan(1e-2);
    });
  });
});


/**
 * Span-ratio values on the reference corpus.
 *
 * NO expectation changes anywhere above: `minSpanRatio` is diagnostic, so window
 * selection — and therefore every asserted PKNCA-parity value — is untouched.
 * These are additive value assertions for the new field.
 *
 * The `< 2` assertions are the empirical claim the contribution rests on: PKNCA
 * returns a valid lambda_z for these profiles and so do we, correctly — but the
 * terminal slope rests on less than two half-lives and nothing in the output said
 * so. `02_indometh` subject 1 is the worst case in the corpus at 0.69, and it
 * carries adj-R² 0.994.
 */
describe('reference suite — LambdaZResult.spanRatio', () => {
  /** PKNCA's own `span.ratio`, as committed by `regen-fixtures.R`. */
  function oracle(fx: FixtureProfile): number {
    const v = fx.provenance?.span_ratio;
    expect(v).toBeDefined();
    expect(v).not.toBeNull();
    return v!;
  }

  /**
   * ALL FOUR committed fixtures, table-driven.
   *
   * Deliberately not two hand-written blocks for theoph/indometh: rule 18 names
   * `rat_simple` explicitly, `04_iv_infusion` is the only fixture exercising the
   * infusion path, and every one of them already carries the oracle. Driving the
   * table means a fifth fixture cannot be silently left uncovered.
   * (Peer review 2026-08-11, finding F2 — the first cut covered 2 of 4.)
   */
  const SPAN_FIXTURES: ReadonlyArray<{
    csv: string; json: string; timeCol: string; route: RouteCode;
    dose: (row: Record<string, number>) => number; tInf?: number;
  }> = [
    {csv: '01_theoph.csv', json: '01_theoph.json', timeCol: 'Time',
      route: ROUTE_PO, dose: (r) => r.Dose * r.Wt},
    {csv: '02_indometh.csv', json: '02_indometh.json', timeCol: 'time',
      route: ROUTE_IV_BOLUS, dose: () => 25},
    {csv: '03_rat_simple.csv', json: '03_rat_simple.json', timeCol: 'Time',
      route: ROUTE_PO, dose: (r) => r.Dose},
    {csv: '04_iv_infusion.csv', json: '04_iv_infusion.json', timeCol: 'time',
      route: ROUTE_IV_INFUSION, dose: () => 100, tInf: 1},
  ];

  it.each(SPAN_FIXTURES.map((f) => [f.json, f] as const))(
    '%s — every subject matches PKNCA span.ratio on PKNCA\'s own window',
    (_json, spec) => {
      const subjects = loadAndGroup(spec.csv);
      const fxBySubject = new Map(loadFixture(spec.json).profiles
        .map((f) => [f.profile_key.subject, f]));
      expect(subjects.length).toBeGreaterThan(0);

      for (const group of subjects) {
        const fx = fxBySubject.get(group.subject)!;
        expect(fx).toBeDefined();
        const inputs = buildInputs(group.rows, spec.timeCol,
          spec.dose(group.rows[0]), spec.route, spec.tInf ?? null);
        const r = computeNca(inputs, PKNCA_RULES);
        expect(r.provenance.lambdaZ).not.toBeNull();
        expect(relErr(r.provenance.lambdaZ!.spanRatio, oracle(fx)))
          .toBeLessThan(TOL.lambdaZ);
        // Same window PKNCA chose — without this the value agreement could be
        // coincidental rather than a parity result.
        expect(r.provenance.lambdaZ!.tStart).toBe(fx.provenance!.lambda_z_time_first);
        expect(r.provenance.lambdaZ!.tEnd).toBe(fx.provenance!.lambda_z_time_last);
      }
    });

  it('the corpus incidence the contribution rests on: 4 of 18 below span 2', () => {
    // The empirical claim, asserted against the committed PKNCA output rather
    // than restated in prose. `02_indometh` subject 1 is the worst case: its
    // slope rests on under ONE half-life while adj-R² reads ~0.994 — "excellent"
    // by every other statistic in the PKNCA output set. That inversion is why
    // spanRatio exists.
    //
    // Scoped to the two HUMAN reference fixtures on purpose — "4 of 18" is a
    // claim about theoph+indometh specifically (12 + 6 profiles). The synthetic
    // rat / infusion fixtures are covered for PARITY by the table above but are
    // not part of this incidence statistic.
    const below: string[] = [];
    for (const file of ['01_theoph.json', '02_indometh.json']) {
      for (const p of loadFixture(file).profiles)
        if (oracle(p) < 2) below.push(`${file}/${p.profile_key.subject}`);
    }

    expect(below.length).toBe(4);

    const ind1 = loadFixture('02_indometh.json').profiles
      .find((f) => f.profile_key.subject === '1')!;
    expect(oracle(ind1)).toBeLessThan(1);
    expect(ind1.provenance!.lambda_z_adj_r_squared).toBeGreaterThan(0.99);
  });
});

/**
 * IV-bolus dose-time gate (GROK-20960 slice U1 / F2) on the real IV corpus.
 *
 * Fx-2 — PKNCA's own `c0`: the committed `c0_extrapolated` was the CORE's c0
 * fed back into PKNCA for the AUC run, so until now nothing independent said
 * the back-extrapolation was right. `c0_pknca` is PKNCA's `c0` PPTESTCD on the
 * raw profile; both the committed value and `provenance.c0.value` are asserted
 * against it (measured equal to ≥ 10 digits on 6/6 subjects).
 *
 * Fx-1 — invariance: a pre-dose `(0, 0)` row — unflagged or BLQ-flagged — must
 * not change a single IV-bolus parameter, because the row is REPLACED by the
 * same `(0, c0)` the no-row profile gets. Before the gate the row was taken as
 * a measured t=0 value and integrated from 0: subject 1 read AUClast 1.7194
 * instead of 2.0099 (−14.5 %, silent). Asserted against the COMMITTED fixture
 * values, not against the no-row run alone, so the test cannot pass by both
 * paths being wrong the same way.
 */
describe('02 Indomethacin — IV-bolus dose-time gate (Fx-1 invariance + Fx-2 c0 oracle)', () => {
  const subjects = loadAndGroup('02_indometh.csv');
  const fxBySubject = new Map(
    loadFixture('02_indometh.json').profiles.map((f) => [f.profile_key.subject, f]));
  const withDoseTimeRow = (rows: Record<string, number>[], flag: 0 | 1): ProfileInputs => {
    const base = buildInputs(rows, 'time', 25, ROUTE_IV_BOLUS);
    const blqMask = new Uint8Array(base.time.length + 1);
    blqMask[0] = flag;
    return {
      ...base,
      time: Float64Array.from([0, ...base.time]),
      conc: Float64Array.from([0, ...base.conc]),
      blqMask,
    };
  };

  it.each(subjects.map((s) => [s.subject, s] as const))(
    'subject %s — c0 matches PKNCA\'s own c0 (fixture AND provenance)', (subjectId, group) => {
      const fx = fxBySubject.get(subjectId)!;
      const oracle = fx.provenance!.c0_pknca!;
      expect(oracle).toBeGreaterThan(0);
      expect(relErr(fx.provenance!.c0_extrapolated!, oracle)).toBeLessThan(TOL.auclast);
      const r = computeNca(buildInputs(group.rows, 'time', 25, ROUTE_IV_BOLUS), PKNCA_RULES);
      expect(r.provenance.c0).not.toBeNull();
      expect(r.provenance.c0!.method).toBe('logslope');
      expect(r.provenance.c0!.replacedDoseTimeRow).toBe(false);
      expect(relErr(r.provenance.c0!.value, oracle)).toBeLessThan(TOL.auclast);
    });

  it.each(subjects.map((s) => [s.subject, s] as const))(
    'subject %s — back-extrapolated share of AUCinf matches the PKNCA-input formula', (subjectId, group) => {
      const fx = fxBySubject.get(subjectId)!;
      const r = computeNca(buildInputs(group.rows, 'time', 25, ROUTE_IV_BOLUS), PKNCA_RULES);
      expect(Math.abs(r.provenance.c0!.pctAucBackExtrap - fx.provenance!.pct_auc_back_extrap!))
        .toBeLessThan(TOL.pctExtrap);
    });

  it('corpus incidence the diagnostic exists for: every indometh subject back-extrapolates > 15 % of AUCinf', () => {
    // A fifth of the exposure on the library's headline IV fixture is in
    // unmeasured territory — with every other statistic reading "fine".
    for (const fx of fxBySubject.values())
      expect(fx.provenance!.pct_auc_back_extrap!).toBeGreaterThan(15);
    expect(fxBySubject.get('1')!.provenance!.pct_auc_back_extrap!).toBeCloseTo(20.6, 1);
  });

  describe.each([['unflagged', 0], ['BLQ-flagged', 1]] as const)(
    'a %s (0, 0) pre-dose row is replaced — every parameter within tolerance of the committed fixture',
    (_label, flag) => {
      it.each(subjects.map((s) => [s.subject, s] as const))('subject %s', (subjectId, group) => {
        const fx = fxBySubject.get(subjectId)!;
        const noRow = computeNca(buildInputs(group.rows, 'time', 25, ROUTE_IV_BOLUS), PKNCA_RULES);
        const inputs = withDoseTimeRow(group.rows, flag);
        // The committed-fixture parity gate, verbatim.
        assertProfile(subjectId, `02_indometh+(0,0,blq=${flag})`, inputs, fx);
        const r = computeNca(inputs, PKNCA_RULES);
        expect(r.provenance.c0!.replacedDoseTimeRow).toBe(true);
        expect(r.provenance.c0!.value).toBe(noRow.provenance.c0!.value);
        expect(r.provenance.c0!.pctAucBackExtrap).toBe(noRow.provenance.c0!.pctAucBackExtrap);
        // Same λz window as the no-row profile (augmented indices coincide: the
        // replaced profile has the no-row profile's length).
        expect(Array.from(r.provenance.lambdaZ!.pointsUsed))
          .toEqual(Array.from(noRow.provenance.lambdaZ!.pointsUsed));
        expect(r.provenance.lambdaZ!.tStart).toBe(fx.provenance!.lambda_z_time_first);
        expect(r.provenance.lambdaZ!.tEnd).toBe(fx.provenance!.lambda_z_time_last);
        expect(r.values.cmax).toBe(fx.parameters.cmax);
        expect(r.values.tmax).toBe(fx.parameters.tmax);
      });
    });
});

/**
 * 06 BLQ rules — per-rule PKNCA parity (GROK-20960 slice U2 / F1, Fx-3).
 *
 * Before this slice every BLQ rule integrated as `exclude`: the integrator and
 * the λz filter read the effective BLQ mask, so a declared `set-zero` /
 * `set-half-lloq` substitution never reached AUC/AUMC/λz. Four hand-authored
 * subjects × the FIVE rule blocks of `06_blq_rules.json`, each block a real
 * PKNCA 0.12.1 run under the `conc.blq` option that means the same thing
 * (see `regen-fixtures.R`). Table-driven from the fixture's own block list, so
 * a sixth block cannot be silently left uncovered.
 *
 * Three cells are DOCUMENTED DIVERGENCES, asserted as such — never skipped:
 *  - R-B (PKNCA `drop` × 3): PKNCA's AUC is NA once the t=0 zero is dropped
 *    (no origin datum; it does not extrapolate to the interval start), while
 *    sci-comp's `exclude` prepends `(0, 0)` by the extravascular convention —
 *    the R-D construction. The fixture's `auc_oracle_block` points R-B's AUC
 *    family at R-D; R-B itself pins Cmax / Tmax / λz / Tlag under PKNCA's own
 *    option.
 *  - R-C with TRAILING substitutes (P2, P3): PKNCA integrates `auclast` through
 *    the LLOQ/2 tail but keeps `tlast` / `clast.obs` at the last above-LOQ
 *    observation and extrapolates AUCinf from THAT (`aucinf = auclast +
 *    clast.obs / λz`, the two limbs on different profiles). sci-comp honours
 *    the substitutes consistently (D1): `cLast = LLOQ/2` at the last sample,
 *    AUCinf = AUClast + cLast/λz. **This is a HOUSE decision, not a cited
 *    vendor behaviour** — it is chosen because it keeps ONE terminal anchor
 *    for both limbs (the same cLast/tLast that ends AUClast starts the tail),
 *    where PKNCA's own treatment is internally mixed. An earlier draft
 *    attributed it to Phoenix; that attribution could not be verified against
 *    a primary Certara source and has been withdrawn (peer review 2026-09-22).
 *    Asserted on both
 *    sides so the divergence is characterised, not hidden.
 *  - R-D Tlag: PKNCA computes `tlag` on the PRE-clean data, so it agrees with
 *    the mask-based definition under keep / drop / numeric; with the BLQ rows
 *    physically removed it cannot see the lag, and reports the predose 0.
 *    sci-comp's Tlag is mask-based under every rule (AD-7).
 */
describe('06 BLQ rules — per-rule PKNCA parity (F1, GROK-20960)', () => {
  interface BlqBlock {
    id: string;
    sci_comp_blq: NcaRules['blq'];
    config: {conc_blq: unknown; auc_oracle_block: string};
    note: string;
    profiles: FixtureProfile[];
  }
  interface BlqFixture {
    config: {lloq: number};
    dataset_meta: {dose: {po: {value: number}; iv_bolus: {value: number}}};
    blocks: BlqBlock[];
  }
  interface BlqRow { time: number; conc: number; blq: 0 | 1; route: string; dose: number }

  const fixture: BlqFixture = JSON.parse(
    readFileSync(join(FIXTURES_DIR, '06_blq_rules.json'), 'utf-8'));
  const LLOQ = fixture.config.lloq;
  const {headers, rows} = parseCsv(readFileSync(join(DATASETS_DIR, '06_blq_rules.csv'), 'utf-8'));
  const col = (name: string) => headers.indexOf(name);
  const bySubject = new Map<string, BlqRow[]>();
  for (const r of rows) {
    const s = r[col('Subject')];
    if (!bySubject.has(s)) bySubject.set(s, []);
    bySubject.get(s)!.push({
      time: parseFloat(r[col('time')]), conc: parseFloat(r[col('conc')]),
      blq: r[col('blq')] === '1' ? 1 : 0, route: r[col('route')], dose: parseFloat(r[col('Dose')]),
    });
  }
  const subjects = Array.from(bySubject.keys());
  const inputsFor = (s: string): ProfileInputs => {
    const rs = bySubject.get(s)!;
    return {
      time: Float64Array.from(rs.map((r) => r.time)),
      conc: Float64Array.from(rs.map((r) => r.conc)),
      blqMask: Uint8Array.from(rs.map((r) => r.blq)),
      lloq: LLOQ, dose: rs[0].dose, doseUnits: 'mg', concentrationUnits: 'mg/L', timeUnits: 'h',
      route: rs[0].route === 'IV-bolus' ? ROUTE_IV_BOLUS : ROUTE_PO,
      infusionDuration: null, bodyWeight: null,
    };
  };
  const rulesFor = (b: BlqBlock): NcaRules => ({...PKNCA_RULES, blq: b.sci_comp_blq});
  const profileOf = (b: BlqBlock, s: string): FixtureProfile =>
    b.profiles.find((p) => p.profile_key.subject === s)!;
  const blockById = (id: string): BlqBlock => fixture.blocks.find((b) => b.id === id)!;
  /** Does the subject's profile END in a BLQ run (a trailing substitute under set-half-lloq)? */
  const hasTrailingBlq = (s: string): boolean => bySubject.get(s)!.slice(-1)[0].blq === 1;
  const isHalfLloq = (b: BlqBlock): boolean => b.sci_comp_blq.afterLast === 'set-half-lloq';
  const isRowsRemoved = (b: BlqBlock): boolean => b.config.conc_blq === 'rows-removed';
  /** sci-comp's documented Tlag: the last BLQ-or-zero sample time before the first
   *  measurable positive one, on the augmented series (0 when none precedes it). */
  const maskBasedTlag = (s: string): number => {
    const rs = bySubject.get(s)!;
    let last = 0;
    for (const r of rs) {
      if (r.blq === 0 && r.conc > 0) return last;
      last = r.time;
    }
    return NaN;
  };

  it('the fixture carries the five rule blocks and every subject in each', () => {
    expect(fixture.blocks.length).toBeGreaterThanOrEqual(5);
    for (const id of ['R-A', 'R-B', 'R-C', 'R-D', 'R-E']) expect(blockById(id)).toBeDefined();
    for (const b of fixture.blocks)
      expect(b.profiles.map((p) => p.profile_key.subject).sort()).toEqual([...subjects].sort());
    expect(subjects).toContain('P1');
    expect(subjects).toContain('I1');
  });

  const cases = fixture.blocks.flatMap((b) => subjects.map((s) => [b.id, s, b] as const));
  it.each(cases)('block %s subject %s — every parameter matches PKNCA (or its documented divergence)',
    (_id, s, b) => {
      const fx = profileOf(b, s);
      const inputs = inputsFor(s);
      const rules = rulesFor(b);
      const r = computeNca(inputs, rules);
      const aucFx = profileOf(blockById(b.config.auc_oracle_block), s);
      const LZ_FAMILY = ['lambda_z', 'half_life'];
      const INF_FAMILY = ['aucinf', 'cl', 'vz', 'pct_aucextrap', 'aumcinf_obs', 'mrt', 'vss', 'pct_aumcextrap'];
      // DOCUMENTED DIVERGENCE 1 — PKNCA reported a fit sci-comp's adj-R² floor
      // rejects. pk.calc.half.life (0.12.1 source) selects by `lambda.z > 0` +
      // the adj-R² tie-break only; `min.hl.r.squared` is not consulted there.
      // sci-comp's `minRSquared` IS a floor (lambda-z.ts: "no PKNCA equivalent").
      // On P2/P3 under set-half-lloq PKNCA fits the flat LLOQ/2 tail at r² 0.76 /
      // 0.68 — a slope through substitutes, not a terminal phase. We report
      // 'partial' (λz not estimated), which is the honest reading.
      const guardRejects = fx.parameters.lambda_z !== null &&
        fx.provenance!.lambda_z_adj_r_squared! < rules.lambdaZ.minRSquared;
      // DOCUMENTED DIVERGENCE 2 — set-half-lloq with a trailing substitute and an
      // ACCEPTED fit: tLast / cLast semantics differ (see the describe header).
      const trailingDiverges = isHalfLloq(b) && hasTrailingBlq(s) && !guardRejects &&
        fx.parameters.lambda_z !== null;
      const skip = new Set<string>(['auclast', 'aumclast', 'tlag']);
      if (guardRejects) for (const k of [...LZ_FAMILY, ...INF_FAMILY]) skip.add(k);
      if (trailingDiverges) for (const k of INF_FAMILY) skip.add(k);

      // The parity gate, verbatim, on everything that is straight PKNCA parity.
      assertProfile(s, `06_blq_rules/${b.id}`, inputs, fx, rules, skip);

      // AUC family through the fixture's oracle indirection (R-B → R-D).
      expect(aucFx.parameters.auclast).not.toBeNull();
      expect(relErr(r.values.aucLast, aucFx.parameters.auclast!)).toBeLessThan(TOL.auclast);
      expect(relErr(r.values.aumcLast, aucFx.parameters.aumclast!)).toBeLessThan(TOL.aumcLast);

      // λz: PKNCA NA ⇒ status 'partial' with the λz family NaN (Fx-5); a value ⇒
      // the SAME window, not just the same number — unless the floor rejects it.
      if (fx.parameters.lambda_z === null || guardRejects) {
        expect(r.status).toBe('partial');
        expect(Number.isNaN(r.values.lambdaZ)).toBe(true);
        expect(Number.isNaN(r.values.aucInf)).toBe(true);
        expect(Number.isFinite(r.values.aucLast)).toBe(true);
        if (guardRejects) {
          // Characterise the divergence from the fixture's own numbers.
          expect(fx.provenance!.lambda_z_adj_r_squared!).toBeLessThan(rules.lambdaZ.minRSquared);
          expect(fx.provenance!.lambda_z_r_squared!).toBeLessThan(rules.lambdaZ.minRSquared);
          expect(isHalfLloq(b) && hasTrailingBlq(s)).toBe(true); // only ever the substituted-tail case
        }
      } else {
        expect(r.status).toBe('ok');
        expect(r.provenance.lambdaZ!.tStart).toBe(fx.provenance!.lambda_z_time_first);
        expect(r.provenance.lambdaZ!.tEnd).toBe(fx.provenance!.lambda_z_time_last);
        expect(r.provenance.lambdaZ!.pointsUsed.length).toBe(fx.provenance!.lambda_z_n_points);
        expect(relErr(r.provenance.lambdaZ!.adjRSquared, fx.provenance!.lambda_z_adj_r_squared!))
          .toBeLessThan(TOL.lambdaZ);
      }

      if (trailingDiverges) {
        // DOCUMENTED DIVERGENCE (set-half-lloq, trailing substitutes). Ours: the
        // substitutes are real points — tLast is the last sample, cLast = LLOQ/2,
        // and AUCinf extrapolates from it (self-consistent — house choice).
        const rs = bySubject.get(s)!;
        const tLastOurs = rs[rs.length - 1].time;
        expect(r.values.aucInf).toBeCloseTo(r.values.aucLast + (LLOQ / 2) / r.values.lambdaZ, 10);
        // PKNCA's: auclast integrated to the last sample (asserted above) but
        // AUCinf extrapolated from clast.obs at the last ABOVE-LOQ time —
        // characterised from the fixture's own numbers, so the divergence is
        // a measured fact, not a skipped cell.
        const p = fx.parameters;
        const prov = fx.provenance!;
        expect(prov.lambda_z_time_last).toBe(tLastOurs); // the substitutes ARE in PKNCA's fit
        expect((prov as {tlast?: number}).tlast).toBeLessThan(tLastOurs); // …but not its tlast
        expect(p.aucinf!).toBeCloseTo(p.auclast! + prov.clast_obs! / p.lambda_z!, 6);
        expect(relErr(r.values.aucInf, p.aucinf!)).toBeGreaterThan(TOL.aucinf); // genuinely different
      }

      // Tlag: mask-based under EVERY rule (AD-7); PKNCA parity where PKNCA can
      // see the BLQ rows (every block except rows-removed).
      if (inputs.route === ROUTE_IV_BOLUS)
        expect(Number.isNaN(r.values.tlag)).toBe(true);
      else {
        expect(r.values.tlag).toBe(maskBasedTlag(s));
        if (!isRowsRemoved(b)) expect(r.values.tlag).toBe(fx.parameters.tlag);
      }
    });

  it('P1: set-zero and exclude are NOT bit-identical — AUClast differs by > 1 % (regression pin)', () => {
    const a = computeNca(inputsFor('P1'), rulesFor(blockById('R-A'))).values.aucLast;
    const bv = computeNca(inputsFor('P1'), rulesFor(blockById('R-B'))).values.aucLast;
    expect(relErr(a, bv)).toBeGreaterThan(0.01);
    // And each is the PKNCA number for its rule, not the other's.
    expect(relErr(a, profileOf(blockById('R-A'), 'P1').parameters.auclast!)).toBeLessThan(TOL.auclast);
    expect(relErr(bv, profileOf(blockById('R-D'), 'P1').parameters.auclast!)).toBeLessThan(TOL.auclast);
  });

  const everyField = (v: ReturnType<typeof computeNca>['values']) =>
    Object.entries(v).map(([k, x]) => [k, Number.isNaN(x) ? 'NaN' : x]);

  it('missing ≡ exclude on every ParameterValues field, for every subject (R-D ≡ R-B)', () => {
    for (const s of subjects) {
      const ex = computeNca(inputsFor(s), rulesFor(blockById('R-B')));
      const mi = computeNca(inputsFor(s), rulesFor(blockById('R-D')));
      expect(everyField(mi.values)).toEqual(everyField(ex.values));
      expect(mi.status).toBe(ex.status);
    }
  });

  it('the nca-studio shipped default ≡ set-zero on every field (R-E ≡ R-A) — GAP-W3', () => {
    for (const s of subjects) {
      const a = computeNca(inputsFor(s), rulesFor(blockById('R-A')));
      const e = computeNca(inputsFor(s), rulesFor(blockById('R-E')));
      expect(everyField(e.values)).toEqual(everyField(a.values));
      expect(e.status).toBe(a.status);
    }
  });

  it('P3 (Fx-5): partial under set-zero / exclude / missing (PKNCA λz NA)', () => {
    for (const id of ['R-A', 'R-B', 'R-D', 'R-E']) {
      expect(profileOf(blockById(id), 'P3').parameters.lambda_z).toBeNull();
      expect(computeNca(inputsFor('P3'), rulesFor(blockById(id))).status).toBe('partial');
    }
  });

  it('P2 / P3 under set-half-lloq: PKNCA fits the flat LLOQ/2 tail below our adj-R² floor → partial', () => {
    // The two documented cells of divergence 1, named explicitly so the case
    // count is pinned: exactly the trailing-substitute subjects whose tail is
    // flat, and no other cell in the table.
    const rc = blockById('R-C');
    const rejected = subjects.filter((s) => {
      const fx = profileOf(rc, s);
      return fx.parameters.lambda_z !== null &&
        fx.provenance!.lambda_z_adj_r_squared! < PKNCA_RULES.lambdaZ.minRSquared;
    });
    expect(rejected.sort()).toEqual(['P2', 'P3']);
    for (const s of rejected) {
      const r = computeNca(inputsFor(s), rulesFor(rc));
      expect(r.status).toBe('partial');
      expect(r.provenance.lambdaZ).toBeNull();
    }
    // …and P4 is the accepted-fit counterpart: its trailing substitute lands on
    // the line, both tools fit it, and only the AUCinf construction differs.
    const p4 = computeNca(inputsFor('P4'), rulesFor(rc));
    expect(p4.status).toBe('ok');
    expect(profileOf(rc, 'P4').provenance!.lambda_z_adj_r_squared!)
      .toBeGreaterThanOrEqual(PKNCA_RULES.lambdaZ.minRSquared);
  });

  it('P4 under set-half-lloq: the ACCEPTED fit contains an unmeasured point, and says so', () => {
    // The case adjusted R² provably cannot catch (peer review, 2026-09-22): the
    // trailing LLOQ/2 substitute at t = 24 is positive, survives the trailing
    // trim, joins the window, and scores 0.999 — high BECAUSE it sits near the
    // line, not because it was measured. It also becomes the terminal anchor
    // the AUCinf tail is extrapolated from. Nothing in the fit statistics can
    // distinguish that from a genuine observation, so the engine reports it.
    const rc = blockById('R-C');
    const inputs = inputsFor('P4');
    const r = computeNca(inputs, rulesFor(rc));
    expect(r.status).toBe('ok');
    expect(r.provenance.lambdaZ!.adjRSquared).toBeGreaterThan(0.99);

    // The window genuinely contains a substituted (BLQ-flagged, not dropped) point.
    const aug = augmentProfile(inputs, rulesFor(rc).blq)!;
    const substituted = Array.from(r.provenance.lambdaZ!.pointsUsed)
      .filter((i) => aug.blqMask[i] !== 0 && aug.dropMask[i] === 0);
    expect(substituted.length).toBe(1);
    expect(aug.time[substituted[0]]).toBe(24);
    expect(aug.conc[substituted[0]]).toBe(LLOQ / 2);

    const w = r.provenance.warnings.find((x) => x.code === 'LAMBDAZ_SUBSTITUTED_BLQ');
    expect(w).toBeDefined();
    expect(w!.severity).toBe('warning');

    // P1 under the same rule fits only measured points → no warning. The signal
    // tracks the fit's contents, not merely "this profile had a BLQ somewhere".
    const p1 = computeNca(inputsFor('P1'), rulesFor(rc));
    expect(p1.provenance.warnings.some((x) => x.code === 'LAMBDAZ_SUBSTITUTED_BLQ')).toBe(false);
  });

  it('I1: the flagged pre-dose row is replaced under EVERY rule; identical results; c0 = PKNCA\'s c0', () => {
    const ref = computeNca(inputsFor('I1'), rulesFor(blockById('R-A')));
    for (const b of fixture.blocks) {
      const r = computeNca(inputsFor('I1'), rulesFor(b));
      expect(r.provenance.c0!.replacedDoseTimeRow).toBe(true);
      expect(r.provenance.c0!.method).toBe('logslope');
      expect(relErr(r.provenance.c0!.value, profileOf(b, 'I1').provenance!.c0_pknca!)).toBeLessThan(TOL.auclast);
      expect(everyField(r.values)).toEqual(everyField(ref.values));
      // Never the substitute: under set-half-lloq the pre-dose row would read 0.025.
      expect(r.provenance.c0!.value).not.toBe(LLOQ / 2);
    }
  });
});

/**
 * Byte-identity pin for the 27 committed dense profiles (01/02/03/04) — a
 * PERMANENT gate, not a one-cycle pin.
 *
 * The tolerance assertions above prove PARITY with PKNCA (≤ 0.1 % / 0.5 % /
 * 1 %); they cannot distinguish "unchanged" from "changed within the gate". A
 * numbers-changing fix needs the exact pin so its delta is attributable: the
 * snapshot was committed from UNCHANGED code (GROK-20960, slice U0) and every
 * later change to any of these profiles fails this test. jest serialises the
 * numbers exactly, so a 1e-12 drift is a failure, not noise.
 *
 * Lifecycle: a future fix that deliberately moves one of these profiles must
 * update the snapshot (`jest -u`) in ITS OWN commit, with the CHANGELOG naming
 * the moved profiles — the same attribution discipline this pin was created
 * for. An unexplained snapshot diff in a review is a defect.
 */
describe('computeNca output snapshot @0.10.0 — byte-identity pin for numbers-changing fixes', () => {
  // The same table the span-ratio block drives: 27 profiles across all four
  // dense fixtures, so a fifth fixture cannot be silently left unpinned.
  const SNAPSHOT_FIXTURES: ReadonlyArray<{
    csv: string; timeCol: string; route: RouteCode;
    dose: (row: Record<string, number>) => number; tInf?: number;
  }> = [
    {csv: '01_theoph.csv', timeCol: 'Time', route: ROUTE_PO, dose: (r) => r.Dose * r.Wt},
    {csv: '02_indometh.csv', timeCol: 'time', route: ROUTE_IV_BOLUS, dose: () => 25},
    {csv: '03_rat_simple.csv', timeCol: 'Time', route: ROUTE_PO, dose: (r) => r.Dose},
    {csv: '04_iv_infusion.csv', timeCol: 'time', route: ROUTE_IV_INFUSION, dose: () => 100, tInf: 1},
  ];

  it.each([false, true])('compensatedSummation=%s — every profile matches the committed snapshot', (compensated) => {
    const rules: NcaRules = {...PKNCA_RULES, compensatedSummation: compensated};
    const out: Record<string, unknown> = {};
    let n = 0;
    for (const spec of SNAPSHOT_FIXTURES) {
      for (const group of loadAndGroup(spec.csv)) {
        const inputs = buildInputs(group.rows, spec.timeCol,
          spec.dose(group.rows[0]), spec.route, spec.tInf ?? null);
        const r = computeNca(inputs, rules);
        const lz = r.provenance.lambdaZ;
        out[`${spec.csv}/${group.subject}`] = {
          status: r.status,
          values: r.values,
          lambdaZ: lz === null ? null : {
            lambdaZ: lz.lambdaZ,
            pointsUsed: Array.from(lz.pointsUsed),
            tStart: lz.tStart,
            tEnd: lz.tEnd,
            spanRatio: lz.spanRatio,
          },
        };
        n++;
      }
    }
    expect(n).toBe(27);
    expect(out).toMatchSnapshot();
  });
});
