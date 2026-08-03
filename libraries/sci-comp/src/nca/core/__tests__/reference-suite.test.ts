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

function assertProfile(
  subject: string, dataset: string,
  inputs: ProfileInputs, fx: FixtureProfile,
) {
  const r = computeNca(inputs, PKNCA_RULES);
  const p = fx.parameters;
  const tag = `${dataset} subject ${subject}`;

  if (p.cmax !== null)
    expect(relErr(r.values.cmax, p.cmax)).toBeLessThan(TOL.cmax);

  if (p.tmax !== null) {
    // Tmax — exact match (per §9.2)
    expect(r.values.tmax).toBe(p.tmax);
  }
  if (p.auclast !== null)
    expect(relErr(r.values.aucLast, p.auclast)).toBeLessThan(TOL.auclast);

  if (p.aucinf !== null)
    expect(relErr(r.values.aucInf, p.aucinf)).toBeLessThan(TOL.aucinf);

  if (p.lambda_z !== null)
    expect(relErr(r.values.lambdaZ, p.lambda_z)).toBeLessThan(TOL.lambdaZ);

  if (p.half_life !== null)
    expect(relErr(r.values.halfLife, p.half_life)).toBeLessThan(TOL.halfLife);

  if (p.cl !== null)
    expect(relErr(r.values.cl, p.cl)).toBeLessThan(TOL.cl);

  if (p.vz !== null)
    expect(relErr(r.values.vz, p.vz)).toBeLessThan(TOL.vz);

  if (p.pct_aucextrap !== null) {
    expect(Math.abs(r.values.pctExtrap - p.pct_aucextrap))
      .toBeLessThan(TOL.pctExtrap);
  }

  // ── FR-200 moment / lag parameters ──────────────────────────────────────
  if (p.aumclast != null)
    expect(relErr(r.values.aumcLast, p.aumclast)).toBeLessThan(TOL.aumcLast);

  if (p.aumcinf_obs != null)
    expect(relErr(r.values.aumcInf, p.aumcinf_obs)).toBeLessThan(TOL.aumcInf);

  if (p.mrt != null)
    expect(relErr(r.values.mrt, p.mrt)).toBeLessThan(TOL.mrt);

  // Vss: gate-first. null (extravascular) ⇒ assert the route gate (NaN);
  // a number (IV) ⇒ assert parity. `== null` also covers a missing key.
  if (p.vss == null)
    expect(Number.isNaN(r.values.vss)).toBe(true);
  else
    expect(relErr(r.values.vss, p.vss)).toBeLessThan(TOL.vss);

  // Tlag: gate-first. null (IV) ⇒ assert the route gate (NaN); a number
  // (extravascular) ⇒ assert the exact observed time.
  if (p.tlag == null)
    expect(Number.isNaN(r.values.tlag)).toBe(true);
  else
    expect(r.values.tlag).toBe(p.tlag);

  if (p.pct_aumcextrap != null) {
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

