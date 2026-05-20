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
import {ROUTE_IV_BOLUS, ROUTE_PO} from '../types';
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
  dose: number, route: RouteCode,
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
    infusionDuration: null,
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
}

const TOL: Tolerances = {
  cmax: 0.001,
  auclast: 0.001,
  aucinf: 0.001,
  lambdaZ: 0.005,
  halfLife: 0.01,
  cl: 0.01,
  vz: 0.01,
  pctExtrap: 0.5,
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
});
