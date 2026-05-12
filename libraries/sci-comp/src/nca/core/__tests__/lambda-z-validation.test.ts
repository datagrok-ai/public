/**
 * Early validation gate for the lambda_z best-fit algorithm (Task 1.7.5).
 *
 * Loads the Phase 0 reference datasets (theoph, indometh) and the matching
 * PKNCA fixtures, runs `lambdaZBestFit` on each subject, and asserts the
 * result agrees with PKNCA within 0.5% (relative).
 *
 * Theoph and Indometh contain no BLQ points, so the BLQ pre-processing
 * pipeline is exercised as a no-op upstream and we run lambda_z directly
 * on the raw concentrations.
 *
 * If this gate fails, do NOT proceed to Task 1.8 — investigate the
 * algorithmic divergence (Risk Register R1).
 */

import {readFileSync} from 'fs';
import {join} from 'path';
import {findCmax} from '../cmax';
import {lambdaZBestFit} from '../lambda-z';
import {insertC0} from '../c0';
import type {LambdaZStrategy} from '../types';

const PKNCA_STRATEGY: LambdaZStrategy = {
  mode: 'auto-best-fit',
  minPoints: 3,
  minRSquared: 0.85,
  excludeCmax: true,
  adjRSquaredFactor: 1e-4, // PKNCA default tie-breaking
};

/** Tolerance for lambda_z comparison: 0.5% (TC-004 in PRD). */
const LAMBDA_Z_REL_TOL = 0.005;

interface FixtureProfile {
  profile_key: {subject: string; route: string};
  parameters: {lambda_z: number | null};
  provenance: {
    lambda_z_n_points: number | null;
    lambda_z_time_first: number | null;
  };
}

interface FixtureFile {
  dataset: string;
  profiles: FixtureProfile[];
}

const DATASETS_DIR = join(__dirname, '..', '..', '__tests__', 'datasets');
const FIXTURES_DIR = join(__dirname, '..', '..', '__tests__', 'fixtures');

/** Parse a write.csv() output: comma-separated, double-quoted strings. */
function parseCsv(text: string): {headers: string[]; rows: string[][]} {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0);
  const stripQuotes = (s: string) => s.replace(/^"|"$/g, '');
  const headers = lines[0].split(',').map(stripQuotes);
  const rows = lines.slice(1).map((l) => l.split(',').map(stripQuotes));
  return {headers, rows};
}

interface SubjectProfile {
  subject: string;
  time: Float64Array;
  conc: Float64Array;
}

/** Group CSV rows by subject and build typed arrays of (time, conc). */
function loadSubjectProfiles(
  csvName: string, timeCol: string,
): SubjectProfile[] {
  const text = readFileSync(join(DATASETS_DIR, csvName), 'utf-8');
  const {headers, rows} = parseCsv(text);
  const subjCol = headers.indexOf('Subject');
  const tIdx = headers.indexOf(timeCol);
  const cIdx = headers.indexOf('conc');
  if (subjCol < 0 || tIdx < 0 || cIdx < 0)
    throw new Error(`bad CSV headers in ${csvName}: ${headers.join(',')}`);

  const grouped = new Map<string, {time: number[]; conc: number[]}>();
  for (const r of rows) {
    const s = r[subjCol];
    if (!grouped.has(s)) grouped.set(s, {time: [], conc: []});
    const g = grouped.get(s)!;
    g.time.push(parseFloat(r[tIdx]));
    g.conc.push(parseFloat(r[cIdx]));
  }
  return Array.from(grouped.entries()).map(([subject, g]) => ({
    subject,
    time: Float64Array.from(g.time),
    conc: Float64Array.from(g.conc),
  }));
}

function loadFixture(fixtureName: string): FixtureFile {
  return JSON.parse(readFileSync(join(FIXTURES_DIR, fixtureName), 'utf-8'));
}

/**
 * Run lambda_z best-fit on one subject. For IV bolus we mirror the PKNCA
 * pipeline by inserting an estimated c0 at t=0 (which becomes the new Cmax),
 * then running the fit on the augmented profile. For oral / extravascular
 * the raw data already contains a t=0 point (or none is needed).
 */
function runLambdaZ(p: SubjectProfile, route: string) {
  let time = p.time;
  let conc = p.conc;
  let blqMask = new Uint8Array(p.conc.length);
  let cmaxIdx: number;

  if (route === 'iv-bolus' && !Array.from(time).includes(0)) {
    const aug = insertC0(time, conc, blqMask);
    if (aug === null) return null;
    time = aug.time;
    conc = aug.conc;
    blqMask = aug.blqMask;
    cmaxIdx = aug.cmaxIdx;
  } else {
    const cmaxR = findCmax(time, conc, blqMask);
    if (cmaxR === null) return null;
    cmaxIdx = cmaxR.cmaxIdx;
  }
  return lambdaZBestFit(time, conc, blqMask, cmaxIdx, PKNCA_STRATEGY);
}

describe('lambda_z validation gate (Task 1.7.5)', () => {
  describe('Theoph (12 subjects, oral, no BLQ)', () => {
    const profiles = loadSubjectProfiles('01_theoph.csv', 'Time');
    const fixture = loadFixture('01_theoph.json');
    const fxBySubject = new Map(
      fixture.profiles.map((f) => [f.profile_key.subject, f]));

    it.each(profiles.map((p) => [p.subject, p] as const))(
      'subject %s — lambda_z within 0.5%% of PKNCA',
      (subjectId, profile) => {
        const fx = fxBySubject.get(subjectId);
        expect(fx).toBeDefined();

        const r = runLambdaZ(profile, 'extravascular');
        const expected = fx!.parameters.lambda_z;

        if (expected === null) {
          expect(r).toBeNull();
          return;
        }

        expect(r).not.toBeNull();
        const relErr = Math.abs(r!.lambdaZ - expected) / expected;
        expect(relErr).toBeLessThan(LAMBDA_Z_REL_TOL);
      });
  });

  describe('Indomethacin (IV bolus, no BLQ — full pipeline with c0 insertion)', () => {
    // PKNCA's `pk.calc.c0` (default chain: c0 → logslope → c1 → cmin → set0)
    // is reproduced in `c0.ts`. For each Indometh subject we estimate c0,
    // insert (0, c0), and run lambda_z on the augmented profile — matching
    // the reference fixture, which was produced the same way.
    const profiles = loadSubjectProfiles('02_indometh.csv', 'time');
    const fixture = loadFixture('02_indometh.json');
    const fxBySubject = new Map(
      fixture.profiles.map((f) => [f.profile_key.subject, f]));

    it.each(profiles.map((p) => [p.subject, p] as const))(
      'subject %s — lambda_z within 0.5%% of PKNCA',
      (subjectId, profile) => {
        const fx = fxBySubject.get(subjectId);
        expect(fx).toBeDefined();

        const r = runLambdaZ(profile, 'iv-bolus');
        const expected = fx!.parameters.lambda_z;

        if (expected === null) {
          expect(r).toBeNull();
          return;
        }

        expect(r).not.toBeNull();
        const relErr = Math.abs(r!.lambdaZ - expected) / expected;
        expect(relErr).toBeLessThan(LAMBDA_Z_REL_TOL);
      });
  });
});
