/**
 * `aucInterval` — PKNCA `aucint.last` parity (rule 18) + contract cases.
 *
 * Oracle: `fixtures/08_aucint.json`, written by `regen-aucint-fixture.R` (a
 * real `pk.nca` run per profile × interval × method — never a transcribed
 * literal). Every row of it is asserted, over every profile of 01 theoph,
 * 02 indometh and 03 rat_simple, in both summation modes, within 0.1 %
 * relative (the AUC class of NFR-04).
 */

import {readFileSync} from 'fs';
import {join} from 'path';
import {aucInterval} from '../auc-interval';
import {computeNca} from '../compute-nca';
import {ROUTE_IV_BOLUS, ROUTE_PO} from '../types';
import type {AucMethod, NcaRules, ProfileInputs, RouteCode} from '../types';

const DATASETS_DIR = join(__dirname, '..', '..', '__tests__', 'datasets');
const FIXTURES_DIR = join(__dirname, '..', '..', '__tests__', 'fixtures');

const AUC_REL_TOL = 0.001;

function rules(aucMethod: AucMethod, compensatedSummation: boolean): NcaRules {
  return {
    aucMethod,
    blq: {
      preFirstMeasurable: 'set-zero', embedded: 'set-zero',
      afterLast: 'set-zero', consecutiveAfterLast: 'set-zero',
    },
    lambdaZ: {mode: 'auto-best-fit', minPoints: 3, minRSquared: 0.85, excludeCmax: true, adjRSquaredFactor: 1e-4},
    extrapWarnPct: 20, extrapErrorPct: 50, extrapWarnPctAumc: 20,
    compensatedSummation,
  };
}

function profile(
  time: number[], conc: number[], route: RouteCode, blq?: number[],
): ProfileInputs {
  return {
    time: Float64Array.from(time), conc: Float64Array.from(conc),
    blqMask: blq ? Uint8Array.from(blq) : new Uint8Array(time.length),
    lloq: 0.05, dose: 1, doseUnits: 'mg', concentrationUnits: 'mg/L', timeUnits: 'h',
    route, infusionDuration: null, bodyWeight: null,
  };
}

/** Per-subject inputs of a committed dataset, built as the reference suite builds them. */
function loadDataset(
  csv: string, timeCol: string, route: RouteCode, dose: (r: Record<string, number>) => number,
): Map<string, ProfileInputs> {
  const lines = readFileSync(join(DATASETS_DIR, csv), 'utf-8')
    .split(/\r?\n/).filter((l) => l.trim().length > 0);
  const unq = (s: string) => s.replace(/^"|"$/g, '');
  const headers = lines[0].split(',').map(unq);
  const grouped = new Map<string, Record<string, number>[]>();
  for (const line of lines.slice(1)) {
    const cells = line.split(',').map(unq);
    const subject = cells[headers.indexOf('Subject')];
    const row: Record<string, number> = {};
    headers.forEach((h, i) => {
      if (h !== 'Subject') row[h] = parseFloat(cells[i]);
    });
    if (!grouped.has(subject)) grouped.set(subject, []);
    grouped.get(subject)!.push(row);
  }
  const out = new Map<string, ProfileInputs>();
  for (const [subject, rows] of grouped) {
    out.set(subject, {
      ...profile(rows.map((r) => r[timeCol]), rows.map((r) => r.conc), route),
      lloq: 0.001, dose: dose(rows[0]),
    });
  }
  return out;
}

const DATASETS: Record<string, Map<string, ProfileInputs>> = {
  '01_theoph': loadDataset('01_theoph.csv', 'Time', ROUTE_PO, (r) => r.Dose * r.Wt),
  '02_indometh': loadDataset('02_indometh.csv', 'time', ROUTE_IV_BOLUS, () => 25),
  '03_rat_simple': loadDataset('03_rat_simple.csv', 'Time', ROUTE_PO, (r) => r.Dose),
};

interface OracleRow {
  dataset: string;
  subject: string;
  auc_method: AucMethod;
  start: number;
  end: number;
  aucint_last: number;
}

const ORACLE: {pknca_version: string; pknca_parameter: string; rows: OracleRow[]} =
  JSON.parse(readFileSync(join(FIXTURES_DIR, '08_aucint.json'), 'utf-8'));

describe('aucInterval vs PKNCA aucint.last (08_aucint.json)', () => {
  it('the oracle is PKNCA 0.12.1 aucint.last and covers every profile of 01/02/03', () => {
    expect(ORACLE.pknca_version).toBe('0.12.1');
    expect(ORACLE.pknca_parameter).toBe('aucint.last');
    for (const [name, subjects] of Object.entries(DATASETS)) {
      const covered = new Set(ORACLE.rows.filter((r) => r.dataset === name).map((r) => r.subject));
      expect([...covered].sort()).toEqual([...subjects.keys()].sort());
    }
    // Both shared methods, and the AC-4.3 windows [0, 12] and [0, 24], for every profile.
    for (const [name, subjects] of Object.entries(DATASETS)) {
      for (const subject of subjects.keys()) {
        for (const m of ['linear-up-log-down', 'linear']) {
          for (const [s, e] of [[0, 12], [0, 24]]) {
            expect(ORACLE.rows.some((r) => r.dataset === name && r.subject === subject &&
              r.auc_method === m && r.start === s && r.end === e)).toBe(true);
          }
        }
      }
    }
  });

  describe.each([false, true])('compensatedSummation = %s', (compensated) => {
    it.each(ORACLE.rows.map((r) => [
      `${r.dataset} subj ${r.subject} ${r.auc_method} [${r.start}, ${r.end}]`, r,
    ] as const))('%s', (_label, row) => {
      const inputs = DATASETS[row.dataset].get(row.subject)!;
      const got = aucInterval(inputs, rules(row.auc_method, compensated), row.start, row.end);
      expect(Number.isFinite(got.value)).toBe(true);
      expect(Math.abs(got.value - row.aucint_last) / Math.abs(row.aucint_last))
        .toBeLessThan(AUC_REL_TOL);
    });
  });
});

describe('aucInterval contract', () => {
  it('[0, Infinity] is AUClast bit for bit, for every fixture profile, method and summation', () => {
    for (const subjects of Object.values(DATASETS)) {
      for (const inputs of subjects.values()) {
        for (const m of ['linear', 'log-linear', 'linear-up-log-down'] as const) {
          for (const comp of [false, true]) {
            const r = rules(m, comp);
            const got = aucInterval(inputs, r, 0, Infinity);
            expect(got.value).toBe(computeNca(inputs, r).values.aucLast);
            expect(got.truncatedAtTlast).toBe(true);
          }
        }
      }
    }
  });

  it('log-interpolates a bound on a descending interval: exact on a mono-exponential', () => {
    // c = 10·e^{-0.1 t}; lin-up/log-down integrates and interpolates it exactly.
    const t = [0, 4, 8, 12];
    const p = profile(t, t.map((x) => 10 * Math.exp(-0.1 * x)), ROUTE_PO);
    const got = aucInterval(p, rules('linear-up-log-down', false), 1, 6);
    const exact = 100 * (Math.exp(-0.1) - Math.exp(-0.6));
    expect(got.value).toBeCloseTo(exact, 12);
    expect(got.startUsed).toBe(1);
    expect(got.endUsed).toBe(6);
    expect(got.truncatedAtTlast).toBe(false);
  });

  it('linearly interpolates a bound on an ascending interval', () => {
    // (0,0) → (2,4): c(1) = 2, area over [0, 1] = 1.
    const p = profile([0, 2, 4], [0, 4, 2], ROUTE_PO);
    expect(aucInterval(p, rules('linear-up-log-down', false), 0, 1).value).toBeCloseTo(1, 14);
  });

  it('does not extrapolate past tlast: end beyond tlast integrates to tlast and says so', () => {
    const p = profile([0, 1, 2, 4], [0, 3, 2, 1], ROUTE_PO);
    const r = rules('linear', false);
    const got = aucInterval(p, r, 0, 24);
    expect(got.value).toBe(computeNca(p, r).values.aucLast);
    expect(got.endUsed).toBe(4);
    expect(got.tLast).toBe(4);
    expect(got.truncatedAtTlast).toBe(true);
  });

  it('a trailing zero is not tlast (trimmed as computeNca does)', () => {
    const p = profile([0, 1, 2, 4, 8], [0, 3, 2, 1, 0], ROUTE_PO);
    const got = aucInterval(p, rules('linear', false), 0, 8);
    expect(got.tLast).toBe(4);
    expect(got.truncatedAtTlast).toBe(true);
  });

  it('an interval wholly past tlast is 0; a degenerate interval is 0', () => {
    const p = profile([0, 1, 2, 4], [0, 3, 2, 1], ROUTE_PO);
    expect(aucInterval(p, rules('linear', false), 6, 12).value).toBe(0);
    expect(aucInterval(p, rules('linear', false), 1.5, 1.5).value).toBe(0);
  });

  it('a start before the dose time adds no area', () => {
    const p = profile([0, 1, 2, 4], [0, 3, 2, 1], ROUTE_PO);
    const r = rules('linear', false);
    const got = aucInterval(p, r, -5, 2);
    expect(got.startUsed).toBe(0);
    expect(got.value).toBe(aucInterval(p, r, 0, 2).value);
  });

  it('integrates the IV-bolus profile from the back-extrapolated c0', () => {
    // No t = 0 sample: (0, c0) is inserted, so [0, 1] includes the dose-time segment.
    const p = profile([1, 2, 4], [8, 4, 1], ROUTE_IV_BOLUS);
    const r = rules('linear-up-log-down', false);
    const res = computeNca(p, r);
    expect(res.provenance.c0!.method).toBe('logslope');
    const c0 = res.provenance.c0!.value; // 16
    const segment = (c0 - 8) * 1 / Math.log(c0 / 8);
    expect(aucInterval(p, r, 0, 1).value).toBeCloseTo(segment, 12);
  });

  it('honours the BLQ rule exactly as computeNca does', () => {
    const t = [0, 1, 2, 3, 4, 6];
    const c = [0, 4, 0, 3, 2, 1];
    const blq = [0, 0, 1, 0, 0, 0];
    const p = profile(t, c, ROUTE_PO, blq);
    for (const embedded of ['set-zero', 'set-half-lloq', 'exclude', 'missing'] as const) {
      const r: NcaRules = {...rules('linear', false), blq: {...rules('linear', false).blq, embedded}};
      expect(aucInterval(p, r, 0, Infinity).value).toBe(computeNca(p, r).values.aucLast);
    }
    const zero = aucInterval(p, rules('linear', false), 1, 3).value; // through the 0 at t = 2
    const excl = aucInterval(p, {...rules('linear', false),
      blq: {...rules('linear', false).blq, embedded: 'exclude'}}, 1, 3).value; // straight 4 → 3
    expect(zero).toBeCloseTo(3.5, 14);
    expect(excl).toBeCloseTo(7, 14);
  });

  it('returns NaN on an all-BLQ profile and on an invalid interval', () => {
    const blqAll = profile([0, 1, 2], [0, 0, 0], ROUTE_PO, [1, 1, 1]);
    expect(Number.isNaN(aucInterval(blqAll, rules('linear', false), 0, 2).value)).toBe(true);
    const p = profile([0, 1, 2], [0, 2, 1], ROUTE_PO);
    expect(Number.isNaN(aucInterval(p, rules('linear', false), 2, 1).value)).toBe(true);
    expect(Number.isNaN(aucInterval(p, rules('linear', false), NaN, 1).value)).toBe(true);
    expect(Number.isNaN(aucInterval(p, rules('linear', false), -Infinity, 1).value)).toBe(true);
    expect(Number.isNaN(aucInterval(p, rules('linear', false), 0, NaN).value)).toBe(true);
  });
});
