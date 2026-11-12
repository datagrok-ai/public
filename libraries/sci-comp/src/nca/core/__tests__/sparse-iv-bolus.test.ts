/**
 * IV-bolus sparse composite: the C0 anchor and its delta-method SE (sparse-iv-bolus-anchor).
 *
 * Oracles (none of them is the code under test re-implemented):
 *  - **Point estimate** (AC-2.2): PKNCA 0.12.1 `pk.calc.c0` + `pk.calc.auciv` on the composite
 *    mean, `fixtures/07_iv_bolus_sparse.json` (`regen-sparse-fixture.R` block 4). Arithmetic
 *    parity on the same Ĉ0 — not evidence of accuracy, nor of Phoenix parity.
 *  - **SE / df** (AC-2.3 a): a finite-difference gradient of the REPORTED `auc` with respect to
 *    each timepoint mean, with the variances recomputed here from the raw data; on the `c1`
 *    composite the gradient is exact (g₁ = w₁ + t₁/2). A small batch design checks
 *    se² = gᵀΣg with Σ (Holder A3) recomputed here.
 *  - **Calibration** (AC-2.3 b): a seeded parametric simulation through the real `sparseAuc`,
 *    with pre-stated tolerances (band [0.85, 1.25]; first segment ±5 / ±15 / ±10 %).
 *  - **Closed form** (AC-2.3 c): a noise-free mono-exponential, where Ĉ0 is exact.
 *  - **Transcription**: the hand-coded delta SE / df of the R script (`delta_se`, `delta_df`).
 */

import {sparseAuc} from '../sparse';
import type {SparseInput, SparseAucResult} from '../sparse';
import {ROUTE_IV_BOLUS, ROUTE_PO} from '../types';
import {mulberry32} from '../prng';
import {studentTInv} from '../../../stats/distributions';
import {ivInputs, mouseInputs} from './sparse-fixture-io';

const IV = {blqRule: 'set-zero' as const, route: ROUTE_IV_BOLUS};
const rel = (a: number, b: number): number => Math.abs(a - b) / Math.abs(b);

/** Observations grouped by nominal time: per-time values, mean, variance (n−1), n. */
function byTime(input: SparseInput): Array<{t: number; idx: number[]; mean: number; s2: number; r: number}> {
  const m = new Map<number, number[]>();
  input.nominalTime.forEach((t, i) => {const l = m.get(t) ?? []; l.push(i); m.set(t, l);});
  return [...m.entries()].sort((a, b) => a[0] - b[0]).map(([t, idx]) => {
    const x = idx.map((i) => input.conc[i]);
    const mean = x.reduce((a, b) => a + b, 0) / x.length;
    const s2 = x.length > 1 ? x.reduce((a, b) => a + (b - mean) ** 2, 0) / (x.length - 1) : 0;
    return {t, idx, mean, s2, r: x.length};
  });
}

/** Finite-difference gradient of the reported AUC w.r.t. each timepoint mean: every observation
 *  at that time is shifted by ±δ (the mean moves by δ, the variance does not). */
function aucGradient(input: SparseInput, opts: object): number[] {
  return byTime(input).map((g) => {
    const d = 1e-6 * Math.max(g.mean, 1e-12);
    const shifted = (s: number): number => {
      const conc = Float64Array.from(input.conc);
      for (const i of g.idx) conc[i] += s;
      return sparseAuc({...input, conc}, opts).auc;
    };
    return (shifted(d) - shifted(-d)) / (2 * d);
  });
}

/** Lever-arm weights over [0, ...times] (the anchor's own weight dropped). */
function leverWeights(times: number[]): number[] {
  const T = [0, ...times];
  return times.map((_, k) => {
    const i = k + 1;
    return (T[i] - T[i - 1]) / 2 + (i < T.length - 1 ? (T[i + 1] - T[i]) / 2 : 0);
  });
}

describe('sparseAuc IV bolus — point estimate vs PKNCA (AC-2.2)', () => {
  for (const {block, input} of ivInputs()) {
    it(`${block.dose_mg_kg} mg/kg: auc = PKNCA auciv on Ĉ0, c0 = pk.calc.c0 (${block.c0_method})`, () => {
      const res = sparseAuc(input, IV);
      expect(rel(res.auc, block.sparse_auclast_c0)).toBeLessThan(1e-3); // NFR-04
      expect(rel(res.auc, block.sparse_auclast_c0)).toBeLessThan(1e-10); // same inputs: round-off
      expect(res.c0).toBeDefined();
      expect(rel(res.c0!.value, block.c0)).toBeLessThan(1e-10);
      expect(res.c0!.method).toBe(block.c0_method);
      // The (0, 0) baseline is what the no-route call still returns.
      expect(rel(sparseAuc(input, {blqRule: 'set-zero'}).auc, block.sparse_auclast)).toBeLessThan(1e-10);
    });
  }

  it('route omitted / extravascular: no c0 on the result', () => {
    const {input} = ivInputs()[0];
    expect(sparseAuc(input).c0).toBeUndefined();
    expect(sparseAuc(input, {route: ROUTE_PO}).c0).toBeUndefined();
  });
});

describe('sparseAuc IV bolus — delta SE / df (AC-2.3, AC-2.4)', () => {
  for (const {block, input} of ivInputs()) {
    it(`${block.dose_mg_kg} mg/kg: se² = Σ g² s²/r and df = Satterthwaite(g), g from the reported auc`, () => {
      const res = sparseAuc(input, {...IV, varianceBorrow: false});
      const g = aucGradient(input, {...IV, varianceBorrow: false});
      const tps = byTime(input);
      const terms = tps.map((tp, i) => g[i] * g[i] * tp.s2 / tp.r);
      const psi = terms.reduce((a, b) => a + b, 0);
      const df = psi * psi / terms.reduce((a, b, i) => a + b * b / Math.max(tps[i].r - 1, 1), 0);
      expect(rel(res.se, Math.sqrt(psi))).toBeLessThan(1e-6);
      expect(rel(res.df, df)).toBeLessThan(1e-6);
    });

    it(`${block.dose_mg_kg} mg/kg: se / df = the R hand-coded delta (transcription check)`, () => {
      const res = sparseAuc(input, IV);
      expect(rel(res.se, block.delta_se)).toBeLessThan(1e-9);
      expect(rel(res.df, block.delta_df)).toBeLessThan(1e-9);
    });

    it(`${block.dose_mg_kg} mg/kg: df within [min(r−1), Σ(r−1)] (AC-2.4)`, () => {
      const res = sparseAuc(input, IV);
      const tps = byTime(input);
      expect(res.df).toBeGreaterThanOrEqual(Math.min(...tps.map((t) => t.r - 1)) - 1e-9);
      expect(res.df).toBeLessThanOrEqual(tps.reduce((a, t) => a + t.r - 1, 0) + 1e-9);
    });
  }

  it('c1 composite: the gradient is exact — g₁ = w₁ + t₁/2, the other weights unchanged', () => {
    const {input} = ivInputs().find(({block}) => block.c0_method === 'c1')!;
    const tps = byTime(input);
    const w = leverWeights(tps.map((t) => t.t));
    const g = aucGradient(input, IV);
    expect(Math.abs(g[0] - (w[0] + tps[0].t / 2))).toBeLessThan(1e-6);
    for (let i = 1; i < g.length; i++) expect(Math.abs(g[i] - w[i])).toBeLessThan(1e-6);
  });

  it('batch IV design: se² = gᵀΣg with Σ (Holder A3) recomputed here; df finite (N-1)', () => {
    // 6 animals, each sampled at 2 of 3 times (r_i = 4, r_ij = 2): a decaying IV profile.
    const plan: Array<[number, string, number]> = [
      [0.25, 'A', 10.0], [0.25, 'B', 12.0], [0.25, 'C', 11.0], [0.25, 'F', 9.5],
      [1, 'A', 6.0], [1, 'B', 7.0], [1, 'D', 8.0], [1, 'E', 9.0],
      [4, 'C', 2.0], [4, 'D', 3.0], [4, 'E', 2.5], [4, 'F', 1.5],
    ];
    const ids = [...new Set(plan.map((p) => p[1]))];
    const input: SparseInput = {
      nominalTime: Float64Array.from(plan.map((p) => p[0])), conc: Float64Array.from(plan.map((p) => p[2])),
      blqMask: new Uint8Array(plan.length), lloq: 0,
      animalId: Int32Array.from(plan.map((p) => ids.indexOf(p[1]) + 1)),
    };
    const res = sparseAuc(input, IV);
    expect(res.topology).toBe('batch');
    const g = aucGradient(input, IV);
    const times = [0.25, 1, 4];
    const at = (t: number): Map<string, number> => new Map(plan.filter((p) => p[0] === t).map((p) => [p[1], p[2]]));
    const sets = times.map(at);
    const mean = sets.map((s) => [...s.values()].reduce((a, b) => a + b, 0) / s.size);
    const r = sets.map((s) => s.size);
    const Sigma = times.map((_, i) => times.map((__, j) => {
      if (i === j) return [...sets[i].values()].reduce((a, x) => a + (x - mean[i]) ** 2, 0) / (r[i] - 1) / r[i];
      const both = [...sets[i].keys()].filter((k) => sets[j].has(k));
      const rij = both.length;
      const sum = both.reduce((a, k) => a + (sets[i].get(k)! - mean[i]) * (sets[j].get(k)! - mean[j]), 0);
      const sij = sum / ((rij - 1) + (1 - rij / r[i]) * (1 - rij / r[j]));
      return rij * sij / (r[i] * r[j]);
    }));
    let psi = 0;
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) psi += g[i] * g[j] * Sigma[i][j];
    expect(rel(res.se, Math.sqrt(psi))).toBeLessThan(1e-6);
    expect(Number.isFinite(res.df) && res.df > 0).toBe(true);
  });
});

describe('sparseAuc IV bolus — noise-free closed form (AC-2.3 c)', () => {
  for (const t1 of [5 / 60, 0.25]) {
    it(`C = e^(−t), t1 = ${(t1 * 60).toFixed(0)} min: Ĉ0 = 1, first segment ½·t1·(1 + e^(−t1))`, () => {
      const times = [t1, 0.5, 1, 2, 4];
      const nominalTime: number[] = []; const conc: number[] = [];
      for (const t of times) for (let k = 0; k < 2; k++) {nominalTime.push(t); conc.push(Math.exp(-t));}
      const input: SparseInput = {
        nominalTime: Float64Array.from(nominalTime), conc: Float64Array.from(conc),
        blqMask: new Uint8Array(conc.length), lloq: 0, animalId: null,
      };
      const p1 = sparseAuc(input, IV);
      const zero = sparseAuc(input);
      expect(p1.c0!.method).toBe('logslope');
      expect(Math.abs(p1.c0!.value - 1)).toBeLessThan(1e-12);
      const c1 = Math.exp(-t1);
      const firstP1 = p1.auc - zero.auc + t1 * c1 / 2;
      expect(Math.abs(firstP1 - t1 * (1 + c1) / 2)).toBeLessThan(1e-12);
      // Why the anchor changed: (0, 0) under-states the first segment by 52.1 % / 56.0 %.
      const truth = 1 - c1;
      const zeroErr = (t1 * c1 / 2) / truth - 1;
      expect(zeroErr).toBeCloseTo(t1 < 0.1 ? -0.521 : -0.560, 3);
    });
  }
});

describe('sparseAuc IV bolus — dose-time rule and BLQ edges (AC-2.3, AC-2.5)', () => {
  const base = (): {t: number[]; c: number[]; b: number[]} => ({
    t: [0.25, 0.25, 0.25, 1, 1, 1, 4, 4, 4, 8, 8, 8],
    c: [10, 12, 11, 6, 7, 6.5, 2, 2.4, 2.2, 0.8, 0.9, 0.7],
    b: new Array(12).fill(0),
  });
  const mk = (d: {t: number[]; c: number[]; b: number[]}): SparseInput => ({
    nominalTime: Float64Array.from(d.t), conc: Float64Array.from(d.c), blqMask: Uint8Array.from(d.b),
    lloq: 0.05, animalId: Int32Array.from(d.t.map((_, i) => i + 1)),
  });
  const numbers = (r: SparseAucResult): number[] => [r.auc, r.se, r.df];

  it('a t = 0 BLQ group and a t < 0 group are not integrated; topology / timepoints keep them', () => {
    const plain = mk(base());
    const d = base();
    d.t.unshift(-0.5, -0.5, 0, 0, 0); d.c.unshift(0, 0, 0, 0, 0); d.b.unshift(1, 1, 1, 1, 1);
    const withPre = mk(d);
    const a = sparseAuc(plain, IV); const b = sparseAuc(withPre, IV);
    expect(numbers(b)).toEqual(numbers(a));
    expect(b.c0).toEqual(a.c0);
    expect(b.composite.timepoints.length).toBe(a.composite.timepoints.length + 2);
    expect(b.topology).toBe(a.topology);
  });

  it('a positive measured t = 0 group is the observed C0: today\'s numbers, method observed', () => {
    const d = base();
    d.t.unshift(0, 0, 0); d.c.unshift(15, 16, 14); d.b.unshift(0, 0, 0);
    const input = mk(d);
    const res = sparseAuc(input, IV);
    expect(numbers(res)).toEqual(numbers(sparseAuc(input, {blqRule: 'set-zero'})));
    expect(res.c0).toEqual({value: 15, method: 'observed'});
  });

  it('first post-dose time BLQ-majority with a positive mean: zeroed, no t1/2 added (R2-1)', () => {
    const d = base();
    d.c[0] = 0; d.b[0] = 1; d.c[1] = 0; d.b[1] = 1; // 2 of 3 BLQ at 0.25 h, mean of the set-zero pool > 0
    const input = mk(d);
    const iv = sparseAuc(input, IV);
    const zero = sparseAuc(input, {blqRule: 'set-zero'});
    expect(iv.c0).toEqual({value: 0, method: 'c1'});
    expect(numbers(iv)).toEqual(numbers(zero));
  });

  it('PO: a BLQ-majority timepoint with a positive mean is zeroed in AUClast (one copy of the rule, N-3)', () => {
    const d = base();
    d.c[9] = 0; d.b[9] = 1; d.c[10] = 0; d.b[10] = 1; // 2 of 3 BLQ at 8 h, the third 0.7
    const res = sparseAuc(mk(d), {blqRule: 'set-zero', route: ROUTE_PO});
    // Hand: linear trapezoid over (0,0),(0.25,11),(1,6.5),(4,2.2); the 8 h mean is zeroed → AUClast stops at 4 h.
    const hand = 0.25 * 11 / 2 + 0.75 * (11 + 6.5) / 2 + 3 * (6.5 + 2.2) / 2;
    expect(Math.abs(res.auc - hand)).toBeLessThan(1e-12);
  });
});

// --- Parametric calibration (AC-2.3 b) -----------------------------------------------------
// Seeded (mulberry32), through the shipped `sparseAuc`. Pre-stated tolerances (plan R1 C-2 /
// R2-1): empSD / mean SE in [0.85, 1.25]; P1 first-segment mean relative error vs the analytic
// first segment within ±5 % / ±15 % / ±10 %. Bias vs AUC(0–tlast), 95 % coverage (P1 and the
// (0, 0) baseline) and the chain-path fractions are REPORTED (console), not gated — METHOD-24
// records them.

interface Model {
  name: string; A: number; alpha: number; B: number; beta: number; times: number[];
  omega: number; n: number; firstSegTol?: number;
}
const MODELS: Model[] = [
  {name: '2-cmt 5 min', A: 8, alpha: 4, B: 2, beta: 0.3, times: [0.083, 0.25, 0.5, 1, 2, 4, 8, 12, 24],
    omega: 0.3, n: 3, firstSegTol: 0.05},
  {name: '2-cmt 15 min', A: 8, alpha: 4, B: 2, beta: 0.3, times: [0.25, 0.5, 1, 2, 4, 8, 24],
    omega: 0.3, n: 3, firstSegTol: 0.15},
  {name: 'chain-flip 1-cmt 15/30 min', A: 0, alpha: 1, B: 1, beta: 0.2, times: [0.25, 0.5, 1, 2, 4, 8, 12],
    omega: 0.3, n: 3, firstSegTol: 0.10},
];
const STEEP: Model = {name: 'steep α (limitation)', A: 20, alpha: 8, B: 1, beta: 0.2,
  times: [0.25, 0.5, 1, 2, 4, 8, 24], omega: 0.7, n: 2};
const SIGMA = 0.1;

function normalGen(rnd: () => number): () => number {
  return () => {
    const u = rnd() || 1e-12; const v = rnd();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  };
}

/** Analytic population means of ∫0^t C over the lognormal macro-constants (Monte Carlo, 4e5). */
function truthArea(m: Model, t: number, rnd: () => number): number {
  const z = normalGen(rnd); const e = (): number => Math.exp(m.omega * z() - m.omega * m.omega / 2);
  let s = 0; const M = 400000;
  for (let k = 0; k < M; k++) {
    const A = m.A * e(); const al = m.alpha * e(); const B = m.B * e(); const be = m.beta * e();
    s += A / al * (1 - Math.exp(-al * t)) + B / be * (1 - Math.exp(-be * t));
  }
  return s / M;
}

function calibrate(m: Model, reps: number, seed: number) {
  const rnd = mulberry32(seed); const z = normalGen(rnd);
  const e = (): number => Math.exp(m.omega * z() - m.omega * m.omega / 2);
  const t1 = m.times[0]; const tl = m.times[m.times.length - 1];
  const truthFirst = truthArea(m, t1, rnd); const truthLast = truthArea(m, tl, rnd);
  const auc: number[] = []; const se: number[] = []; const first: number[] = [];
  const paths = new Map<string, number>(); let covP1 = 0; let cov00 = 0; const auc00: number[] = [];
  for (let r = 0; r < reps; r++) {
    const nominalTime: number[] = []; const conc: number[] = [];
    for (const t of m.times) {
      for (let k = 0; k < m.n; k++) {
        const A = m.A * e(); const al = m.alpha * e(); const B = m.B * e(); const be = m.beta * e();
        nominalTime.push(t);
        conc.push((A * Math.exp(-al * t) + B * Math.exp(-be * t)) * Math.exp(SIGMA * z() - SIGMA * SIGMA / 2));
      }
    }
    const input: SparseInput = {nominalTime: Float64Array.from(nominalTime), conc: Float64Array.from(conc),
      blqMask: new Uint8Array(conc.length), lloq: 0, animalId: null};
    const p1 = sparseAuc(input, IV); const z0 = sparseAuc(input);
    const c1 = conc.slice(0, m.n).reduce((a, b) => a + b, 0) / m.n;
    auc.push(p1.auc); se.push(p1.se); auc00.push(z0.auc);
    first.push(p1.auc - z0.auc + t1 * c1 / 2);
    paths.set(p1.c0!.method, (paths.get(p1.c0!.method) ?? 0) + 1);
    if (Math.abs(p1.auc - truthLast) <= studentTInv(0.975, p1.df) * p1.se) covP1++;
    if (Math.abs(z0.auc - truthLast) <= studentTInv(0.975, z0.df) * z0.se) cov00++;
  }
  const mean = (x: number[]): number => x.reduce((a, b) => a + b, 0) / x.length;
  const sd = (x: number[]): number => {
    const mu = mean(x);
    return Math.sqrt(x.reduce((a, b) => a + (b - mu) ** 2, 0) / (x.length - 1));
  };
  return {
    ratio: sd(auc) / mean(se),
    firstSegErr: mean(first) / truthFirst - 1,
    biasP1: mean(auc) / truthLast - 1, bias00: mean(auc00) / truthLast - 1,
    coverageP1: covP1 / reps, coverage00: cov00 / reps,
    paths: Object.fromEntries([...paths].map(([k, v]) => [k, v / reps])),
  };
}

describe('sparseAuc IV bolus — parametric calibration (AC-2.3 b)', () => {
  const report: Record<string, unknown> = {};
  afterAll(() => console.log(`IV-bolus calibration (2000 replicates):\n${JSON.stringify(report, null, 1)}`));

  MODELS.forEach((m, k) => {
    it(`${m.name}: empSD / SE in [0.85, 1.25]; first segment within ±${(m.firstSegTol! * 100).toFixed(0)} %`, () => {
      const c = calibrate(m, 2000, 1000 + k);
      report[m.name] = c;
      expect(c.ratio).toBeGreaterThanOrEqual(0.85);
      expect(c.ratio).toBeLessThanOrEqual(1.25);
      expect(Math.abs(c.firstSegErr)).toBeLessThanOrEqual(m.firstSegTol!);
    });
  });

  it('chain-flip model exercises both chain paths', () => {
    const c = calibrate(MODELS[2], 2000, 1002);
    expect((c.paths['logslope'] ?? 0) > 0.05 && (c.paths['c1'] ?? 0) > 0.05).toBe(true);
  });

  it('documented limitation: steep α, ω 0.7, n 2 — the delta SE is small (empSD / SE > 1.15)', () => {
    const c = calibrate(STEEP, 2000, 1009);
    report[STEEP.name] = c;
    expect(c.ratio).toBeGreaterThan(1.15);
  });
});

describe('SPARSE_TERMINAL_OVEREST — the estimated excess (AC-5.5)', () => {
  const mk = (t: number[], c: number[]): SparseInput => ({
    nominalTime: Float64Array.from(t.flatMap((x) => [x, x])), conc: Float64Array.from(c.flatMap((x) => [x, x])),
    blqMask: new Uint8Array(2 * t.length), lloq: 0, animalId: null,
  });
  const excess = (dt: number, c1: number, c2: number): number =>
    dt * (c1 + c2) / 2 - dt * (c1 - c2) / Math.log(c1 / c2);

  it('a two-point exponential tail: % = (linear − log) / AUClast, exact', () => {
    const c1 = 10; const c2 = 10 * Math.exp(-2); // last mean < 20 % Cmax; no t½ → concentration-only rule
    const res = sparseAuc(mk([1, 3], [c1, c2]));
    expect(res.warnings.some((w) => w.code === 'SPARSE_TERMINAL_OVEREST')).toBe(true);
    const want = 100 * excess(2, c1, c2) / res.auc;
    expect(Math.abs(res.terminalExcessPct! - want)).toBeLessThan(1e-12);
    const msg = res.warnings.find((w) => w.code === 'SPARSE_TERMINAL_OVEREST')!.message;
    expect(msg).toContain(`estimated excess ≈ ${want.toFixed(1)} % of AUClast`);
    expect(msg).toContain('the AUClast CI does not include this excess');
  });

  it('summed over every descending interval after Cmax; a rising interval adds nothing', () => {
    const t = [0.5, 1, 2, 4, 8]; const c = [6, 10, 8, 3, 0.5];
    const res = sparseAuc(mk(t, c));
    const want = 100 * (excess(1, 10, 8) + excess(2, 8, 3) + excess(4, 3, 0.5)) / res.auc;
    expect(Math.abs(res.terminalExcessPct! - want)).toBeLessThan(1e-12);
  });

  it('an embedded zero mean adds nothing but does not end the sum (AUClast goes on through it)', () => {
    // 2 h: two of the pair BLQ → the mean is zeroed; the later descending intervals still count.
    const input: SparseInput = {
      nominalTime: Float64Array.from([1, 1, 2, 2, 4, 4, 8, 8]),
      conc: Float64Array.from([10, 10, 0, 0, 3, 3, 0.5, 0.5]),
      blqMask: Uint8Array.from([0, 0, 1, 1, 0, 0, 0, 0]), lloq: 0.1, animalId: null,
    };
    const res = sparseAuc(input);
    const want = 100 * excess(4, 3, 0.5) / res.auc; // 1→2 and 2→4 touch the zero: nothing
    expect(Math.abs(res.terminalExcessPct! - want)).toBeLessThan(1e-12);
  });

  it('a flat interval (c₁ = c₂) adds nothing (no 0/0); silent note → no excess field', () => {
    const res = sparseAuc(mk([1, 2, 3, 6], [10, 10, 5, 1]));
    expect(Math.abs(res.terminalExcessPct! - 100 * (excess(1, 10, 5) + excess(3, 5, 1)) / res.auc)).toBeLessThan(1e-12);
    const silent = sparseAuc(mk([1, 2, 3], [10, 8, 6])); // last mean ≥ 20 % Cmax
    expect(silent.terminalExcessPct).toBeUndefined();
  });
});

describe('sparseAuc — every other route unchanged (AC-2.1 guard)', () => {
  for (const {dose, input} of mouseInputs()) {
    it(`PO mouse ${dose} mg/kg: extravascular route = no route`, () => {
      expect(sparseAuc(input, {route: ROUTE_PO})).toEqual(sparseAuc(input));
    });
  }
});
