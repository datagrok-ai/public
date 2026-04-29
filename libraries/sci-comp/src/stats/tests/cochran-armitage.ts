/**
 * Cochran-Armitage trend test for incidence data.
 *
 * - `cochranArmitageBasic` — original `statistics_fixed.py` API: returns
 *   `{statistic, p_value}` with `null` on degenerate input. Always uses
 *   default scores `0..k-1`, binomial variance, two-sided alternative.
 * - `cochranArmitage` — modified API (Buonaccorsi 2014, Zhou 2017): scores,
 *   alternative, variance method, optional modified statistic. Raises on
 *   invalid input, returns `Z=0, p=1` on degenerate input.
 */

import {normalCdf, normalSf} from '../distributions';
import {sum, toFloat64} from '../internal/normalize';
import {Alternative, NumericInput, TestResult} from '../types';

/** Original (basic) Cochran-Armitage trend test, two-sided. */
export function cochranArmitageBasic(
  counts: NumericInput,
  totals: NumericInput,
): TestResult {
  const c = toFloat64(counts);
  const t = toFloat64(totals);
  const k = c.length;
  if (k < 2 || sum(t) === 0) return {statistic: null, pValue: null};

  const n = sum(t);
  const pBar = sum(c) / n;
  if (pBar === 0 || pBar === 1) return {statistic: null, pValue: null};

  const d = new Float64Array(k);
  for (let i = 0; i < k; i++) d[i] = i;

  let num = 0; let dT = 0;
  for (let i = 0; i < k; i++) {num += d[i] * c[i]; dT += d[i] * t[i];}
  num -= pBar * dT;

  const dBar = dT / n;
  let Sxx = 0;
  for (let i = 0; i < k; i++) {
    const dev = d[i] - dBar;
    Sxx += t[i] * dev * dev;
  }
  const denomSq = pBar * (1 - pBar) * Sxx;
  if (denomSq <= 0) return {statistic: null, pValue: null};

  const z = num / Math.sqrt(denomSq);
  const p = 2 * (1 - normalCdf(Math.abs(z)));
  return {statistic: z, pValue: p};
}

// ── Modified API ─────────────────────────────────────────────────

export type Variance = 'binomial' | 'hypergeometric';

export interface CASettings {
  scores?: NumericInput;
  alternative?: Alternative;
  variance?: Variance;
  modified?: boolean;
}

export interface CAResult {
  zStatistic: number;
  chi2Statistic: number;
  pValue: number;
  alternative: Alternative;
  varianceMethod: Variance;
  scores: number[];
  pBar: number;
  nGroups: number;
  zModified?: number;
  pValueModified?: number;
}

/**
 * Modified Cochran-Armitage trend test.
 *
 * Throws `Error` on invalid inputs (count > total, negative values, length
 * mismatch, k < 2, all totals zero, unknown alternative or variance).
 * Returns `{z=0, p=1, …}` for degenerate cases (p̄ = 0, p̄ = 1, Sxx ≤ 0).
 */
export function cochranArmitage(
  counts: NumericInput,
  totals: NumericInput,
  settings: CASettings = {},
): CAResult {
  const c = toFloat64(counts);
  const t = toFloat64(totals);
  if (c.length !== t.length)
    throw new Error(`cochranArmitage: counts/totals length mismatch (${c.length} vs ${t.length})`);
  const k = c.length;
  if (k < 2)
    throw new Error(`cochranArmitage: need at least 2 groups, got ${k}.`);
  for (let i = 0; i < k; i++) {
    if (c[i] < 0 || t[i] < 0)
      throw new Error('cochranArmitage: counts and totals must be non-negative.');
    if (c[i] > t[i])
      throw new Error('cochranArmitage: each count must be ≤ corresponding total.');
  }
  const n = sum(t);
  if (n === 0)
    throw new Error('cochranArmitage: total sample size is 0.');

  const alternative: Alternative = settings.alternative ?? 'two-sided';
  const variance: Variance = settings.variance ?? 'binomial';
  if (alternative !== 'two-sided' && alternative !== 'increasing' && alternative !== 'decreasing')
    throw new Error(`cochranArmitage: invalid alternative '${alternative}'.`);
  if (variance !== 'binomial' && variance !== 'hypergeometric')
    throw new Error(`cochranArmitage: invalid variance '${variance}'.`);

  const d = new Float64Array(k);
  if (settings.scores !== undefined) {
    const s = toFloat64(settings.scores);
    if (s.length !== k)
      throw new Error(`cochranArmitage: scores length (${s.length}) ≠ number of groups (${k}).`);
    d.set(s);
  } else
    for (let i = 0; i < k; i++) d[i] = i;


  const pBar = sum(c) / n;

  if (pBar === 0 || pBar === 1)
    return degenerate(d, pBar, k, alternative, variance);

  let dotDC = 0; let dotDT = 0;
  for (let i = 0; i < k; i++) {dotDC += d[i] * c[i]; dotDT += d[i] * t[i];}
  const num = dotDC - pBar * dotDT;
  const dBar = dotDT / n;

  let Sxx = 0;
  for (let i = 0; i < k; i++) {
    const dev = d[i] - dBar;
    Sxx += t[i] * dev * dev;
  }
  if (Sxx <= 0)
    return degenerate(d, pBar, k, alternative, variance);

  const denomSqBase = pBar * (1 - pBar) * Sxx;
  const denomSq = variance === 'binomial' ? denomSqBase : denomSqBase * n / (n - 1);

  const z = num / Math.sqrt(denomSq);
  const chi2 = z * z;
  const p = pFromZ(z, alternative);

  const result: CAResult = {
    zStatistic: z,
    chi2Statistic: chi2,
    pValue: p,
    alternative,
    varianceMethod: variance,
    scores: Array.from(d),
    pBar,
    nGroups: k,
  };

  if (settings.modified) {
    let sigma2m = 0;
    for (let i = 0; i < k; i++) {
      const pHat = t[i] > 0 ? c[i] / t[i] : 0;
      const dev = d[i] - dBar;
      sigma2m += t[i] * dev * dev * pHat * (1 - pHat);
    }
    if (sigma2m <= 0) {
      result.zModified = 0;
      result.pValueModified = 1;
    } else {
      const zMod = num / Math.sqrt(sigma2m);
      result.zModified = zMod;
      result.pValueModified = pFromZ(zMod, alternative);
    }
  }

  return result;
}

function degenerate(
  d: Float64Array, pBar: number, k: number,
  alt: Alternative, variance: Variance,
): CAResult {
  return {
    zStatistic: 0,
    chi2Statistic: 0,
    pValue: 1,
    alternative: alt,
    varianceMethod: variance,
    scores: Array.from(d),
    pBar,
    nGroups: k,
  };
}

function pFromZ(z: number, alt: Alternative): number {
  if (alt === 'two-sided') return 2 * normalSf(Math.abs(z));
  if (alt === 'increasing') return normalSf(z);
  return normalCdf(z);
}

// ── Threshold (Williams-style sequential) test ───────────────────

export interface ThresholdSettings {
  alpha?: number;
  adjustAlpha?: boolean;
}

export interface ThresholdStep {
  test: string;
  controlCount: number;
  controlTotal: number;
  controlPct: number;
  treatedCount: number;
  treatedTotal: number;
  treatedPct: number;
  z: number;
  p: number;
  alphaAdj: number;
  significant: boolean;
  effectGroup?: number | null;
  noelGroups?: number[];
}

/**
 * Williams-type sequential threshold test for proportions (Young 1985).
 *
 * Walks from the lowest dose, pooling each non-significant group into a
 * cumulative "control". Stops at the first significant comparison, marking
 * that group as the Effect Level (EL) and the prior groups as NOELs.
 */
export function thresholdTest(
  counts: NumericInput,
  totals: NumericInput,
  settings: ThresholdSettings = {},
): ThresholdStep[] {
  const c = Array.from(toFloat64(counts));
  const t = Array.from(toFloat64(totals));
  const k = c.length;
  if (k < 2)
    throw new Error('thresholdTest: need at least 2 groups (control + 1 dose).');

  const alpha = settings.alpha ?? 0.05;
  const adjust = settings.adjustAlpha ?? true;
  const nComparisons = k - 1;
  const alphaAdj = adjust ? 1.0 - Math.pow(1.0 - alpha, 1.0 / nComparisons) : alpha;

  const results: ThresholdStep[] = [];
  let poolCount = c[0];
  let poolTotal = t[0];

  for (let i = 1; i < k; i++) {
    const r = cochranArmitage([poolCount, c[i]], [poolTotal, t[i]], {
      scores: [0, 1],
      alternative: 'increasing',
    });
    const sig = r.pValue <= alphaAdj;
    const entry: ThresholdStep = {
      test: `groups [${range(i).join(', ')}] vs group ${i}`,
      controlCount: poolCount,
      controlTotal: poolTotal,
      controlPct: poolTotal > 0 ? poolCount / poolTotal * 100 : 0,
      treatedCount: c[i],
      treatedTotal: t[i],
      treatedPct: t[i] > 0 ? c[i] / t[i] * 100 : 0,
      z: r.zStatistic,
      p: r.pValue,
      alphaAdj,
      significant: sig,
    };
    if (sig) {
      entry.effectGroup = i;
      entry.noelGroups = range(i);
      results.push(entry);
      return results;
    }
    poolCount += c[i];
    poolTotal += t[i];
    results.push(entry);
  }

  if (results.length > 0) {
    results[results.length - 1].noelGroups = range(k);
    results[results.length - 1].effectGroup = null;
  }
  return results;
}

function range(n: number): number[] {
  const out: number[] = [];
  for (let i = 0; i < n; i++) out.push(i);
  return out;
}
