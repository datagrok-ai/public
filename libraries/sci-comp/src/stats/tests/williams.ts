/**
 * Williams' step-down test for dose-response with PAVA isotonic regression.
 */

import {studentTCdf, studentTInv} from '../distributions';
import {sum, toFloat64} from '../internal/normalize';
import {NumericInput} from '../types';
import {lookup1971, lookup1972} from './williams-tables';

export type WilliamsDirection = 'increase' | 'decrease' | 'auto';

export interface WilliamsSettings {
  direction?: WilliamsDirection;
  alpha?: number;
}

export interface WilliamsStepResult {
  doseLabel: string;
  doseIndex: number;
  constrainedMean: number;
  controlMean: number;
  testStatistic: number;
  criticalValue: number;
  pValue: number;
  significant: boolean;
  alpha: number;
}

export interface WilliamsResult {
  direction: 'increase' | 'decrease';
  pooledVariance: number;
  pooledDf: number;
  constrainedMeans: number[];
  stepDownResults: WilliamsStepResult[];
  minimumEffectiveDose: string | null;
  minimumEffectiveIndex: number | null;
  allGroupsTested: boolean;
}

// ── PAVA — Pool-Adjacent-Violators ───────────────────────────────

/**
 * Isotonic regression under a non-decreasing constraint via the stack-based
 * pool-adjacent-violators algorithm.
 */
export function pavaIncreasing(
  values: NumericInput, weights: NumericInput,
): Float64Array {
  const v = toFloat64(values);
  const w = toFloat64(weights);
  const n = v.length;
  if (n <= 1) return v.slice();

  const blockValues: number[] = [];
  const blockWeights: number[] = [];
  const blockSizes: number[] = [];

  for (let i = 0; i < n; i++) {
    blockValues.push(v[i]);
    blockWeights.push(w[i]);
    blockSizes.push(1);
    while (blockValues.length >= 2 && blockValues[blockValues.length - 2] > blockValues[blockValues.length - 1]) {
      const w1 = blockWeights[blockWeights.length - 2];
      const w2 = blockWeights[blockWeights.length - 1];
      const v1 = blockValues[blockValues.length - 2];
      const v2 = blockValues[blockValues.length - 1];
      const pooled = (w1 * v1 + w2 * v2) / (w1 + w2);
      const s = blockSizes[blockSizes.length - 2] + blockSizes[blockSizes.length - 1];
      blockValues.pop(); blockWeights.pop(); blockSizes.pop();
      blockValues[blockValues.length - 1] = pooled;
      blockWeights[blockWeights.length - 1] = w1 + w2;
      blockSizes[blockSizes.length - 1] = s;
    }
  }

  const out = new Float64Array(n);
  let idx = 0;
  for (let b = 0; b < blockValues.length; b++)
    for (let k = 0; k < blockSizes[b]; k++) out[idx++] = blockValues[b];

  return out;
}

/** Isotonic regression under a non-increasing constraint. */
export function pavaDecreasing(
  values: NumericInput, weights: NumericInput,
): Float64Array {
  const v = toFloat64(values);
  const negated = new Float64Array(v.length);
  for (let i = 0; i < v.length; i++) negated[i] = -v[i];
  const result = pavaIncreasing(negated, weights);
  for (let i = 0; i < result.length; i++) result[i] = -result[i];
  return result;
}

// ── Critical-value selection ─────────────────────────────────────

const ALPHA_1971 = new Set([0.05, 0.01]);
const DOSE_LEVELS_1972 = new Set([2, 3, 4, 5, 6, 8, 10]);

function getCriticalValue(
  doseIndex: number, df: number, ns: Float64Array, alpha: number,
): number {
  const i = doseIndex;
  if (i === 1) return studentTInv(1 - alpha, df);

  const w = ns[0] / ns[i];

  if (w <= 1.0 && ALPHA_1971.has(alpha) && i <= 10)
    return lookup1971(i, df, alpha);

  const supported1972 = new Set([0.050, 0.025, 0.010, 0.005]);
  if (supported1972.has(alpha) && i <= 10) {
    if (DOSE_LEVELS_1972.has(i))
      return lookup1972(i, df, alpha, w);
    if (i === 7) {
      const cv6 = lookup1972(6, df, alpha, w);
      const cv8 = lookup1972(8, df, alpha, w);
      return Math.round(((cv6 + cv8) / 2) * 1000) / 1000;
    }
    if (i === 9) {
      const cv8 = lookup1972(8, df, alpha, w);
      const cv10 = lookup1972(10, df, alpha, w);
      return Math.round(((cv8 + cv10) / 2) * 1000) / 1000;
    }
    return lookup1972(10, df, alpha, w);
  }

  return studentTInv(1 - alpha, df);
}

// ── Williams' test ───────────────────────────────────────────────

/**
 * Williams' step-down test for dose-response trend.
 *
 * Inputs: per-group means, standard deviations, sample sizes, and dose
 * labels. Index 0 is the control. Direction `'auto'` infers from highest
 * dose vs control.
 */
export function williamsTest(
  means: NumericInput,
  sds: NumericInput,
  ns: NumericInput,
  doseLabels: readonly string[],
  settings: WilliamsSettings = {},
): WilliamsResult {
  const meanArr = toFloat64(means);
  const sdArr = toFloat64(sds);
  const nArr = toFloat64(ns);
  let direction = settings.direction ?? 'auto';
  const alpha = settings.alpha ?? 0.05;

  const k = meanArr.length - 1;
  const N = sum(nArr);
  const dfPooled = N - k - 1;

  if (k < 1 || dfPooled < 1) {
    return {
      direction: direction === 'auto' ? 'increase' : direction,
      pooledVariance: 0,
      pooledDf: Math.max(dfPooled, 0),
      constrainedMeans: Array.from(meanArr),
      stepDownResults: [],
      minimumEffectiveDose: null,
      minimumEffectiveIndex: null,
      allGroupsTested: false,
    };
  }

  let ssWithin = 0;
  for (let i = 0; i < nArr.length; i++)
    ssWithin += (nArr[i] - 1) * sdArr[i] * sdArr[i];
  const s2Pooled = ssWithin / dfPooled;
  const sPooled = Math.sqrt(s2Pooled);

  if (direction === 'auto')
    direction = meanArr[meanArr.length - 1] > meanArr[0] ? 'increase' : 'decrease';

  const doseMeans = meanArr.slice(1);
  const doseNs = nArr.slice(1);
  const doseConstrained = direction === 'increase' ?
    pavaIncreasing(doseMeans, doseNs) :
    pavaDecreasing(doseMeans, doseNs);

  const constrained = new Float64Array(meanArr.length);
  constrained[0] = meanArr[0];
  for (let i = 0; i < doseConstrained.length; i++) constrained[i + 1] = doseConstrained[i];

  const results: WilliamsStepResult[] = [];
  for (let i = k; i >= 1; i--) {
    const se = sPooled * Math.sqrt(1 / nArr[0] + 1 / nArr[i]);
    if (se <= 0) break;
    const t = direction === 'increase' ?
      (constrained[i] - meanArr[0]) / se :
      (meanArr[0] - constrained[i]) / se;
    const cv = getCriticalValue(i, dfPooled, nArr, alpha);
    const sig = t > cv;
    const pApprox = Math.max(1 - studentTCdf(t, dfPooled), 0);

    results.push({
      doseLabel: doseLabels[i],
      doseIndex: i,
      constrainedMean: constrained[i],
      controlMean: meanArr[0],
      testStatistic: round(t, 4),
      criticalValue: round(cv, 4),
      pValue: round(pApprox, 6),
      significant: sig,
      alpha,
    });

    if (!sig) break;
  }

  const sigResults = results.filter((r) => r.significant);
  let medLabel: string | null = null;
  let medIndex: number | null = null;
  if (sigResults.length > 0) {
    const med = sigResults[sigResults.length - 1];
    medLabel = med.doseLabel;
    medIndex = med.doseIndex;
  }

  return {
    direction,
    pooledVariance: round(s2Pooled, 6),
    pooledDf: dfPooled,
    constrainedMeans: Array.from(constrained, (v) => round(v, 6)),
    stepDownResults: results,
    minimumEffectiveDose: medLabel,
    minimumEffectiveIndex: medIndex,
    allGroupsTested: results.length === k,
  };
}

function round(x: number, digits: number): number {
  const f = Math.pow(10, digits);
  return Math.round(x * f) / f;
}
