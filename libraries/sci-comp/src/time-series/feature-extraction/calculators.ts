/**
 * Multi-pass feature calculators.
 *
 * Computations are grouped into three passes plus median and linear trend
 * to minimize iterations over each sample slice.
 */

import {NumericArray} from './types';

/* ================================================================== */
/*  Pass 1: basic statistics                                           */
/* ================================================================== */

/** Output of {@link pass1}: basic aggregates collected in a single loop. */
export interface Pass1Result {
  /** Number of data points. */
  n: number;
  /** Sum of values. */
  sum: number;
  /** Sum of squared values (abs_energy). */
  sumSq: number;
  /** Minimum value. */
  min: number;
  /** Maximum value. */
  max: number;
  /** Index of the first occurrence of the minimum. */
  firstMinIdx: number;
  /** Index of the first occurrence of the maximum. */
  firstMaxIdx: number;
  /** Index of the last occurrence of the minimum. */
  lastMinIdx: number;
  /** Index of the last occurrence of the maximum. */
  lastMaxIdx: number;
  /** Number of elements equal to the minimum. */
  minCount: number;
  /** Number of elements equal to the maximum. */
  maxCount: number;
  /** Arithmetic mean (sum / n). */
  mean: number;
  /** Population variance (ddof=0). */
  variance: number;
  /** Population standard deviation. */
  std: number;
  /** Range (max − min). */
  range: number;
}

/**
 * Pass 1: single loop collecting sum, sumSq, min/max locations and counts.
 *
 * Derived fields (mean, variance, std, range) are computed at the end.
 */
export function pass1(data: NumericArray, n: number): Pass1Result {
  let sum = 0;
  let sumSq = 0;
  let min = Infinity;
  let max = -Infinity;
  let firstMinIdx = 0;
  let firstMaxIdx = 0;
  let lastMinIdx = 0;
  let lastMaxIdx = 0;
  let minCount = 0;
  let maxCount = 0;

  for (let i = 0; i < n; i++) {
    const v = data[i];
    sum += v;
    sumSq += v * v;

    if (v < min) {
      min = v; firstMinIdx = i; lastMinIdx = i; minCount = 1;
    } else if (v === min) {
      lastMinIdx = i; minCount++;
    }

    if (v > max) {
      max = v; firstMaxIdx = i; lastMaxIdx = i; maxCount = 1;
    } else if (v === max) {
      lastMaxIdx = i; maxCount++;
    }
  }

  const mean = sum / n;
  let variance = sumSq / n - mean * mean;
  if (variance < 0) variance = 0;
  const std = Math.sqrt(variance);
  const range = max - min;

  return {
    n, sum, sumSq, min, max,
    firstMinIdx, firstMaxIdx, lastMinIdx, lastMaxIdx,
    minCount, maxCount,
    mean, variance, std, range,
  };
}

/* ================================================================== */
/*  Pass 2: differences, uniqueness, crossings, strikes                */
/* ================================================================== */

/** Output of {@link pass2}: consecutive-difference and uniqueness metrics. */
export interface Pass2Result {
  /** Sum of absolute consecutive differences. */
  sumAbsDiff: number;
  /** Sum of squared consecutive differences. */
  sumSqDiff: number;
  /** Number of unique values. */
  uniqueCount: number;
  /** Whether any value appears more than once. */
  hasDuplicate: boolean;
  /** Total count of values that appear more than once. */
  reoccurringDatapointCount: number;
  /** Number of times the series crosses `mValue`. */
  crossingCount: number;
  /** Longest consecutive run of values strictly above the mean. */
  longestStrikeAbove: number;
  /** Longest consecutive run of values strictly below the mean. */
  longestStrikeBelow: number;
}

/**
 * Pass 2: single loop collecting differences, uniqueness, crossings, and strikes.
 *
 * @param data   Sample slice.
 * @param n      Number of points.
 * @param mean   Mean from pass 1 (used for strike computation).
 * @param mValue Crossing threshold (typically 0).
 */
export function pass2(
  data: NumericArray, n: number, mean: number, mValue: number,
): Pass2Result {
  let sumAbsDiff = 0;
  let sumSqDiff = 0;
  let crossingCount = 0;

  const valueCounts = new Map<number, number>();

  let curAbove = 0;
  let maxAbove = 0;
  let curBelow = 0;
  let maxBelow = 0;

  const first = data[0];
  valueCounts.set(first, 1);
  if (first > mean) curAbove = 1;
  if (first < mean) curBelow = 1;

  for (let i = 1; i < n; i++) {
    const v = data[i];
    const prev = data[i - 1];
    const d = v - prev;

    sumAbsDiff += Math.abs(d);
    sumSqDiff += d * d;

    if ((prev - mValue) * (v - mValue) < 0) crossingCount++;

    valueCounts.set(v, (valueCounts.get(v) ?? 0) + 1);

    if (v > mean) {
      curAbove++;
      if (curAbove > maxAbove) maxAbove = curAbove;
      curBelow = 0;
    } else if (v < mean) {
      curBelow++;
      if (curBelow > maxBelow) maxBelow = curBelow;
      curAbove = 0;
    } else {
      curAbove = 0;
      curBelow = 0;
    }
  }

  if (curAbove > maxAbove) maxAbove = curAbove;
  if (curBelow > maxBelow) maxBelow = curBelow;

  const uniqueCount = valueCounts.size;
  const hasDuplicate = uniqueCount !== n;
  let reoccurringDatapointCount = 0;
  for (const count of valueCounts.values())
    if (count > 1) reoccurringDatapointCount += count;

  return {
    sumAbsDiff, sumSqDiff,
    uniqueCount, hasDuplicate, reoccurringDatapointCount,
    crossingCount,
    longestStrikeAbove: maxAbove,
    longestStrikeBelow: maxBelow,
  };
}

/* ================================================================== */
/*  Pass 3: threshold-based features                                   */
/* ================================================================== */

/** Output of {@link pass3}: count above/below mean and ratio-beyond-r-sigma counts. */
export interface Pass3Result {
  /** Number of values strictly above the mean. */
  countAboveMean: number;
  /** Number of values strictly below the mean. */
  countBelowMean: number;
  /** Per-threshold count of values beyond r × σ from the mean. */
  beyondRSigmaCounts: number[];
}

/**
 * Pass 3: single loop collecting count above/below mean and ratio_beyond_r_sigma.
 *
 * @param data         Sample slice.
 * @param n            Number of points.
 * @param mean         Mean from pass 1.
 * @param std          Standard deviation from pass 1.
 * @param rSigmaValues Array of r thresholds (e.g. [1, 1.5, 2]).
 */
export function pass3(
  data: NumericArray, n: number, mean: number, std: number, rSigmaValues: number[],
): Pass3Result {
  let countAbove = 0;
  let countBelow = 0;
  const beyondCounts = new Array<number>(rSigmaValues.length).fill(0);

  for (let i = 0; i < n; i++) {
    const v = data[i];
    if (v > mean) countAbove++;
    if (v < mean) countBelow++;

    if (std > 0) {
      const absDeviation = Math.abs(v - mean);
      for (let j = 0; j < rSigmaValues.length; j++) {
        if (absDeviation > rSigmaValues[j] * std)
          beyondCounts[j]++;
      }
    }
  }

  return {
    countAboveMean: countAbove,
    countBelowMean: countBelow,
    beyondRSigmaCounts: beyondCounts,
  };
}

/* ================================================================== */
/*  Median                                                             */
/* ================================================================== */

/**
 * Compute the median of a sample slice.
 *
 * Copies data into `buf`, sorts the copy, and returns the middle element.
 * The original `data` array is never mutated.
 *
 * @param data Sample slice (not modified).
 * @param n    Number of points.
 * @param buf  Pre-allocated scratch buffer (length ≥ n).
 */
export function computeMedian(data: NumericArray, n: number, buf: Float64Array): number {
  for (let i = 0; i < n; i++)
    buf[i] = data[i];

  const slice = buf.subarray(0, n);
  slice.sort();

  if (n % 2 === 0)
    return (slice[n / 2 - 1] + slice[n / 2]) / 2;
  else
    return slice[(n - 1) / 2];
}
