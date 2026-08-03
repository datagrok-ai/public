/**
 * Welch's t-test (unequal variances).
 */

import {studentTTwoTail} from '../distributions';
import {mean, stripNaN, toFloat64, variance} from '../internal/normalize';
import {NumericInput, TestResult} from '../types';

/**
 * Welch's two-sample t-test. Strips NaN before computation.
 *
 * Returns `{statistic: null, pValue: null}` when either group has fewer than
 * 2 finite observations. Returns `{statistic: NaN, pValue: NaN}` (matching
 * scipy) when both groups have zero variance.
 */
export function welchTTest(x: NumericInput, y: NumericInput): TestResult {
  const a = stripNaN(toFloat64(x));
  const b = stripNaN(toFloat64(y));
  if (a.length < 2 || b.length < 2)
    return {statistic: null, pValue: null};

  const m1 = mean(a);
  const m2 = mean(b);
  const v1 = variance(a);
  const v2 = variance(b);
  const n1 = a.length;
  const n2 = b.length;

  const seSq = v1 / n1 + v2 / n2;
  if (seSq === 0)
    return {statistic: NaN, pValue: NaN};

  const t = (m1 - m2) / Math.sqrt(seSq);
  // Welch–Satterthwaite degrees of freedom
  const df = (seSq * seSq) / (
    (v1 * v1) / (n1 * n1 * (n1 - 1)) +
    (v2 * v2) / (n2 * n2 * (n2 - 1))
  );
  const p = studentTTwoTail(t, df);

  return {statistic: t, pValue: p};
}
