/**
 * Hedges' g — bias-corrected Cohen's d for small samples.
 */

import {mean, stripNaN, toFloat64, variance} from '../internal/normalize';
import {NumericInput} from '../types';

/**
 * Hedges' g effect size with the small-sample J-correction.
 *
 * `g = d · J`, where `d` is Cohen's d with weighted pooled SD and
 * `J ≈ 1 − 3/(4·df − 1)` (Hedges 1981; df = n1 + n2 − 2).
 *
 * Returns `null` when either group has fewer than 2 finite observations or
 * when the pooled SD is zero.
 */
export function hedgesG(g1: NumericInput, g2: NumericInput): number | null {
  const a = stripNaN(toFloat64(g1));
  const b = stripNaN(toFloat64(g2));
  if (a.length < 2 || b.length < 2) return null;

  const n1 = a.length;
  const n2 = b.length;
  const df = n1 + n2 - 2;
  const pooledStd = Math.sqrt(
    ((n1 - 1) * variance(a) + (n2 - 1) * variance(b)) / df,
  );
  if (pooledStd === 0) return null;

  const d = (mean(a) - mean(b)) / pooledStd;
  const j = 1 - 3 / (4 * df - 1);
  return d * j;
}
