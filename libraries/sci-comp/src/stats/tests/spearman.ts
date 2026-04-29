/**
 * Spearman rank correlation and severity-trend (a thin wrapper).
 */

import {studentTTwoTail} from '../distributions';
import {mean, pairwiseNonNaN, toFloat64} from '../internal/normalize';
import {rankdata} from '../internal/rank';
import {NumericInput, SpearmanResult} from '../types';

/**
 * Spearman's rank correlation. Strips NaN pairwise, requires n ≥ 3.
 *
 * The p-value uses the t-distribution approximation:
 *   t = ρ · √((n−2) / (1 − ρ²)),  df = n − 2.
 */
export function spearman(x: NumericInput, y: NumericInput): SpearmanResult {
  const ax = toFloat64(x);
  const ay = toFloat64(y);
  if (ax.length !== ay.length)
    throw new Error(`spearman: length mismatch (${ax.length} vs ${ay.length})`);
  const {x: a, y: b} = pairwiseNonNaN(ax, ay);
  const n = a.length;
  if (n < 3) return {rho: null, pValue: null};

  const rho = pearson(rankdata(a), rankdata(b));
  if (Number.isNaN(rho)) return {rho: null, pValue: null};

  if (rho === 1 || rho === -1) return {rho, pValue: 0};
  const t = rho * Math.sqrt((n - 2) / (1 - rho * rho));
  return {rho, pValue: studentTTwoTail(t, n - 2)};
}

/**
 * Spearman correlation between dose levels and average severities.
 *
 * Guards against constant severity (returns null) and otherwise delegates to
 * `spearman`. Matches the Python `severity_trend` semantics exactly.
 */
export function severityTrend(
  doseLevels: NumericInput,
  avgSeverities: NumericInput,
): SpearmanResult {
  const dl = toFloat64(doseLevels);
  const sev = toFloat64(avgSeverities);
  if (dl.length !== sev.length)
    throw new Error(`severityTrend: length mismatch`);
  const {x: a, y: b} = pairwiseNonNaN(dl, sev);
  if (a.length < 3) return {rho: null, pValue: null};
  // Constant severity → correlation undefined
  let constant = true;
  for (let i = 1; i < b.length; i++)
    if (b[i] !== b[0]) {constant = false; break;}

  if (constant) return {rho: null, pValue: null};
  return spearman(a, b);
}

function pearson(x: Float64Array, y: Float64Array): number {
  const n = x.length;
  const mx = mean(x);
  const my = mean(y);
  let sxy = 0; let sxx = 0; let syy = 0;
  for (let i = 0; i < n; i++) {
    const dx = x[i] - mx;
    const dy = y[i] - my;
    sxy += dx * dy;
    sxx += dx * dx;
    syy += dy * dy;
  }
  if (sxx === 0 || syy === 0) return NaN;
  return sxy / Math.sqrt(sxx * syy);
}
