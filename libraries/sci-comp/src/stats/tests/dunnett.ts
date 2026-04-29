/**
 * Dunnett's many-to-one comparison test (treated groups vs control).
 *
 * Adjusted p-values come from the joint distribution of the contrasts under
 * H₀, computed via 2D numerical integration:
 *
 *   t_i = (α_i Z_i − β_i Z_0) / s,
 *
 * where (Z_i, Z_0) are i.i.d. standard normals, s = S/σ ~ chi(df)/√df,
 * and α_i = √(n_0/(n_0+n_i)), β_i = √(n_i/(n_0+n_i)) (Dunnett 1955).
 */

import {gammaln, normalCdf} from '../distributions';
import {mean, stripNaN, toFloat64} from '../internal/normalize';
import {NumericInput} from '../types';
import {hedgesG} from './hedges-g';
import {welchTTest} from './welch-t';

export interface TreatedGroup {
  doseLevel: number;
  values: NumericInput;
}

export interface DunnettPairwise {
  doseLevel: number;
  pValue: number | null;
  pValueAdj: number | null;
  statistic: number | null;
  effectSize: number | null;
}

/**
 * Run Dunnett's test for each treated group vs the control.
 *
 * Returns an empty list when the control has fewer than 2 finite values or
 * no treated groups are provided. Treated groups with fewer than 2 finite
 * values are returned with `null` statistic/p_value but a computed
 * `effectSize` (Hedges' g).
 *
 * Falls back to Welch + Bonferroni if Dunnett's pooled-variance computation
 * is not feasible (e.g. residual df < 1) — matching the Python fallback.
 */
export function dunnettPairwise(
  control: NumericInput,
  treated: readonly TreatedGroup[],
): DunnettPairwise[] {
  const ctrl = stripNaN(toFloat64(control));
  if (ctrl.length < 2 || treated.length === 0) return [];

  const doseLevels: number[] = [];
  const effectSizes: (number | null)[] = [];
  const valid: Float64Array[] = [];
  const validIdx: number[] = [];

  for (let i = 0; i < treated.length; i++) {
    const arr = stripNaN(toFloat64(treated[i].values));
    doseLevels.push(treated[i].doseLevel);
    effectSizes.push(hedgesG(arr, ctrl));
    if (arr.length >= 2) {valid.push(arr); validIdx.push(i);}
  }

  // Pooled within-group variance across control + valid treated groups.
  const k = valid.length;
  let nTotal = ctrl.length;
  let ss = sumSqDev(ctrl);
  for (const v of valid) {nTotal += v.length; ss += sumSqDev(v);}
  const dfPooled = nTotal - (k + 1);

  // Build fallback: Welch + Bonferroni when Dunnett is infeasible.
  if (dfPooled < 1 || k === 0)
    return assembleFallback(control, treated, doseLevels, effectSizes, validIdx);


  const mse = ss / dfPooled;
  const s = Math.sqrt(mse);
  const meanCtrl = mean(ctrl);
  const ts: number[] = [];
  const alphas: number[] = [];
  const betas: number[] = [];
  const n0 = ctrl.length;
  for (const arr of valid) {
    const ni = arr.length;
    const seFactor = Math.sqrt(1 / ni + 1 / n0);
    ts.push((mean(arr) - meanCtrl) / (s * seFactor));
    alphas.push(Math.sqrt(n0 / (n0 + ni)));
    betas.push(Math.sqrt(ni / (n0 + ni)));
  }

  // For each contrast i, p_adj_i = 1 − P(|T_max| < |t_i| | df, α, β).
  const pAdj = ts.map((t) => 1 - dunnettMaxCdf(Math.abs(t), dfPooled, alphas, betas));

  return assembleDunnett(treated, doseLevels, effectSizes, validIdx, ts, pAdj);
}

// ── Result assembly ──────────────────────────────────────────────

function assembleDunnett(
  treated: readonly TreatedGroup[],
  doseLevels: number[],
  effectSizes: (number | null)[],
  validIdx: number[],
  ts: number[],
  pAdj: number[],
): DunnettPairwise[] {
  const out: DunnettPairwise[] = [];
  for (let i = 0; i < treated.length; i++) {
    const j = validIdx.indexOf(i);
    if (j < 0) {
      out.push({
        doseLevel: doseLevels[i],
        pValue: null, pValueAdj: null, statistic: null,
        effectSize: round(effectSizes[i], 4),
      });
    } else {
      const p = pAdj[j];
      const pRounded = round(p, 6);
      out.push({
        doseLevel: doseLevels[i],
        pValue: pRounded,
        pValueAdj: pRounded,
        statistic: ts[j],
        effectSize: round(effectSizes[i], 4),
      });
    }
  }
  return out;
}

function assembleFallback(
  control: NumericInput,
  treated: readonly TreatedGroup[],
  doseLevels: number[],
  effectSizes: (number | null)[],
  validIdx: number[],
): DunnettPairwise[] {
  const nValid = validIdx.length;
  const out: DunnettPairwise[] = [];
  for (let i = 0; i < treated.length; i++) {
    if (validIdx.indexOf(i) < 0) {
      out.push({
        doseLevel: doseLevels[i],
        pValue: null, pValueAdj: null, statistic: null,
        effectSize: round(effectSizes[i], 4),
      });
      continue;
    }
    const r = welchTTest(treated[i].values, control);
    const rawP = r.pValue;
    const adjP = rawP === null || Number.isNaN(rawP) ?
      null :
      Math.min(rawP * nValid, 1.0);
    out.push({
      doseLevel: doseLevels[i],
      pValue: round(adjP, 6),
      pValueAdj: round(adjP, 6),
      statistic: r.statistic === null || Number.isNaN(r.statistic) ? null : r.statistic,
      effectSize: round(effectSizes[i], 4),
    });
  }
  return out;
}

// ── 2D numerical integration ──────────────────────────────────────

/**
 * P(|T_1| < c, ..., |T_k| < c) under H₀ with Dunnett structure.
 *
 * Outer Simpson over s ∈ (0, sMax], inner Simpson over Z_0 ∈ [-zMax, zMax].
 * The integrand is the product of the conditional rectangle probabilities.
 */
function dunnettMaxCdf(
  c: number, df: number, alphas: number[], betas: number[],
): number {
  const sMax = 6.0;
  const sMin = 1e-3;
  const zMax = 8.0;
  const nOuter = 80; // Simpson panels for s
  const nInner = 80; // Simpson panels for z

  const outerH = (sMax - sMin) / nOuter;
  const innerH = (2 * zMax) / nInner;

  let total = 0;
  for (let i = 0; i <= nOuter; i++) {
    const s = sMin + i * outerH;
    const fs = chiDensity(s, df);
    if (fs === 0) continue;
    const wOuter = (i === 0 || i === nOuter) ? 1 : (i % 2 === 1 ? 4 : 2);

    let inner = 0;
    for (let j = 0; j <= nInner; j++) {
      const z = -zMax + j * innerH;
      const phi = Math.exp(-0.5 * z * z) * 0.39894228040143268; // 1/√(2π)
      let prod = 1;
      for (let m = 0; m < alphas.length; m++) {
        const upper = (c * s + betas[m] * z) / alphas[m];
        const lower = (-c * s + betas[m] * z) / alphas[m];
        prod *= normalCdf(upper) - normalCdf(lower);
        if (prod === 0) break;
      }
      const wInner = (j === 0 || j === nInner) ? 1 : (j % 2 === 1 ? 4 : 2);
      inner += wInner * phi * prod;
    }
    inner *= innerH / 3;

    total += wOuter * fs * inner;
  }
  return total * outerH / 3;
}

/** Density of s = S/σ where S² · df / σ² ~ χ²_df. */
function chiDensity(s: number, df: number): number {
  if (s <= 0) return 0;
  const logF =
    (df / 2) * Math.log(df) -
    (df / 2 - 1) * Math.LN2 -
    gammaln(df / 2) +
    (df - 1) * Math.log(s) -
    s * s * df / 2;
  return Math.exp(logF);
}

// ── Helpers ──────────────────────────────────────────────────────

function sumSqDev(arr: Float64Array): number {
  const m = mean(arr);
  let s = 0;
  for (let i = 0; i < arr.length; i++) {
    const d = arr[i] - m;
    s += d * d;
  }
  return s;
}

function round(x: number | null, digits: number): number | null {
  if (x === null) return null;
  const f = Math.pow(10, digits);
  return Math.round(x * f) / f;
}
