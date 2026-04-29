/**
 * ANCOVA — analysis of covariance for organ-weight normalisation.
 *
 * Fits `y ~ C(group) + covariate` via OLS, reports adjusted (LS) means,
 * pairwise contrasts vs control, slope homogeneity F-test, and an effect
 * decomposition (total / direct / indirect with Hedges' g).
 */

import {fCdf, studentTCdf} from '../distributions';
import {mean, toFloat64} from '../internal/normalize';
import {inverse, matmul, Matrix, matvec, quadForm, transpose} from '../internal/matrix';
import {NumericInput} from '../types';

export interface AncovaSettings {
  controlGroup?: number;
  useOrganFreeBw?: boolean;
  alpha?: number;
}

export interface AdjustedMean {
  group: number;
  rawMean: number;
  adjustedMean: number;
  n: number;
  se: number;
}

export interface AncovaPairwise {
  group: number;
  difference: number;
  se: number;
  tStatistic: number;
  pValue: number;
  significant: boolean;
}

export interface AncovaSlope {
  estimate: number;
  se: number;
  tStatistic: number;
  pValue: number;
}

export interface AncovaSlopeHomogeneity {
  fStatistic: number | null;
  pValue: number | null;
  homogeneous: boolean;
}

export interface AncovaEffectDecomposition {
  group: number;
  totalEffect: number;
  directEffect: number;
  indirectEffect: number;
  proportionDirect: number;
  directG: number;
  directP: number;
}

export interface AncovaResult {
  adjustedMeans: AdjustedMean[];
  pairwise: AncovaPairwise[];
  slope: AncovaSlope;
  slopeHomogeneity: AncovaSlopeHomogeneity;
  effectDecomposition: AncovaEffectDecomposition[];
  modelRSquared: number;
  mse: number;
  useOrganFreeBw: boolean;
  covariateMean: number;
}

/**
 * One-way ANCOVA: `organ_weight ~ C(dose_group) + body_weight`.
 *
 * Returns `null` when there are insufficient observations
 * (`n < k + 2` or `k < 2`).
 */
export function runAncova(
  organValues: NumericInput,
  bodyWeights: NumericInput,
  groups: NumericInput,
  settings: AncovaSettings = {},
): AncovaResult | null {
  const ovRaw = toFloat64(organValues);
  const bwRaw = toFloat64(bodyWeights);
  const gpRaw = toFloat64(groups);
  if (ovRaw.length !== bwRaw.length || ovRaw.length !== gpRaw.length)
    throw new Error('runAncova: input length mismatch.');

  const ov: number[] = [];
  const bw: number[] = [];
  const gp: number[] = [];
  for (let i = 0; i < ovRaw.length; i++) {
    if (Number.isNaN(ovRaw[i]) || Number.isNaN(bwRaw[i])) continue;
    ov.push(ovRaw[i]);
    bw.push(bwRaw[i]);
    gp.push(gpRaw[i] | 0);
  }

  const uniqueGroups = Array.from(new Set(gp)).sort((a, b) => a - b);
  const k = uniqueGroups.length;
  const n = ov.length;
  if (n < k + 2 || k < 2) return null;

  const controlGroup = settings.controlGroup ?? 0;
  const useOrganFreeBw = settings.useOrganFreeBw ?? false;
  const alpha = settings.alpha ?? 0.05;

  const cov: number[] = useOrganFreeBw ?
    bw.map((b, i) => b - ov[i]) :
    bw.slice();
  const covMean = mean(Float64Array.from(cov));

  const treated = uniqueGroups.filter((g) => g !== controlGroup);
  const pA = 1 + treated.length + 1;

  const Xa: Matrix = [];
  for (let i = 0; i < n; i++) {
    const row = new Array<number>(pA).fill(0);
    row[0] = 1;
    for (let j = 0; j < treated.length; j++)
      row[1 + j] = gp[i] === treated[j] ? 1 : 0;
    row[pA - 1] = cov[i];
    Xa.push(row);
  }

  const fit = fitOls(Xa, ov);
  const dfA = fit.df;
  if (dfA <= 0) return null;
  const mse = fit.rss / dfA;

  // Interaction model for slope-homogeneity F-test
  const pInt = pA + treated.length;
  const Xint: Matrix = [];
  for (let i = 0; i < n; i++) {
    const row = new Array<number>(pInt).fill(0);
    for (let j = 0; j < pA; j++) row[j] = Xa[i][j];
    for (let j = 0; j < treated.length; j++)
      row[pA + j] = Xa[i][1 + j] * cov[i];
    Xint.push(row);
  }
  let fHom: number | null = null;
  let pHom: number | null = null;
  try {
    const fitInt = fitOls(Xint, ov);
    const dfDiff = dfA - fitInt.df;
    if (dfDiff > 0 && fitInt.df > 0 && fitInt.rss > 0) {
      fHom = ((fit.rss - fitInt.rss) / dfDiff) / (fitInt.rss / fitInt.df);
      pHom = 1 - fCdf(fHom, dfDiff, fitInt.df);
    }
  } catch {
    // Singular interaction model — leave fHom/pHom as null.
  }
  const homogeneous = pHom === null || pHom >= alpha;

  const slopeEst = fit.beta[pA - 1];
  const slopeSe = Math.sqrt(Math.max(0, fit.vcov[pA - 1][pA - 1]));
  const slopeT = slopeSe > 0 ? slopeEst / slopeSe : 0;
  const slopeP = 2 * (1 - studentTCdf(Math.abs(slopeT), dfA));

  const ovMean = mean(Float64Array.from(ov));
  let tss = 0;
  for (let i = 0; i < n; i++) tss += (ov[i] - ovMean) ** 2;
  const rSq = tss > 0 ? 1 - fit.rss / tss : 0;

  const adjustedMeans: AdjustedMean[] = [];
  const rawMeans = new Map<number, number>();
  for (const g of uniqueGroups) {
    const grpVals: number[] = [];
    for (let i = 0; i < n; i++) if (gp[i] === g) grpVals.push(ov[i]);
    const rm = mean(Float64Array.from(grpVals));
    rawMeans.set(g, rm);
    const xPred = new Array<number>(pA).fill(0);
    xPred[0] = 1;
    if (treated.includes(g)) xPred[1 + treated.indexOf(g)] = 1;
    xPred[pA - 1] = covMean;
    const adjMean = dot(xPred, fit.beta);
    const adjSe = Math.sqrt(Math.max(0, quadForm(xPred, fit.vcov)));
    adjustedMeans.push({
      group: g,
      rawMean: round4(rm),
      adjustedMean: round4(adjMean),
      n: grpVals.length,
      se: round4(adjSe),
    });
  }

  const pairwise: AncovaPairwise[] = [];
  for (let j = 0; j < treated.length; j++) {
    const cVec = new Array<number>(pA).fill(0);
    cVec[1 + j] = 1;
    const diff = dot(cVec, fit.beta);
    const seDiff = Math.sqrt(Math.max(0, quadForm(cVec, fit.vcov)));
    const tStat = seDiff > 0 ? diff / seDiff : 0;
    const pVal = 2 * (1 - studentTCdf(Math.abs(tStat), dfA));
    pairwise.push({
      group: treated[j],
      difference: round4(diff),
      se: round4(seDiff),
      tStatistic: round4(tStat),
      pValue: round6(pVal),
      significant: pVal < alpha,
    });
  }

  const sqrtMse = mse > 0 ? Math.sqrt(mse) : 1;
  const ctrlRaw = rawMeans.get(controlGroup) ?? 0;
  const effectDecomposition: AncovaEffectDecomposition[] = pairwise.map((pw) => {
    const total = (rawMeans.get(pw.group) ?? 0) - ctrlRaw;
    const direct = pw.difference;
    const indirect = total - direct;
    const propDirect = Math.abs(total) > 1e-10 ? direct / total : 1;
    const directD = sqrtMse > 0 ? direct / sqrtMse : 0;
    const jCorr = dfA > 1 ? 1 - 3 / (4 * dfA - 1) : 1;
    const directG = Math.abs(directD * jCorr);
    return {
      group: pw.group,
      totalEffect: round4(total),
      directEffect: round4(direct),
      indirectEffect: round4(indirect),
      proportionDirect: round4(propDirect),
      directG: round4(directG),
      directP: pw.pValue,
    };
  });

  return {
    adjustedMeans,
    pairwise,
    slope: {
      estimate: round6(slopeEst),
      se: round6(slopeSe),
      tStatistic: round4(slopeT),
      pValue: round6(slopeP),
    },
    slopeHomogeneity: {
      fStatistic: fHom !== null ? round4(fHom) : null,
      pValue: pHom !== null ? round6(pHom) : null,
      homogeneous,
    },
    effectDecomposition,
    modelRSquared: round4(rSq),
    mse: round6(mse),
    useOrganFreeBw,
    covariateMean: round4(covMean),
  };
}

// ── OLS via normal equations ─────────────────────────────────────

function fitOls(
  X: Matrix, y: number[],
): {beta: number[]; rss: number; df: number; vcov: Matrix} {
  const n = X.length;
  const p = X[0].length;
  const Xt = transpose(X);
  const XtX = matmul(Xt, X);
  const Xty = matvec(Xt, y);
  const XtXinv = inverse(XtX);
  const beta = matvec(XtXinv, Xty);
  const yPred = matvec(X, beta);
  let rss = 0;
  for (let i = 0; i < n; i++) {
    const r = y[i] - yPred[i];
    rss += r * r;
  }
  const df = n - p;
  const mse = df > 0 ? rss / df : Infinity;
  const vcov: Matrix = [];
  for (let i = 0; i < p; i++) {
    const row = new Array<number>(p);
    for (let j = 0; j < p; j++) row[j] = XtXinv[i][j] * mse;
    vcov.push(row);
  }
  return {beta, rss, df, vcov};
}

function dot(a: number[], b: number[]): number {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}

function round4(x: number): number {return Math.round(x * 1e4) / 1e4;}
function round6(x: number): number {return Math.round(x * 1e6) / 1e6;}
