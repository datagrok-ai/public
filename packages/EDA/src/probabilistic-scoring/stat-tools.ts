// Probabilistic scoring (pMPO) statistical tools
// Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//@ts-ignore: no types
import * as jStat from 'jstat';

import {Cutoff, DescriptorStatistics, SigmoidParams} from './pmpo-defs';

/** Splits the dataframe into desired and non-desired tables based on the desirability column */
export function getDesiredTables(df: DG.DataFrame, desirability: DG.Column) {
  const groups = df.groupBy([desirability.name]).getGroups() as any;
  let desired: DG.DataFrame;
  let nonDesired: DG.DataFrame;

  for (const name in groups) {
    if (name.toLowerCase().includes('true'))
      desired = groups[name];
    else
      nonDesired = groups[name];
  }

  //@ts-ignore
  return {desired, nonDesired};
} // getDesiredTables

/* Welch two-sample t-test (two-sided) */
export function getDescriptorStatistics(des: DG.Column, nonDes: DG.Column): DescriptorStatistics {
  const desLen = des.length;
  const nonDesLen = nonDes.length;
  if (desLen < 2 || nonDesLen < 2) {
    throw new Error(`Failed to compute the "${des.name}" descriptor statistics: 
      both samples must have at least two observations.`);
  }

  const desAvg = des.stats.avg;
  const nonDesAvg = nonDes.stats.avg;
  const desVar = des.stats.variance;
  const nonDesVar = nonDes.stats.variance;
  const desStd = des.stats.stdev;
  const nonDesStd = nonDes.stats.stdev;

  const se = Math.sqrt(desVar / desLen + nonDesVar / nonDesLen);
  if (se === 0) {
    throw new Error(`Failed to compute the "${des.name}" descriptor statistics: 
      zero variance.`);
  }

  const t = (desAvg - nonDesAvg) / se;

  // Welch–Satterthwaite degrees of freedom
  const numerator = (desVar / desLen + nonDesVar / nonDesLen) ** 2;
  const denom = (desVar * desVar) / (desLen * desLen * (desLen - 1)) +
    (nonDesVar * nonDesVar) / (nonDesLen * nonDesLen * (nonDesLen - 1));
  const df = numerator / denom;

  // two-sided p-value
  const cdf = jStat.studentt.cdf(Math.abs(t), df);
  const pValue = 2 * (1 - cdf);

  return {
    desAvg: desAvg,
    desStd: desStd,
    desLen: desLen,
    nonDesAvg: nonDesAvg,
    nonDesStd: nonDesStd,
    nonSesLen: nonDesLen,
    min: Math.min(des.stats.min, nonDes.stats.min),
    max: Math.max(des.stats.max, nonDes.stats.max),
    tstat: t,
    pValue: pValue,
  };
} // getDescriptorStatistics

/** Compute cutoffs for the pMPO method */
export function getCutoffs(muDesired: number, stdDesired: number, muNotDesired: number,
  stdNotDesired: number): Cutoff {
  if (muDesired < muNotDesired) {
    return {
      cutoff: ((muNotDesired - muDesired) / (stdDesired + stdNotDesired)) * stdDesired + muDesired,
      cutoffDesired: Math.max(muDesired, muNotDesired - stdNotDesired),
      cutoffNotDesired: Math.max(muDesired + stdDesired, muNotDesired),
    };
  } else {
    return {
      cutoff: ((muDesired - muNotDesired) / (stdDesired + stdNotDesired)) * stdNotDesired + muNotDesired,
      cutoffDesired: Math.min(muNotDesired + stdNotDesired, muDesired),
      cutoffNotDesired: Math.max(muNotDesired, muDesired - stdDesired),
    };
  }
} // getCutoffs

/** Solve normal intersection for the pMPO method */
export function solveNormalIntersection(mu1: number, s1: number, mu2: number, s2: number): number[] {
  const a = 1 / (2 * s1 ** 2) - 1 / (2 * s2 ** 2);
  const b = mu2 / (s2 ** 2) - mu1 / (s1 ** 2);
  const c = (mu1 ** 2) / (2 * s1 ** 2) - (mu2 ** 2) / (2 * s2 ** 2) - Math.log(s2 / s1);

  // If a is nearly zero, solve linear equation
  if (Math.abs(a) < 1e-12) {
    if (Math.abs(b) < 1e-12) return [];
    return [-c / b];
  }

  const disc = b * b - 4 * a * c;
  if (disc < 0) return [];

  const sqrtDisc = Math.sqrt(disc);
  const x1 = (-b + sqrtDisc) / (2 * a);
  const x2 = (-b - sqrtDisc) / (2 * a);

  return [x1, x2];
} // solveNormalIntersection

/** Compute sigmoid parameters for the pMPO method */
export function computeSigmoidParamsFromX0(muDes: number, sigmaDes: number, x0: number, xBound: number,
  qCutoff: number = 0.05): SigmoidParams {
  let pX0: number;

  if (sigmaDes <= 0)
    pX0 = x0 === muDes ? 1.0 : 0.0;
  else {
    // normal pdf
    const coef = 1 / (sigmaDes * Math.sqrt(2 * Math.PI));
    const exponent = -0.5 * ((x0 - muDes) / sigmaDes) ** 2;
    pX0 = coef * Math.exp(exponent);
  }

  const eps = 1e-12;
  pX0 = Math.max(pX0, eps);

  const b = Math.max(1.0 / pX0 - 1.0, eps);
  const n = 1.0 / qCutoff - 1.0;
  const dx = xBound - x0;

  let c: number;
  if (Math.abs(dx) < 1e-12)
    c = 1.0;
  else {
    const ratio = n / b;
    if (ratio <= 0)
      c = 1.0;
    else {
      try {
        c = Math.exp(-Math.log(ratio) / dx);
      } catch {
        c = 1.0;
      }
    }
  }

  return {pX0: pX0, b: b, c: c};
} // computeSigmoidParamsFromX0

/** Generalized sigmoid function */
export function sigmoidS(x: number, x0: number, b: number, c: number): number {
  if (c > 0)
    return 1.0 / (1.0 + b * (c ** (-(x - x0))));
  return 1.0/(1.0 + b);
}

/** Normal probability density function */
export function gaussDesirabilityFunc(x: number, mu: number, sigma: number): number {
  return Math.exp(-((x - mu)**2) / (2 * sigma**2));
}

