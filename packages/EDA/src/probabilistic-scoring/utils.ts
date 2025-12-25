import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export type DescriptorStatistics = {
  desAvg: number,
  desStd: number,
  desLen: number,
  nonDesAvg: number,
  nonDesStd: number,
  nonSesLen: number,
  tstat: number,
  pValue: number,
  zScore: number,
};

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

/* Lanczos approximation for log gamma */
function logGamma(z: number): number {
  const g = 7;
  const p = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
  ];
  if (z < 0.5) {
    // Reflection formula
    return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - logGamma(1 - z);
  }
  z -= 1;
  let x = p[0];
  for (let i = 1; i < p.length; i++) x += p[i] / (z + i);
  const t = z + g + 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
} // logGamma

/* Continued fraction for incomplete beta (betacf) from Numerical Recipes */
function betacf(a: number, b: number, x: number): number {
  const MAXIT = 200;
  const EPS = 3e-14;
  const FPMIN = 1e-300;

  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;
  let c = 1;
  let d = 1 - (qab * x) / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;

  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;

    // even step
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h = h * d * c;

    // odd step
    aa = -((a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h = h * del;

    if (Math.abs(del - 1.0) < EPS) break;
  }

  return h;
} // betacf

/* Regularized incomplete beta I_x(a,b) */
function regularizedIncompleteBeta(x: number, a: number, b: number): number {
  if (x <= 0) return 0;
  if (x >= 1) return 1;

  // Use symmetry to improve convergence
  const bt =
    Math.exp(
      logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x),
    );

  let result: number;
  if (x < (a + 1) / (a + b + 2))
    result = (bt * betacf(a, b, x)) / a;
  else {
    // Use symmetry relation
    result = 1 - (bt * betacf(b, a, 1 - x)) / b;
  }

  // Clamp numerical errors
  if (result < 0) result = 0;
  if (result > 1) result = 1;
  return result;
} // regularizedIncompleteBeta

/* CDF of Student's t-distribution with df degrees of freedom */
function studentTCDF(t: number, df: number): number {
  // Uses relationship to incomplete beta
  const x = df / (df + t * t);
  const a = df / 2;
  const b = 0.5;
  const ib = regularizedIncompleteBeta(x, a, b);
  let cdf: number;
  if (t >= 0)
    cdf = 1 - 0.5 * ib;
  else
    cdf = 0.5 * ib;

  // Clamp
  if (cdf < 0) cdf = 0;
  if (cdf > 1) cdf = 1;
  return cdf;
} // studentTCDF

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
  const cdf = studentTCDF(Math.abs(t), df);
  const pValue = 2 * (1 - cdf);

  return {
    desAvg: desAvg,
    desStd: desStd,
    desLen: desLen,
    nonDesAvg: nonDesAvg,
    nonDesStd: nonDesStd,
    nonSesLen: nonDesLen,
    tstat: t,
    pValue: pValue,
    zScore: Math.abs(desAvg - nonDesAvg) / (desStd + nonDesStd),
  };
} // getDescriptorStatistics
