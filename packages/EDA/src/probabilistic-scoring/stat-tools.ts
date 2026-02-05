// Probabilistic scoring (pMPO) statistical tools
// Link: https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//@ts-ignore: no types
import * as jStat from 'jstat';

import {ConfusionMatrix, Cutoff, DescriptorStatistics, ModelEvaluationResult,
  ROC_TRESHOLDS, ROC_TRESHOLDS_COUNT, SigmoidParams} from './pmpo-defs';

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

/** Computes the confusion matrix given desirability (labels) and prediction columns
 * @param desirability - desirability column (boolean)
 * @param prediction - prediction column (numeric)
 * @param threshold - threshold to convert prediction scores to binary labels
 * @return ConfusionMatrix object with TP, TN, FP, FN counts
 */
export function getConfusionMatrix(desirability: DG.Column, prediction: DG.Column, threshold: number): ConfusionMatrix {
  if (desirability.length !== prediction.length)
    throw new Error('Failed to compute confusion matrix: columns have different lengths.');

  if (desirability.type !== DG.COLUMN_TYPE.BOOL)
    throw new Error('Failed to compute confusion matrix: desirability column must be boolean.');

  if (!prediction.isNumerical)
    throw new Error('Failed to compute confusion matrix: prediction column must be numerical.');

  let TP = 0;
  let TN = 0;
  let FP = 0;
  let FN = 0;

  const desRaw = desirability.getRawData();
  const predRaw = prediction.getRawData();

  let desIdx = 0;
  let curPos = 0;
  let desElem = desRaw[0];

  // Here, we extract bits from the desirability boolean column in chunks of 32 bits
  for (let predIdx = 0; predIdx < prediction.length; ++predIdx) {
    // console.log(predIdx + 1, ': ',
    //   desirability.get(predIdx), '<-->', (desElem >>> curPos) & 1, ' vs ', predRaw[predIdx] >= threshold);

    if (((desElem >>> curPos) & 1) == 1) { // True actual
      if (predRaw[predIdx] >= threshold) { // True predicted
        ++TP;
      } else { // False predicted
        ++FN;
      }
    } else { // False actual
      if (predRaw[predIdx] >= threshold) { // True predicted
        ++FP;
      } else { // False predicted
        ++TN;
      }
    }

    ++curPos;

    // Move to the next desirability element if we have processed 32 bits
    if (curPos >= 32) {
      curPos = 0;
      ++desIdx;
      desElem = desRaw[desIdx];
    }
  } // for predIdx

  return {TP: TP, TN: TN, FP: FP, FN: FN};
} // getConfusionMatrix

/** Computes Area Under Curve (AUC) given TPR and FPR arrays
 * @param tpr - True Positive Rate array
 * @param fpr - False Positive Rate array
 * @return AUC value
 */
export function getAuc(tpr: Float32Array, fpr: Float32Array): number {
  if (tpr.length !== fpr.length)
    throw new Error('Failed to compute AUC: TPR and FPR arrays have different lengths.');

  let auc = 0.0;

  for (let i = 1; i < tpr.length; ++i) {
    const xDiff = Math.abs(fpr[i] - fpr[i - 1]);
    const yAvg = (tpr[i] + tpr[i - 1]) / 2.0;
    auc += xDiff * yAvg;
  }

  return auc;
} // getAuc

/** Converts numeric prediction column to boolean based on the given threshold
 * @param numericPrediction - numeric prediction column
 * @param threshold - threshold to convert prediction scores to binary labels
 * @param name - name for the resulting boolean column
 * @return Boolean prediction column
 */
export function getBoolPredictionColumn(numericPrediction: DG.Column, threshold: number, name: string): DG.Column {
  if (!numericPrediction.isNumerical)
    throw new Error('Failed to compute confusion matrix: prediction column must be numerical.');

  const size = numericPrediction.length;
  const boolPredData = new Array<boolean>(size);
  const predRaw = numericPrediction.getRawData();

  for (let i = 0; i < size; ++i)
    boolPredData[i] = (predRaw[i] >= threshold);

  return DG.Column.fromList(DG.COLUMN_TYPE.BOOL, name, boolPredData);
} // getBoolPredictionColumn

/** Computes pMPO model evaluation metrics: AUC, optimal threshold, TPR and FPR arrays
 * @param desirability - desirability column (boolean)
 * @param prediction - prediction column (numeric)
 * @return ModelEvaluationResult object with AUC, optimal threshold, TPR and FPR arrays
 */
export function getPmpoEvaluation(desirability: DG.Column, prediction: DG.Column): ModelEvaluationResult {
  const tpr = new Float32Array(ROC_TRESHOLDS_COUNT);
  const fpr = new Float32Array(ROC_TRESHOLDS_COUNT);

  let bestJ = -1;
  let currentJ = -1;
  let bestThreshold = ROC_TRESHOLDS[0];

  // Compute TPR and FPR for each threshold
  for (let i = 0; i < ROC_TRESHOLDS_COUNT; ++i) {
    const confusion = getConfusionMatrix(desirability, prediction, ROC_TRESHOLDS[i]);
    tpr[i] = (confusion.TP + confusion.FN) > 0 ? confusion.TP / (confusion.TP + confusion.FN) : 0;
    fpr[i] = (confusion.FP + confusion.TN) > 0 ? confusion.FP / (confusion.FP + confusion.TN) : 0;
    currentJ = tpr[i] - fpr[i];

    if (currentJ > bestJ) {
      bestJ = currentJ;
      bestThreshold = ROC_TRESHOLDS[i];
    }
  }

  return {
    auc: getAuc(tpr, fpr),
    threshold: bestThreshold,
    tpr: tpr,
    fpr: fpr,
  };
} // getPmpoEvaluation
