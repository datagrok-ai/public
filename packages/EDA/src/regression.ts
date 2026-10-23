// Regression tools

import * as DG from 'datagrok-api/dg';

import {cruxFit} from './crux/crux-model';
import {isNumeric, isRegressionTarget} from './crux/methods';

/** Work units per design-matrix cell of a fit */
const WORK_PER_CELL = 1000;

/** Applicability and interactivity thresholds for linear regression. */
enum LIN_REG_LIMITS {
  // isApplicable: hard bounds — fit must complete without exhausting browser memory
  MAX_FEATURES = 1000,
  MAX_FEATURES_X_SAMPLES = 1e7,
  // isInteractive: responsive (~0.5 s) budget
  INTERACTIVE_TIME_BUDGET = 2.5e8,
}

/** Check whether linear regression can be applied to the given data */
export function isLinearRegressionApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  if (!isNumeric(features) || !isRegressionTarget(target))
    return false;

  const featuresCount = features.length;
  const samplesCount = target.length;

  if (featuresCount > LIN_REG_LIMITS.MAX_FEATURES)
    return false;

  if (featuresCount * samplesCount > LIN_REG_LIMITS.MAX_FEATURES_X_SAMPLES)
    return false;

  return true;
}

/** Check whether linear regression training is fast enough for interactive use */
export function isLinearRegressionInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length * features.length * WORK_PER_CELL <= LIN_REG_LIMITS.INTERACTIVE_TIME_BUDGET;
}

/** Ordinary least squares, ridge (L2 only) or elastic net, on standardised features */
export function trainLinearRegression(df: DG.DataFrame, target: DG.Column, alpha: number, lambda: number,
  iterations: number): Promise<Uint8Array> {
  const {model, estimator}: {model: string, estimator: Record<string, number>} = alpha !== 0 ? {model: 'elasticnet',
    estimator: {alpha: alpha + lambda, l1_ratio: alpha / (alpha + lambda), max_iter: iterations}} :
    lambda !== 0 ? {model: 'ridge', estimator: {alpha: lambda}} : {model: 'linear', estimator: {}};
  return cruxFit('regression', model, df, target, {estimator, preprocess: 'standard'});
}
