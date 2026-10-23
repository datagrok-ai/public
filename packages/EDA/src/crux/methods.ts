import * as DG from 'datagrok-api/dg';

import {cruxFit} from './crux-model';

const SVM_KERNELS: Record<string, string> = {Linear: 'linear', Polynomial: 'poly', RBF: 'rbf', Sigmoid: 'sigmoid'};

/** INT or FLOAT: the only feature types the Crux bridge accepts */
export function isNumericColumn(col: DG.Column): boolean {
  return col.type === DG.COLUMN_TYPE.INT || col.type === DG.COLUMN_TYPE.FLOAT;
}

export function isNumeric(features: DG.ColumnList): boolean {
  return features.toList().every(isNumericColumn);
}

export function isClassTarget(target: DG.Column): boolean {
  return target.type === DG.COLUMN_TYPE.STRING || target.type === DG.COLUMN_TYPE.BOOL;
}

export function isRegressionTarget(target: DG.Column): boolean {
  return isNumericColumn(target);
}

/** Classes of a STRING / BOOL target; missing values are not a class (Crux skips their rows) */
function classCount(target: DG.Column): number {
  return target.type === DG.COLUMN_TYPE.BOOL ? new Set(target.toList()).size :
    target.categories.filter((c) => c !== '').length;
}

/** Rejects the inputs the pre-Crux methods rejected, with their messages */
function validate(method: string, rules: [boolean, string][]): void {
  for (const [ok, message] of rules) {
    if (!ok)
      throw new Error(`${method}: ${message}`);
  }
}

export function isSoftmaxApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && isClassTarget(target);
}

export function isSoftmaxInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return features.length <= 100 && target.length <= 50000;
}

/** Multinomial logistic regression; the old gradient-descent penalty maps to scikit-learn's C */
export function trainSoftmax(df: DG.DataFrame, target: DG.Column, iterations: number, penalty: number,
  tolerance: number): Promise<Uint8Array> {
  validate('Softmax', [[penalty > 0, 'penalty must be positive']]);
  const C = (classCount(target) <= 2 ? 2 : 1) / penalty;
  return cruxFit('classification', 'logistic', df, target,
    {estimator: {C, max_iter: iterations, tol: tolerance}, preprocess: 'standard'});
}

/** XGBoost and SVM: numeric features, a class or a regression target */
export function isXgboostApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && (isClassTarget(target) || isRegressionTarget(target));
}

export const isSvmApplicable = isXgboostApplicable;

export function isXgboostInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  const samples = target.length;
  const count = features.length;
  if (samples <= 10000)
    return count <= 100;
  if (samples <= 50000)
    return count <= 20;
  return samples <= 100000 && count <= 10;
}

export function trainXgboost(df: DG.DataFrame, target: DG.Column, iterations: number, eta: number, maxDepth: number,
  lambda: number, alpha: number): Promise<Uint8Array> {
  const classification = isClassTarget(target);
  validate('XGBoost', [
    [Number.isInteger(iterations) && iterations >= 1, 'iterations must be a positive integer'],
    [eta > 0 && eta <= 1, 'eta must be in (0, 1]'],
    [Number.isInteger(maxDepth) && maxDepth >= 1, 'max depth must be a positive integer'],
    [lambda >= 0, 'lambda must be non-negative'],
    [alpha >= 0, 'alpha must be non-negative'],
    [!classification || classCount(target) >= 2, 'classification target must have at least 2 categories'],
  ]);
  return cruxFit(classification ? 'classification' : 'regression', 'xgboost', df, target,
    {estimator: {n_estimators: iterations, learning_rate: eta, max_depth: maxDepth, reg_lambda: lambda,
      reg_alpha: alpha}});
}

export function isSvmInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length <= 10000 && features.length <= 50;
}

/** libsvm's C-SVC / epsilon-SVR; `gamma = 0` keeps libsvm's 1 / n_features */
export function trainSvm(df: DG.DataFrame, target: DG.Column, kernel: string, cost: number, gamma: number,
  degree: number, coef0: number, epsilon: number): Promise<Uint8Array> {
  validate('SVM', [[Number.isInteger(degree) && degree >= 1, 'degree must be a positive integer']]);
  const classification = isClassTarget(target);
  const estimator = {kernel: SVM_KERNELS[kernel] ?? 'rbf', C: cost, gamma: gamma === 0 ? 1 / df.columns.length : gamma,
    degree, coef0, ...(classification ? {} : {svr_epsilon: epsilon})};
  return cruxFit(classification ? 'classification' : 'regression', classification ? 'svc' : 'svr', df, target,
    {estimator, preprocess: 'standard'});
}
