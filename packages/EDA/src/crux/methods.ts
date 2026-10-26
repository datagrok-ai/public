import * as DG from 'datagrok-api/dg';

import {cruxFit, CruxTask} from './crux-model';

const SVM_KERNELS: Record<string, string> = {Linear: 'linear', Polynomial: 'poly', RBF: 'rbf', Sigmoid: 'sigmoid'};
const KNN_WEIGHTS: Record<string, string> = {Uniform: 'uniform', Distance: 'distance'};
const MLP_ACTIVATIONS: Record<string, string> =
  {ReLU: 'relu', Tanh: 'tanh', Logistic: 'logistic', Identity: 'identity'};
const REALMLP_VARIANTS: Record<string, string> = {'TD': 'td', 'TD-S': 'td_s'};
const TABM_EMBEDDINGS: Record<string, string> =
  {'Piecewise linear': 'piecewise-linear', 'Periodic': 'periodic', 'Linear-ReLU': 'linear-relu', 'None': 'none'};
const HIDDEN_LAYERS = /^[1-9]\d*(,[1-9]\d*)*$/;

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

const isAnyTarget = (t: DG.Column) => isClassTarget(t) || isRegressionTarget(t);

function taskOf(target: DG.Column): CruxTask {
  return isClassTarget(target) ? 'classification' : 'regression';
}

function cells(features: DG.ColumnList, target: DG.Column): number {
  return features.length * target.length;
}

function targetRows(target: DG.Column): number {
  return target.length - target.stats.missingValueCount;
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

function hiddenLayers(method: string, text: string, example: string): string {
  const layers = text.replace(/\s/g, '');
  validate(method, [[HIDDEN_LAYERS.test(layers),
    `hidden layers must be comma-separated positive integers, e.g. ${example}`]]);
  return layers;
}

const positiveInt = (v: number, name: string): [boolean, string] =>
  [Number.isInteger(v) && v >= 1, `${name} must be a positive integer`];

const validFraction = (f: number): [boolean, string] =>
  [f >= 0.05 && f <= 0.5, 'validation fraction must be from 0.05 to 0.5'];

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
  return isNumeric(features) && isAnyTarget(target);
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
    positiveInt(iterations, 'iterations'),
    [eta > 0 && eta <= 1, 'eta must be in (0, 1]'],
    positiveInt(maxDepth, 'max depth'),
    [lambda >= 0, 'lambda must be non-negative'],
    [alpha >= 0, 'alpha must be non-negative'],
    [!classification || classCount(target) >= 2, 'classification target must have at least 2 categories'],
  ]);
  return cruxFit(taskOf(target), 'xgboost', df, target,
    {estimator: {n_estimators: iterations, learning_rate: eta, max_depth: maxDepth, reg_lambda: lambda,
      reg_alpha: alpha}});
}

export function isSvmInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length <= 10000 && features.length <= 50;
}

/** libsvm's C-SVC / epsilon-SVR; `gamma = 0` keeps libsvm's 1 / n_features */
export function trainSvm(df: DG.DataFrame, target: DG.Column, kernel: string, cost: number, gamma: number,
  degree: number, coef0: number, epsilon: number): Promise<Uint8Array> {
  validate('SVM', [positiveInt(degree, 'degree')]);
  const classification = isClassTarget(target);
  const estimator = {kernel: SVM_KERNELS[kernel] ?? 'rbf', C: cost, gamma: gamma === 0 ? 1 / df.columns.length : gamma,
    degree, coef0, ...(classification ? {} : {svr_epsilon: epsilon})};
  return cruxFit(taskOf(target), classification ? 'svc' : 'svr', df, target,
    {estimator, preprocess: 'standard'});
}

export function isRandomForestApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && isAnyTarget(target) && cells(features, target) <= 1e6;
}

export function isRandomForestInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return cells(features, target) <= 2e4;
}

export function isExtraTreesApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && isAnyTarget(target) && cells(features, target) <= 2e6;
}

export const isExtraTreesInteractive = isRandomForestInteractive;

/** `maxFeatures = 'auto'` is scikit-learn's default: sqrt for classification, all features for regression */
function trainForest(model: string, df: DG.DataFrame, target: DG.Column, trees: number, maxDepth: number,
  minSamplesLeaf: number, maxFeatures: string, seed: number): Promise<Uint8Array> {
  const features = maxFeatures === 'auto' ? (isClassTarget(target) ? 'sqrt' : 'all') : maxFeatures;
  return cruxFit(taskOf(target), model, df, target, {estimator: {n_estimators: trees, max_depth: maxDepth,
    min_samples_leaf: minSamplesLeaf, max_features: features}, seed});
}

export function trainRandomForest(df: DG.DataFrame, target: DG.Column, trees: number, maxDepth: number,
  minSamplesLeaf: number, maxFeatures: string, seed: number): Promise<Uint8Array> {
  return trainForest('random_forest', df, target, trees, maxDepth, minSamplesLeaf, maxFeatures, seed);
}

export function trainExtraTrees(df: DG.DataFrame, target: DG.Column, trees: number, maxDepth: number,
  minSamplesLeaf: number, maxFeatures: string, seed: number): Promise<Uint8Array> {
  return trainForest('extra_trees', df, target, trees, maxDepth, minSamplesLeaf, maxFeatures, seed);
}

/** Hist Gradient Boosting, CatBoost, Decision Tree and Neural Network */
export function isHistGradientBoostingApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && isAnyTarget(target) && cells(features, target) <= 1e7;
}

export function isHistGradientBoostingInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length <= 20000 && cells(features, target) <= 2e5;
}

/** Early stopping is off: scikit-learn's `auto` silently holds out 10% of the rows above 10k */
export function trainHistGradientBoosting(df: DG.DataFrame, target: DG.Column, iterations: number,
  learningRate: number, maxLeafNodes: number, maxDepth: number, minSamplesLeaf: number, l2: number,
  seed: number): Promise<Uint8Array> {
  const estimator = {n_estimators: iterations, learning_rate: learningRate, max_leaf_nodes: maxLeafNodes,
    min_samples_leaf: minSamplesLeaf, l2_regularization: l2, early_stopping: 'off',
    ...(maxDepth > 0 ? {max_depth: maxDepth} : {})};
  return cruxFit(taskOf(target), 'hist_gradient_boosting', df, target, {estimator, seed});
}

export const isCatBoostApplicable = isHistGradientBoostingApplicable;
export function isCatBoostInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length <= 20000 && cells(features, target) <= 1e5;
}

/** CatBoost's own parameter names; Crux accepts iterations 0 and depth above 16, so they are checked here */
export function trainCatBoost(df: DG.DataFrame, target: DG.Column, iterations: number, learningRate: number,
  depth: number, l2LeafReg: number, seed: number): Promise<Uint8Array> {
  validate('CatBoost', [
    positiveInt(iterations, 'iterations'),
    [learningRate > 0 && learningRate <= 1, 'learning rate must be in (0, 1]'],
    [Number.isInteger(depth) && depth >= 1 && depth <= 16, 'depth must be an integer from 1 to 16'],
    [l2LeafReg >= 0, 'L2 must be non-negative'],
  ]);
  return cruxFit(taskOf(target), 'catboost', df, target,
    {estimator: {iterations, learning_rate: learningRate, depth, l2_leaf_reg: l2LeafReg}, seed});
}

export const isDecisionTreeApplicable = isHistGradientBoostingApplicable;

export function isDecisionTreeInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return cells(features, target) <= 2e5;
}

export function trainDecisionTree(df: DG.DataFrame, target: DG.Column, maxDepth: number, minSamplesLeaf: number,
  seed: number): Promise<Uint8Array> {
  return cruxFit(taskOf(target), 'decision_tree', df, target,
    {estimator: {max_depth: maxDepth, min_samples_leaf: minSamplesLeaf}, seed});
}

/** The model holds the training matrix: rows and cells are capped */
export function isKnnApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && isAnyTarget(target) && target.length <= 50000 && cells(features, target) <= 2.5e6;
}

export function isKnnInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length <= 10000 && cells(features, target) <= 1e5;
}

/** Crux fits `neighbors` above the training rows but fails at predict; Forge trains per fold, so the check runs on
 * the rows of the frame it gets */
export function trainKnn(df: DG.DataFrame, target: DG.Column, neighbors: number, weights: string):
  Promise<Uint8Array> {
  const rows = targetRows(target);
  validate('k-Nearest Neighbors', [[Number.isInteger(neighbors) && neighbors >= 1 && neighbors <= rows,
    `neighbors (${neighbors}) must be a positive integer not above the training rows (${rows})`]]);
  return cruxFit(taskOf(target), 'knn', df, target,
    {estimator: {n_neighbors: neighbors, weights: KNN_WEIGHTS[weights] ?? 'uniform'}, preprocess: 'standard'});
}

export function isNaiveBayesApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && isClassTarget(target) && cells(features, target) <= 1e7;
}

export function isNaiveBayesInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length <= 100000;
}

export function trainNaiveBayes(df: DG.DataFrame, target: DG.Column): Promise<Uint8Array> {
  return cruxFit('classification', 'gaussian_nb', df, target);
}

export const isLdaApplicable = isNaiveBayesApplicable;

export function isLdaInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length <= 100000 && features.length <= 100;
}

/** `shrinkage = 0` sends none (scikit-learn's svd solver) */
export function trainLda(df: DG.DataFrame, target: DG.Column, shrinkage: number): Promise<Uint8Array> {
  return cruxFit('classification', 'lda', df, target, {estimator: shrinkage > 0 ? {shrinkage} : {}});
}

export const isMlpApplicable = isHistGradientBoostingApplicable;

export function isMlpInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length <= 2000 && features.length <= 50;
}

/** `alpha` is always sent: Crux's classifier default is 1.0, not scikit-learn's 1e-4 */
export function trainMlp(df: DG.DataFrame, target: DG.Column, hiddenLayerSizes: string, activation: string,
  iterations: number, learningRate: number, alpha: number, seed: number): Promise<Uint8Array> {
  const layers = hiddenLayers('Neural Network', hiddenLayerSizes, '100 or 64,32');
  validate('Neural Network', [[alpha >= 0, 'L2 must be non-negative']]);
  return cruxFit(taskOf(target), 'mlp', df, target, {estimator: {hidden_layer_sizes: layers,
    activation: MLP_ACTIVATIONS[activation] ?? 'relu', max_iter: iterations, learning_rate: learningRate, alpha,
    tol: 1e-4}, preprocess: 'standard', seed});
}

export function isRealMlpApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && isAnyTarget(target) && target.length <= 10000 && cells(features, target) <= 1e5;
}

/** One fit takes 0.7-1 s on iris and Forge runs six per retrain */
export function isRealMlpInteractive(_features: DG.ColumnList, _target: DG.Column): boolean {
  return false;
}

/** `learningRate = 0` sends no `lr` (RealMLP's default: 0.04 for classification, 0.2 for regression) */
export function trainRealMlp(df: DG.DataFrame, target: DG.Column, variant: string, epochs: number,
  hiddenLayerSizes: string, learningRate: number, batchSize: number, validationFraction: number, seed: number):
  Promise<Uint8Array> {
  const layers = hiddenLayers('RealMLP', hiddenLayerSizes, '256,256,256 or 128,128');
  validate('RealMLP', [
    positiveInt(epochs, 'epochs'),
    positiveInt(batchSize, 'batch size'),
    [learningRate >= 0, 'learning rate must be non-negative (0 keeps the RealMLP default)'],
    validFraction(validationFraction),
  ]);
  return cruxFit(taskOf(target), 'realmlp', df, target, {estimator: {variant: REALMLP_VARIANTS[variant] ?? 'td',
    n_epochs: epochs, hidden_sizes: layers.split(',').map(Number), batch_size: batchSize,
    val_fraction: validationFraction, ...(learningRate > 0 ? {lr: learningRate} : {})}, seed});
}

export function isTabmApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && isAnyTarget(target) && target.length <= 20000 && cells(features, target) <= 2e5;
}

export function isTabmInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return target.length <= 200 && cells(features, target) <= 2000;
}

/** TabM ignores unknown estimator keys silently, so the names below must match crux's exactly */
export function trainTabm(df: DG.DataFrame, target: DG.Column, ensembleSize: number, width: number,
  embeddings: string, epochs: number, patience: number, learningRate: number, dropout: number,
  validationFraction: number, seed: number): Promise<Uint8Array> {
  validate('TabM', [
    positiveInt(ensembleSize, 'ensemble size'),
    positiveInt(width, 'width'),
    positiveInt(epochs, 'epochs'),
    positiveInt(patience, 'patience'),
    [learningRate > 0, 'learning rate must be positive'],
    [dropout >= 0 && dropout < 1, 'dropout must be from 0 to below 1'],
    validFraction(validationFraction),
  ]);
  return cruxFit(taskOf(target), 'tabm', df, target, {estimator: {k: ensembleSize, d_block: width,
    embeddings: TABM_EMBEDDINGS[embeddings] ?? 'piecewise-linear', max_epochs: epochs, patience, lr: learningRate,
    dropout, val_fraction: validationFraction}, seed});
}
