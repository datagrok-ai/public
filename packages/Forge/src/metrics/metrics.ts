import * as DG from 'datagrok-api/dg';
import {ForgeError} from '../forge-error';
import {ModelTask} from '../generated/db';

export const METRIC_IDS = ['mse', 'rmse', 'mae', 'r2', 'accuracy', 'f1', 'sensitivity', 'specificity', 'precision',
  'npv', 'auc', 'balancedAccuracy', 'mcc', 'averagePrecision', 'logLoss', 'brier'] as const;
export type MetricId = typeof METRIC_IDS[number];
export type MetricValues = Partial<Record<MetricId, number>>;

/** The per-class ids of a classifier's `metrics.classes`: precision, sensitivity (recall), F1, AUC-ROC and average
 * precision of the class against the rest. */
export const CLASS_METRIC_IDS =
  ['precision', 'sensitivity', 'f1', 'auc', 'averagePrecision'] as const satisfies readonly MetricId[];

/** One class of a training: its rows among the measured ones and its values per split. */
export interface ClassMetrics { label: string; rows: number; train: MetricValues; validation: MetricValues }

/** One class of one measurement (an evaluation). */
export interface ClassValues { label: string; rows: number; values: MetricValues }

/** The ROC and precision-recall curves of the positive class, as Crux computes them: one point per distinct
 * probability; `roc.threshold[0]` is Infinity and the PR's last point (precision 1, recall 0) carries Infinity too. */
export interface Curves {
  roc: {fpr: Float32Array; tpr: Float32Array; threshold: Float32Array};
  pr: {precision: Float32Array; recall: Float32Array; threshold: Float32Array};
}

/** A classifier's class probabilities: one column per class, as the method's `proba_<class>` columns. */
export interface ClassProbabilities {
  classes: string[];
  columns: DG.Column[];
}

export function probabilityOf(probabilities: ClassProbabilities, cls: string): DG.Column | undefined {
  const k = probabilities.classes.indexOf(cls);
  return k < 0 ? undefined : probabilities.columns[k];
}

export const METRIC_LABELS: Record<MetricId, string> = {
  mse: 'MSE',
  rmse: 'RMSE',
  mae: 'MAE',
  r2: 'R2',
  accuracy: 'Accuracy',
  f1: 'F1 score',
  sensitivity: 'Recall (Sensitivity)',
  specificity: 'Specificity',
  precision: 'Precision',
  npv: 'Negative Predictive Value',
  auc: 'AUC-ROC',
  balancedAccuracy: 'Balanced accuracy',
  mcc: 'MCC',
  averagePrecision: 'Average precision',
  logLoss: 'Log loss',
  brier: 'Brier score',
};

export const METRIC_DESCRIPTIONS: Record<MetricId, string> = {
  mse: 'Mean Squared Error of values in the column. It is the average of the squared differences between ' +
    'predicted and actual values. Lower values indicate better model accuracy.',
  rmse: 'Root Mean Squared Error of values in the column. It represents the square root of the average squared ' +
    'differences between predicted and actual values. Lower values indicate better model performance.',
  mae: 'Mean Absolute Error of values in the column. It is the average of the absolute differences between ' +
    'predicted and actual values. Lower values indicate better model accuracy.',
  r2: 'R2 value of the column. It represents the proportion of the variance for a dependent variable that ' +
    'is explained by an independent variable. Values close to 1 indicate a good fit.',
  accuracy: 'Accuracy of the model. It measures the proportion of true results (both true positives and true ' +
    'negatives) among the total number of cases examined. Higher values indicate better overall performance.',
  f1: 'F1 score of the model. It is the harmonic mean of precision and recall; for more than two classes, ' +
    'the average of the per-class F1 scores. Higher values indicate better overall performance.',
  sensitivity: 'Recall (sensitivity) of the model. It measures the proportion of actual positives that are ' +
    'correctly identified. Higher values indicate fewer false negatives.',
  specificity: 'Specificity of the model. It measures the proportion of actual negatives that are correctly ' +
    'identified. Higher values indicate fewer false positives.',
  precision: 'Precision of the model. It measures the proportion of true positives among the total number of ' +
    'positive predictions. Higher values indicate fewer false positives.',
  npv: 'Negative Predictive Value of the model. It measures the proportion of true negatives among the total number ' +
    'of negative predictions. Higher values indicate fewer false negatives.',
  auc: 'How well the model separates the classes: 1 is perfect, 0.5 is random guessing. With more than two classes, ' +
    'the average over the classes.',
  balancedAccuracy: 'The average over the classes of the share of each class\'s rows predicted correctly: 1 is ' +
    'perfect. Unlike Accuracy, a large class does not inflate it.',
  mcc: 'Matthews correlation coefficient: how well the predicted classes agree with the actual ones over all ' +
    'classes. 1 is perfect, 0 is random guessing, -1 is total disagreement.',
  averagePrecision: 'The area under the precision-recall curve of the class probabilities: 1 is perfect, the ' +
    'positive class\'s share of the rows is random guessing. With more than two classes, the average over the classes.',
  logLoss: 'Log loss: the average of -ln(probability of the actual class): 0 is perfect, lower is better; a ' +
    'confident wrong prediction costs most.',
  brier: 'Brier score: the mean squared difference between the class probabilities and the actual class (0 or 1): ' +
    '0 is perfect, lower is better; 0.25 is guessing 0.5 for two classes. With more than two classes, summed over ' +
    'the classes.',
};

/** The metrics the classification threshold can be optimized for (**Best threshold by**, `thresholdMetric`). */
export const BEST_THRESHOLD_METRICS =
  ['f1', 'balancedAccuracy', 'mcc', 'accuracy'] as const satisfies readonly MetricId[];
export type BestThresholdMetric = typeof BEST_THRESHOLD_METRICS[number];

/** A training's two measurements: the final model on its training rows, and the out-of-fold predictions. */
export type Split = 'train' | 'validation';

/** A training's train and validation values and, since schema 0.1.25, the trivial model's values on its training
 * rows. */
export interface SplitMetrics { train: MetricValues; validation: MetricValues; baseline?: MetricValues }

/** The metrics judged for underfitting by another one: the shares of one class by Balanced accuracy, the errors of a
 * regressor by R2. */
const JUDGED_BY: Partial<Record<MetricId, MetricId>> = {sensitivity: 'balancedAccuracy',
  specificity: 'balancedAccuracy', precision: 'balancedAccuracy', npv: 'balancedAccuracy', mse: 'r2', rmse: 'r2',
  mae: 'r2'};
/** The share of the way from the trivial model to a perfect one below which a train value is underfitting. */
const MIN_SKILL = 0.1;

/** The metrics for which lower is better: the errors. */
export const MINIMIZED_METRICS: ReadonlySet<MetricId> = new Set<MetricId>(['mse', 'rmse', 'mae', 'logLoss', 'brier']);

const SCORE_GAP_LIMIT = 0.05;
const CORRELATION_GAP_LIMITS: Partial<Record<MetricId, number>> = {mcc: 0.1, r2: 0.1};
const ERROR_RATIO_LIMIT = 1.5;

export interface Overfitting {
  /** `gap` for a score, `ratio` for an error. */
  kind: 'gap' | 'ratio';
  /** A score's train − validation: positive is worse; undefined for an error. */
  gap?: number;
  /** An error's validation / train, 2 decimals; undefined for a score and for a train value of 0. */
  ratio?: number;
  isOverfitting: boolean;
}

/** The user's overfitting rule, on the values rounded as shown; it replaces the built-in tool's areValuesSimilar. */
export function overfittingOf(id: MetricId, train: number, validation: number): Overfitting {
  const [t, v] = [roundTo3(train), roundTo3(validation)];
  if (!MINIMIZED_METRICS.has(id)) {
    const gap = roundTo3(t - v);
    return {kind: 'gap', gap, isOverfitting: gap > (CORRELATION_GAP_LIMITS[id] ?? SCORE_GAP_LIMIT)};
  }
  const ratio = t === 0 ? undefined : Number((v / t).toFixed(2));
  return {kind: 'ratio', ratio, isOverfitting: ratio === undefined ? v > 0 : ratio > ERROR_RATIO_LIMIT};
}

/** [x] rounded to the 3 decimals the values are shown with. */
export function roundTo3(x: number): number {
  return Number(x.toFixed(3));
}

/** The underfitting of a metric: [id] is the metric judged ({@link JUDGED_BY}), [train] and [baseline] its values
 * rounded as shown, [skill] the share of the way from the trivial model to a perfect one. */
export interface Underfitting { id: MetricId; train: number; baseline: number; skill: number; isUnderfitting: boolean }

/** The trivial model's values on the training rows of [target]: for a classifier from the shares p_k of the K classes
 * with rows (π the share of [positiveClass] for two classes), for a regressor R2 0. */
export function baselineMetrics(task: ModelTask, target: DG.Column, positiveClass?: string): MetricValues {
  if (task === 'regression')
    return {r2: 0};
  const counts = [...classCounts(target)];
  const total = counts.reduce((sum, [, n]) => sum + n, 0);
  const k = counts.length;
  if (k === 0)
    return {};
  const shares = counts.map(([, n]) => n / total);
  const positive = positiveClass === undefined || k !== 2 ? undefined :
    (counts.find(([label]) => label === positiveClass)?.[1] ?? 0) / total;
  const share = positive ?? 1 / k;
  return {
    balancedAccuracy: 1 / k,
    accuracy: Math.max(...shares),
    mcc: 0,
    auc: 0.5,
    averagePrecision: share,
    f1: share,
    logLoss: -shares.reduce((sum, p) => sum + p * Math.log(p), 0),
    brier: positive === undefined ? shares.reduce((sum, p) => sum + p * (1 - p), 0) : positive * (1 - positive),
  };
}

/** The underfitting of [id] by the metric it is judged by: its train value against the [baseline] one; undefined
 * without either, or when the baseline is already perfect. */
export function underfittingOf(id: MetricId, train: MetricValues, baseline: MetricValues): Underfitting | undefined {
  const judged = JUDGED_BY[id] ?? id;
  const [trainValue, baseValue] = [train[judged], baseline[judged]];
  if (trainValue === undefined || baseValue === undefined)
    return undefined;
  const [t, b] = [roundTo3(trainValue), roundTo3(baseValue)];
  const isError = MINIMIZED_METRICS.has(judged);
  if (isError ? b <= 0 : b >= 1)
    return undefined;
  const skill = roundTo3(isError ? 1 - t / b : (t - b) / (1 - b));
  return {id: judged, train: t, baseline: b, skill, isUnderfitting: skill < MIN_SKILL};
}

export type VerdictStatus = 'underfitting' | 'overfitting' | 'noOverfitting';
export interface Verdict { status: VerdictStatus; metric: MetricId; underfitting?: Underfitting;
  overfitting?: Overfitting }

/** Underfitting first, then overfitting by [metric] itself, then no overfitting; undefined when neither can be
 * judged. */
export function verdictOf(metrics: SplitMetrics, metric: MetricId): Verdict | undefined {
  const {train, validation, baseline} = metrics;
  const underfitting = baseline === undefined ? undefined : underfittingOf(metric, train, baseline);
  const [t, v] = [train[metric], validation[metric]];
  const overfitting = t === undefined || v === undefined ? undefined : overfittingOf(metric, t, v);
  if (underfitting?.isUnderfitting)
    return {status: 'underfitting', metric, underfitting, overfitting};
  if (overfitting === undefined)
    return undefined;
  return {status: overfitting.isOverfitting ? 'overfitting' : 'noOverfitting', metric, underfitting, overfitting};
}

/** The metric a model is judged and compared by until the user picks one: [searchMetric] (the one its search
 * optimized), else Balanced accuracy, else R2 (the tasks' default search metrics), else the first of [ids]; only one
 * of [ids]. */
export function defaultMetricOf(ids: readonly MetricId[], searchMetric?: MetricId): MetricId | undefined {
  const candidates: (MetricId | undefined)[] = [searchMetric, 'balancedAccuracy', 'r2'];
  return candidates.find((id): id is MetricId => id !== undefined && ids.includes(id)) ?? ids[0];
}

/** The verdict metric of one model until the user picks one: the {@link defaultMetricOf} the metrics with a
 * validation or test value, as the model cards pick it. */
export function verdictMetricOf({validation, test}: {validation: MetricValues; test?: MetricValues},
  searchMetric?: MetricId): MetricId | undefined {
  return defaultMetricOf(METRIC_IDS.filter((id) => validation[id] !== undefined || test?.[id] !== undefined),
    searchMetric);
}

/** `guessing`, or `predicting the mean` when [id] is judged by R2. */
export function trivialModelOf(id: MetricId): string {
  return (JUDGED_BY[id] ?? id) === 'r2' ? 'predicting the mean' : 'guessing';
}

/** The rows of each class of [col] that has any, empty cells aside, in the column's category order. */
export function classCounts(col: DG.Column): Map<string, number> {
  const counts = new Map<string, number>();
  if (col.type === DG.COLUMN_TYPE.STRING) {
    // A text column's raw data are the indexes of its categories.
    const categories = col.categories;
    const indexes = col.getRawData();
    const tally = new Int32Array(categories.length);
    for (let i = 0; i < col.length; i++)
      tally[indexes[i]]++;
    for (let k = 0; k < categories.length; k++) {
      if (categories[k] !== '' && tally[k] > 0)
        counts.set(categories[k], tally[k]);
    }
    return counts;
  }
  for (let i = 0; i < col.length; i++) {
    if (!col.isNone(i)) {
      const value = String(col.get(i));
      counts.set(value, (counts.get(value) ?? 0) + 1);
    }
  }
  return counts;
}

export function metricOfLabel(label: string | null): MetricId | undefined {
  return METRIC_IDS.find((id) => METRIC_LABELS[id] === label);
}

/** A metric value worth a word: a train value barely better than the trivial model's ({@link underfittingOf}), or a
 * validation value worse than the train one by the overfitting rule ({@link overfittingOf}). */
export interface MetricFinding {
  id: MetricId;
  kind: 'underfitting' | 'mismatch';
  /** The splits whose values the finding marks. */
  splits: Split[];
  text: string;
  detail: string;
  /** The underfitting an `underfitting` finding reports. */
  underfitting?: Underfitting;
}

interface ClassCounts { tp: number; fp: number; fn: number }

export function regressionMetrics(actual: DG.Column, predicted: DG.Column): MetricValues {
  const rows = pairedRows(actual, predicted);
  const n = rows.length;
  let sum = 0;
  for (const i of rows)
    sum += actual.getNumber(i);
  const mean = sum / n;
  let ssRes = 0;
  let ssTot = 0;
  let absSum = 0;
  for (const i of rows) {
    const y = actual.getNumber(i);
    const error = y - predicted.getNumber(i);
    ssRes += error * error;
    absSum += Math.abs(error);
    ssTot += (y - mean) * (y - mean);
  }
  const mse = ssRes / n;
  // A constant target has no variance: scikit-learn's convention.
  const r2 = ssTot === 0 ? (ssRes === 0 ? 1 : 0) : 1 - ssRes / ssTot;
  return {mse, rmse: Math.sqrt(mse), mae: absSum / n, r2};
}

/** The rows where both have a value, and each label's counts over them (compared as text). */
function labelCounts(actual: DG.Column, predicted: DG.Column): {n: number; counts: Map<string, ClassCounts>} {
  const rows = pairedRows(actual, predicted);
  const counts = new Map<string, ClassCounts>();
  const countsOf = (label: string): ClassCounts => {
    let c = counts.get(label);
    if (c === undefined) {
      c = {tp: 0, fp: 0, fn: 0};
      counts.set(label, c);
    }
    return c;
  };
  for (const i of rows) {
    const a = String(actual.get(i));
    const p = String(predicted.get(i));
    if (a === p)
      countsOf(a).tp++;
    else {
      countsOf(a).fn++;
      countsOf(p).fp++;
    }
  }
  return {n: rows.length, counts};
}

/** Accuracy, balanced accuracy (the mean recall over the actual classes), MCC (multiclass, as scikit-learn), F1 and,
 * for at most two classes, the shares of [positiveClass]. */
export function classificationMetrics(actual: DG.Column, predicted: DG.Column,
  positiveClass: string): MetricValues {
  const {n, counts} = labelCounts(actual, predicted);
  const metrics = countMetrics([...counts.values()], n);
  if (counts.size > 2) {
    let f1Sum = 0;
    for (const c of counts.values())
      f1Sum += f1Of(c);
    metrics.f1 = f1Sum / counts.size;
    return metrics;
  }

  const positive = counts.get(positiveClass) ?? {tp: 0, fp: 0, fn: 0};
  const {tp, fp, fn} = positive;
  const tn = n - tp - fp - fn;
  const shares: [MetricId, number, number][] = [
    ['sensitivity', tp, tp + fn],
    ['specificity', tn, tn + fp],
    ['precision', tp, tp + fp],
    ['npv', tn, tn + fn],
  ];
  for (const [id, numerator, denominator] of shares) {
    if (denominator > 0)
      metrics[id] = numerator / denominator;
  }
  metrics.f1 = f1Of(positive);
  return metrics;
}

/** Accuracy, balanced accuracy (the mean recall over the actual classes) and MCC (multiclass, as scikit-learn) of each
 * class's [counts] over [n] rows. */
function countMetrics(counts: ClassCounts[], n: number): MetricValues {
  let correct = 0;
  let recallSum = 0;
  let actualClasses = 0;
  let trueTimesPredicted = 0;
  let trueSquares = 0;
  let predictedSquares = 0;
  for (const c of counts) {
    correct += c.tp;
    const [t, p] = [c.tp + c.fn, c.tp + c.fp];
    if (t > 0) {
      recallSum += c.tp / t;
      actualClasses++;
    }
    trueTimesPredicted += t * p;
    trueSquares += t * t;
    predictedSquares += p * p;
  }
  const covariance = correct * n - trueTimesPredicted;
  const spread = (n * n - trueSquares) * (n * n - predictedSquares);
  return {accuracy: correct / n, balancedAccuracy: recallSum / actualClasses,
    mcc: spread === 0 ? 0 : covariance / Math.sqrt(spread)};
}

/** Each label's precision, sensitivity (recall) and F1 against the rest, as EDA's report gives them per class (a
 * share without rows left out). */
export function classLabelMetrics(actual: DG.Column, predicted: DG.Column): Map<string, MetricValues> {
  const values = new Map<string, MetricValues>();
  for (const [label, c] of labelCounts(actual, predicted).counts) {
    const v: MetricValues = {};
    if (c.tp + c.fp > 0)
      v.precision = c.tp / (c.tp + c.fp);
    if (c.tp + c.fn > 0)
      v.sensitivity = c.tp / (c.tp + c.fn);
    v.f1 = f1Of(c);
    values.set(label, v);
  }
  return values;
}

/** The classes of [target], its empty category aside. */
export function classesOf(target: DG.Column): string[] {
  return target.categories.filter((c) => c !== '');
}

/** The metrics of [predicted]; with a classifier's [probabilities], also `auc` ({@link aucOfProbabilities}), `logLoss`
 * and `brier` ({@link probabilityMetrics}). */
export function metricsOf(task: ModelTask, actual: DG.Column, predicted: DG.Column, positiveClass?: string,
  probabilities?: ClassProbabilities): MetricValues {
  if (task === 'regression')
    return regressionMetrics(actual, predicted);
  const values = classificationMetrics(actual, predicted, positiveClass ?? actual.categories[0]);
  if (probabilities === undefined)
    return values;
  const auc = aucOfProbabilities(actual, probabilities, positiveClass);
  return {...values, ...(auc === undefined ? {} : {auc}), ...probabilityMetrics(actual, probabilities, positiveClass)};
}

const MIN_PROBABILITY = 1e-15;

/** Log loss and the Brier score of [probabilities] against [actual], over the rows where both have values (classes
 * compared as text): log loss is the mean of -ln of the actual class's probability (clipped to [1e-15, 1 - 1e-15]; 0
 * for a class without a column); Brier is, for two classes with a [positiveClass], the mean squared difference of its
 * probability and 1 / 0, otherwise the mean over the rows of that difference summed over the classes. */
export function probabilityMetrics(actual: DG.Column, probabilities: ClassProbabilities,
  positiveClass?: string): MetricValues {
  const {classes, columns} = probabilities;
  if (columns.length === 0)
    return {};
  const positive = positiveClass === undefined || classes.length !== 2 ? -1 : classes.indexOf(positiveClass);
  let rows = 0;
  let logLoss = 0;
  let brier = 0;
  for (let i = 0; i < actual.length; i++) {
    if (actual.isNone(i) || columns.some((c) => c.isNone(i)))
      continue;
    const label = String(actual.get(i));
    const k = classes.indexOf(label);
    const p = k < 0 ? 0 : columns[k].getNumber(i);
    logLoss -= Math.log(Math.min(Math.max(p, MIN_PROBABILITY), 1 - MIN_PROBABILITY));
    if (positive >= 0)
      brier += (columns[positive].getNumber(i) - (k === positive ? 1 : 0)) ** 2;
    else {
      for (const [j, col] of columns.entries())
        brier += (col.getNumber(i) - (j === k ? 1 : 0)) ** 2;
    }
    rows++;
  }
  return rows === 0 ? {} : {logLoss: logLoss / rows, brier: brier / rows};
}

/** The threshold on [probability] (of [positiveClass]) that maximizes [metric] on the rows where [actual] and
 * [probability] have values: each distinct probability is a candidate (`p >= t` is positive; ties: the candidate
 * nearest 0.5); the value is the midpoint between the best candidate and the next lower probability (the candidate
 * itself when it is the lowest), rounded to 3 decimals when the rounding stays in that gap; the gap is checked in
 * single precision, as the class-probabilities cut compares; undefined when the rows hold one class only. */
export function bestThreshold(actual: DG.Column, probability: DG.Column, positiveClass: string,
  metric: BestThresholdMetric): number | undefined {
  const rows: {p: number; isPositive: boolean}[] = [];
  for (let i = 0; i < actual.length; i++) {
    if (!actual.isNone(i) && !probability.isNone(i))
      rows.push({p: probability.getNumber(i), isPositive: String(actual.get(i)) === positiveClass});
  }
  rows.sort((a, b) => b.p - a.p);
  const positives = rows.filter((r) => r.isPositive).length;
  const negatives = rows.length - positives;
  if (positives === 0 || negatives === 0)
    return undefined;
  const candidates: number[] = [];
  let best = -1;
  let bestValue = -Infinity;
  let tp = 0;
  let fp = 0;
  for (let k = 0; k < rows.length;) {
    const t = rows[k].p;
    for (; k < rows.length && rows[k].p === t; k++) {
      if (rows[k].isPositive)
        tp++;
      else
        fp++;
    }
    candidates.push(t);
    const positive = {tp, fp, fn: positives - tp};
    const value = metric === 'f1' ? f1Of(positive) : countMetrics([positive,
      {tp: negatives - fp, fp: positives - tp, fn: fp}], rows.length)[metric]!;
    const isTie = value === bestValue && Math.abs(t - 0.5) < Math.abs(candidates[best] - 0.5);
    if (value > bestValue || isTie) {
      bestValue = value;
      best = candidates.length - 1;
    }
  }
  const t = candidates[best];
  if (best === candidates.length - 1)
    return t;
  const next = candidates[best + 1];
  const mid = (t + next) / 2;
  const rounded = Number(mid.toFixed(3));
  const cuts = (value: number) => isAtThreshold(t, value) && !isAtThreshold(next, value);
  return cuts(rounded) ? rounded : cuts(mid) ? mid : t;
}

/** Probability [p] reaches [threshold], compared in single precision, as the probabilities are stored. */
export function isAtThreshold(p: number, threshold: number): boolean {
  return p >= Math.fround(threshold);
}

/** A numerical target makes a regression, any other a classification. */
export function taskOf(target: DG.Column): ModelTask {
  return target.matches(DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME) ? 'regression' : 'classification';
}

/** The area under the ROC curve of [score] for [positiveClass] of [actual]: the trapezoid over the rows sorted by
 * score, highest first, with the rows of one score as one step (a diagonal segment, as scikit-learn; the built-in
 * tool stepped row by row, so ties counted in their row order); undefined when the rows hold only one of the
 * classes. */
export function aucOf(actual: DG.Column, score: DG.Column, positiveClass: string): number | undefined {
  const paired = pairedRows(actual, score);
  const values = Float64Array.from(paired, (i) => score.getNumber(i));
  const order = Array.from(paired.keys()).sort((a, b) => values[b] - values[a]);
  const scores = order.map((k) => values[k]);
  const isPositive = order.map((k) => String(actual.get(paired[k])) === positiveClass);
  const positives = isPositive.filter((p) => p).length;
  const negatives = order.length - positives;
  if (positives === 0 || negatives === 0)
    return undefined;
  let auc = 0;
  let tp = 0;
  let fp = 0;
  for (let k = 0; k < order.length;) {
    const [tpr, fpr] = [tp / positives, fp / negatives];
    for (const tied = scores[k]; k < order.length && scores[k] === tied; k++) {
      if (isPositive[k])
        tp++;
      else
        fp++;
    }
    auc += (fp / negatives - fpr) * (tp / positives + tpr) / 2;
  }
  return auc;
}

/** AUC-ROC of [probabilities] against [actual]: with a [positiveClass], the AUC of its probability; otherwise the
 * macro average of the one-vs-rest AUCs over the classes whose rows hold both classes; undefined when none does. */
export function aucOfProbabilities(actual: DG.Column, probabilities: ClassProbabilities,
  positiveClass?: string): number | undefined {
  if (positiveClass !== undefined) {
    const column = probabilityOf(probabilities, positiveClass);
    return column === undefined ? undefined : aucOf(actual, column, positiveClass);
  }
  const aucs = probabilities.classes.map((cls, k) => aucOf(actual, probabilities.columns[k], cls))
    .filter((auc): auc is number => auc !== undefined);
  return aucs.length === 0 ? undefined : aucs.reduce((sum, auc) => sum + auc, 0) / aucs.length;
}

/** The built-in tool's closeness of two values: within 0.1 when both are small (their mean magnitude below 1),
 * otherwise within 90% of that mean, so that only values several times apart differ (the statistics insight; the
 * metric findings use overfittingOf and underfittingOf). */
export function areValuesSimilar(a: number, b: number): boolean {
  const mean = (Math.abs(a) + Math.abs(b)) / 2;
  return mean < 1 ? Math.abs(a - b) <= 0.1 : Math.abs(a - b) / mean <= 0.9;
}

/** For each metric: a train value barely better than the trivial model's (the metrics with a baseline of their own,
 * {@link underfittingOf}), then a validation value worse than the train one by the overfitting rule
 * ({@link overfittingOf}). */
export function metricFindings(metrics: SplitMetrics): MetricFinding[] {
  const findings: MetricFinding[] = [];
  const {baseline} = metrics;
  for (const id of METRIC_IDS) {
    const train = metrics.train[id];
    const validation = metrics.validation[id];
    const label = METRIC_LABELS[id];
    const underfitting = baseline === undefined || JUDGED_BY[id] !== undefined ? undefined :
      underfittingOf(id, metrics.train, baseline);
    if (underfitting?.isUnderfitting) {
      const trivial = trivialModelOf(id);
      findings.push({id, kind: 'underfitting', splits: ['train'],
        text: `${label} is barely better than ${trivial} on the training rows.`,
        detail: `Train ${fixed(underfitting.train)}, ${trivial} ${fixed(underfitting.baseline)}.`, underfitting});
    }
    if (train !== undefined && validation !== undefined && overfittingOf(id, train, validation).isOverfitting) {
      findings.push({id, kind: 'mismatch', splits: ['validation'],
        detail: `Train ${fixed(train)}, validation ${fixed(validation)}.`,
        text: `The validation ${label} is far from the training one.`});
    }
  }
  return findings;
}

/** [value] with 3 decimals, `-` when there is none. */
export function fixed(value: number | undefined): string {
  return value === undefined ? '-' : value.toFixed(3);
}

function pairedRows(actual: DG.Column, predicted: DG.Column): number[] {
  const rows: number[] = [];
  for (let i = 0; i < actual.length; i++) {
    if (!actual.isNone(i) && !predicted.isNone(i))
      rows.push(i);
  }
  if (rows.length === 0)
    throw new ForgeError('No row has both an actual and a predicted value, so the quality cannot be measured.');
  return rows;
}

function f1Of(c: ClassCounts): number {
  const precision = c.tp + c.fp > 0 ? c.tp / (c.tp + c.fp) : 0;
  const recall = c.tp + c.fn > 0 ? c.tp / (c.tp + c.fn) : 0;
  return precision + recall > 0 ? 2 * precision * recall / (precision + recall) : 0;
}
