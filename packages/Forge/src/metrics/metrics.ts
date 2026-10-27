import * as DG from 'datagrok-api/dg';
import {ForgeError} from '../forge-error';
import {ModelTask} from '../generated/db';

export const METRIC_IDS =
  ['mse', 'rmse', 'mae', 'r2', 'accuracy', 'f1', 'sensitivity', 'specificity', 'precision', 'npv', 'auc'] as const;
export type MetricId = typeof METRIC_IDS[number];
export type MetricValues = Partial<Record<MetricId, number>>;

/** A classifier's class probabilities: one column per class, as the method's `proba_<class>` columns. */
export interface ClassProbabilities {
  classes: string[];
  columns: DG.Column[];
}

export function probabilityOf(probabilities: ClassProbabilities, cls: string): DG.Column | undefined {
  const k = probabilities.classes.indexOf(cls);
  return k < 0 ? undefined : probabilities.columns[k];
}

/** The column a two-class result's ROC curve is drawn from: [positiveClass]'s probability when [probabilities] has
 * exactly two columns; undefined otherwise. */
export function rocProbabilityOf(probabilities: ClassProbabilities | undefined,
  positiveClass: string): DG.Column | undefined {
  return probabilities === undefined || probabilities.classes.length !== 2 ? undefined :
    probabilityOf(probabilities, positiveClass);
}

export const METRIC_LABELS: Record<MetricId, string> = {
  mse: 'MSE',
  rmse: 'RMSE',
  mae: 'MAE',
  r2: 'R2',
  accuracy: 'Accuracy',
  f1: 'F1',
  sensitivity: 'Sensitivity',
  specificity: 'Specificity',
  precision: 'Precision',
  npv: 'Negative Predicted Value',
  auc: 'AUC-ROC',
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
  f1: 'F1 score of the model. It is the harmonic mean of precision and sensitivity; for more than two classes, ' +
    'the average of the per-class F1 scores. Higher values indicate better overall performance.',
  sensitivity: 'Sensitivity (or recall) of the model. It measures the proportion of actual positives that are ' +
    'correctly identified. Higher values indicate fewer false negatives.',
  specificity: 'Specificity of the model. It measures the proportion of actual negatives that are correctly ' +
    'identified. Higher values indicate fewer false positives.',
  precision: 'Precision of the model. It measures the proportion of true positives among the total number of ' +
    'positive predictions. Higher values indicate fewer false positives.',
  npv: 'Negative Predicted Value of the model. It measures the proportion of true negatives among the total number ' +
    'of negative predictions. Higher values indicate fewer false negatives.',
  auc: 'How well the model separates the classes: 1 is perfect, 0.5 is random guessing. With more than two classes, ' +
    'the average over the classes.',
};

/** A training's two measurements: the final model on its training rows, and the out-of-fold predictions. */
export type Split = 'train' | 'validation';
const SPLITS: Split[] = ['train', 'validation'];

/** The values below which a metric is marked (the built-in tool's thresholds). */
const METRIC_THRESHOLDS: Partial<Record<MetricId, number>> = {r2: 0.3, auc: 0.5};

/** A metric value worth a word: below its threshold, or a validation value far from the train one. */
export interface MetricFinding {
  id: MetricId;
  kind: 'threshold' | 'mismatch';
  /** The splits whose values the finding marks. */
  splits: Split[];
  text: string;
  detail: string;
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

export function classificationMetrics(actual: DG.Column, predicted: DG.Column,
  positiveClass: string): MetricValues {
  const rows = pairedRows(actual, predicted);
  const n = rows.length;
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

  let correct = 0;
  let f1Sum = 0;
  for (const c of counts.values()) {
    correct += c.tp;
    f1Sum += f1Of(c);
  }
  const metrics: MetricValues = {accuracy: correct / n};
  if (counts.size > 2) {
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

/** The metrics of [predicted]; with a classifier's [probabilities], also `auc` ({@link aucOfProbabilities}). */
export function metricsOf(task: ModelTask, actual: DG.Column, predicted: DG.Column, positiveClass?: string,
  probabilities?: ClassProbabilities): MetricValues {
  if (task === 'regression')
    return regressionMetrics(actual, predicted);
  const values = classificationMetrics(actual, predicted, positiveClass ?? actual.categories[0]);
  const auc = probabilities === undefined ? undefined : aucOfProbabilities(actual, probabilities, positiveClass);
  return auc === undefined ? values : {...values, auc};
}

/** A numerical target makes a regression, any other a classification. */
export function taskOf(target: DG.Column): ModelTask {
  return target.matches(DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME) ? 'regression' : 'classification';
}

/** The metrics of the [labels] a Predict probability model cut its scores into, with the scores' [auc] when there is
 * one. */
export function probabilityValues(actual: DG.Column, labels: DG.Column, positiveClass: string,
  auc: number | undefined): MetricValues {
  const values = classificationMetrics(actual, labels, positiveClass);
  return auc === undefined ? values : {...values, auc};
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

/** The cutoffs the ROC curve is drawn at: 0, 0.01, ..., 1. */
export const ROC_THRESHOLDS: readonly number[] = Array.from({length: 101}, (_, i) => i / 100);

/** The points of a ROC curve, one per threshold of {@link ROC_THRESHOLDS}. */
interface RocCurve {
  /** The share of the negative rows scored at or above the threshold. */
  fpr: Float32Array;
  /** The share of the positive rows scored at or above the threshold. */
  tpr: Float32Array;
}

/** The ROC curve of [score] for [positiveClass] of [actual]: at each of {@link ROC_THRESHOLDS}, a row is predicted
 * positive when its score is at least the threshold (as the cutoff cuts the scores); a rate whose class has no row is
 * 0. */
export function rocCurveOf(actual: DG.Column, score: DG.Column, positiveClass: string): RocCurve {
  const paired = pairedRows(actual, score);
  const values = Float64Array.from(paired, (i) => score.getNumber(i));
  const isPositive = paired.map((i) => String(actual.get(i)) === positiveClass);
  const positives = isPositive.filter((p) => p).length;
  const negatives = paired.length - positives;
  const fpr = new Float32Array(ROC_THRESHOLDS.length);
  const tpr = new Float32Array(ROC_THRESHOLDS.length);
  for (const [t, threshold] of ROC_THRESHOLDS.entries()) {
    let tp = 0;
    let fp = 0;
    for (let k = 0; k < values.length; k++) {
      if (values[k] >= threshold) {
        if (isPositive[k])
          tp++;
        else
          fp++;
      }
    }
    tpr[t] = positives > 0 ? tp / positives : 0;
    fpr[t] = negatives > 0 ? fp / negatives : 0;
  }
  return {fpr, tpr};
}

/** The built-in tool's closeness of two values: within 0.1 when both are small (their mean magnitude below 1),
 * otherwise within 90% of that mean, so that only values several times apart differ. */
export function areValuesSimilar(a: number, b: number): boolean {
  const mean = (Math.abs(a) + Math.abs(b)) / 2;
  return mean < 1 ? Math.abs(a - b) <= 0.1 : Math.abs(a - b) / mean <= 0.9;
}

/** For each metric: the splits whose value is below its threshold, then a validation value far from the train one. */
export function metricFindings(metrics: Record<Split, MetricValues>): MetricFinding[] {
  const findings: MetricFinding[] = [];
  for (const id of METRIC_IDS) {
    const train = metrics.train[id];
    const validation = metrics.validation[id];
    const label = METRIC_LABELS[id];
    const detail = `Train ${fixed(train)}, validation ${fixed(validation)}.`;
    const threshold = METRIC_THRESHOLDS[id];
    const below = threshold === undefined ? [] :
      SPLITS.filter((split) => (metrics[split][id] ?? threshold) < threshold);
    if (below.length > 0)
      findings.push({id, kind: 'threshold', splits: below, text: `${label} is below ${threshold}.`, detail});
    if (train !== undefined && validation !== undefined && !areValuesSimilar(train, validation)) {
      findings.push({id, kind: 'mismatch', splits: ['validation'], detail,
        text: `The validation ${label} is far from the training one: the model may be overfitting.`});
    }
  }
  return findings;
}

/** [value] with 3 decimals, `-` when there is none. */
function fixed(value: number | undefined): string {
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
