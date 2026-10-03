import * as DG from 'datagrok-api/dg';
import {ForgeError} from '../forge-error';
import {ModelTask} from '../generated/db';

export const METRIC_IDS =
  ['mse', 'rmse', 'mae', 'r2', 'accuracy', 'f1', 'sensitivity', 'specificity', 'precision', 'npv'] as const;
export type MetricId = typeof METRIC_IDS[number];
export type MetricValues = Partial<Record<MetricId, number>>;

export const METRIC_LABELS: Record<MetricId, string> = {
  mse: 'MSE',
  rmse: 'RMSE',
  mae: 'MAE',
  r2: 'R squared',
  accuracy: 'Accuracy',
  f1: 'F1',
  sensitivity: 'Sensitivity',
  specificity: 'Specificity',
  precision: 'Precision',
  npv: 'Negative Predicted Value',
};

export const METRIC_DESCRIPTIONS: Record<MetricId, string> = {
  mse: 'Mean Squared Error of values in the column. It is the average of the squared differences between ' +
    'predicted and actual values. Lower values indicate better model accuracy.',
  rmse: 'Root Mean Squared Error of values in the column. It represents the square root of the average squared ' +
    'differences between predicted and actual values. Lower values indicate better model performance.',
  mae: 'Mean Absolute Error of values in the column. It is the average of the absolute differences between ' +
    'predicted and actual values. Lower values indicate better model accuracy.',
  r2: 'R-squared value of the column. It represents the proportion of the variance for a dependent variable that ' +
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
};

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

export function metricsOf(task: ModelTask, actual: DG.Column, predicted: DG.Column,
  positiveClass?: string): MetricValues {
  return task === 'regression' ? regressionMetrics(actual, predicted) :
    classificationMetrics(actual, predicted, positiveClass ?? actual.categories[0]);
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
