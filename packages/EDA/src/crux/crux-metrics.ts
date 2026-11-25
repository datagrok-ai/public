import * as DG from 'datagrok-api/dg';

import {makeVector, Table} from 'apache-arrow';

import {runCrux} from './crux-model';
import {columnTexts, nullSentinel, numbers, scalar, stringValues, utf8Vector, vectorsIpc} from './crux-bridge';

export interface ReportProbabilities {
  classes: string[];
  columns: DG.Column[];
}

export interface ClassReport {
  label: string;
  rows: number;
  precision?: number;
  recall?: number;
  f1: number;
  auc?: number;
  averagePrecision?: number;
}

export interface ClassificationReport {
  n: number;
  labels: string[];
  accuracy: number;
  balancedAccuracy: number;
  mcc: number;
  cohenKappa?: number;
  precision?: number;
  recall?: number;
  f1: number;
  auc?: number;
  averagePrecision?: number;
  classes: ClassReport[];
  confusion: number[][];
  roc?: {fpr: number[]; tpr: number[]; threshold: number[]};
  pr?: {precision: number[]; recall: number[]; threshold: (number | null)[]};
}

function texts(col: DG.Column): (string | null)[] {
  switch (col.type) {
  case DG.COLUMN_TYPE.STRING:
    return columnTexts(col);
  case DG.COLUMN_TYPE.BOOL:
    return col.toList().map((v) => v ? 'true' : 'false');
  case DG.COLUMN_TYPE.INT:
  case DG.COLUMN_TYPE.FLOAT: {
    const raw = col.getRawData();
    const none = nullSentinel(col);
    return Array.from({length: col.length}, (_, i) => raw[i] === none ? null : String(raw[i]));
  }
  default:
    throw new Error(`Crux metrics: "${col.name}" is ${col.type}; a text, boolean or numeric column expected`);
  }
}

function required(table: Table, name: string): number {
  const v = scalar(table, name);
  if (v === undefined)
    throw new Error(`Crux metrics: ml.classification_metrics returned no ${name}`);
  return v;
}

function mean(values: (number | undefined)[]): number | undefined {
  const given = values.filter((v): v is number => v !== undefined);
  return given.length === 0 ? undefined : given.reduce((s, v) => s + v, 0) / given.length;
}

/** The Crux (scikit-learn) classification scores of `predicted` against `actual` over the rows where both have a
 * value (classes compared as text: a boolean is `true` / `false`): accuracy, balanced accuracy, MCC, Cohen's kappa,
 * precision / recall / F1 of `positiveClass` when the rows hold exactly two classes and it is one of them, else macro;
 * the confusion matrix (`labels` sorted, true-major) and each class's rows, precision, recall and F1 from it; with
 * `probabilities` (one column per class) each class's one-vs-rest AUC-ROC and average precision over the rows with a
 * value, the overall ones (the positive class's, else the macro average over the classes with both positive and
 * negative rows) and, when `curves`, the positive class's ROC and precision-recall curves (one point per distinct
 * probability; the ROC's first threshold is +Infinity, the PR's last is null) */
export async function classificationReport(actual: DG.Column, predicted: DG.Column, positiveClass?: string,
  probabilities?: ReportProbabilities, curves = false): Promise<ClassificationReport> {
  const actualTexts = texts(actual);
  const predictedTexts = texts(predicted);
  const rows: number[] = [];
  for (let i = 0; i < actualTexts.length; ++i) {
    if (actualTexts[i] !== null && predictedTexts[i] !== null)
      rows.push(i);
  }
  if (rows.length === 0)
    throw new Error('Crux metrics: no row has both an actual and a predicted class');
  const yTrueTexts = rows.map((i) => actualTexts[i]);
  const yPredTexts = rows.map((i) => predictedTexts[i]);

  const matrix = await runCrux('ml.confusion_matrix',
    vectorsIpc({y_true: utf8Vector(yTrueTexts), y_pred: utf8Vector(yPredTexts)}), {});
  const k = Math.round(Math.sqrt(matrix.numRows));
  const trueLabels = stringValues(matrix, 'true_label');
  const counts = numbers(matrix, 'count');
  const labels = Array.from({length: k}, (_, i) => trueLabels[i * k]!);
  const confusion = labels.map((_, i) => labels.map((_, j) => counts[i * k + j] ?? 0));

  const index = new Map(labels.map((l, i) => [l, i]));
  const yTrue = Int32Array.from(yTrueTexts, (t) => index.get(t!)!);
  const yPred = Int32Array.from(yPredTexts, (t) => index.get(t!)!);
  const pos = positiveClass === undefined ? -1 : labels.indexOf(positiveClass);
  const binary = labels.length === 2 && pos >= 0;
  const scores = await runCrux('ml.classification_metrics',
    vectorsIpc({y_true: makeVector(yTrue), y_pred: makeVector(yPred)}),
    binary ? {average: 'binary', pos_label: pos} : {average: 'macro'});

  const classes: ClassReport[] = labels.map((label, i) => {
    let rowSum = 0;
    let colSum = 0;
    for (let j = 0; j < k; ++j) {
      rowSum += confusion[i][j];
      colSum += confusion[j][i];
    }
    const precision = colSum > 0 ? confusion[i][i] / colSum : undefined;
    const recall = rowSum > 0 ? confusion[i][i] / rowSum : undefined;
    const [p, r] = [precision ?? 0, recall ?? 0];
    return {label, rows: rowSum, precision, recall, f1: p + r > 0 ? 2 * p * r / (p + r) : 0};
  });

  const report: ClassificationReport = {
    n: rows.length,
    labels,
    accuracy: required(scores, 'accuracy'),
    balancedAccuracy: required(scores, 'balanced_accuracy'),
    mcc: required(scores, 'mcc'),
    cohenKappa: scalar(scores, 'cohen_kappa'),
    precision: scalar(scores, 'precision'),
    recall: scalar(scores, 'recall'),
    f1: scalar(scores, 'f1') ?? 0,
    classes,
    confusion,
  };
  if (probabilities === undefined)
    return report;

  let positive: {y01: Int32Array, values: Float64Array} | undefined;
  for (const [c, label] of probabilities.classes.entries()) {
    const i = index.get(label);
    if (i === undefined)
      continue;
    const column = probabilities.columns[c];
    const raw = column.getRawData();
    const none = nullSentinel(column);
    const scored = rows.filter((r) => raw[r] !== none && Number.isFinite(raw[r]));
    const y01 = Int32Array.from(scored, (r) => actualTexts[r] === label ? 1 : 0);
    const values = Float64Array.from(scored, (r) => raw[r]);
    const positives = y01.reduce((s, y) => s + y, 0);
    if (positives === 0 || positives === y01.length)
      continue;
    const ranking = await runCrux('ml.roc_auc', vectorsIpc({y_true: makeVector(y01), y_score: makeVector(values)}),
      {});
    classes[i].auc = scalar(ranking, 'roc_auc');
    classes[i].averagePrecision = scalar(ranking, 'average_precision');
    if (i === pos)
      positive = {y01, values};
  }
  report.auc = binary ? classes[pos].auc : mean(classes.map((c) => c.auc));
  report.averagePrecision = binary ? classes[pos].averagePrecision : mean(classes.map((c) => c.averagePrecision));

  if (curves && binary && positive !== undefined) {
    const ipc = vectorsIpc({y_true: makeVector(positive.y01), y_score: makeVector(positive.values)});
    // runCrux transfers the buffer to the worker: the first call gets a copy.
    const roc = await runCrux('ml.roc_curve', ipc.slice(), {drop_intermediate: false});
    report.roc = {fpr: numbers(roc, 'fpr') as number[], tpr: numbers(roc, 'tpr') as number[],
      threshold: numbers(roc, 'threshold') as number[]};
    const pr = await runCrux('ml.precision_recall_curve', ipc, {});
    report.pr = {precision: numbers(pr, 'precision') as number[], recall: numbers(pr, 'recall') as number[],
      threshold: numbers(pr, 'threshold')};
  }
  return report;
}
