import * as DG from 'datagrok-api/dg';
import {ModelTask} from '../generated/db';
import {onFrame} from '../preparation/shared-frame';
import {CLASS_METRIC_IDS, ClassMetrics, classesOf, ClassProbabilities, ClassValues, Curves, MetricValues, metricsOf}
  from './metrics';

/** Mirrors `ClassReport` of EDA's `src/crux/crux-metrics.ts`: one class against the rest. */
export interface ClassReport {
  label: string;
  rows: number;
  precision?: number;
  recall?: number;
  f1: number;
  auc?: number;
  averagePrecision?: number;
}

/** Mirrors `ClassificationReport` of EDA's `src/crux/crux-metrics.ts`: the object `Eda:classificationReport`
 * returns. */
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

const REPORT_FUNCTION = {package: 'Eda', name: 'classificationReport'};

function reportFunction(): DG.Func | undefined {
  return DG.Func.find(REPORT_FUNCTION)[0];
}

/** The Crux classification report of [predicted] against [actual] (EDA's `classificationReport`), with
 * [probabilities] and, when [curves], the positive class's curves; undefined when the EDA package has no such
 * function or the function fails (Crux not loaded, a Crux error): the metrics go on without it. */
export async function classificationReport(actual: DG.Column, predicted: DG.Column, positiveClass?: string,
  probabilities?: ClassProbabilities, curves = false): Promise<ClassificationReport | undefined> {
  const func = reportFunction();
  if (func === undefined)
    return undefined;
  const call = (frame: DG.DataFrame | null) => func.apply({actual, predicted, positiveClass: positiveClass ?? '',
    classes: probabilities?.classes ?? [], probabilities: frame, curves});
  try {
    const result: unknown = probabilities === undefined ? await call(null) : await onFrame(probabilities.columns, call);
    return result as ClassificationReport;
  } catch {
    return undefined;
  }
}

export interface Measured { values: MetricValues; report?: ClassificationReport }

/** The metrics of the class probabilities (or scores), which no cutoff changes. */
export const SCORE_METRIC_IDS = ['auc', 'averagePrecision'] as const;

/** The metrics of [predicted] ({@link metricsOf}); for a classifier with the report, AUC-ROC and average precision
 * are the report's (with [probabilities]), and the report itself gives the per-class values and the curves. Without
 * the report, AUC-ROC is {@link metricsOf}'s. */
export async function measure(task: ModelTask, actual: DG.Column, predicted: DG.Column, positiveClass?: string,
  probabilities?: ClassProbabilities, curves = false): Promise<Measured> {
  if (task === 'regression')
    return {values: metricsOf(task, actual, predicted)};
  const report = await classificationReport(actual, predicted, positiveClass ?? classesOf(actual)[0], probabilities,
    curves);
  if (report === undefined)
    return {values: metricsOf(task, actual, predicted, positiveClass, probabilities)};
  const values = metricsOf(task, actual, predicted, positiveClass);
  for (const id of SCORE_METRIC_IDS) {
    const v = report[id];
    if (v !== undefined)
      values[id] = v;
  }
  return {values, report};
}

function classMetricValues(c: ClassReport): MetricValues {
  const all: MetricValues = {precision: c.precision, sensitivity: c.recall, f1: c.f1, auc: c.auc,
    averagePrecision: c.averagePrecision};
  const values: MetricValues = {};
  for (const id of CLASS_METRIC_IDS) {
    const v = all[id];
    if (typeof v === 'number')
      values[id] = v;
  }
  return values;
}

/** The per-class values of [report] ({@link CLASS_METRIC_IDS}: the report's precision, recall as sensitivity, f1,
 * auc, averagePrecision; a missing one left out). */
export function classValuesOf(report: ClassificationReport): ClassValues[] {
  return report.classes.map((c) => ({label: c.label, rows: c.rows, values: classMetricValues(c)}));
}

/** The per-class record of a training from both splits' reports, by label (a label of one report only has empty
 * values for the other); undefined without a report. */
export function classMetricsOf(train: ClassificationReport | undefined,
  validation: ClassificationReport | undefined): ClassMetrics[] | undefined {
  if (train === undefined && validation === undefined)
    return undefined;
  const trainClasses = train?.classes ?? [];
  const validationClasses = validation?.classes ?? [];
  const labels = [...new Set([...trainClasses, ...validationClasses].map((c) => c.label))];
  return labels.map((label) => {
    const t = trainClasses.find((c) => c.label === label);
    const v = validationClasses.find((c) => c.label === label);
    return {label, rows: (t ?? v)!.rows, train: t === undefined ? {} : classMetricValues(t),
      validation: v === undefined ? {} : classMetricValues(v)};
  });
}

/** The curves of [report] as float arrays; the PR curve's null threshold becomes Infinity; undefined without
 * curves. */
export function curvesOf(report: ClassificationReport | undefined): Curves | undefined {
  if (report?.roc === undefined || report.pr === undefined)
    return undefined;
  const {roc, pr} = report;
  return {
    roc: {fpr: Float32Array.from(roc.fpr), tpr: Float32Array.from(roc.tpr),
      threshold: Float32Array.from(roc.threshold)},
    pr: {precision: Float32Array.from(pr.precision), recall: Float32Array.from(pr.recall),
      threshold: Float32Array.from(pr.threshold, (t) => t ?? Infinity)},
  };
}
