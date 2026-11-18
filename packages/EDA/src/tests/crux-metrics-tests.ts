import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectArray, expectExceptionAsync, expectFloat,
  test} from '@datagrok-libraries/test/src/test';

import {ClassificationReport, classificationReport} from '../crux/crux-metrics';
import {cruxFit, cruxPredict} from '../crux/crux-model';
import {IRIS_SPECIES, accuracy, copyCol, mulberry32, numericFeatures} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';
const FIRST_TIMEOUT = 120000;
const TIMEOUT = 60000;
const BENCH_TIMEOUT = 600000;
const TOLERANCE = 1e-9;

const labels = (name: string, v: (string | null)[]) => DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, v);
const floats = (name: string, v: (number | null)[]) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, v);

function none(value: unknown, what: string): void {
  expect(value === undefined, true, `${what} must be undefined, got ${JSON.stringify(value)}`);
}

function near(actual: number | undefined, expected: number, what: string): void {
  expect(actual !== undefined, true, `${what} is undefined`);
  expectFloat(actual!, expected, TOLERANCE, what);
}

/** `tolerance` 1e-6 for the thresholds: they are the FLOAT (float32) column's values */
function nearAll(actual: (number | null)[] | undefined, expected: (number | null)[], what: string,
  tolerance = TOLERANCE): void {
  expect(actual?.length, expected.length, `${what} length`);
  for (const [i, e] of expected.entries()) {
    if (e === null || !Number.isFinite(e))
      expect(actual![i], e, `${what}[${i}]`);
    else
      expectFloat(actual![i]!, e, tolerance, `${what}[${i}]`);
  }
}

/** The trapezoid AUC-ROC of `scores` for the rows where `positive` is true, the rows of one score as one step */
function aucJs(positive: boolean[], scores: number[]): number {
  const order = scores.map((_, i) => i).sort((a, b) => scores[b] - scores[a]);
  const positives = positive.filter((p) => p).length;
  const negatives = positive.length - positives;
  let auc = 0;
  let tp = 0;
  let fp = 0;
  for (let k = 0; k < order.length;) {
    const [tpr, fpr] = [tp / positives, fp / negatives];
    for (const tied = scores[order[k]]; k < order.length && scores[order[k]] === tied; k++) {
      if (positive[order[k]])
        tp++;
      else
        fp++;
    }
    auc += (fp / negatives - fpr) * (tp / positives + tpr) / 2;
  }
  return auc;
}

function twoClasses(): {actual: DG.Column, predicted: DG.Column, probaA: DG.Column, probaB: DG.Column} {
  return {
    actual: labels('actual', ['a', 'a', 'b', 'b']),
    predicted: labels('predicted', ['a', 'b', 'b', 'b']),
    probaA: floats('proba_a', [0.9, 0.3, 0.8, 0.1]),
    probaB: floats('proba_b', [0.1, 0.7, 0.2, 0.9]),
  };
}

function expectTwoClassLabels(report: ClassificationReport): void {
  expect(report.n, 4);
  expectArray(report.labels, ['a', 'b']);
  near(report.accuracy, 0.75, 'accuracy');
  near(report.balancedAccuracy, 0.75, 'balancedAccuracy');
  near(report.mcc, 2 / Math.sqrt(12), 'mcc');
  near(report.precision, 1, 'precision');
  near(report.recall, 0.5, 'recall');
  near(report.f1, 2 / 3, 'f1');
  expectArray(report.confusion, [[1, 1], [0, 2]]);
  expect(report.classes.length, 2);
  const [a, b] = report.classes;
  expect(a.label, 'a');
  expect(a.rows, 2);
  near(a.precision, 1, 'a precision');
  near(a.recall, 0.5, 'a recall');
  near(a.f1, 2 / 3, 'a f1');
  expect(b.label, 'b');
  expect(b.rows, 2);
  near(b.precision, 2 / 3, 'b precision');
  near(b.recall, 1, 'b recall');
  near(b.f1, 0.8, 'b f1');
  expect(typeof report.cohenKappa, 'number');
}

async function irisModel(): Promise<{features: DG.Column[], target: DG.Column, bytes: Uint8Array}> {
  const df = await grok.dapi.files.readCsv(IRIS);
  const features = numericFeatures(df, 'Species', ['col 1']).toList();
  const target = df.col('Species')!;
  const bytes = await cruxFit('classification', 'logistic', DG.DataFrame.fromColumns(features.map(copyCol)), target,
    {preprocess: 'standard'});
  return {features, target, bytes};
}

category('Crux metrics', () => {
  test('Two classes: scores, the confusion matrix and the per-class values', async () => {
    const {actual, predicted} = twoClasses();
    const report = await classificationReport(actual, predicted, 'a');
    expectTwoClassLabels(report);
    for (const [name, value] of Object.entries({auc: report.auc, averagePrecision: report.averagePrecision,
      roc: report.roc, pr: report.pr}))
      none(value, name);
  }, {timeout: FIRST_TIMEOUT});

  test('Probabilities: AUC-ROC, average precision and the curves of the positive class', async () => {
    const {actual, predicted, probaA, probaB} = twoClasses();
    const probabilities = {classes: ['a', 'b'], columns: [probaA, probaB]};
    const report = await classificationReport(actual, predicted, 'a', probabilities, true);
    near(report.auc, 0.75, 'auc');
    near(report.classes[0].auc, 0.75, 'a auc');
    near(report.averagePrecision, 5 / 6, 'averagePrecision');
    near(report.classes[1].auc, 0.75, 'b auc');
    nearAll(report.roc?.fpr, [0, 0, 0.5, 0.5, 1], 'roc.fpr');
    nearAll(report.roc?.tpr, [0, 0.5, 0.5, 1, 1], 'roc.tpr');
    expect(report.roc!.threshold[0], Infinity);
    nearAll(report.roc!.threshold.slice(1), [0.9, 0.8, 0.3, 0.1], 'roc.threshold', 1e-6);
    nearAll(report.pr?.recall, [1, 1, 0.5, 0.5, 0], 'pr.recall');
    nearAll(report.pr?.precision, [0.5, 2 / 3, 0.5, 1, 1], 'pr.precision');
    nearAll(report.pr?.threshold, [0.1, 0.3, 0.8, 0.9, null], 'pr.threshold', 1e-6);

    const plain = await classificationReport(actual, predicted, 'a', probabilities, false);
    none(plain.roc, 'roc');
    none(plain.pr, 'pr');
  }, {timeout: TIMEOUT});

  test('Three classes: macro scores and the one-vs-rest averages', async () => {
    const actual = labels('actual', ['a', 'a', 'b', 'b', 'c']);
    const predicted = labels('predicted', ['a', 'b', 'b', 'b', 'c']);
    const values = [[0.7, 0.2, 0.1], [0.4, 0.5, 0.1], [0.2, 0.6, 0.2], [0.1, 0.8, 0.1], [0.1, 0.1, 0.8]];
    const classes = ['a', 'b', 'c'];
    const columns = classes.map((c, k) => floats(`proba_${c}`, values.map((row) => row[k])));
    const report = await classificationReport(actual, predicted, undefined, {classes, columns}, true);

    near(report.accuracy, 0.8, 'accuracy');
    near(report.f1, (2 / 3 + 0.8 + 1) / 3, 'f1');
    const meanOf = (xs: number[]) => xs.reduce((s, x) => s + x, 0) / xs.length;
    near(report.precision, meanOf(report.classes.map((c) => c.precision!)), 'precision');
    near(report.recall, meanOf(report.classes.map((c) => c.recall!)), 'recall');
    near(report.balancedAccuracy, meanOf(report.classes.map((c) => c.recall!)), 'balancedAccuracy');
    const actualValues = actual.toList();
    for (const [k, c] of classes.entries()) {
      near(report.classes[k].auc, aucJs(actualValues.map((v) => v === c), values.map((row) => row[k])),
        `${c} auc`);
    }
    near(report.auc, meanOf(report.classes.map((c) => c.auc!)), 'auc');
    none(report.roc, 'roc');
    none(report.pr, 'pr');
  }, {timeout: TIMEOUT});

  test('A class with no positive rows and an unknown class get no AUC', async () => {
    const {actual, predicted, probaA, probaB} = twoClasses();
    const probaC = floats('proba_c', [0, 0, 0, 0]);
    const report = await classificationReport(actual, predicted, undefined,
      {classes: ['a', 'b', 'c'], columns: [probaA, probaB, probaC]});
    expect(report.classes.length, 2);
    near(report.auc, (report.classes[0].auc! + report.classes[1].auc!) / 2, 'auc');
  }, {timeout: TIMEOUT});

  test('Boolean columns', async () => {
    const actual = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'actual', [true, false, true, false]);
    const predicted = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'predicted', [true, true, true, false]);
    const report = await classificationReport(actual, predicted, 'true');
    expectArray(report.labels, ['false', 'true']);
    near(report.precision, 2 / 3, 'precision');
    near(report.recall, 1, 'recall');

    const text = labels('predicted', ['true', 'true', 'true', 'false']);
    const mixed = await classificationReport(actual, text, 'true');
    expectArray(mixed.labels, report.labels);
    expectArray(mixed.confusion, report.confusion);
    near(mixed.mcc, report.mcc, 'mcc');
  }, {timeout: TIMEOUT});

  test('Rows without a value are left out', async () => {
    const actual = labels('actual', ['a', null, 'b', 'b']);
    const predicted = labels('predicted', ['a', 'b', null, 'b']);
    const probaA = floats('proba_a', [0.9, 0.2, 0.95, 0.1]);
    const probaB = floats('proba_b', [null, 0.5, 0.5, NaN]);
    const report = await classificationReport(actual, predicted, 'a', {classes: ['a', 'b'], columns: [probaA, probaB]});
    expect(report.n, 2);
    expectArray(report.confusion, [[1, 0], [0, 1]]);
    near(report.auc, 1, 'auc');
    none(report.classes[1].auc, 'b auc');

    const empty = labels('empty', [null, null, null, null]);
    await expectExceptionAsync(async () => {
      await classificationReport(empty, predicted);
    }, (e) => String(e?.message ?? e).includes('no row has both'));
  }, {timeout: TIMEOUT});

  test('A positive class absent from the rows', async () => {
    const actual = labels('actual', ['b', 'b', 'b', 'c']);
    const predicted = labels('predicted', ['b', 'b', 'c', 'c']);
    const probabilities = {classes: ['b', 'c'], columns: [floats('proba_b', [0.9, 0.8, 0.4, 0.1]),
      floats('proba_c', [0.1, 0.2, 0.6, 0.9])]};
    const report = await classificationReport(actual, predicted, 'a', probabilities, true);
    near(report.precision, 0.75, 'precision');
    none(report.roc, 'roc');
    none(report.pr, 'pr');
  }, {timeout: TIMEOUT});

  test('The package function', async () => {
    const {actual, predicted, probaA, probaB} = twoClasses();
    const probabilities = DG.DataFrame.fromColumns([probaA, probaB]);
    const full: ClassificationReport = await grok.functions.call('Eda:classificationReport', {actual, predicted,
      positiveClass: 'a', classes: ['a', 'b'], probabilities, curves: true});
    near(full.accuracy, 0.75, 'accuracy');
    expect(full.roc?.fpr.length, 5);

    const plain: ClassificationReport = await grok.functions.call('Eda:classificationReport', {actual, predicted,
      positiveClass: '', classes: [], probabilities: null, curves: false});
    const direct = await classificationReport(actual, predicted);
    for (const key of ['accuracy', 'balancedAccuracy', 'mcc', 'precision', 'recall', 'f1'] as const)
      near(plain[key], direct[key]!, key);
    expectArray(plain.confusion, direct.confusion);
    none(plain.auc, 'auc');
    none(plain.roc, 'roc');
  }, {timeout: TIMEOUT});

  test('Parity with cruxPredict probabilities on iris', async () => {
    const {features, target, bytes} = await irisModel();
    const out = await cruxPredict(DG.DataFrame.fromColumns(features.map(copyCol)), bytes);
    const prediction = out.col('Prediction')!;
    const columns = IRIS_SPECIES.map((s) => out.col(`proba_${s}`)!);
    const report = await classificationReport(target, prediction, undefined, {classes: IRIS_SPECIES, columns});
    near(report.accuracy, accuracy(target, prediction), 'accuracy');
    expect(report.classes.length, 3);
    for (const c of report.classes)
      expect(c.auc !== undefined && c.auc > 0.5 && c.auc <= 1, true, `${c.label} auc ${c.auc}`);
    expectArray(report.confusion.map((row) => row.reduce((s, v) => s + v, 0)), [50, 50, 50]);
  }, {timeout: TIMEOUT});

  test('BENCH 100k', async () => {
    const n = 100000;
    const random = mulberry32(7);
    const actualValues = Array.from({length: n}, () => random() < 0.4 ? 'F' : 'M');
    const score = new Float32Array(n);
    for (let i = 0; i < n; ++i)
      score[i] = Math.min(1, Math.max(0, (actualValues[i] === 'F' ? 0.65 : 0.35) + (random() - 0.5) * 0.8));
    const actual = DG.Column.fromStrings('actual', actualValues);
    const predicted = DG.Column.fromStrings('predicted', Array.from(score, (s) => s >= 0.5 ? 'F' : 'M'));
    const probaF = DG.Column.fromFloat32Array('proba_F', score);
    const probaM = DG.Column.fromFloat32Array('proba_M', score.map((s) => 1 - s));
    const start = performance.now();
    const report = await classificationReport(actual, predicted, 'F', {classes: ['F', 'M'], columns: [probaF, probaM]},
      true);
    const ms = performance.now() - start;
    console.log(`BENCH classificationReport 100k: ${ms.toFixed(0)} ms`);
    expect(report.roc !== undefined && report.pr !== undefined, true, 'the curves');
    if (DG.Test.isInBenchmark)
      expect(ms < 3000, true, `${ms.toFixed(0)} ms`);
  }, {timeout: BENCH_TIMEOUT, benchmark: true});
});
