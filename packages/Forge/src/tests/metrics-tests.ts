import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {ForgeError} from '../forge-error';
import {areValuesSimilar, aucOf, aucOfProbabilities, CLASS_METRIC_IDS, ClassProbabilities, classificationMetrics,
  METRIC_DESCRIPTIONS, METRIC_IDS, METRIC_LABELS, metricFindings, MetricId, MetricValues, probabilityOf,
  regressionMetrics} from '../metrics/metrics';
import {ClassificationReport, classificationReport, curvesOf} from '../metrics/report';
import {metricsRecordOf} from '../training/train-model';
import {residualsOf, STATISTIC_IDS, statisticsOf} from '../metrics/statistics';

const numbers = (values: (number | null)[]) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'n', values);
const labels = (values: string[]) => DG.Column.fromList(DG.COLUMN_TYPE.STRING, 's', values);

function expectMetric(metrics: MetricValues, id: MetricId, expected: number): void {
  expectFloat(metrics[id] ?? NaN, expected, 1e-3, id);
}

function expectNoMetric(metrics: MetricValues, id: MetricId): void {
  expect(metrics[id] === undefined, true, `Unexpected ${id} ${metrics[id]}`);
}

category('Metrics', () => {
  test('regression metrics match hand-computed values', async () => {
    const metrics = regressionMetrics(numbers([1, 2, 3, 4]), numbers([1.5, 2, 2.5, 5]));
    expectMetric(metrics, 'mse', 0.375);
    expectMetric(metrics, 'rmse', 0.6124);
    expectMetric(metrics, 'mae', 0.5);
    expectMetric(metrics, 'r2', 0.7);
  });

  test('r2 on a constant target follows the scikit-learn convention', async () => {
    expect(regressionMetrics(numbers([2, 2, 2]), numbers([2, 2, 2])).r2, 1);
    expect(regressionMetrics(numbers([2, 2, 2]), numbers([1, 2, 3])).r2, 0);
  });

  test('missing pairs are skipped', async () => {
    const metrics = regressionMetrics(numbers([1, null, 3, 4]), numbers([1.5, 2, 2.5, 5]));
    expectMetric(metrics, 'mse', (0.25 + 0.25 + 1) / 3);
  });

  test('binary metrics match the confusion-matrix definitions', async () => {
    const metrics = classificationMetrics(labels(['a', 'a', 'b', 'b']), labels(['a', 'b', 'b', 'b']), 'a');
    expectMetric(metrics, 'accuracy', 0.75);
    expectMetric(metrics, 'sensitivity', 0.5);
    expectMetric(metrics, 'specificity', 1);
    expectMetric(metrics, 'precision', 1);
    expectMetric(metrics, 'npv', 0.6667);
    expectMetric(metrics, 'f1', 0.6667);
  });

  test('a share with a zero denominator is omitted', async () => {
    const metrics = classificationMetrics(labels(['a', 'a']), labels(['a', 'a']), 'a');
    expectMetric(metrics, 'sensitivity', 1);
    expectMetric(metrics, 'precision', 1);
    expectMetric(metrics, 'accuracy', 1);
    expectMetric(metrics, 'f1', 1);
    expectNoMetric(metrics, 'specificity');
    expectNoMetric(metrics, 'npv');
  });

  test('multiclass metrics match hand-computed values', async () => {
    const metrics = classificationMetrics(labels(['a', 'a', 'b', 'b', 'c']), labels(['a', 'b', 'b', 'b', 'c']), 'a');
    expectMetric(metrics, 'accuracy', 0.8);
    expectMetric(metrics, 'f1', (2 / 3 + 0.8 + 1) / 3);
    for (const id of ['sensitivity', 'specificity', 'precision', 'npv'] as const)
      expectNoMetric(metrics, id);
  });

  test('balanced accuracy and MCC match the hand-computed values and the Crux report', async () => {
    const metrics = classificationMetrics(labels(['a', 'a', 'b', 'b', 'c']), labels(['a', 'b', 'b', 'b', 'c']), 'a');
    expectMetric(metrics, 'balancedAccuracy', 2.5 / 3);
    expectMetric(metrics, 'mcc', 11 / Math.sqrt(16 * 14));
    const cases: [string[], string[]][] = [[['a', 'a', 'b', 'b', 'a', 'b'], ['a', 'b', 'b', 'a', 'a', 'b']],
      [['a', 'a', 'b', 'b', 'c', 'c', 'a'], ['a', 'b', 'b', 'c', 'c', 'a', 'a']],
      [['a', 'a', 'b', 'b'], ['a', 'c', 'b', 'b']], [['a', 'a', 'b'], ['a', 'a', 'a']]];
    for (const [actual, predicted] of cases) {
      const ts = classificationMetrics(labels(actual), labels(predicted), 'a');
      const report = (await classificationReport(labels(actual), labels(predicted)))!;
      expectFloat(ts.balancedAccuracy ?? NaN, report.balancedAccuracy, 1e-9, `balanced accuracy ${predicted}`);
      expectFloat(ts.mcc ?? NaN, report.mcc, 1e-9, `mcc ${predicted}`);
    }
  });

  test('bool and string labels compare as text', async () => {
    const actual = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'b', [true, false, true]);
    const metrics = classificationMetrics(actual, labels(['true', 'true', 'true']), 'true');
    expectMetric(metrics, 'accuracy', 0.6667);
    expectMetric(metrics, 'sensitivity', 1);
    expectMetric(metrics, 'specificity', 0);
  });

  test('empty overlap throws ForgeError', async () => {
    const isForgeError = (e: unknown) => e instanceof ForgeError;
    await expectExceptionAsync(async () => {
      regressionMetrics(numbers([null, null]), numbers([1, 2]));
    }, isForgeError);
    await expectExceptionAsync(async () => {
      classificationMetrics(labels(['a', 'b']), DG.Column.fromList(DG.COLUMN_TYPE.STRING, 's', [null, null]), 'a');
    }, isForgeError);
  });

  test('areValuesSimilar: an absolute tolerance for small values, a relative one otherwise', async () => {
    expect(areValuesSimilar(0.5, 0.55), true);
    expect(areValuesSimilar(0.5, 0.7), false);
    expect(areValuesSimilar(10, 20), true);
    expect(areValuesSimilar(10, 30), false);
  });

  test('metricFindings: values below the threshold and a validation value far from the train one', async () => {
    const low = metricFindings({train: {r2: 0.2, mse: 1}, validation: {r2: 0.2, mse: 1.1}});
    expect(low.length, 1, JSON.stringify(low));
    expect(low[0].kind, 'threshold');
    expectArray(low[0].splits, ['train', 'validation']);
    expect(low[0].text, 'R2 is below 0.3.');
    expect(low[0].detail, 'Train 0.200, validation 0.200.');

    const far = metricFindings({train: {accuracy: 0.9}, validation: {accuracy: 0.4}});
    expect(far.length, 1, JSON.stringify(far));
    expect(far[0].kind, 'mismatch');
    expectArray(far[0].splits, ['validation']);
    expect(far[0].text, 'The validation Accuracy is far from the training one: the model may be overfitting.');
    expect(far[0].detail, 'Train 0.900, validation 0.400.');

    expect(metricFindings({train: {accuracy: 0.97, f1: 0.97, auc: 0.99}, validation: {accuracy: 0.95, f1: 0.95,
      auc: 0.98}}).length, 0);
    const oneSplit = metricFindings({train: {auc: 0.6}, validation: {auc: 0.55}});
    expect(oneSplit.length, 0, JSON.stringify(oneSplit));
    const belowOnValidation = metricFindings({train: {auc: 0.52}, validation: {auc: 0.45}});
    expectArray(belowOnValidation.map((f) => `${f.kind}: ${f.splits.join()}`), ['threshold: validation']);
  });

  test('METRIC_IDS ends with the Crux ids and every id has a label and a description', async () => {
    expectArray(METRIC_IDS.slice(-3), ['balancedAccuracy', 'mcc', 'averagePrecision']);
    expectArray(CLASS_METRIC_IDS, ['precision', 'sensitivity', 'f1', 'auc', 'averagePrecision']);
    for (const id of METRIC_IDS)
      expect(METRIC_LABELS[id].length > 0 && METRIC_DESCRIPTIONS[id].length > 0, true, `${id} has no text`);
  });

  test('metricsRecordOf keeps the classes and drops a malformed one', async () => {
    const record = metricsRecordOf({train: {accuracy: 1}, validation: {accuracy: 0.5},
      classes: [{label: 'a', rows: 2, train: {precision: 1}, validation: {precision: 0.5}}, {label: 3}]});
    expect(record?.classes?.length, 1);
    const [a] = record!.classes!;
    expect(a.label, 'a');
    expect(a.rows, 2);
    expect(a.train.precision, 1);
    expect(a.validation.precision, 0.5);
    expect(metricsRecordOf({train: {}, validation: {}})?.classes === undefined, true, 'classes without any');
  });

  test('curvesOf turns the PR curve\'s null threshold into Infinity and gives none without curves', async () => {
    const report: ClassificationReport = {n: 4, labels: ['a', 'b'], accuracy: 0.75, balancedAccuracy: 0.75, mcc: 0.5,
      f1: 0.5, classes: [], confusion: [[1, 1], [0, 2]],
      roc: {fpr: [0, 0, 1], tpr: [0, 1, 1], threshold: [Infinity, 0.9, 0.1]},
      pr: {precision: [0.5, 1, 1], recall: [1, 0.5, 0], threshold: [0.1, 0.9, null]}};
    const curves = curvesOf(report)!;
    expect(curves.roc.threshold[0], Infinity);
    expectArray(Array.from(curves.roc.tpr), [0, 1, 1]);
    expect(curves.pr.threshold[2], Infinity);
    expectFloat(curves.pr.threshold[1], 0.9, 1e-6);
    expect(curvesOf({...report, roc: undefined, pr: undefined}) === undefined, true, 'curves without any');
    expect(curvesOf(undefined) === undefined, true, 'curves without a report');
  });

  test('classificationReport gives no report when EDA\'s function fails', async () => {
    const actual = labels(['a', 'b', 'a', 'b']);
    const predicted = labels(['a', 'b', 'b', 'b']);
    expect((await classificationReport(actual, predicted, 'a'))?.accuracy, 0.75);
    const broken = {classes: ['a', 'b'], columns: [numbers([0.9, 0.1, 0.4, 0.2])]};
    expect(await classificationReport(actual, predicted, 'a', broken) === undefined, true,
      'A failed call gave a report');
  });

  test('metricFindings: a validation MCC far from the train one is a mismatch', async () => {
    const findings = metricFindings({train: {mcc: 0.9}, validation: {mcc: 0.2}});
    expect(findings.length, 1, JSON.stringify(findings));
    expect(findings[0].kind, 'mismatch');
    expect(findings[0].text, 'The validation MCC is far from the training one: the model may be overfitting.');
  });

  test('aucOfProbabilities: the positive class for two classes, the macro one-vs-rest average for more', async () => {
    const two = labels(['a', 'a', 'b', 'b']);
    // a's positives score 0.9 and 0.3, its negatives 0.8 and 0.1: 3 of 4 pairs ordered.
    const probaA = numbers([0.9, 0.3, 0.8, 0.1]);
    const flat = numbers([0.5, 0.5, 0.5, 0.5]);
    const binary: ClassProbabilities = {classes: ['a', 'b'], columns: [probaA, flat]};
    expectFloat(aucOfProbabilities(two, binary, 'a')!, 0.75, 1e-9);
    expectFloat(aucOfProbabilities(two, binary, 'b')!, 0.5, 1e-9);
    expectFloat(aucOfProbabilities(two, binary)!, (0.75 + 0.5) / 2, 1e-9);
    expect(probabilityOf(binary, 'z') === undefined, true, 'an unknown class has a column');
    expect(aucOfProbabilities(two, binary, 'z') === undefined, true, 'an unknown positive class has an AUC');

    const three = labels(['a', 'a', 'b', 'c']);
    const columns = [numbers([0.7, 0.2, 0.4, 0.1]), numbers([0.2, 0.5, 0.5, 0.3]), numbers([0.1, 0.3, 0.1, 0.6])];
    const macro = aucOfProbabilities(three, {classes: ['a', 'b', 'c'], columns})!;
    const each = ['a', 'b', 'c'].map((cls, k) => aucOf(three, columns[k], cls)!);
    expectFloat(macro, (each[0] + each[1] + each[2]) / 3, 1e-9);

    const absent: ClassProbabilities = {classes: ['a', 'b', 'c'], columns: [probaA, flat, numbers([0, 0, 0, 0])]};
    expectFloat(aucOfProbabilities(two, absent)!, (0.75 + 0.5) / 2, 1e-9);
  });

  test('statisticsOf lists the eleven statistics, q2 being the median', async () => {
    const rows = statisticsOf(numbers([1, 2, 3, 4, 10]), numbers([1, 1, 2, 2, 3]));
    expectArray(rows.map((r) => r.id), [...STATISTIC_IDS]);
    const valueOf = (id: string, side: 'actual' | 'predicted') => rows.find((r) => r.id === id)?.[side] ?? NaN;
    for (const side of ['actual', 'predicted'] as const)
      expect(valueOf('q2', side), valueOf('med', side));
    expect(valueOf('min', 'actual'), 1);
    expect(valueOf('max', 'actual'), 10);
    expectFloat(valueOf('avg', 'predicted'), 1.8, 1e-6);
  });

  test('residualsOf: the prediction minus the target, empty where either is', async () => {
    const residuals = residualsOf(numbers([1, 2, null, 4]), numbers([1.5, 1, 3, null]));
    expect(residuals.name, 'residuals');
    expect(residuals.type, DG.COLUMN_TYPE.FLOAT);
    expectFloat(residuals.get(0), 0.5, 1e-6);
    expectFloat(residuals.get(1), -1, 1e-6);
    expect(residuals.isNone(2) && residuals.isNone(3), true, 'A residual without both values');
  });
});
