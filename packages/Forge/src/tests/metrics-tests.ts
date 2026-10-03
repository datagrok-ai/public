import * as DG from 'datagrok-api/dg';
import {category, expect, expectExceptionAsync, expectFloat, test} from '@datagrok-libraries/test/src/test';
import {ForgeError} from '../forge-error';
import {classificationMetrics, MetricId, MetricValues, regressionMetrics} from '../metrics/metrics';

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
});
