import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {ForgeError} from '../forge-error';
import {areValuesSimilar, aucOf, aucOfProbabilities, baselineMetrics, BEST_THRESHOLD_METRICS, bestThreshold,
  CLASS_METRIC_IDS, ClassProbabilities, classificationMetrics, defaultMetricOf, METRIC_DESCRIPTIONS, METRIC_IDS,
  METRIC_LABELS,
  metricFindings, MetricId, MetricValues, overfittingOf, probabilityMetrics, probabilityOf, regressionMetrics,
  underfittingOf, verdictOf} from '../metrics/metrics';
import {ClassificationReport, classificationReport, curvesOf} from '../metrics/report';
import {defaultMetric, directionOf, SEARCH_METRICS} from '../training/search';
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

  test('missing pairs are left out', async () => {
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

  test('metricFindings: a train value barely better than the baseline and a validation value far apart', async () => {
    const low = metricFindings({train: {r2: 0.05, mse: 1}, validation: {r2: 0.05, mse: 1.1}, baseline: {r2: 0}});
    expect(low.length, 1, JSON.stringify(low));
    expect(low[0].kind, 'underfitting');
    expectArray(low[0].splits, ['train']);
    expect(low[0].text, 'R2 is barely better than predicting the mean on the training rows.');
    expect(low[0].detail, 'Train 0.050, predicting the mean 0.000.');
    expect(metricFindings({train: {r2: 0.05}, validation: {r2: 0.05}}).length, 0, 'underfitting without a baseline');

    const classifier = metricFindings({train: {balancedAccuracy: 0.3, sensitivity: 0.3, accuracy: 0.9},
      validation: {balancedAccuracy: 0.3, sensitivity: 0.3, accuracy: 0.9},
      baseline: {balancedAccuracy: 0.25, accuracy: 0.9}});
    expectArray(classifier.map((f) => `${f.id} ${f.kind}`), ['accuracy underfitting', 'balancedAccuracy underfitting']);
    expect(classifier[1].text, 'Balanced accuracy is barely better than guessing on the training rows.');
    expect(classifier[1].detail, 'Train 0.300, guessing 0.250.');

    const far = metricFindings({train: {accuracy: 0.9}, validation: {accuracy: 0.4}});
    expect(far.length, 1, JSON.stringify(far));
    expect(far[0].kind, 'mismatch');
    expectArray(far[0].splits, ['validation']);
    expect(far[0].text, 'The validation Accuracy is far from the training one.');
    expect(far[0].detail, 'Train 0.900, validation 0.400.');

    expect(metricFindings({train: {accuracy: 0.97, f1: 0.97, auc: 0.99}, validation: {accuracy: 0.95, f1: 0.95,
      auc: 0.98}, baseline: {accuracy: 0.5, f1: 0.5, auc: 0.5}}).length, 0);
    expect(metricFindings({train: {auc: 0.52}, validation: {auc: 0.49}}).length, 0, 'no threshold finding');
    expect(metricFindings({train: {accuracy: 0.4}, validation: {accuracy: 0.9}}).length, 0);
    const error = metricFindings({train: {rmse: 1}, validation: {rmse: 1.6}});
    expectArray(error.map((f) => `${f.id} ${f.kind}`), ['rmse mismatch']);
  });

  test('baselineMetrics: the trivial model from the class shares, R2 0 for a regressor', async () => {
    const shares = (counts: [string, number][]) => DG.Column.fromStrings('t',
      counts.flatMap(([label, n]) => Array<string>(n).fill(label)));
    const four = baselineMetrics('classification', shares([['a', 90], ['b', 3], ['c', 3], ['d', 4]]));
    expectFloat(four.balancedAccuracy ?? NaN, 0.25, 1e-9, 'balancedAccuracy');
    expectFloat(four.accuracy ?? NaN, 0.9, 1e-9, 'accuracy');
    expect(four.mcc, 0);
    expect(four.auc, 0.5);
    expectFloat(four.averagePrecision ?? NaN, 0.25, 1e-9, 'averagePrecision');
    expectFloat(four.f1 ?? NaN, 0.25, 1e-9, 'f1');
    expectFloat(four.logLoss ?? NaN, 0.433973, 1e-6, 'logLoss');
    expectFloat(four.brier ?? NaN, 0.1866, 1e-9, 'brier');

    const two = baselineMetrics('classification', shares([['yes', 70], ['no', 30]]), 'yes');
    expectFloat(two.averagePrecision ?? NaN, 0.7, 1e-9, 'two-class averagePrecision');
    expectFloat(two.f1 ?? NaN, 0.7, 1e-9, 'two-class f1');
    expectFloat(two.brier ?? NaN, 0.21, 1e-9, 'two-class brier');
    expectFloat(two.logLoss ?? NaN, 0.610864, 1e-6, 'two-class logLoss');

    const gapped = baselineMetrics('classification', shares([['a', 3], ['', 4], ['b', 1]]));
    expectFloat(gapped.balancedAccuracy ?? NaN, 0.5, 1e-9, 'the empty category is not a class');
    expectFloat(gapped.accuracy ?? NaN, 0.75, 1e-9, 'the empty rows are not counted');
    expect(JSON.stringify(baselineMetrics('regression', numbers([1, 2, 3]))), JSON.stringify({r2: 0}));
  });

  test('underfittingOf: the skill below 0.1 on the shown values; the delegated metrics', async () => {
    const ba = (train: number) => underfittingOf('balancedAccuracy', {balancedAccuracy: train},
      {balancedAccuracy: 0.25})?.isUnderfitting;
    expect(ba(0.32), true);
    expect(ba(0.33), false);
    expect(ba(0.325), false, 'a skill of exactly 0.100');
    const logLoss = (train: number) => underfittingOf('logLoss', {logLoss: train}, {logLoss: 0.434})?.isUnderfitting;
    expect(logLoss(0.40), true);
    expect(logLoss(0.38), false);
    const r2 = (train: number) => underfittingOf('r2', {r2: train}, {r2: 0})?.isUnderfitting;
    expect(r2(0.099), true);
    expect(r2(0.1), false);
    const rmse = underfittingOf('rmse', {rmse: 1, r2: 0.05}, {r2: 0});
    expect(`${rmse?.id} ${rmse?.isUnderfitting}`, 'r2 true');
    const recall = underfittingOf('sensitivity', {sensitivity: 0.9, balancedAccuracy: 0.6}, {balancedAccuracy: 0.5});
    expect(`${recall?.id} ${recall?.isUnderfitting}`, 'balancedAccuracy false');
    expect(underfittingOf('auc', {auc: 0.6}, {}) === undefined, true, 'without a base');
    expect(underfittingOf('auc', {}, {auc: 0.5}) === undefined, true, 'without a train value');
  });

  test('verdictOf: underfitting, then overfitting, then no overfitting', async () => {
    const baseline = {balancedAccuracy: 0.25};
    const statusOf = (train: number, validation: number, base?: MetricValues) => verdictOf({
      train: {balancedAccuracy: train}, validation: {balancedAccuracy: validation},
      ...(base === undefined ? {} : {baseline: base})}, 'balancedAccuracy')?.status;
    expect(statusOf(0.3, 0.1, baseline), 'underfitting', 'an underfitting record that also overfits');
    expect(statusOf(0.9, 0.6, baseline), 'overfitting');
    expect(statusOf(0.9, 0.88, baseline), 'noOverfitting');
    expect(statusOf(0.3, 0.29), 'noOverfitting', 'never underfitting without a baseline');
    expect(verdictOf({train: {balancedAccuracy: 0.3}, validation: {}}, 'balancedAccuracy') === undefined, true,
      'no validation value and no baseline');
    expect(verdictOf({train: {balancedAccuracy: 0.3}, validation: {}, baseline}, 'balancedAccuracy')?.status,
      'underfitting', 'underfitting without a validation value');
  });

  test('overfittingOf: 0.05 for the scores, 0.1 for MCC and R2, 1.5 times for the errors', async () => {
    const accuracy = overfittingOf('accuracy', 0.96, 0.90);
    expect(accuracy.isOverfitting, true);
    expectFloat(accuracy.gap ?? NaN, 0.06, 1e-9);
    expect(overfittingOf('accuracy', 0.95, 0.91).isOverfitting, false);
    expect(overfittingOf('auc', 0.90, 0.95).isOverfitting, false);
    expect(overfittingOf('mcc', 0.9, 0.85).isOverfitting, false);
    expect(overfittingOf('mcc', 0.9, 0.79).isOverfitting, true);
    expect(overfittingOf('r2', 0.8, 0.65).isOverfitting, true);
    const logLoss = overfittingOf('logLoss', 0.1, 0.16);
    expect(logLoss.isOverfitting, true);
    expect(`${logLoss.ratio} ${logLoss.gap}`, '1.6 undefined');
    expect(overfittingOf('logLoss', 0.1, 0.14).isOverfitting, false);
    expect(overfittingOf('brier', 0, 0).isOverfitting, false);
    expect(overfittingOf('mse', 0, 0.1).isOverfitting, true);
  });

  test('overfittingOf decides on the values rounded to 3 decimals, as shown', async () => {
    const shown = overfittingOf('accuracy', 0.9734, 0.9466);
    expect(shown.gap, 0.026);
    const atLimit = overfittingOf('accuracy', 0.96, 0.91);
    expect(`${atLimit.gap} ${atLimit.isOverfitting}`, '0.05 false');
    const roundedDown = overfittingOf('accuracy', 0.95049, 0.9);
    expect(`${roundedDown.gap} ${roundedDown.isOverfitting}`, '0.05 false');
    expect(overfittingOf('accuracy', 0.9506, 0.9).isOverfitting, true);
    expect(overfittingOf('mcc', 0.9, 0.8).isOverfitting, false);
    const ratioOf = (train: number, validation: number) => {
      const {kind, ratio, isOverfitting} = overfittingOf('logLoss', train, validation);
      return `${kind} ${ratio} ${isOverfitting}`;
    };
    expect(ratioOf(0.1, 0.15), 'ratio 1.5 false');
    expect(ratioOf(0.0996, 0.1497), 'ratio 1.5 false', 'The raw values overfit, the shown ones not');
    expect(ratioOf(0.1, 0.151), 'ratio 1.51 true');
    expect(ratioOf(0.123, 0.185), 'ratio 1.5 false', 'The shown ratio is 1.50');
    expect(ratioOf(0, 0.2), 'ratio undefined true');
    expect(ratioOf(0, 0), 'ratio undefined false');
    const {kind, ratio} = overfittingOf('accuracy', 0.9, 0.8);
    expect(`${kind} ${ratio}`, 'gap undefined');
  });

  test('METRIC_IDS ends with the Crux and the probability ids; every id has a label and a description', async () => {
    expectArray(METRIC_IDS.slice(-5), ['balancedAccuracy', 'mcc', 'averagePrecision', 'logLoss', 'brier']);
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

  test('metricsRecordOf keeps test and baseline', async () => {
    const record = metricsRecordOf({train: {r2: 0.9}, validation: {r2: 0.7}, test: {r2: 0.5}, baseline: {r2: 0}});
    expect(record?.test?.r2, 0.5);
    expect(record?.baseline?.r2, 0);
    expect(metricsRecordOf({train: {}, validation: {}})?.test === undefined, true, 'test without one');
    expect(metricsRecordOf({train: {}, validation: {}})?.baseline === undefined, true, 'baseline without one');
  });

  test('SEARCH_METRICS and the directions', async () => {
    const ids: readonly string[] = METRIC_IDS;
    const all = [...SEARCH_METRICS.regression, ...SEARCH_METRICS.classification];
    expectArray(all.filter((id) => !ids.includes(id)), []);
    expect(directionOf('rmse'), -1);
    expect(directionOf('mae'), -1);
    expect(directionOf('r2'), 1);
    expect(directionOf('logLoss'), -1);
    expect(directionOf('brier'), -1);
    expectArray(SEARCH_METRICS.classification.slice(-2), ['logLoss', 'brier']);
    expect(defaultMetric('classification'), 'balancedAccuracy');
    expect(defaultMetric('regression'), 'r2');
    expect(defaultMetricOf(SEARCH_METRICS.classification), defaultMetric('classification'));
    expect(defaultMetricOf(SEARCH_METRICS.regression), defaultMetric('regression'));
    expect(defaultMetricOf(METRIC_IDS, 'f1'), 'f1');
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
    expect(findings[0].text, 'The validation MCC is far from the training one.');
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

  test('probabilityMetrics: log loss and Brier score match hand-computed values', async () => {
    const doubles = (values: number[]) => DG.Column.fromFloat64Array('p', Float64Array.from(values));
    const two: ClassProbabilities = {classes: ['a', 'b'],
      columns: [doubles([0.8, 0.3, 0.6]), doubles([0.2, 0.7, 0.4])]};
    const binary = probabilityMetrics(labels(['a', 'b', 'a']), two, 'a');
    expectFloat(binary.logLoss ?? NaN, -(Math.log(0.8) + Math.log(0.7) + Math.log(0.6)) / 3, 1e-12, 'logLoss');
    expectFloat(binary.brier ?? NaN, (0.2 ** 2 + 0.3 ** 2 + 0.4 ** 2) / 3, 1e-12, 'brier');

    const classes = ['a', 'b', 'c'];
    const oneHot: ClassProbabilities = {classes,
      columns: [doubles([1, 0, 0]), doubles([0, 1, 0]), doubles([0, 0, 1])]};
    const perfect = probabilityMetrics(labels(['a', 'b', 'c']), oneHot);
    expectFloat(perfect.logLoss ?? NaN, 0, 1e-12, 'one-hot logLoss');
    expectFloat(perfect.brier ?? NaN, 0, 1e-12, 'one-hot brier');
    const spread: ClassProbabilities = {classes,
      columns: [doubles([0.5, 0.2]), doubles([0.3, 0.6]), doubles([0.2, 0.2])]};
    const summed = probabilityMetrics(labels(['a', 'b']), spread);
    expectFloat(summed.logLoss ?? NaN, -(Math.log(0.5) + Math.log(0.6)) / 2, 1e-12, 'three-class logLoss');
    expectFloat(summed.brier ?? NaN, (0.5 ** 2 + 0.3 ** 2 + 0.2 ** 2 + 0.2 ** 2 + 0.4 ** 2 + 0.2 ** 2) / 2, 1e-12,
      'three-class brier');

    const wrong = probabilityMetrics(labels(['a']), {classes: ['a', 'b'], columns: [doubles([0]), doubles([1])]}, 'a');
    expectFloat(wrong.logLoss ?? NaN, -Math.log(1e-15), 1e-9, 'clipped logLoss');
    expect(Number.isFinite(wrong.logLoss), true, `logLoss ${wrong.logLoss}`);
    const unseen = probabilityMetrics(labels(['c', 'a']), {classes: ['a', 'b'],
      columns: [doubles([0.3, 0.8]), doubles([0.7, 0.2])]});
    expectFloat(unseen.logLoss ?? NaN, (-Math.log(1e-15) - Math.log(0.8)) / 2, 1e-9, 'logLoss of an unseen class');
    expectFloat(unseen.brier ?? NaN, (0.3 ** 2 + 0.7 ** 2 + 0.2 ** 2 + 0.2 ** 2) / 2, 1e-12, 'unseen class brier');

    const gapped = probabilityMetrics(labels(['a', 'b', 'a']), {classes: ['a', 'b'],
      columns: [numbers([0.8, null, 0.6]), numbers([0.2, null, 0.4])]}, 'a');
    expectFloat(gapped.logLoss ?? NaN, -(Math.log(0.8) + Math.log(0.6)) / 2, 1e-6, 'logLoss without the gap');
    expectFloat(gapped.brier ?? NaN, (0.2 ** 2 + 0.4 ** 2) / 2, 1e-6, 'brier without the gap');
  });

  test('bestThreshold: the midpoint below the best candidate, ties to the one nearest 0.5', async () => {
    const separated = numbers([0.9, 0.6, 0.4, 0.1]);
    const overlapping = numbers([0.9, 0.8, 0.3, 0.1]);
    for (const metric of BEST_THRESHOLD_METRICS) {
      expect(bestThreshold(labels(['a', 'a', 'b', 'b']), separated, 'a', metric), 0.5, `${metric} separated`);
      expect(bestThreshold(labels(['a', 'b', 'a', 'b']), overlapping, 'a', metric), 0.2, `${metric} overlapping`);
    }
    expect(bestThreshold(labels(['a', 'a', 'a', 'a']), separated, 'a', 'f1') === undefined, true,
      'A threshold of one class');
    // 0.7 is stored as 0.699999988: the rounded value is checked in single precision, as the step cuts.
    expect(bestThreshold(labels(['a', 'b']), numbers([0.7, 0.6995]), 'a', 'f1'), 0.7, 'A single-precision candidate');
    expectArray([...BEST_THRESHOLD_METRICS], ['f1', 'balancedAccuracy', 'mcc', 'accuracy']);
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
