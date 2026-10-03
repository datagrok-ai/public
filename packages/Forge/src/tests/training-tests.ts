import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {ForgeError} from '../forge-error';
import {METRIC_IDS} from '../metrics/metrics';
import {defaultFeatures} from '../training/default-features';
import {kFold} from '../training/k-fold';
import {checkTrainable, TrainingRequest, trainingProblems, trainModel} from '../training/train-model';
import {IRIS, MEASUREMENTS, requestOf} from './test-data';

const TIMEOUT = 60000;

async function irisRequest(features: string[], target: string): Promise<TrainingRequest> {
  const iris = await grok.data.files.openTable(IRIS);
  return requestOf(iris.clone(null, features), iris.getCol(target));
}

category('Training', () => {
  test('kFold is deterministic and balanced', async () => {
    const folds = kFold(103, 5, 42);
    expectArray(kFold(103, 5, 42), folds);
    const sizes = [0, 0, 0, 0, 0];
    for (const f of folds)
      sizes[f]++;
    expectArray(sizes, [21, 21, 21, 20, 20]);
    expect(kFold(103, 5, 43).some((f, i) => f !== folds[i]), true, 'Another seed gives the same folds');
  });

  test('defaultFeatures skips the target and all-unique integer columns', async () => {
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'row', [1, 2, 3, 4, 5, 6]),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'id', [101, 205, 303, 404, 550, 606]),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'count', [1, 1, 2, 3, 3, 4]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'label', ['a', 'b', 'a', 'b', 'a', 'b']),
    ]);
    expectArray(defaultFeatures(table, table.getCol('label')).map((c) => c.name), ['count', 'x']);
    expectArray(defaultFeatures(table, table.getCol('x')).map((c) => c.name), ['count']);
  });

  test('trainingProblems reports per input', async () => {
    const values = Array.from({length: 12}, (_, i) => i + 1);
    const n = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', values);
    const cls = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'cls', values.map((v) => v % 2 === 0 ? 'a' : 'b'));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', values.map((v) => v * 1.5 + (v % 3)));
    const yWithGap = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', values.map((v) => v === 5 ? null : v * 1.5));
    const constant = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c', values.map(() => 'a'));
    const firstSix = DG.BitSet.create(12, (i) => i < 6);

    let problems = await trainingProblems(requestOf([cls], y));
    expect(problems.features.length, 1);
    expect(problems.features[0].includes('cls') && problems.features[0].includes('numerical'), true,
      problems.features[0]);
    expect(problems.target.length, 0);

    problems = await trainingProblems(requestOf([y], y));
    expect(problems.features.some((m) => m.includes('also a feature')), true, problems.features.join(' '));

    problems = await trainingProblems(requestOf([n.clone(firstSix)], y.clone(firstSix)));
    expect(problems.target.some((m) => m.includes('at least 10 rows')), true, problems.target.join(' '));

    problems = await trainingProblems(requestOf([n], yWithGap));
    expect(problems.target.some((m) => m.includes('empty values')), true, problems.target.join(' '));

    problems = await trainingProblems(requestOf([n], constant));
    expect(problems.target.some((m) => m.includes('only one value')), true, problems.target.join(' '));

    const valid = requestOf([n], y);
    problems = await trainingProblems(valid);
    let messages = [...problems.target, ...problems.features];
    expect(messages.length, 0, messages.join(' '));
    await checkTrainable(valid);

    const invalid = requestOf([cls], yWithGap);
    problems = await trainingProblems(invalid);
    messages = [...problems.target, ...problems.features];
    expect(messages.length, 2, messages.join(' '));
    await expectExceptionAsync(() => checkTrainable(invalid),
      (e) => e instanceof ForgeError && messages.every((m) => e.message.includes(m)));
  });

  test('trains a classifier on iris', async () => {
    const result = await trainModel(await irisRequest(MEASUREMENTS, 'Species'));
    const {train, validation} = result.metrics;
    expect(result.task, 'classification');
    expect((validation.accuracy ?? 0) > 0.85, true, `Validation accuracy ${validation.accuracy}`);
    expect((train.accuracy ?? 0) >= (validation.accuracy ?? 0) - 0.05, true,
      `Train accuracy ${train.accuracy}, validation ${validation.accuracy}`);
    expect(validation.f1 !== undefined, true, 'F1 is missing');
    expect(validation.sensitivity === undefined, true, 'Sensitivity on three classes');
    expect(result.blob.length > 0, true, 'Empty blob');
    expect(JSON.stringify(result.splitting), JSON.stringify({scheme: 'kfold', folds: 5, isStratified: false}));
    expect(JSON.stringify(result.options), JSON.stringify({preprocessingInfo: [], postprocessingInfo: []}));
    expect(result.target.categories?.length, 3);
    expectArray(result.features.columns.map((c) => c.name), MEASUREMENTS);
  }, {timeout: TIMEOUT});

  test('trains a regressor on iris', async () => {
    const features = MEASUREMENTS.filter((name) => name !== 'Petal.Length');
    const result = await trainModel(await irisRequest(features, 'Petal.Length'));
    const {validation} = result.metrics;
    expect(result.task, 'regression');
    expect((validation.r2 ?? 0) > 0.8, true, `Validation r2 ${validation.r2}`);
    for (const id of ['mse', 'rmse', 'mae'] as const)
      expect(validation[id] !== undefined, true, `${id} is missing`);
    expect(validation.accuracy === undefined, true, 'Accuracy on a regression');
  }, {timeout: TIMEOUT});

  test('the same seed gives the same validation metrics', async () => {
    const first = (await trainModel(await irisRequest(MEASUREMENTS, 'Species'))).metrics.validation;
    const second = (await trainModel(await irisRequest(MEASUREMENTS, 'Species'))).metrics.validation;
    for (const id of METRIC_IDS) {
      const value = first[id];
      if (value === undefined)
        expect(second[id] === undefined, true, id);
      else
        expectFloat(second[id] ?? NaN, value, 1e-9, id);
    }
  }, {timeout: TIMEOUT});

  test('trainingProblems on iris: defaultFeatures excludes the row-number column and a text feature is reported',
    async () => {
      const iris = await grok.data.files.openTable(IRIS);
      expectArray(defaultFeatures(iris, iris.getCol('Species')).map((c) => c.name), MEASUREMENTS);
      const features = ['Sepal.Length', 'Sepal.Width', 'Petal.Width', 'Species'];
      const problems = await trainingProblems(requestOf(iris.clone(null, features), iris.getCol('Petal.Length')));
      expect(problems.features.some((m) => m.includes('Species')), true, problems.features.join(' '));
    }, {timeout: TIMEOUT});

  test('binary classification records the positive class', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    const species = iris.getCol('Species');
    const excluded = species.categories[2];
    const twoSpecies = iris.clone(DG.BitSet.create(iris.rowCount, (i) => species.get(i) !== excluded));
    const result = await trainModel(requestOf(twoSpecies.clone(null, MEASUREMENTS), twoSpecies.getCol('Species')));
    expect(result.target.categories?.length, 2);
    expect(result.metrics.positiveClass, result.target.categories?.[0]);
    for (const id of ['sensitivity', 'specificity', 'precision', 'npv'] as const)
      expect(result.metrics.validation[id] !== undefined, true, `${id} is missing`);
  }, {timeout: TIMEOUT});

  test('cancellation stops training', async () => {
    const request = await irisRequest(MEASUREMENTS, 'Species');
    await expectExceptionAsync(async () => {
      await trainModel(request, {canceled: true, update: () => {}});
    }, (e) => e instanceof ForgeError && e.message.includes('cancelled'));
  }, {timeout: TIMEOUT});
});
