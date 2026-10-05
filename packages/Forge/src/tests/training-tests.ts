import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';
import {ForgeError} from '../forge-error';
import {defaultFeatures} from '../training/default-features';
import {kFold} from '../training/k-fold';
import {releaseFrame} from '../preparation/shared-frame';
import {datasetFingerprint} from '../storage/dataset-fingerprint';
import {checkTrainable, prepareTraining, TrainingRequest, trainingProblems, trainModel}
  from '../training/train-model';
import {columnsOf, expectMetrics, expectReleased, framesSharing, IMPUTE, IRIS, MEASUREMENTS, openIris, requestOf,
  selectionOf, valuesOf} from './test-data';

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

  test('defaultFeatures skips the target, all-unique integer, date and very large integer columns', async () => {
    const table = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'row', [1, 2, 3, 4, 5, 6]),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'id', [101, 205, 303, 404, 550, 606]),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, 'count', [1, 1, 2, 3, 3, 4]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'label', ['a', 'b', 'a', 'b', 'a', 'b']),
      DG.Column.dateTime('date', 6),
      DG.Column.fromBigInt64Array('big', new BigInt64Array([1n, 1n, 2n, 3n, 3n, 4n])),
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
    const date = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'date', values.map((v) => new Date(2020, 0, v)));
    const big = DG.Column.fromBigInt64Array('big', new BigInt64Array(values.map((v) => BigInt(v))));
    const firstSix = DG.BitSet.create(12, (i) => i < 6);

    let problems = await trainingProblems(selectionOf([cls], y));
    expect(problems.features.length, 1);
    expect(problems.features[0].includes('cls') && problems.features[0].includes('numerical'), true,
      problems.features[0]);
    expect(problems.target.length, 0);

    problems = await trainingProblems(selectionOf([y], y));
    expect(problems.features.some((m) => m.includes('also a feature')), true, problems.features.join(' '));

    problems = await trainingProblems(selectionOf([n.clone(firstSix)], y.clone(firstSix)));
    expect(problems.target.some((m) => m.includes('at least 10 rows')), true, problems.target.join(' '));

    problems = await trainingProblems(selectionOf([n], yWithGap));
    expect(problems.target.length, 0, problems.target.join(' '));
    const request = await prepareTraining(selectionOf([n], yWithGap));
    expect(request.options.missingValues?.skippedRows, 1);
    expect(request.target.length, 11);

    problems = await trainingProblems(selectionOf([n], constant));
    expect(problems.target.some((m) => m.includes('only one value')), true, problems.target.join(' '));

    problems = await trainingProblems(selectionOf([n, date], y));
    expect(problems.features.some((m) => m.includes('needs numerical features. Uncheck: date.')), true,
      problems.features.join(' '));

    problems = await trainingProblems(selectionOf([n, big], y));
    expect(problems.features.includes('\'big\' holds very large whole numbers, which XGBoost cannot read. ' +
      'Convert the column to a decimal type or uncheck it.'), true, problems.features.join(' '));

    problems = await trainingProblems(selectionOf([n], big));
    expect(problems.target.includes('\'big\' holds very large whole numbers, which XGBoost cannot read. ' +
      'Convert the column to a decimal type or choose another target.'), true, problems.target.join(' '));

    const nWithGap = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', values.map((v) => v === 3 ? null : v));
    problems = await trainingProblems(selectionOf([nWithGap], y, IMPUTE));
    expectArray(problems.missingValues,
      ['Impute needs at least two features. Choose Skip rows or check more features.']);
    expect((await trainingProblems(selectionOf([nWithGap], y))).missingValues.length, 0);

    const valid = selectionOf([n], y);
    const [validProblems, frames] = await framesSharing([n], () => trainingProblems(valid));
    let messages = [...validProblems.target, ...validProblems.features, ...validProblems.missingValues];
    expect(messages.length, 0, messages.join(' '));
    expectReleased(frames);
    await checkTrainable(valid);

    const invalid = selectionOf([cls], constant);
    problems = await trainingProblems(invalid);
    messages = [...problems.target, ...problems.features];
    expect(messages.length, 2, messages.join(' '));
    await expectExceptionAsync(() => checkTrainable(invalid),
      (e) => e instanceof ForgeError && messages.every((m) => e.message.includes(m)));
  });

  test('datetime target is refused', async () => {
    const values = Array.from({length: 12}, (_, i) => i + 1);
    const n = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', values);
    const date = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'date', values.map((v) => new Date(2020, 0, v)));
    const problems = await trainingProblems(selectionOf([n], date));
    expect(problems.target.some((m) => m.includes('Choose a numerical, text or boolean target.')), true,
      problems.target.join(' '));
  });

  test('prepareTraining skips rows with a missing target in both modes', async () => {
    const rows = 12;
    const a = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', valuesOf(rows, (i) => i, [2]));
    const b = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => 2 * i));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(rows, (i) => 3 * i, [9]));
    const features = DG.DataFrame.fromColumns([a, b]);

    const [[skipped, imputed], frames] = await framesSharing([a, b], async () =>
      [await prepareTraining(selectionOf(features, y)), await prepareTraining(selectionOf(features, y, IMPUTE))]);
    expect(skipped.options.missingValues?.skippedRows, 2);
    expect(skipped.target.length, 10);
    expect(JSON.stringify(skipped.options), JSON.stringify({preprocessingInfo: ['ignore-missing'],
      postprocessingInfo: [], missingValues: {mode: 'skip', skippedRows: 2}}));

    expect(imputed.options.missingValues?.skippedRows, 1);
    expect(imputed.target.length, 11);
    expect(imputed.features.getCol('a').isNone(2), false, 'The gap in a is not imputed');
    expect(JSON.stringify(imputed.options), JSON.stringify({preprocessingInfo: ['impute-missing'],
      postprocessingInfo: [], missingValues: {mode: 'impute', neighbors: 4, distance: 'Euclidean', skippedRows: 1}}));

    expect(a.isNone(2) && y.isNone(9), true, 'The input columns changed');
    expect(features.rowCount, rows);
    releaseFrame(skipped.features);
    releaseFrame(imputed.features);
    expectReleased(frames);
  }, {timeout: TIMEOUT});

  test('Skip rows leaves a text target only the classes its rows have', async () => {
    const rows = 20;
    const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', valuesOf(rows, (i) => i));
    const z = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'z', valuesOf(rows, (i) => (i * 7) % 5));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'y',
      Array.from({length: rows}, (_, i) => i === 4 ? null : i < 10 ? 'a' : 'b'));
    const request = await prepareTraining(selectionOf([x, z], y));
    expectArray(request.target.categories, ['a', 'b']);
    const result = await trainModel(request);
    releaseFrame(request.features);
    expect(result.metrics.positiveClass, 'a');
    for (const id of ['sensitivity', 'specificity', 'precision', 'npv'] as const)
      expect(result.metrics.validation[id] !== undefined, true, `${id} is missing`);
    expect(y.categories.includes(''), true, 'The user\'s column lost its empty category');

    const iris = await openIris();
    iris.getCol('Species').set(20, null);
    const irisRequest = await prepareTraining(selectionOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species')));
    expect(irisRequest.target.categories.length, 3, irisRequest.target.categories.join(', '));
    const fingerprint = datasetFingerprint(irisRequest.features, irisRequest.target);
    const irisResult = await trainModel(irisRequest);
    releaseFrame(irisRequest.features);
    expect(irisResult.target.categories?.length, 3);
    expect(irisResult.target.categories?.includes(''), false, 'The stored target lists an empty class');
    expect(fingerprint.columns[4].categories?.includes(''), false, 'The data summary lists an empty class');
  }, {timeout: TIMEOUT});

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
    expect(JSON.stringify(result.options), JSON.stringify({preprocessingInfo: [], postprocessingInfo: [],
      missingValues: {mode: 'skip', skippedRows: 0}}));
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
    expectMetrics(second, first, 1e-9);
  }, {timeout: TIMEOUT});

  test('trainingProblems on iris: defaultFeatures excludes the row-number column and a text feature is reported',
    async () => {
      const iris = await grok.data.files.openTable(IRIS);
      expectArray(defaultFeatures(iris, iris.getCol('Species')).map((c) => c.name), MEASUREMENTS);
      const features = ['Sepal.Length', 'Sepal.Width', 'Petal.Width', 'Species'];
      const problems = await trainingProblems(selectionOf(iris.clone(null, features), iris.getCol('Petal.Length')));
      expect(problems.features.some((m) => m.includes('Species')), true, problems.features.join(' '));
    }, {timeout: TIMEOUT});

  test('binary classification records the positive class', async () => {
    const iris = await grok.data.files.openTable(IRIS);
    const species = iris.getCol('Species');
    const excluded = species.categories[2];
    const twoSpecies = iris.clone(DG.BitSet.create(iris.rowCount, (i) => species.get(i) !== excluded));
    const result = await trainModel(await requestOf(twoSpecies.clone(null, MEASUREMENTS),
      twoSpecies.getCol('Species')));
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
