import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {defaultHyperparameters, Engine} from '../engines/engine';
import {applicableEngines} from '../engines/applicable-engines';
import {apply, LoopProgress, yieldToEventLoop} from '../engines/engine-calls';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {defaultFeatures} from '../training/default-features';
import {kFold} from '../training/k-fold';
import {aucOf} from '../metrics/metrics';
import {MissingValuesSettings, prepareMissingValues} from '../preparation/missing-values';
import {prepareFeatures, replayPreprocessing} from '../preparation/pipeline';
import {PreparationOptions} from '../preparation/preparation-options';
import {datasetFingerprint} from '../storage/dataset-fingerprint';
import {checkSelection, prepareTraining, recutResult, TrainingProblems, TrainingRequest, TrainingResult,
  TrainingSelection, trainModel} from '../training/train-model';
import {QueuedTraining, TrainingQueue} from '../training/training-queue';
import {columnNamed, columnsOf, engineByName, expectMetrics, expectNoFrames, expectReleased, framesSharing, IMPUTE,
  inDiscoveryOrder, IRIS, MEASUREMENTS, names, NO_OPTIONS, ONE_HOT_STEPS, openIris,
  PREDICT_PROBABILITY, requestOf, selectionOf, SKIP_UNIQUE_STEPS, twoSpeciesIris, valuesOf} from './test-data';

const TIMEOUT = 60000;
const ENGINE_TIMEOUT = 120000;
const SKIP: MissingValuesSettings = {mode: 'skip'};
/** The targets each EDA method can learn on iris: Species (classification), Petal.Length (regression). */
const ENGINE_TARGETS: [string, string[]][] = [['XGBoost', ['Species', 'Petal.Length']],
  ['SVM', ['Species', 'Petal.Length']], ['Softmax', ['Species']], ['Linear Regression', ['Petal.Length']],
  ['PLS Regression', ['Petal.Length']]];
const NOT_SOFTMAX = 'Softmax cannot learn from this selection. ' +
  'It needs numerical features and a numerical, text or boolean target.';
const NO_METHOD = 'No method can learn from this selection. Check the features and the target.';
const EDA_METHODS = ['Linear Regression', 'Softmax', 'PLS Regression', 'XGBoost', 'SVM'];

/** A selection of [engine] with its default hyperparameters. */
function engineSelection(engine: Engine, features: DG.Column[], target: DG.Column,
  missingValues: MissingValuesSettings = SKIP): TrainingSelection {
  return {...selectionOf(features, target, missingValues), engine, hyperparameters: defaultHyperparameters(engine)};
}

/** [engine] whose train and apply calls are listed in [calls] and record, in [gaps], every table or column passed
 * with a missing value. */
function gapChecked(engine: Engine, calls: string[], gaps: string[]): Engine {
  const functions = {...engine.functions};
  for (const role of ['train', 'apply'] as const) {
    const func = engine.functions[role];
    if (func === undefined)
      continue;
    const spy: DG.Func = Object.create(func);
    spy.apply = async (parameters: {[name: string]: unknown} = {}) => {
      calls.push(role);
      for (const value of Object.values(parameters)) {
        const columns = value instanceof DG.DataFrame ? value.columns.toList() : value instanceof DG.Column ?
          [value] : [];
        for (const col of columns.filter((c) => c.stats.missingValueCount > 0))
          gaps.push(`${role}: ${col.name}`);
      }
      return func.apply(parameters);
    };
    functions[role] = spy;
  }
  return {...engine, functions};
}

/** A stub training of [fits] steps that checks the progress before each one, as `trainModel` does. */
function stubTraining(fits: number, done: number[]): (progress: LoopProgress) => Promise<number> {
  return async (progress) => {
    for (let i = 0; i < fits; i++) {
      await yieldToEventLoop();
      if (progress.canceled)
        throw new ForgeError('Training was cancelled.');
      done.push(i);
    }
    return fits;
  };
}

/** The problems [checkSelection] finds in [selection], with only the selection's method to ask. */
async function problemsOf(selection: TrainingSelection): Promise<TrainingProblems> {
  return (await checkSelection(selection, [selection.engine])).problems;
}

/** A task-bar indicator stand-in whose cancel is clicked [ms] after it is created. */
function cancelledAfter(ms: number): () => LoopProgress {
  return () => {
    const indicator = {canceled: false, update: () => {}};
    setTimeout(() => indicator.canceled = true, ms);
    return indicator;
  };
}

const errorOf = (queued: QueuedTraining<unknown>) =>
  queued.outcome === 'failed' || queued.outcome === 'cancelled' ? queued.error : undefined;

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

  test('checkSelection reports the rules per input', async () => {
    const values = Array.from({length: 12}, (_, i) => i + 1);
    const n = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', values);
    const cls = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'cls', values.map((v) => v % 2 === 0 ? 'a' : 'b'));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', values.map((v) => v * 1.5 + (v % 3)));
    const yWithGap = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', values.map((v) => v === 5 ? null : v * 1.5));
    const constant = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c', values.map(() => 'a'));
    const date = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'date', values.map((v) => new Date(2020, 0, v)));
    const big = DG.Column.fromBigInt64Array('big', new BigInt64Array(values.map((v) => BigInt(v))));
    const firstSix = DG.BitSet.create(12, (i) => i < 6);

    let problems = await problemsOf(selectionOf([cls], y));
    expect(problems.features.length, 1);
    expect(problems.features[0].includes('cls') && problems.features[0].includes('numerical'), true,
      problems.features[0]);
    expect(problems.target.length, 0);

    problems = await problemsOf(selectionOf([y], y));
    expect(problems.features.some((m) => m.includes('also a feature')), true, problems.features.join(' '));

    problems = await problemsOf(selectionOf([n.clone(firstSix)], y.clone(firstSix)));
    expect(problems.target.some((m) => m.includes('at least 10 rows')), true, problems.target.join(' '));

    problems = await problemsOf(selectionOf([n], yWithGap));
    expect(problems.target.length, 0, problems.target.join(' '));
    const request = await prepareTraining(selectionOf([n], yWithGap));
    expect(request.prepared.options.missingValues?.skippedRows, 1);
    expect(request.target.length, 11);

    problems = await problemsOf(selectionOf([n], constant));
    expect(problems.target.some((m) => m.includes('only one value')), true, problems.target.join(' '));

    problems = await problemsOf(selectionOf([n, date], y));
    expect(problems.features.some((m) => m.includes('needs numerical features. Uncheck: date.')), true,
      problems.features.join(' '));

    problems = await problemsOf(selectionOf([n, big], y));
    expect(problems.features.includes('\'big\' holds very large whole numbers, which XGBoost cannot read. ' +
      'Convert the column to a decimal type or uncheck it.'), true, problems.features.join(' '));

    problems = await problemsOf(selectionOf([n], big));
    expect(problems.target.includes('\'big\' holds very large whole numbers, which XGBoost cannot read. ' +
      'Convert the column to a decimal type or choose another target.'), true, problems.target.join(' '));

    const nWithGap = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', values.map((v) => v === 3 ? null : v));
    problems = await problemsOf(selectionOf([nWithGap], y, IMPUTE));
    expectArray(problems.missingValues,
      ['Impute needs at least two features. Choose Skip rows or check more features.']);
    expect((await problemsOf(selectionOf([nWithGap], y))).missingValues.length, 0);

    const valid = selectionOf([n], y);
    const [validProblems, frames] = await framesSharing([n], () => problemsOf(valid));
    let messages = [...validProblems.target, ...validProblems.features, ...validProblems.missingValues,
      ...validProblems.method];
    expect(messages.length, 0, messages.join(' '));
    expectReleased(frames);

    problems = await problemsOf(selectionOf([cls], constant));
    messages = [...problems.target, ...problems.features];
    expect(messages.length, 2, messages.join(' '));
  });

  test('datetime target is refused', async () => {
    const values = Array.from({length: 12}, (_, i) => i + 1);
    const n = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', values);
    const date = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'date', values.map((v) => new Date(2020, 0, v)));
    const problems = await problemsOf(selectionOf([n], date));
    expect(problems.target.some((m) => m.includes('Choose a numerical, text or boolean target.')), true,
      problems.target.join(' '));
  });

  test('prepareTraining skips rows with a missing target in both modes', async () => {
    const rows = 12;
    const a = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'a', valuesOf(rows, (i) => i, [2]));
    const b = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'b', valuesOf(rows, (i) => 2 * i));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(rows, (i) => 3 * i, [9]));
    const features = [a, b];

    const [[skipped, imputed], frames] = await framesSharing(features, async () =>
      [await prepareTraining(selectionOf(features, y)), await prepareTraining(selectionOf(features, y, IMPUTE))]);
    expect(skipped.prepared.options.missingValues?.skippedRows, 2);
    expect(skipped.target.length, 10);
    expect(JSON.stringify(skipped.prepared.options), JSON.stringify({preprocessingInfo: ['ignore-missing'],
      postprocessingInfo: [], missingValues: {mode: 'skip', skippedRows: 2}}));

    expect(imputed.prepared.options.missingValues?.skippedRows, 1);
    expect(imputed.target.length, 11);
    expect(columnNamed(imputed.features, 'a').isNone(2), false, 'The gap in a is not imputed');
    expect(JSON.stringify(imputed.prepared.options), JSON.stringify({preprocessingInfo: ['impute-missing'],
      postprocessingInfo: [], missingValues: {mode: 'impute', neighbors: 4, distance: 'Euclidean', skippedRows: 1}}));

    expect(a.isNone(2) && y.isNone(9), true, 'The input columns changed');
    expect(a.length, rows);
    // A skipped row makes every column a copy, so the imputer's frame holds no column of the user's.
    expectNoFrames(frames);
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

  test('checkSelection on iris: defaultFeatures excludes the row-number column and a text feature is reported',
    async () => {
      const iris = await grok.data.files.openTable(IRIS);
      expectArray(defaultFeatures(iris, iris.getCol('Species')).map((c) => c.name), MEASUREMENTS);
      const features = ['Sepal.Length', 'Sepal.Width', 'Petal.Width', 'Species'];
      const problems = await problemsOf(selectionOf(iris.clone(null, features), iris.getCol('Petal.Length')));
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

  test('TrainingQueue: a newer training supersedes the running one before its next fit', async () => {
    const queue = new TrainingQueue();
    const first: number[] = [];
    const second: number[] = [];
    const firstRun = queue.run(stubTraining(5, first));
    while (first.length === 0)
      await yieldToEventLoop();
    // A fit follows its check in the same task, so no fit can start after this point.
    const fitsBefore = first.length;
    const secondRun = queue.run(stubTraining(3, second));
    const [older, newer] = await Promise.all([firstRun, secondRun]);
    expect(older.outcome, 'superseded');
    expect(first.length, fitsBefore, 'The superseded training made another fit');
    expect(fitsBefore < 5, true, 'The first training ended before it was superseded');
    expect(newer.outcome === 'completed' ? newer.result : undefined, 3);
    expectArray(second, [0, 1, 2]);
  });

  test('TrainingQueue: supersede() stops the running training without a newer one', async () => {
    const queue = new TrainingQueue();
    const fits: number[] = [];
    const run = queue.run(stubTraining(5, fits));
    while (fits.length === 0)
      await yieldToEventLoop();
    queue.supersede();
    expect((await run).outcome, 'superseded');
    expect(fits.length < 5, true, 'The superseded training made every fit');
  });

  test('TrainingQueue: a waiting training superseded in turn never starts and shows no indicator', async () => {
    const queue = new TrainingQueue();
    const fits: number[][] = [[], [], []];
    const indicators: number[] = [];
    const outcomes = await Promise.all(fits.map((done, i) => queue.run(stubTraining(3, done), () => {
      indicators.push(i);
      return {canceled: false, update: () => {}};
    })));
    expectArray(outcomes.map((o) => o.outcome), ['superseded', 'superseded', 'completed']);
    expect(fits[1].length, 0, 'The waiting training started');
    expect(indicators.includes(1), false, 'The waiting training created its indicator');
  });

  test('TrainingQueue: the indicator\'s cancel and update, and a failure', async () => {
    const queue = new TrainingQueue();
    const done: number[] = [];
    const long = async (progress: LoopProgress) => {
      const start = Date.now();
      while (Date.now() - start < 10000) {
        await yieldToEventLoop();
        if (progress.canceled)
          throw new ForgeError('Training was cancelled.');
        done.push(0);
      }
      return 0;
    };
    const cancelled = await queue.run(long, cancelledAfter(100));
    expect(cancelled.outcome, 'cancelled');
    expect(errorOf(cancelled) instanceof ForgeError, true, String(errorOf(cancelled)));
    expect(done.length > 0, true, 'The training was cancelled before it started');

    const updates: string[] = [];
    const reported = await queue.run(async (progress) => {
      progress.update(50, 'Fold 1 of 5');
      return 1;
    }, () => ({canceled: false, update: (_: number, description: string) => updates.push(description)}));
    expect(reported.outcome, 'completed');
    expectArray(updates, ['Fold 1 of 5']);

    const failed = await queue.run(async () => {
      throw new Error('forge-test: the method failed');
    });
    expect(failed.outcome, 'failed');
    expect(String(errorOf(failed)).includes('the method failed'), true, String(errorOf(failed)));
  });

  test('checkSelection lists the methods, suggests one and says whether it retrains live', async () => {
    const engines = EngineRegistry.discover();
    const edaNames = (list: Engine[]) => list.map((e) => e.name).filter((name) => ENGINE_TARGETS.some(([n]) =>
      n === name));
    const iris = await openIris();
    const features = columnsOf(iris, MEASUREMENTS);
    const [classifier, frames] = await framesSharing(features,
      () => checkSelection(selectionOf(features, iris.getCol('Species')), engines));
    expectArray(edaNames(classifier.engines), inDiscoveryOrder(engines, ['XGBoost', 'SVM', 'Softmax']));
    expect(classifier.best?.name, 'XGBoost');
    expect(classifier.isInteractive, true);
    const {target, features: featureProblems, missingValues, method} = classifier.problems;
    expect([...target, ...featureProblems, ...missingValues, ...method].length, 0);
    expect(classifier.failed.length, 0);
    expectReleased(frames);

    const petalLength = iris.getCol('Petal.Length');
    const three = columnsOf(iris, MEASUREMENTS.filter((name) => name !== 'Petal.Length'));
    const regressor = await checkSelection(selectionOf(three, petalLength), engines);
    expectArray(edaNames(regressor.engines),
      inDiscoveryOrder(engines, ['XGBoost', 'SVM', 'Linear Regression', 'PLS Regression']));
    expect(regressor.best?.name, 'Linear Regression');
    expect(regressor.problems.method.length, 0);

    const softmax = engineByName(engines, 'Softmax');
    const notListed = await checkSelection(engineSelection(softmax, three, petalLength), engines);
    expectArray(notListed.problems.method, [NOT_SOFTMAX]);
    expect(notListed.isInteractive, false);
    const none = await checkSelection(engineSelection(softmax, three, petalLength), [softmax]);
    expectArray(none.problems.method, [NO_METHOD]);
    expect(none.engines.length, 0);
    expect(none.best === undefined, true, 'A method is suggested from an empty list');

    const ruled = await checkSelection(selectionOf([...three, petalLength], petalLength), engines);
    expect(ruled.problems.features.length > 0, true, 'The target as a feature is not reported');
    expect(ruled.engines.length, 0);
    expect(ruled.problems.method.length, 0);
  }, {timeout: TIMEOUT});

  test('checkSelection: SVM does not retrain live on 20000 rows', async () => {
    const engines = EngineRegistry.discover();
    const svm = engineByName(engines, 'SVM');
    const demog = grok.data.demo.demog(20000);
    const features = columnsOf(demog, ['age', 'height', 'weight']);
    const check = await checkSelection(engineSelection(svm, features, demog.getCol('sex')), engines);
    expect(check.engines.some((e) => e.name === 'SVM'), true, 'SVM is not listed');
    expect(check.problems.method.length, 0, check.problems.method.join(' '));
    expect(check.isInteractive, false);
  }, {timeout: TIMEOUT});
});

// The preparation steps and Predict probability, and every method end to end (one category, split for its length).
category('Training', () => {
  test('checkSelection: text and yes/no features need One-hot encoding, dates are refused', async () => {
    const values = Array.from({length: 12}, (_, i) => i + 1);
    const n = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', values);
    const cls = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'cls', values.map((v) => v % 2 === 0 ? 'a' : 'b'));
    const flag = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'flag', values.map((v) => v % 3 === 0));
    const date = DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'date', values.map((v) => new Date(2020, 0, v)));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', values.map((v) => v * 1.5 + (v % 3)));
    let problems = await problemsOf(selectionOf([n, cls, flag], y));
    expectArray(problems.features, ['XGBoost needs numerical features. Check One-hot encoding in Preparation, ' +
      'or uncheck: cls, flag.']);
    problems = await problemsOf({...selectionOf([n, cls, flag], y), steps: ONE_HOT_STEPS});
    expectArray([...problems.features, ...problems.method], []);
    problems = await problemsOf({...selectionOf([n, cls, date], y), steps: ONE_HOT_STEPS});
    expectArray(problems.features, ['XGBoost needs numerical features. Uncheck: date.']);
  }, {timeout: TIMEOUT});

  test('prepareFeatures: one-hot records the categories and builds their columns in that order', async () => {
    const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [0.5, 1.5, 2.5, 3.5]);
    const color = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'color', ['red', 'blue', 'green', 'blue']);
    const flag = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'flag', [true, false, false, true]);
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [1, 2, 3, 4]);
    const prepared = prepareFeatures([color, x, flag], y, ONE_HOT_STEPS, NO_OPTIONS);
    expectArray(prepared.options.preprocessingInfo, ['one-hot']);
    expect(JSON.stringify(prepared.options.oneHotCategories),
      JSON.stringify({color: ['blue', 'green', 'red'], flag: ['false', 'true']}));
    expectArray(names(prepared.columns),
      ['x', 'color=blue', 'color=green', 'color=red', 'flag=false', 'flag=true']);
    expectArray(columnNamed(prepared.columns, 'color=blue').toList(), [0, 1, 0, 1]);
    expectArray(columnNamed(prepared.columns, 'flag=true').toList(), [1, 0, 0, 1]);
    expect(prepared.columns[0] === x && prepared.target === y, true, 'A numerical column or the target was copied');
    expect(color.type === DG.COLUMN_TYPE.STRING && color.length === 4, true, 'The text column changed');
  });

  test('applicableEngines lists the methods for a one-hot encoded demog-like frame, in discovery order', async () => {
    const rows = 40;
    const age = DG.Column.fromList(DG.COLUMN_TYPE.INT, 'age', valuesOf(rows, (i) => 20 + i));
    const height = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'height', valuesOf(rows, (i) => 150 + (i * 7) % 40));
    const weight = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'weight', valuesOf(rows, (i) => 50 + (i * 5) % 45));
    const sex = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'sex', Array.from({length: rows}, (_, i) =>
      i % 3 === 0 ? 'M' : 'F'));
    const race = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'race', Array.from({length: rows}, (_, i) =>
      ['Asian', 'Black', 'Caucasian', 'Other'][i % 4]));
    const engines = EngineRegistry.discover();
    const listed = async (features: DG.Column[], target: DG.Column) => {
      const {columns, target: prepared} = prepareFeatures(features, target, ONE_HOT_STEPS, NO_OPTIONS);
      const {applicable, failed} = await applicableEngines(engines, columns, prepared);
      const edaFailed = failed.map((f) => f.engine.name).filter((name) => EDA_METHODS.includes(name));
      expect(edaFailed.length, 0, edaFailed.join(', '));
      return applicable.map((e) => e.name).filter((name) => EDA_METHODS.includes(name));
    };
    expectArray(await listed([height, weight, sex, race], age),
      inDiscoveryOrder(engines, ['XGBoost', 'SVM', 'Linear Regression', 'PLS Regression']));
    expectArray(await listed([age, height, weight, race], sex),
      inDiscoveryOrder(engines, ['XGBoost', 'SVM', 'Softmax']));
  }, {timeout: TIMEOUT});

  test('replayPreprocessing builds the recorded one-hot columns for a table with a missing and an extra category',
    async () => {
      const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [1, 2, 3]);
      const color = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'color', ['red', 'purple', 'red']);
      const options: PreparationOptions = {...NO_OPTIONS, preprocessingInfo: ['one-hot'],
        oneHotCategories: {color: ['blue', 'green', 'red']}};
      const replayed = replayPreprocessing([x, color], options);
      expectArray(names(replayed), ['x', 'color=blue', 'color=green', 'color=red']);
      expectArray(columnNamed(replayed, 'color=blue').toList(), [0, 0, 0]);
      expectArray(columnNamed(replayed, 'color=green').toList(), [0, 0, 0]);
      expectArray(columnNamed(replayed, 'color=red').toList(), [1, 0, 1]);
      // A model without the record (the built-in tool's) uses the applied data's categories.
      expectArray(names(replayPreprocessing([x, color], {...options, oneHotCategories: undefined})),
        ['x', 'color=purple', 'color=red']);
    });

  test('prepareFeatures: skip unique categories drops the id column before one-hot', async () => {
    const id = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'id', ['s1', 's2', 's3', 's4']);
    const group = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'group', ['a', 'b', 'a', 'b']);
    const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [1, 2, 3, 4]);
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [2, 4, 6, 8]);
    const skipped = prepareFeatures([id, group, x], y, SKIP_UNIQUE_STEPS, NO_OPTIONS);
    expectArray(names(skipped.columns), ['group', 'x']);
    expectArray(skipped.options.preprocessingInfo, ['skip-unique-categories']);
    expectArray(skipped.options.skippedColumns ?? [], ['id']);
    const nothingSkipped = prepareFeatures([group, x], y, SKIP_UNIQUE_STEPS, NO_OPTIONS).options;
    expectArray(nothingSkipped.preprocessingInfo, []);
    expect(nothingSkipped.skippedColumns === undefined, true, 'skippedColumns is recorded without a skipped column');

    const both = prepareFeatures([id, group, x], y, {...SKIP_UNIQUE_STEPS, oneHot: true},
      {...NO_OPTIONS, preprocessingInfo: ['ignore-missing']});
    expectArray(both.options.preprocessingInfo, ['ignore-missing', 'skip-unique-categories', 'one-hot']);
    expectArray(names(both.columns), ['x', 'group=a', 'group=b']);
    expect(JSON.stringify(both.options.oneHotCategories), JSON.stringify({group: ['a', 'b']}));
    // The replay leaves the id out even where its values repeat.
    const repeated = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'id', ['s1', 's1', 's2']);
    const appliedGroup = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'group', ['b', 'a', 'b']);
    const appliedX = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [1, 2, 3]);
    expectArray(names(replayPreprocessing([repeated, appliedGroup, appliedX], both.options)),
      ['x', 'group=a', 'group=b']);
  });

  test('checkSelection: Skip unique categories that leaves no feature is a Features problem', async () => {
    const id = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'id', Array.from({length: 12}, (_, i) => `s${i}`));
    const y = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(12, (i) => i * 1.5));
    const problems = await problemsOf({...selectionOf([id], y),
      steps: {...SKIP_UNIQUE_STEPS, oneHot: true}});
    expectArray(problems.features, ['Skip unique categories leaves no feature. Check more features.']);
    expectArray(problems.method, []);
  }, {timeout: TIMEOUT});

  test('prepareFeatures: predict probability makes a two-class target a float 0/1 target', async () => {
    const x = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [1, 2, 3, 4, 5]);
    const sex = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'sex', ['M', 'F', 'F', 'M', null]);
    const prepared = prepareFeatures([x], sex, {...PREDICT_PROBABILITY, cutoff: 0.4}, NO_OPTIONS);
    expect(prepared.target.type, DG.COLUMN_TYPE.FLOAT);
    expect(prepared.target.name, 'sex');
    expectArray(prepared.target.toList(), [0, 1, 1, 0, null]);
    expect(JSON.stringify(prepared.options), JSON.stringify({preprocessingInfo: [],
      postprocessingInfo: ['binary-classification'], positiveClass: 'F', negativeClass: 'M',
      binaryClassificationThreshold: 0.4, targetType: 'string'}));
    expect(sex.type, DG.COLUMN_TYPE.STRING);

    const three = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'kind', ['a', 'b', 'c', 'a', 'b']);
    const ignored = prepareFeatures([x], three, PREDICT_PROBABILITY, NO_OPTIONS);
    expect(ignored.target === three, true, 'A three-class target was changed');
    expectArray(ignored.options.postprocessingInfo, []);
  });

  test('aucOf: the trapezoid over the rows sorted by score', async () => {
    const labels = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'y', ['a', 'b', 'a', 'b']);
    const scores = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 's', [0.9, 0.8, 0.3, 0.1]);
    expectFloat(aucOf(labels, scores, 'a') ?? NaN, 0.75, 1e-9);
    expectFloat(aucOf(labels, scores, 'b') ?? NaN, 0.25, 1e-9);
    const oneClass = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'y', ['a', 'a', 'a', 'a']);
    expect(aucOf(oneClass, scores, 'a') === undefined, true, 'AUC of one class');
  });

  test('aucOf: tied scores form one diagonal segment', async () => {
    const tiedLabels = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'y', ['a', 'b']);
    const tied = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 's', [0.5, 0.5]);
    expectFloat(aucOf(tiedLabels, tied, 'a') ?? NaN, 0.5, 1e-9);
    expectFloat(aucOf(tiedLabels, tied, 'b') ?? NaN, 0.5, 1e-9);
    // Of the 9 positive-negative pairs, 0.9 beats 3, the tied 0.7 beats 1 and ties 2, 0.2 beats 1: (3 + 2 + 1) / 9.
    const labels = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'y', ['a', 'b', 'a', 'b', 'a', 'b']);
    const mixed = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 's', [0.9, 0.7, 0.7, 0.7, 0.2, 0.1]);
    expectFloat(aucOf(labels, mixed, 'a') ?? NaN, 2 / 3, 1e-9);
  });

  test('recutResult cuts the kept scores at another cutoff without retraining', async () => {
    const actual = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'y', ['a', 'b', 'a', 'b']);
    const scores = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', [0.9, 0.8, 0.3, 0.1]);
    const options: PreparationOptions = {...NO_OPTIONS, postprocessingInfo: ['binary-classification'],
      positiveClass: 'a', negativeClass: 'b', binaryClassificationThreshold: 0.5, targetType: 'string'};
    const result: TrainingResult = {task: 'regression', target: {name: 'y', type: 'string', categories: ['a', 'b']},
      features: {columns: []}, options, splitting: {scheme: 'kfold', folds: 5}, blob: new Uint8Array(),
      metrics: {train: {}, validation: {}}, seed: 1, hyperparameters: {}, rowCount: 4,
      scores: {actual, positiveClass: 'a', train: scores, validation: scores,
        auc: {train: aucOf(actual, scores, 'a'), validation: aucOf(actual, scores, 'a')}}};
    const allPositive = recutResult(result, 0);
    expect(allPositive.options.binaryClassificationThreshold, 0);
    for (const values of [allPositive.metrics.train, allPositive.metrics.validation]) {
      expect(values.sensitivity, 1);
      expect(values.specificity, 0);
      expectFloat(values.auc ?? NaN, 0.75, 1e-9);
    }
    const allNegative = recutResult(result, 1);
    expect(allNegative.metrics.validation.sensitivity, 0);
    expect(allNegative.metrics.validation.specificity, 1);
    expectFloat(allNegative.metrics.validation.auc ?? NaN, 0.75, 1e-9);
    expect(allNegative.metrics.positiveClass, 'a');
    expect(result.options.binaryClassificationThreshold, 0.5);
    const plain: TrainingResult = {...result, scores: undefined};
    expect(recutResult(plain, 0) === plain, true, 'A result without scores was cut');
  });

  test('Predict probability: a regressor on the classes, AUC-ROC, metrics of the labels, non-integer scores',
    async () => {
      const iris = await twoSpeciesIris();
      const features = columnsOf(iris, MEASUREMENTS);
      const species = iris.getCol('Species');
      const engines = EngineRegistry.discover();
      const check = await checkSelection({...selectionOf(features, species), steps: PREDICT_PROBABILITY}, engines);
      const listed = check.engines.map((e) => e.name);
      expect(listed.includes('Linear Regression') && !listed.includes('Softmax'), true, listed.join(', '));

      const request = await requestOf(features, species, PREDICT_PROBABILITY);
      expect(request.prepared.target.type, DG.COLUMN_TYPE.FLOAT);
      expect(request.target === species, true, 'The selection\'s target was replaced');
      const result = await trainModel(request);
      expect(result.task, 'regression');
      expect(JSON.stringify(result.target), JSON.stringify({name: 'Species', type: 'string',
        categories: ['versicolor', 'virginica']}));
      const {options, metrics, scores} = result;
      expect(JSON.stringify([options.positiveClass, options.negativeClass, options.binaryClassificationThreshold,
        options.targetType]), JSON.stringify(['versicolor', 'virginica', 0.5, 'string']));
      expectArray(options.postprocessingInfo, ['binary-classification']);
      expect(metrics.positiveClass, 'versicolor');
      for (const values of [metrics.train, metrics.validation]) {
        const auc = values.auc ?? NaN;
        expect(auc > 0.5 && auc <= 1, true, `AUC-ROC ${auc}`);
        expect(values.accuracy !== undefined && values.sensitivity !== undefined, true, JSON.stringify(values));
        expect(values.mse === undefined, true, 'A regression metric of the scores');
      }
      expect(scores?.validation.toList().some((v) => typeof v === 'number' && !Number.isInteger(v)), true,
        'XGBoost returned 0/1 scores');
      expect(scores?.actual === species, true, 'The scores are measured against another target');
    }, {timeout: TIMEOUT});

  for (const [name, targets] of ENGINE_TARGETS) {
    test(`${name} trains and applies with missing values under Skip rows and Impute`, async () => {
      const calls: string[] = [];
      const gaps: string[] = [];
      const engine = gapChecked(engineByName(EngineRegistry.discover(), name), calls, gaps);
      for (const targetName of targets) {
        for (const missingValues of [SKIP, IMPUTE]) {
          const iris = await openIris();
          const featureNames = MEASUREMENTS.filter((n) => n !== targetName);
          const features = columnsOf(iris, featureNames);
          features[0].set(3, null);
          features[1].set(70, null);
          const expectedRows = missingValues.mode === 'skip' ? 148 : 150;
          const what = `${name}, ${targetName}, ${missingValues.mode}`;

          const request = await prepareTraining(engineSelection(engine, features, iris.getCol(targetName),
            missingValues));
          const result = await trainModel(request);
          expect(result.rowCount, expectedRows, what);
          expect(result.task, targetName === 'Species' ? 'classification' : 'regression', what);

          const prepared = await prepareMissingValues(features, undefined, missingValues);
          const prediction = await apply(engine, prepared.features, result.blob);
          expect(prediction.length, expectedRows, what);
          expect(prediction.stats.missingValueCount, 0, `${what}: empty predictions`);
          expect(features[0].isNone(3) && features[1].isNone(70), true, `${what}: the gaps were filled in place`);
        }
      }
      // Per target and mode: a train and an apply per fold and for the final model, then the application.
      const runs = 2 * targets.length;
      expect(calls.filter((role) => role === 'train').length, 6 * runs, 'The method was called around the spy');
      expect(calls.filter((role) => role === 'apply').length, 7 * runs, 'The method was called around the spy');
      expectArray(gaps, []);
    }, {timeout: ENGINE_TIMEOUT});
  }
});
