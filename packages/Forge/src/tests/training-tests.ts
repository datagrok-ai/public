import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, expect, expectArray, expectExceptionAsync, expectFloat, test}
  from '@datagrok-libraries/test/src/test';
import {defaultHyperparameters, Engine, hyperparametersOf} from '../engines/engine';
import {applicableEngines} from '../engines/applicable-engines';
import {runGetModel} from '../catalog/get-model';
import {apply, EnginePrediction, LoopProgress, yieldToEventLoop} from '../engines/engine-calls';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {defaultFeatures} from '../training/default-features';
import {DEFAULT_VALIDATION, heldOutCount, holdOutFold, kFold, minimumRows} from '../training/k-fold';
import {aucOf, aucOfProbabilities, classesOf, METRIC_DESCRIPTIONS, metricsOf, MetricValues, probabilityMetrics,
  probabilityOf} from '../metrics/metrics';
import {classificationReport} from '../metrics/report';
import {MissingValuesSettings, prepareMissingValues} from '../preparation/missing-values';
import {prepareFeatures, replayPostprocessing, replayPreprocessing} from '../preparation/pipeline';
import {isRecord, PreparationOptions, preparationOptionsOf, stringsOf} from '../preparation/preparation-options';
import {forgeDb} from '../generated/db';
import {deleteModel} from '../storage/model-store';
import {datasetFingerprint} from '../storage/dataset-fingerprint';
import {dataInsights, metricInsights, resultInsights, summaryInsight} from '../training/insights';
import {searchRecordOf} from '../training/search';
import {runTrainModel, TrainModelInput} from '../training/train-api';
import {bestThresholdOf, checkSelection, foldMasks, foldPrediction, needsRowsWithoutGaps, outOfFoldProbabilities,
  prepareTraining, recutResult, TrainingProblems, TrainingRequest, TrainingResult, TrainingSelection, trainingEngines,
  trainModel, trainPredictionOf, validationSplit} from '../training/train-model';
import {QueuedTraining, TrainingQueue} from '../training/training-queue';
import {cancelTaskBarEntry, CLASSIFIERS, columnNamed, columnsOf, deleteModelsNamed, deleteRunsOf, EDA_METHODS,
  engineByName, expectMetrics, expectNoFrames, expectReleased, failsWith, framesSharing, GAP_ROWS, IMPUTE,
  inDiscoveryOrder, IRIS, MEASUREMENTS, names, NO_OPTIONS, NO_PREPARATION_STEPS, ONE_HOT_STEPS, openDemoTable, openIris,
  PREDICT_PROBABILITY, REGRESSORS, requestOf, runsOf, selectionOf, SKIP_UNIQUE_STEPS, taskBarEntry, twoSpeciesIris,
  valuesOf, withGaps} from './test-data';

const TIMEOUT = 60000;
const ENGINE_TIMEOUT = 120000;
const SKIP: MissingValuesSettings = {mode: 'skip'};
/** The targets each EDA method without `mlnativenans` learns on iris: Species (classification), Petal.Length
 * (regression). */
const SKIP_IMPUTE_TARGETS: [string, string[]][] = [['SVM', ['Species', 'Petal.Length']], ['Softmax', ['Species']],
  ['Linear Regression', ['Petal.Length']], ['PLS Regression', ['Petal.Length']],
  ['k-Nearest Neighbors', ['Species', 'Petal.Length']], ['Naive Bayes', ['Species']], ['LDA', ['Species']],
  ['Neural Network', ['Species', 'Petal.Length']], ['RealMLP', ['Species']]];
/** The EDA methods with `mlnativenans` and the iris targets they are checked on. */
const KEEP_TARGETS: [string, string[]][] = [['XGBoost', ['Species', 'Petal.Length']], ['Random Forest', ['Species']],
  ['Extra Trees', ['Species']], ['Hist Gradient Boosting', ['Species']], ['CatBoost', ['Species']],
  ['Decision Tree', ['Species']], ['TabM', ['Species']]];
const NOT_SOFTMAX = 'Softmax cannot learn from this selection. ' +
  'It needs numerical features and a numerical, text or boolean target.';
const NO_METHOD = 'No method can learn from this selection. Check the features and the target.';
const needsRows = (name: string) => needsRowsWithoutGaps(engineByName(EngineRegistry.discover(), name));

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

/** [value] as JSON with every object's keys sorted: a stored json value comes back with its keys reordered. */
function sortedJson(value: unknown): string {
  return JSON.stringify(value, (_, v) => isRecord(v) ?
    Object.fromEntries(Object.keys(v).sort().map((k) => [k, v[k]])) : v);
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

  test('kFold stratifies by the target', async () => {
    const classes = DG.Column.fromStrings('class',
      Array.from({length: 103}, (_, i) => i < 50 ? 'a' : i < 83 ? 'b' : 'c'));
    const folds = kFold(103, 5, 42, classes);
    expectArray(kFold(103, 5, 42, classes), folds);
    const spreads = ['a', 'b', 'c'].map((cls) => {
      const counts = [0, 0, 0, 0, 0];
      for (let r = 0; r < 103; r++) {
        if (classes.get(r) === cls)
          counts[folds[r]]++;
      }
      return Math.max(...counts) - Math.min(...counts);
    });
    expectArray(spreads, [0, 1, 0]);
    const sizes = [0, 0, 0, 0, 0];
    for (const f of folds)
      sizes[f]++;
    expect(Math.max(...sizes) - Math.min(...sizes) <= 1, true, `Fold sizes ${sizes}`);
  });

  test('kFold without shuffle: contiguous folds, classes dealt in row order, the seed unused', async () => {
    const plain = Array.from(kFold(10, 5, 42, undefined, false));
    expectArray(plain, [0, 0, 1, 1, 2, 2, 3, 3, 4, 4]);
    expectArray(Array.from(kFold(10, 5, 7, undefined, false)), plain);
    const classes = DG.Column.fromStrings('class', ['a', 'b', 'a', 'b', 'a', 'a']);
    expectArray(Array.from(kFold(6, 2, 42, classes, false)), [0, 0, 1, 1, 0, 1]);
  });

  test('holdOutFold holds out each class\'s share; without shuffle the last rows', async () => {
    const classes = DG.Column.fromStrings('class', Array.from({length: 100}, (_, i) => i < 60 ? 'a' : 'b'));
    const fold = holdOutFold(100, 0.2, 42, classes);
    expectArray(Array.from(holdOutFold(100, 0.2, 42, classes)), Array.from(fold));
    const held = (cls: string) => Array.from(fold).filter((f, r) => f === 0 && classes.get(r) === cls).length;
    expectArray([held('a'), held('b')], [12, 8]);
    expect(Array.from(fold).filter((f) => f === -1).length, 80);
    const last = Array.from(holdOutFold(10, 0.3, 42, undefined, false));
    expectArray(last, [-1, -1, -1, -1, -1, -1, -1, 0, 0, 0]);
    const perClass = Array.from(holdOutFold(100, 0.2, 1, classes, false));
    expectArray(perClass.flatMap((f, r) => f === 0 ? [r] : []),
      [...Array.from({length: 12}, (_, i) => 48 + i), ...Array.from({length: 8}, (_, i) => 92 + i)]);
    expectArray([minimumRows({...DEFAULT_VALIDATION, scheme: 'holdout', testFraction: 0.2}),
      minimumRows(DEFAULT_VALIDATION), heldOutCount(150, 0.2)], [10, 10, 30]);
  });

  test('trainModel with a hold out measures Validation on the held out rows, then fits all rows', async () => {
    const request = {...await irisRequest(MEASUREMENTS, 'Species'),
      validation: {...DEFAULT_VALIDATION, scheme: 'holdout' as const}};
    const result = await trainModel(request);
    expect(sortedJson(result.splitting), sortedJson({scheme: 'holdout', trainFraction: 0.8, isStratified: true,
      isShuffled: true}));
    expect(result.rowCount, 150);
    const fold = validationSplit(request);
    expect(Array.from(fold).filter((f) => f === 0).length, 30);
    const {columns, target} = request.prepared;
    const held = DG.BitSet.create(150, (r) => fold[r] === 0);
    const fit = await foldPrediction(request.engine, columns, target, request.hyperparameters, foldMasks(fold, 0));
    expectFloat(result.metrics.validation.accuracy ?? NaN,
      metricsOf('classification', target.clone(held), fit.prediction).accuracy ?? NaN, 1e-9);
    expect(result.prediction.length, 150);
    expect(typeof result.metrics.validation.auc, 'number');
  }, {timeout: TIMEOUT});

  test('the validation choice: K-fold unstratified, unshuffled, and Predict probability on a hold out', async () => {
    const request = await irisRequest(MEASUREMENTS, 'Species');
    const plain = await trainModel({...request, validation: {...DEFAULT_VALIDATION, stratified: false}});
    expect(sortedJson(plain.splitting), sortedJson({scheme: 'kfold', folds: 5, isStratified: false,
      isShuffled: true}));
    expectArray(Array.from(validationSplit({...request, validation: {...DEFAULT_VALIDATION, stratified: false}})),
      Array.from(kFold(150, 5, request.seed)));
    const ordered = await trainModel({...request, validation: {...DEFAULT_VALIDATION, shuffle: false}});
    expect(ordered.splitting.isShuffled, false);
    const otherSeed = await trainModel({...request, seed: 7, validation: {...DEFAULT_VALIDATION, shuffle: false}});
    expectMetrics(otherSeed.metrics.validation, ordered.metrics.validation, 1e-9);

    const iris = await twoSpeciesIris();
    const probability = await trainModel({...await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'),
      PREDICT_PROBABILITY), validation: {...DEFAULT_VALIDATION, scheme: 'holdout'}});
    expect(probability.splitting.isStratified, true);
    expect(probability.outOfFold.probabilities?.columns[0].stats.missingValueCount, 80);
    const recut = recutResult(probability, {predictProbability: true, threshold: 0.3});
    expect(typeof recut.metrics.validation.accuracy, 'number');
    expect(typeof recut.metrics.validation.auc, 'number');
  }, {timeout: TIMEOUT});

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
    expect(JSON.stringify(result.splitting), JSON.stringify({scheme: 'kfold', folds: 5, isStratified: true,
      isShuffled: true}));
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
    const edaNames = (list: Engine[]) => list.map((e) => e.name).filter((name) => EDA_METHODS.includes(name));
    const iris = await openIris();
    const features = columnsOf(iris, MEASUREMENTS);
    const [classifier, frames] = await framesSharing(features,
      () => checkSelection(selectionOf(features, iris.getCol('Species')), engines));
    expectArray(edaNames(classifier.engines), inDiscoveryOrder(engines, CLASSIFIERS));
    expect(classifier.best?.name, 'XGBoost');
    expect(classifier.isInteractive, true);
    const {target, features: featureProblems, missingValues, method} = classifier.problems;
    expect([...target, ...featureProblems, ...missingValues, ...method].length, 0);
    expect(classifier.failed.length, 0);
    expectReleased(frames);

    const twoSpecies = await twoSpeciesIris();
    const base = selectionOf(columnsOf(twoSpecies, MEASUREMENTS), twoSpecies.getCol('Species'));
    const probability = await checkSelection({...base, steps: {...base.steps, predictProbability: true}}, engines);
    expect(probability.best?.name, 'XGBoost', 'Predict probability suggests another method');
    expectArray(edaNames(probability.engines), inDiscoveryOrder(engines, CLASSIFIERS));

    const petalLength = iris.getCol('Petal.Length');
    const three = columnsOf(iris, MEASUREMENTS.filter((name) => name !== 'Petal.Length'));
    const regressor = await checkSelection(selectionOf(three, petalLength), engines);
    expectArray(edaNames(regressor.engines), inDiscoveryOrder(engines, REGRESSORS));
    expect(regressor.best?.name, 'Linear Regression');
    expect(regressor.problems.method.length, 0);

    const demog = grok.data.demo.demog(100);
    const encoded = selectionOf(columnsOf(demog, ['HEIGHT', 'SEX', 'RACE']), demog.getCol('AGE'));
    const oneHot = await checkSelection({...encoded, steps: {...encoded.steps, oneHot: true}}, engines);
    expect(oneHot.best?.name, 'PLS Regression', 'The one-hot columns are not counted as numerical features');

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
    expectArray(await listed([height, weight, sex, race], age), inDiscoveryOrder(engines, REGRESSORS));
    expectArray(await listed([age, height, weight, race], sex), inDiscoveryOrder(engines, CLASSIFIERS));
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

  test('recutResult re-cuts the kept class probabilities without retraining', async () => {
    const actual = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'y', ['a', 'b', 'a', 'b']);
    const prediction = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Prediction', ['a', 'b', 'b', 'b']);
    const probabilities = {classes: ['a', 'b'], columns: [
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'proba_a', [0.9, 0.8, 0.3, 0.1]),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'proba_b', [0.1, 0.2, 0.7, 0.9])]};
    const scoreValues = {auc: 0.75, logLoss: 0.6};
    const result: TrainingResult = {task: 'classification',
      target: {name: 'y', type: 'string', categories: ['a', 'b']}, features: {columns: []}, options: NO_OPTIONS,
      splitting: {scheme: 'kfold', folds: 5}, blob: new Uint8Array(), prediction, probabilities,
      metrics: {train: scoreValues, validation: scoreValues, positiveClass: 'a'}, seed: 1, hyperparameters: {},
      rowCount: 4, actual, outOfFold: {prediction, probabilities}};
    const allPositive = recutResult(result, {predictProbability: true, threshold: 0});
    expect(JSON.stringify(allPositive.options), JSON.stringify({preprocessingInfo: [],
      postprocessingInfo: ['class-probabilities'], positiveClass: 'a', negativeClass: 'b',
      binaryClassificationThreshold: 0}));
    expectArray(trainPredictionOf(allPositive).toList(), ['a', 'a', 'a', 'a']);
    for (const values of [allPositive.metrics.train, allPositive.metrics.validation]) {
      expect(values.sensitivity, 1);
      expect(values.specificity, 0);
      expectFloat(values.auc ?? NaN, 0.75, 1e-9);
      expectFloat(values.logLoss ?? NaN, 0.6, 1e-9);
    }
    const allNegative = recutResult(result, {predictProbability: true, threshold: 1});
    expectArray(trainPredictionOf(allNegative).toList(), ['b', 'b', 'b', 'b']);
    expect(allNegative.metrics.validation.sensitivity, 0);
    expect(allNegative.metrics.validation.specificity, 1);
    expectFloat(allNegative.metrics.validation.auc ?? NaN, 0.75, 1e-9);
    expect(allNegative.metrics.positiveClass, 'a');
    const off = recutResult(allNegative, {predictProbability: false, threshold: 1});
    expect(JSON.stringify(off.options), JSON.stringify(NO_OPTIONS));
    expectFloat(off.metrics.validation.accuracy ?? NaN, 0.75, 1e-9);

    const optimized = recutResult(result, {predictProbability: true, threshold: 0, thresholdMetric: 'f1'});
    expect(optimized.options.binaryClassificationThreshold, 0.2);
    expect(optimized.options.thresholdMetric, 'f1');
    expectArray(trainPredictionOf(optimized).toList(), ['a', 'a', 'a', 'b']);
    const noBest = recutResult({...result, outOfFold: {prediction}},
      {predictProbability: true, threshold: 0.4, thresholdMetric: 'f1'});
    expect(noBest.options.binaryClassificationThreshold, 0.4, 'No best keeps the threshold');
    expect(noBest.options.thresholdMetric === undefined, true, 'No best records the metric');
    expect(recutResult(optimized, {predictProbability: true, threshold: 0.3}).options.thresholdMetric === undefined,
      true, 'A threshold without a metric keeps the metric');

    const plain: TrainingResult = {...result, probabilities: undefined};
    expect(recutResult(plain, {predictProbability: true, threshold: 0}) === plain, true, 'Cut without probabilities');
    const stale: TrainingResult = {...plain, options: allPositive.options};
    expect(JSON.stringify(recutResult(stale, {predictProbability: true, threshold: 0}).options),
      JSON.stringify(NO_OPTIONS));
  }, {timeout: TIMEOUT});

  test('Predict probability on a boolean target: the first category positive, boolean labels', async () => {
    const iris = await twoSpeciesIris();
    const species = iris.getCol('Species');
    const target = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'isVersicolor',
      Array.from({length: iris.rowCount}, (_, i) => species.get(i) === 'versicolor'));
    const result = await trainModel(await requestOf(columnsOf(iris, MEASUREMENTS), target,
      {...PREDICT_PROBABILITY, threshold: 0.7}));
    const [positive, negative] = classesOf(target);
    expect(target.categories[0], positive);
    expect(result.metrics.positiveClass, positive);
    const {options} = result;
    expect(JSON.stringify([options.postprocessingInfo, options.positiveClass, options.negativeClass,
      options.binaryClassificationThreshold]), JSON.stringify([['class-probabilities'], positive, negative, 0.7]));
    const probability = probabilityOf(result.probabilities!, positive);
    expect(probability?.name, `proba_${positive}`);
    const labels = trainPredictionOf(result);
    expect(labels.type, DG.COLUMN_TYPE.BOOL);
    const expected = Array.from({length: labels.length}, (_, i) =>
      String(probability!.getNumber(i) >= Math.fround(0.7) ? positive : negative) === 'true');
    expectArray(labels.toList(), expected);
    expect(expected.includes(true) && expected.includes(false), true, 'One class at 0.7');
    const recut = trainPredictionOf(recutResult(result, {predictProbability: true, threshold: 0}));
    expect(recut.toList().every((v) => String(v) === positive), true, 'A threshold of 0 leaves a negative label');
  }, {timeout: TIMEOUT});

  test('Predict probability keeps the classifier and cuts the positive class\'s probability', async () => {
    const iris = await twoSpeciesIris();
    const features = columnsOf(iris, MEASUREMENTS);
    const species = iris.getCol('Species');
    const request = await requestOf(features, species, {...PREDICT_PROBABILITY, threshold: 0.9});
    expect(request.prepared.target === species, true, 'The prepared target is not the selection\'s');
    const result = await trainModel(request);
    expect(result.task, 'classification');
    const {options, metrics} = result;
    expectArray(options.postprocessingInfo, ['class-probabilities']);
    expect(JSON.stringify([options.positiveClass, options.negativeClass, options.binaryClassificationThreshold]),
      JSON.stringify(['versicolor', 'virginica', 0.9]));
    expect(metrics.positiveClass, 'versicolor');
    const labels = trainPredictionOf(result);
    expectArray(labels.toList(), replayPostprocessing(result.prediction, options, result.probabilities).toList());
    expect(labels.toList().some((v, i) => v !== result.prediction.get(i)), true,
      'The threshold 0.9 gives the method\'s own labels');
    for (const values of [metrics.train, metrics.validation]) {
      for (const id of ['accuracy', 'sensitivity', 'auc', 'logLoss', 'brier'] as const)
        expect(typeof values[id], 'number', id);
    }
    expect(result.actual === species, true, 'The predictions are measured against another target');

    const full = await openIris();
    const three = await trainModel(await requestOf(columnsOf(full, MEASUREMENTS), full.getCol('Species'),
      PREDICT_PROBABILITY));
    expectArray(three.options.postprocessingInfo, ['class-probabilities']);
    expect([three.options.positiveClass, three.options.negativeClass, three.options.binaryClassificationThreshold]
      .every((v) => v === undefined), true, JSON.stringify(three.options));
    expectArray(trainPredictionOf(three).toList(), three.prediction.toList());
  }, {timeout: TIMEOUT});

  test('dataInsights: unbalanced classes, correlated columns, almost unique values', async () => {
    const rows = 10;
    const float = (name: string, value: (i: number) => number) =>
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, valuesOf(rows, value));
    const text = (name: string, value: (i: number) => string) =>
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, Array.from({length: rows}, (_, i) => value(i)));
    const x = float('x', (i) => i);
    const z = float('z', (i) => (i * 7) % 5);
    const balanced = text('t', (i) => i % 2 === 0 ? 'a' : 'b');

    const unbalanced = dataInsights([text('f', (i) => i < 1 ? 'c' : 'd')], text('t', (i) => i < 8 ? 'a' : 'b'));
    expect(JSON.stringify(unbalanced), JSON.stringify([
      {text: 'The classes are unbalanced.', detail: 'a: 80.0%\nb: 20.0%', panel: 'classBalance'},
      {text: 'The classes of f are unbalanced.', detail: 'c: 10.0%\nd: 90.0%'}]));
    // 1 : 1 : 4 is 16.67% : 16.67% : 66.67%; the two largest remainders round up so the shares sum to 100%.
    const sixths = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 't', ['a', 'b', 'c', 'c', 'c', 'c']);
    const shares = dataInsights([], sixths)[0].detail ?? '';
    expect(shares, 'a: 16.7%\nb: 16.7%\nc: 66.6%');
    const percents = shares.split('\n').map((line) => Number(line.slice(line.indexOf(': ') + 2, -1)));
    expectFloat(percents.reduce((sum, p) => sum + p, 0), 100, 1e-9);
    const wine = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'type',
      Array.from({length: 6497}, (_, i) => i < 1599 ? 'red' : 'white'));
    expect(dataInsights([], wine)[0].detail, 'red: 24.6%\nwhite: 75.4%');
    // One short line for every correlated pair, the pairs in its list.
    const correlated = dataInsights([x, float('y', (i) => 2 * i), float('w', (i) => -i), float('v', (i) => 3 * i)],
      balanced);
    expectArray(correlated.map((i) => i.text), ['Features are highly correlated.']);
    expectArray(correlated[0].items ?? [], ['x and y: 1.00', 'x and v: 1.00', 'y and v: 1.00']);
    expect(correlated[0].detail === undefined, true, 'The correlation has a detail line');
    expect(correlated[0].panel, 'correlation');
    expect(dataInsights([x], float('y', (i) => 2 * i))[0].panel === undefined, true, 'A pair with the target points ' +
      'to Features correlation');
    const unique = dataInsights([z], text('t', (i) => `v${Math.min(i, 8)}`));
    const uniqueInsight = unique.find((i) =>
      i.text === 'Almost every value of t is different; it will not help the model.');
    expect(uniqueInsight?.detail, '90% of its values are unique (above 80%).');
    expectArray(dataInsights([x, z], balanced), []);
  });

  test('checkSelection gives the data insights of iris', async () => {
    const iris = await openIris();
    const check = await checkSelection(selectionOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species')),
      EngineRegistry.discover());
    const correlated = check.insights.find((i) => i.text === 'Features are highly correlated.');
    expect(correlated?.items?.includes('Petal.Length and Petal.Width: 0.96'), true, JSON.stringify(check.insights));
    const failed = await checkSelection(selectionOf([], iris.getCol('Species')), EngineRegistry.discover());
    expectArray(failed.insights, []);
  }, {timeout: TIMEOUT});

  test('resultInsights: a classifier\'s accuracy and class balance, a regressor\'s residuals and statistics',
    async () => {
      const rows = 20;
      const labels = (value: (i: number) => string) =>
        DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'c', Array.from({length: rows}, (_, i) => value(i)));
      const classes = labels((i) => i % 2 === 0 ? 'a' : 'b');
      const texts = (values: MetricValues, positiveClass?: string) =>
        resultInsights({values, target: classes, prediction: classes, positiveClass}).map((i) => i.text);
      expectArray(texts({accuracy: 0.1}), ['Only 10% of the rows are predicted correctly.']);
      expectArray(texts({sensitivity: 0.9, specificity: 0.3}, 'a'), ['The model finds a much better than b.']);
      expectArray(texts({sensitivity: 0.3, specificity: 0.9}, 'a'), ['The model finds b much better than a.']);
      expectArray(texts({accuracy: 1, sensitivity: 1, specificity: 1}, 'a'), []);

      const target = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'y', valuesOf(rows, (i) => i));
      const predicted = (value: (i: number) => number) =>
        DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'p', valuesOf(rows, value));
      const near = (i: number) => i + (i % 2 === 0 ? 0.1 : -0.1);
      const insightsOf = (prediction: DG.Column) => resultInsights({values: {}, target, prediction});
      const regression = (prediction: DG.Column) => insightsOf(prediction).map((i) => i.text);
      expect(regression(predicted((i) => i === 5 ? 100 : near(i))).includes(
        'Large errors are mostly on one side: the model over-predicts.'), true);
      expect(regression(predicted((i) => i === 5 ? -100 : near(i))).includes(
        'Large errors are mostly on one side: the model under-predicts.'), true);
      const constant = insightsOf(predicted(() => 9.5));
      const deviations = constant.filter((i) => i.panel === 'statistics');
      expectArray(deviations.map((i) => i.text),
        ['The predictions\' statistics deviate significantly from the actual ones.']);
      const items = deviations[0].items ?? [];
      expect(items.length > 1 && items.some((item) => item.startsWith('variance: predicted 0.000, actual ')), true,
        items.join(' | '));
      expectArray(regression(target), []);
    });

  test('metricInsights: a value below its threshold as is, the overfitting metrics in one insight', async () => {
    const insights = metricInsights({train: {r2: 0.2, mse: 1, mae: 1}, validation: {r2: 0.9, mse: 4, mae: 1.05}});
    expectArray(insights.map((i) => i.text), ['R2 is below 0.3.', 'The model may be overfitting.']);
    expect(insights[0].detail, 'Train 0.200, validation 0.900.');
    expectArray(insights[1].items ?? [], ['MSE: train 1.000, validation 4.000', 'R2: train 0.200, validation 0.900']);
    expect(insights.every((i) => i.panel === 'metrics'), true, JSON.stringify(insights));
    expectArray(metricInsights({train: {accuracy: 0.97}, validation: {accuracy: 0.95}}), []);
  });

  test('summaryInsight: Balanced accuracy or R2 on the evaluated, validation or training rows', async () => {
    const classifier = summaryInsight({train: {balancedAccuracy: 0.9}, validation: {balancedAccuracy: 0.85}});
    expect(classifier?.text,
      'Balanced accuracy 0.850 on the validation rows: on average 85% of each class is predicted correctly.');
    expect(classifier?.panel, 'metrics');
    expect(classifier?.detail, METRIC_DESCRIPTIONS.balancedAccuracy);
    const regressor = summaryInsight({train: {r2: 0.8}, validation: {}});
    expect(regressor?.text, 'R2 0.800 on the training rows: the model explains 80% of the target\'s variation.');
    expect(regressor?.detail, METRIC_DESCRIPTIONS.r2);
    expect(summaryInsight({train: {r2: 0.5}, validation: {r2: 0.5}}, {r2: -0.2})?.text,
      'R2 -0.200 on the evaluated rows: the model predicts worse than the target\'s mean.');
    expect(summaryInsight({train: {balancedAccuracy: 0.9}, validation: {balancedAccuracy: 0.85}},
      {balancedAccuracy: 0.7})?.text,
    'Balanced accuracy 0.700 on the evaluated rows: on average 70% of each class is predicted correctly.');
    expect(summaryInsight({train: {}, validation: {}}) === undefined, true, 'A line without either metric');
  });

  test('trainModel keeps the prediction of the training rows; prepareTraining the kept rows', async () => {
    const iris = await openIris();
    const result = await trainModel(await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species')));
    expect(result.prediction.length, 150);
    const labels = trainPredictionOf(result);
    expect(labels === result.prediction && labels.type === DG.COLUMN_TYPE.STRING, true, labels.type);

    const twoSpecies = await twoSpeciesIris();
    const probability = await trainModel(await requestOf(columnsOf(twoSpecies, MEASUREMENTS),
      twoSpecies.getCol('Species'), PREDICT_PROBABILITY));
    expect(probability.probabilities !== undefined, true, 'The probability model kept no probabilities');
    expect(trainPredictionOf(probability).toList().every((v) => v === 'versicolor' || v === 'virginica'), true,
      trainPredictionOf(probability).categories.join());
    expectArray(trainPredictionOf(probability).categories, ['versicolor', 'virginica']);

    expect((await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'))).keptRows === null, true,
      'Rows were skipped without a missing value');
    iris.getCol('Sepal.Width').set(3, null);
    const skipped = await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'));
    expect(skipped.keptRows?.trueCount, 149);
    expect(skipped.keptRows?.get(3), false);
  }, {timeout: TIMEOUT});

  test('trainModel measures AUC-ROC from the class probabilities: three classes macro, two classes the positive class',
    async () => {
      const iris = await openIris();
      const species = iris.getCol('Species');
      const result = await trainModel(await requestOf(columnsOf(iris, MEASUREMENTS), species));
      const probabilities = result.probabilities;
      expectArray(probabilities?.classes ?? [], ['setosa', 'versicolor', 'virginica']);
      for (const col of probabilities?.columns ?? [])
        expect(`${col.type} ${col.length}`, `${DG.COLUMN_TYPE.FLOAT} 150`, col.name);
      for (const split of ['train', 'validation'] as const) {
        const auc = result.metrics[split].auc ?? NaN;
        expect(auc > 0.5 && auc <= 1, true, `${split} AUC-ROC ${auc}`);
      }
      expectFloat(result.metrics.train.auc ?? NaN, aucOfProbabilities(species, probabilities!) ?? NaN, 1e-9);
      expect(result.task, 'classification');
      expect(result.outOfFold.probabilities !== undefined, true, 'No out-of-fold probabilities');
      const {logLoss, brier} = probabilityMetrics(species, probabilities!);
      expectFloat(result.metrics.train.logLoss ?? NaN, logLoss ?? NaN, 1e-9);
      expectFloat(result.metrics.train.brier ?? NaN, brier ?? NaN, 1e-9);
      expect(typeof result.metrics.train.mcc, 'number');
      for (const id of ['balancedAccuracy', 'averagePrecision'] as const) {
        const value = result.metrics.validation[id] ?? NaN;
        expect(value > 0 && value <= 1, true, `validation ${id} ${value}`);
      }
      expectArray(result.metrics.classes?.map((c) => c.label) ?? [], ['setosa', 'versicolor', 'virginica']);
      for (const c of result.metrics.classes ?? [])
        expect(`${c.rows} ${typeof c.train.auc} ${typeof c.validation.sensitivity}`, '50 number number', c.label);
      expect(result.curves === undefined, true, 'Three classes have curves');

      const twoSpecies = await twoSpeciesIris();
      const two = twoSpecies.getCol('Species');
      const binary = await trainModel(await requestOf(columnsOf(twoSpecies, MEASUREMENTS), two));
      expect(binary.metrics.positiveClass, 'versicolor');
      const positive = probabilityOf(binary.probabilities!, 'versicolor')!;
      expectFloat(binary.metrics.train.auc ?? NaN, aucOf(two, positive, 'versicolor') ?? NaN, 1e-9);
      const {roc, pr} = binary.curves!;
      expectArray([roc.threshold[0], roc.fpr.length, pr.threshold[pr.threshold.length - 1],
        pr.recall[pr.recall.length - 1]], [Infinity, roc.tpr.length, Infinity, 0]);
    }, {timeout: ENGINE_TIMEOUT});

  test('a regression gets no probability metrics', async () => {
    const iris = await openIris();
    const regression = await trainModel(await requestOf(columnsOf(iris, ['Sepal.Length', 'Sepal.Width']),
      iris.getCol('Petal.Length')));
    expect(regression.probabilities === undefined && regression.outOfFold.probabilities === undefined, true,
      'A regression has class probabilities');
    for (const values of [regression.metrics.train, regression.metrics.validation]) {
      for (const id of ['auc', 'logLoss', 'brier', 'mcc'] as const)
        expect(values[id] === undefined, true, `A regression has ${id}`);
    }
    expect(regression.curves === undefined, true, 'A regression has curves');

    const twoSpecies = await twoSpeciesIris();
    const probability = await trainModel(await requestOf(columnsOf(twoSpecies, MEASUREMENTS),
      twoSpecies.getCol('Species'), PREDICT_PROBABILITY));
    expect(probability.probabilities !== undefined && probability.outOfFold.probabilities !== undefined, true,
      'Predict probability has no class probabilities');
    expect(typeof probability.metrics.validation.logLoss, 'number');
    expect(probability.curves !== undefined, true, 'Predict probability has no curves');
  }, {timeout: ENGINE_TIMEOUT});

  test('outOfFoldProbabilities gathers by class name and gives 0 to a class a fold never saw', async () => {
    const floats = (name: string, values: number[]) => DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values);
    const prediction = DG.Column.fromStrings('Prediction', ['a', 'b']);
    const fold0: EnginePrediction = {prediction, probabilities: {classes: ['a', 'b'],
      columns: [floats('proba_a', [0.9, 0.2]), floats('proba_b', [0.1, 0.8])]}};
    const fold1: EnginePrediction = {prediction, probabilities: {classes: ['a', 'b', 'c'],
      columns: [floats('proba_a', [0.5, 0.1]), floats('proba_b', [0.3, 0.2]), floats('proba_c', [0.2, 0.7])]}};
    const fold = [0, 1, 0, 1];
    const position = [0, 0, 1, 1];
    const gathered = outOfFoldProbabilities([fold0, fold1], fold1.probabilities!, fold, position)!;
    expectArray(gathered.columns.map((c) => c.name), ['proba_a', 'proba_b', 'proba_c']);
    const rounded = (col: DG.Column) => col.toList().map((v: number) => Math.round(v * 100) / 100);
    expectArray(rounded(gathered.columns[0]), [0.9, 0.5, 0.2, 0.1]);
    expectArray(rounded(gathered.columns[1]), [0.1, 0.3, 0.8, 0.2]);
    expectArray(rounded(gathered.columns[2]), [0, 0.2, 0, 0.7]);
    expect(outOfFoldProbabilities([fold0, {prediction}], fold0.probabilities!, fold, position) === undefined, true,
      'A fold without probabilities gave probabilities');
  });

  test('a boolean target gets AUC-ROC of its positive class', async () => {
    const demog = await openDemoTable('System:DemoFiles/demog.csv', 'demog');
    const control = demog.getCol('CONTROL');
    const result = await trainModel(await requestOf(columnsOf(demog, ['AGE', 'HEIGHT', 'WEIGHT']), control,
      NO_PREPARATION_STEPS, undefined, {mode: 'keep'}));
    expectArray(result.probabilities?.classes ?? [], ['false', 'true']);
    expect(result.metrics.positiveClass, control.categories[0]);
    expect(probabilityOf(result.probabilities!, result.metrics.positiveClass!) !== undefined, true,
      `No probability column of ${result.metrics.positiveClass}`);
    expect(typeof result.metrics.validation.auc, 'number');
    expectArray(result.metrics.classes?.map((c) => c.label) ?? [], ['false', 'true']);
  }, {timeout: ENGINE_TIMEOUT});

  for (const [name, targets] of SKIP_IMPUTE_TARGETS) {
    test(`${name} trains and applies with missing values under Skip rows and Impute`, async () => {
      const calls: string[] = [];
      const gaps: string[] = [];
      const engine = gapChecked(engineByName(EngineRegistry.discover(), name), calls, gaps);
      for (const targetName of targets) {
        for (const missingValues of [SKIP, IMPUTE]) {
          const iris = await openIris();
          const featureNames = MEASUREMENTS.filter((n) => n !== targetName);
          const features = withGaps(columnsOf(iris, featureNames));
          const expectedRows = missingValues.mode === 'skip' ? 148 : 150;
          const what = `${name}, ${targetName}, ${missingValues.mode}`;

          const request = await prepareTraining(engineSelection(engine, features, iris.getCol(targetName),
            missingValues));
          const result = await trainModel(request);
          expect(result.rowCount, expectedRows, what);
          expect(result.task, targetName === 'Species' ? 'classification' : 'regression', what);
          if (targetName === 'Species') {
            expect(name === 'SVM' ? result.metrics.validation.auc === undefined :
              typeof result.metrics.validation.auc === 'number', true, `${what}: AUC-ROC`);
            if (name === 'SVM') {
              expect(result.metrics.validation.averagePrecision === undefined &&
                typeof result.metrics.validation.mcc === 'number', true, `${what}: average precision or MCC`);
            }
          }

          const prepared = await prepareMissingValues(features, undefined, missingValues);
          const prediction = (await apply(engine, prepared.features, result.blob)).prediction;
          expect(prediction.length, expectedRows, what);
          expect(prediction.stats.missingValueCount, 0, `${what}: empty predictions`);
          expect(features[0].isNone(GAP_ROWS[0]) && features[1].isNone(GAP_ROWS[1]), true,
            `${what}: the gaps were filled in place`);
        }
      }
      // Per target and mode: a train and an apply per fold and for the final model, then the application.
      const runs = 2 * targets.length;
      expect(calls.filter((role) => role === 'train').length, 6 * runs, 'The method was called around the spy');
      expect(calls.filter((role) => role === 'apply').length, 7 * runs, 'The method was called around the spy');
      expectArray(gaps, []);
    }, {timeout: ENGINE_TIMEOUT});
  }

  for (const [name, targets] of KEEP_TARGETS) {
    test(`${name} trains and applies on the rows with missing values under keep`, async () => {
      const calls: string[] = [];
      const gaps: string[] = [];
      const engine = gapChecked(engineByName(EngineRegistry.discover(), name), calls, gaps);
      for (const targetName of targets) {
        const iris = await openIris();
        const features = withGaps(columnsOf(iris, MEASUREMENTS.filter((n) => n !== targetName)));
        const what = `${name}, ${targetName}`;

        const request = await prepareTraining(engineSelection(engine, features, iris.getCol(targetName),
          {mode: 'keep'}));
        expect(request.keptRows === null, true, `${what}: rows were skipped`);
        const result = await trainModel(request);
        expect(result.rowCount, 150, what);
        expect(JSON.stringify(result.options.missingValues), JSON.stringify({mode: 'keep', skippedRows: 0}), what);
        expectArray(result.options.preprocessingInfo, []);
        if (targetName === 'Species') {
          expect(typeof result.metrics.validation.auc, 'number', `${what}: no AUC-ROC`);
          expect(typeof result.metrics.validation.balancedAccuracy, 'number', `${what}: no balanced accuracy`);
          expect(Array.isArray(result.metrics.classes), true, `${what}: no classes`);
        }

        const prediction = (await apply(engine, features, result.blob)).prediction;
        expect(prediction.length, 150, what);
        expect(prediction.stats.missingValueCount, 0, `${what}: empty predictions`);
      }
      expect(gaps.some((g) => g === 'train: Sepal.Length' || g === 'train: Sepal.Width'), true,
        `The gaps did not reach the method: ${gaps.join(', ')}`);
    }, {timeout: ENGINE_TIMEOUT});
  }

  test('an explicit skip or impute is honoured for a method that handles missing values', async () => {
    const calls: string[] = [];
    const gaps: string[] = [];
    const engines = EngineRegistry.discover();
    const xgboost = gapChecked(engineByName(engines, 'XGBoost'), calls, gaps);
    for (const [missingValues, rows] of [[SKIP, 148], [IMPUTE, 150]] as const) {
      const iris = await openIris();
      const features = withGaps(columnsOf(iris, MEASUREMENTS.filter((n) => n !== 'Species')));
      const request = await prepareTraining(engineSelection(xgboost, features, iris.getCol('Species'), missingValues));
      expect((await trainModel(request)).rowCount, rows, missingValues.mode);
    }
    expectArray(gaps, []);

    const iris = await openIris();
    const features = withGaps(columnsOf(iris, MEASUREMENTS.filter((n) => n !== 'Species')));
    const check = await checkSelection(engineSelection(engineByName(engines, 'SVM'), features, iris.getCol('Species'),
      {mode: 'keep'}), engines);
    expect(check.problems.missingValues[0], needsRows('SVM'));
  }, {timeout: ENGINE_TIMEOUT});
});

// The Crux metrics (EDA's classificationReport) next to the TS ones.
category('Training', () => {
  test('Crux report parity: the TS metrics equal the report\'s', async () => {
    const twoSpecies = await twoSpeciesIris();
    const binaryRequest = await requestOf(columnsOf(twoSpecies, MEASUREMENTS), twoSpecies.getCol('Species'));
    const binary = await trainModel(binaryRequest);
    const report = (await classificationReport(binaryRequest.target, binary.prediction, binary.metrics.positiveClass,
      binary.probabilities, true))!;
    const train = binary.metrics.train;
    const pairs: [string, number | undefined, number | undefined][] = [['accuracy', report.accuracy, train.accuracy],
      ['precision', report.precision, train.precision], ['recall', report.recall, train.sensitivity],
      ['f1', report.f1, train.f1], ['auc', report.auc, train.auc],
      ['balancedAccuracy', report.balancedAccuracy, train.balancedAccuracy], ['mcc', report.mcc, train.mcc]];
    for (const [id, crux, ts] of pairs)
      expectFloat(crux ?? NaN, ts ?? NaN, 1e-9, id);
    const diagonal = report.confusion.reduce((sum, row, i) => sum + row[i], 0);
    expectFloat(diagonal / report.n, report.accuracy, 1e-9, 'confusion diagonal');

    const iris = await openIris();
    const request = await requestOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'));
    const multi = await trainModel(request);
    const multiReport = (await classificationReport(request.target, multi.prediction, undefined, multi.probabilities))!;
    expectFloat(multiReport.f1, multi.metrics.train.f1 ?? NaN, 1e-9, 'macro f1');
    expectFloat(multiReport.auc ?? NaN, multi.metrics.train.auc ?? NaN, 1e-9, 'macro auc');
    expectFloat(multiReport.balancedAccuracy, multi.metrics.train.balancedAccuracy ?? NaN, 1e-9, 'balanced accuracy');
    expectFloat(multiReport.mcc, multi.metrics.train.mcc ?? NaN, 1e-9, 'multiclass mcc');
    const aucs = multiReport.classes.map((c) => c.auc ?? NaN);
    expectFloat(aucs.reduce((s, a) => s + a, 0) / aucs.length, multiReport.auc ?? NaN, 1e-9, 'mean of the class AUCs');
  }, {timeout: ENGINE_TIMEOUT});

  test('a two-class text target with empty values is binary: the positive class\'s AUC-ROC and curves', async () => {
    const twoSpecies = await twoSpeciesIris();
    const species = twoSpecies.getCol('Species');
    species.set(0, null);
    species.set(60, null);
    const request = await requestOf(columnsOf(twoSpecies, MEASUREMENTS), species);
    const result = await trainModel(request);
    expect(result.metrics.positiveClass, 'versicolor');
    expect(result.curves !== undefined, true, 'No curves');
    const positive = probabilityOf(result.probabilities!, 'versicolor')!;
    expectFloat(result.metrics.train.auc ?? NaN, aucOf(request.target, positive, 'versicolor') ?? NaN, 1e-9);
  }, {timeout: ENGINE_TIMEOUT});

  test('recutResult re-cuts the label metrics, keeps AUC-ROC, AP, log loss, Brier and the curves', async () => {
    const twoSpecies = await twoSpeciesIris();
    const probability = await trainModel(await requestOf(columnsOf(twoSpecies, MEASUREMENTS),
      twoSpecies.getCol('Species'), PREDICT_PROBABILITY));
    const recut = recutResult(probability, {predictProbability: true, threshold: 0.9});
    const [before, after] = [probability.metrics.train, recut.metrics.train];
    const isRecut = (['train', 'validation'] as const).some((split) =>
      recut.metrics[split].sensitivity !== probability.metrics[split].sensitivity ||
      recut.metrics[split].specificity !== probability.metrics[split].specificity);
    expect(isRecut, true, 'The threshold changed no label metric');
    expect(after.auc, before.auc);
    expect(after.averagePrecision, before.averagePrecision);
    expect(recut.metrics.validation.auc, probability.metrics.validation.auc);
    expect(recut.metrics.validation.averagePrecision, probability.metrics.validation.averagePrecision);
    for (const split of ['train', 'validation'] as const) {
      for (const id of ['logLoss', 'brier'] as const) {
        const value = probability.metrics[split][id];
        expect(typeof value === 'number' && recut.metrics[split][id] === value, true, `${split} ${id} ${value}`);
      }
    }
    const classes = probability.metrics.classes ?? [];
    expect(classes.some((c) => c.train.auc !== undefined), true, 'No class has AUC-ROC');
    for (const c of classes) {
      const recutClass = recut.metrics.classes?.find((r) => r.label === c.label);
      for (const split of ['train', 'validation'] as const) {
        // expect() would read an undefined expected value as its default, true.
        expect(recutClass?.[split].auc === c[split].auc, true, `${c.label} ${split} AUC-ROC`);
        expect(recutClass?.[split].averagePrecision === c[split].averagePrecision, true,
          `${c.label} ${split} average precision`);
      }
    }
    expect(recut.curves === probability.curves, true, 'The curves were recomputed');
    const {actual, metrics: {positiveClass}} = probability;
    const report = (await classificationReport(actual, trainPredictionOf(recut), positiveClass))!;
    expectFloat(report.balancedAccuracy, after.balancedAccuracy ?? NaN, 1e-9, 'recut balanced accuracy');
    expectFloat(report.mcc, after.mcc ?? NaN, 1e-9, 'recut mcc');
    for (const c of report.classes) {
      const values = recut.metrics.classes?.find((r) => r.label === c.label)?.train;
      const pairs: [string, number | undefined, number | undefined][] = [['precision', c.precision, values?.precision],
        ['recall', c.recall, values?.sensitivity], ['f1', c.f1, values?.f1]];
      for (const [id, crux, ts] of pairs)
        expectFloat(crux ?? NaN, ts ?? NaN, 1e-9, `recut ${c.label} ${id}`);
    }
    const iris = await openIris();
    const regression = await trainModel(await requestOf(columnsOf(iris, ['Sepal.Length', 'Sepal.Width']),
      iris.getCol('Petal.Length')));
    expect(recutResult(regression, {predictProbability: true, threshold: 0.5}) === regression, true,
      'A regression was cut');
  }, {timeout: ENGINE_TIMEOUT});
});

/** `forge-test-api-<time>` from the demo file [path] (no creation script). */
const apiTable = (path: string = IRIS) => openDemoTable(path, 'api');

// Forge:trainModel: the API training, its checks, its run records and its registration.
category('Training', () => {
  test('runTrainModel with a target only trains as the view would and records an api run', async () => {
    const iris = await apiTable();
    try {
      const result = await runTrainModel({table: iris, target: 'Species', showProgress: false});
      expect(result.method, 'XGBoost');
      expect(result.task, 'classification');
      expectArray(result.features, MEASUREMENTS);
      expect(result.rowCount, 150);
      expect(result.skippedRows, 0);
      expect(typeof result.metrics.validation.accuracy, 'number');
      expect(typeof result.metrics.validation.auc, 'number');
      expect(typeof result.metrics.validation.balancedAccuracy, 'number');
      expect(Array.isArray(result.metrics.classes), true, 'No classes');
      expect(result.modelId === null && result.name === null, true, JSON.stringify(result));
      const runs = await runsOf([iris]);
      expectArray(runs.map((r) => `${r.id} ${r.status} ${r.source} ${r.model_id ?? ''}`),
        [`${result.runId} completed api `]);
      expect(runs[0].seed, result.seed);
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: TIMEOUT});

  test('runTrainModel: the validation options and the split seed 42 by default', async () => {
    const iris = await apiTable();
    try {
      const plain = await runTrainModel({table: iris, target: 'Species', showProgress: false});
      expect(plain.seed, 42);
      const held = await runTrainModel({table: iris, target: 'Species', options: {scheme: 'holdout',
        testFraction: 0.3, shuffle: false}, showProgress: false});
      const regression = await runTrainModel({table: iris, target: 'Petal.Length', features: ['Sepal.Length',
        'Sepal.Width'], method: 'Linear Regression', options: {folds: 3, stratified: true}, showProgress: false});
      const runs = await runsOf([iris]);
      const splittingOf = (runId: string) => sortedJson(runs.find((r) => r.id === runId)?.splitting);
      expect(splittingOf(plain.runId), sortedJson({scheme: 'kfold', folds: 5, isStratified: true, isShuffled: true}));
      expect(splittingOf(held.runId), sortedJson({scheme: 'holdout', trainFraction: 0.7, isStratified: true,
        isShuffled: false}));
      expect(splittingOf(regression.runId), sortedJson({scheme: 'kfold', folds: 3, isStratified: false,
        isShuffled: true}));
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: TIMEOUT});

  test('runTrainModel: a regression with the given features, method and settings', async () => {
    const iris = await apiTable();
    try {
      const result = await runTrainModel({table: iris, target: 'Petal.Length', features: ['Sepal.Length',
        'Sepal.Width'], method: 'Linear Regression', hyperparameters: {iterations: 50}, showProgress: false});
      expect(result.task, 'regression');
      expect(result.method, 'Linear Regression');
      expect(typeof result.metrics.validation.r2, 'number');
      const [run] = await runsOf([iris]);
      expect(run.hyperparameters?.iterations, 50);
      expectArray(run.features?.columns?.map((c: {name: string}) => c.name) ?? [], ['Sepal.Length', 'Sepal.Width']);
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: TIMEOUT});

  test('runTrainModel: Predict probability on two classes records the step and the threshold', async () => {
    const iris = await twoSpeciesIris();
    try {
      for (const name of ['threshold', 'cutoff']) {
        const result = await runTrainModel({table: iris, target: 'Species', options: {predictProbability: true,
          [name]: 0.4, seed: 7}, showProgress: false});
        expect(typeof result.metrics.validation.auc, 'number');
        expect(typeof result.metrics.validation.logLoss, 'number');
        expect(result.seed, 7);
        const run = (await runsOf([iris])).find((r) => r.id === result.runId);
        expectArray(run?.options?.postprocessingInfo ?? [], ['class-probabilities']);
        expect(run?.options?.binaryClassificationThreshold, 0.4, name);
        expect(run?.options?.thresholdMetric === undefined, true, `${name} stored a metric`);
      }
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: TIMEOUT});

  test('runTrainModel: without a threshold, Predict probability optimizes it for thresholdMetric', async () => {
    const iris = await twoSpeciesIris();
    const name = `forge-test-api-threshold-${Date.now()}`;
    let modelId: string | null = null;
    try {
      for (const [given, metric] of [[undefined, 'f1'], ['mcc', 'mcc']] as const) {
        const options = {predictProbability: true, seed: 7, ...(given === undefined ? {} : {thresholdMetric: given})};
        const result = await runTrainModel({table: iris, target: 'Species', options, showProgress: false,
          ...(given === undefined ? {name} : {})});
        modelId ??= result.modelId;
        const engine = engineByName(trainingEngines().engines, result.method);
        const request = await prepareTraining({...selectionOf(columnsOf(iris, MEASUREMENTS), iris.getCol('Species'),
          SKIP, engine), seed: 7, steps: {...PREDICT_PROBABILITY, thresholdMetric: metric}});
        const expected = bestThresholdOf(await trainModel(request), metric);
        expect(expected !== undefined, true, `No best ${metric} threshold`);
        const run = (await runsOf([iris])).find((r) => r.id === result.runId);
        expect(run?.options?.binaryClassificationThreshold, expected, metric);
        expect(run?.options?.thresholdMetric, metric);
      }
      const saved = preparationOptionsOf((await forgeDb.models.get(modelId!)).options);
      expect(saved.thresholdMetric, 'f1');
    } finally {
      if (modelId !== null)
        await deleteModel(modelId);
      await deleteRunsOf([iris]);
    }
  }, {timeout: TIMEOUT});

  test('runTrainModel SVM with Predict probability fits probabilities', async () => {
    const iris = await apiTable();
    try {
      const probability = await runTrainModel({table: iris, target: 'Species', method: 'SVM',
        options: {predictProbability: true}, showProgress: false});
      expect(typeof probability.metrics.validation.auc, 'number');
      expect(typeof probability.metrics.validation.logLoss, 'number');
      const plain = await runTrainModel({table: iris, target: 'Species', method: 'SVM', showProgress: false});
      expect(plain.metrics.validation.auc === undefined, true, 'SVM without probabilities has AUC-ROC');
      const runs = await runsOf([iris]);
      const settingsOf = (runId: string) => runs.find((r) => r.id === runId)?.hyperparameters;
      expect(settingsOf(probability.runId)?.probability, true);
      expect(settingsOf(plain.runId)?.probability, false);
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: ENGINE_TIMEOUT});

  test('runTrainModel: Impute skips only the rows without a target', async () => {
    const demog = await apiTable('System:DemoFiles/demog.csv');
    const weightGaps = demog.getCol('WEIGHT').stats.missingValueCount;
    try {
      expect(weightGaps > 0 && demog.getCol('HEIGHT').stats.missingValueCount > 0, true, 'demog has no gaps');
      const result = await runTrainModel({table: demog, target: 'WEIGHT', features: ['AGE', 'HEIGHT'],
        options: {missingValues: 'impute', neighbors: 3}, showProgress: false});
      expect(result.skippedRows, weightGaps);
      const [run] = await runsOf([demog]);
      const {mode, neighbors, distance, skippedRows} = run.options?.missingValues ?? {};
      expect(`${mode} ${neighbors} ${distance} ${skippedRows}`, `impute 3 Euclidean ${weightGaps}`);
    } finally {
      await deleteRunsOf([demog]);
    }
  }, {timeout: ENGINE_TIMEOUT});

  test('runTrainModel: missing values are kept by default only for a method that handles them', async () => {
    const iris = await apiTable();
    withGaps(columnsOf(iris, ['Sepal.Length', 'Sepal.Width']));
    try {
      const suggested = await runTrainModel({table: iris, target: 'Species', showProgress: false});
      expect(`${suggested.method} ${suggested.rowCount} ${suggested.skippedRows}`, 'XGBoost 150 0');
      const svm = await runTrainModel({table: iris, target: 'Species', method: 'SVM', showProgress: false});
      expect(`${svm.rowCount} ${svm.skippedRows}`, '148 2');
      await expectExceptionAsync(() => runTrainModel({table: iris, target: 'Species', method: 'SVM',
        options: {missingValues: 'keep'}, showProgress: false}).then(() => {}), failsWith(needsRows('SVM')));
      await expectExceptionAsync(() => runTrainModel({table: iris, target: 'Petal.Length', features: ['Sepal.Length',
        'Sepal.Width'], options: {missingValues: 'keep'}, showProgress: false}).then(() => {}),
      failsWith(needsRows('Linear Regression')));
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: ENGINE_TIMEOUT});

  test('runTrainModel refuses a wrong input in one sentence and records nothing', async () => {
    const iris = await apiTable();
    const {engines} = trainingEngines();
    const xgboostSettings = hyperparametersOf(engineByName(engines, 'XGBoost')).map((p) => p.name).join(', ');
    const refusals: [Partial<TrainModelInput>, string][] = [
      [{target: 'Specie'}, 'The table has no column \'Specie\'.'],
      [{features: ['Sepal.Length', 'a', 'b']}, 'The table has no columns \'a\', \'b\'.'],
      [{method: 'Boost'}, `No method 'Boost' is installed. The methods are: ${engines.map((e) => e.name).join(', ')}.`],
      [{method: 'XGBoost', hyperparameters: {speed: 1}},
        `The method XGBoost has no setting 'speed'. Its settings are: ${xgboostSettings}.`],
      [{method: 'Softmax', hyperparameters: {penalty: 'fast'}}, 'The setting \'penalty\' of Softmax needs a number.'],
      [{hyperparameters: {iterations: 500}}, 'The setting \'iterations\' of XGBoost must be at most 100.'],
      [{options: {cutof: 0.5}}, 'Unknown option \'cutof\'. The options are: missingValues, neighbors, distance, ' +
        'oneHot, skipUniqueCategories, predictProbability, threshold, thresholdMetric, seed, scheme, folds, ' +
        'testFraction, stratified, shuffle, search.'],
      [{options: {predictProbability: true, thresholdMetric: 'auc'}}, 'The option \'thresholdMetric\' needs one of: ' +
        'f1, balancedAccuracy, mcc, accuracy.'],
      [{options: {scheme: 'loo'}}, 'The option \'scheme\' needs kfold or holdout.'],
      [{options: {scheme: 'holdout', testFraction: 0.6}}, 'The option \'testFraction\' needs a number from 0.05 to ' +
        '0.5.'],
      [{options: {shuffle: 'no'}}, 'The option \'shuffle\' needs true or false.'],
      [{target: 'Petal.Length', options: {predictProbability: true}}, 'Predict probability needs a text or boolean ' +
        'target.'],
      [{method: 'SVM', hyperparameters: {probability: true}}, 'The method SVM has no setting \'probability\'. Its ' +
        'settings are: kernel, cost, gamma, degree, coef0, epsilon.'],
      [{options: {missingValues: 'fill'}}, 'The option \'missingValues\' needs skip, impute, keep ' +
        '(keep: only for a method that handles missing values).'],
      [{name: 'forge-test-api-ref', storage: 'reference'}, 'The table has no known source, so the model cannot ' +
        'keep a reference to it. Use storage none or copy.'],
      [{name: '  '}, 'Enter a name for the model.'],
      [{features: []}, 'Choose at least one feature.'],
    ];
    try {
      const wrong: string[] = [];
      for (const [input, message] of refusals) {
        const error = await runTrainModel({table: iris, target: 'Species', showProgress: false, ...input})
          .then(() => null, (e: unknown) => e);
        if (!(error instanceof ForgeError) || error.message !== message)
          wrong.push(`${JSON.stringify(input)}: ${String(error)}`);
      }
      expectArray(wrong, []);
      expect((await runsOf([iris])).length, 0, 'A refused training was recorded');
    } finally {
      await deleteModelsNamed('forge-test-api-ref');
      await deleteRunsOf([iris]);
    }
  }, {timeout: TIMEOUT});

  test('runTrainModel: a cancel in the task bar records a cancelled run', async () => {
    const demog = grok.data.demo.demog(100000);
    demog.name = `forge-test-api-${Date.now()}`;
    try {
      const training = runTrainModel({table: demog, target: 'sex', features: ['age', 'height', 'weight'],
        method: 'XGBoost', hyperparameters: {iterations: 100, maxDepth: 20}});
      const outcome = training.then(() => null, (e: unknown) => e);
      await awaitCheck(() => taskBarEntry('Training XGBoost model') !== undefined, 'No task-bar entry', 10000);
      cancelTaskBarEntry(taskBarEntry('Training XGBoost model'));
      const error = await outcome;
      expect(error instanceof ForgeError && error.message === 'Training was cancelled.', true, String(error));
      const runs = await runsOf([demog]);
      expectArray(runs.map((r) => `${r.status} ${r.source}`), ['cancelled api']);
    } finally {
      await deleteRunsOf([demog]);
    }
  }, {timeout: ENGINE_TIMEOUT});

  test('Forge:trainModel through grok.functions.call', async () => {
    const iris = await apiTable();
    try {
      const plain: unknown = await grok.functions.call('Forge:trainModel', {table: iris, target: 'Species'});
      expect(isRecord(plain), true, `The result is ${String(plain)}`);
      const result = isRecord(plain) ? plain : {};
      expect(result.method, 'XGBoost');
      expect(result.modelId === null, true, JSON.stringify(result));
      expectArray(stringsOf(result.features), MEASUREMENTS);
      const metrics = isRecord(result.metrics) && isRecord(result.metrics.validation) ? result.metrics.validation : {};
      expect(typeof metrics.accuracy, 'number');
      expect(isRecord(result.metrics) && Array.isArray(result.metrics.classes), true, 'No metrics.classes');

      const given: unknown = await grok.functions.call('Forge:trainModel', {table: iris, target: 'Petal.Length',
        features: ['Sepal.Length', 'Sepal.Width'], method: 'Linear Regression', hyperparameters: {iterations: 50},
        options: {seed: 11}, showProgress: false});
      expect(isRecord(given) && given.seed === 11 && given.task === 'regression', true, JSON.stringify(given));
      expectArray(stringsOf(isRecord(given) ? given.features : null), ['Sepal.Length', 'Sepal.Width']);
      const runs = await runsOf([iris]);
      expect(runs.find((r) => r.seed === 11)?.hyperparameters?.iterations, 50);
      expect(runs.every((r) => r.source === 'api'), true);
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: TIMEOUT});

  test('Forge:trainModel with a method and search tunes it', async () => {
    const iris = await apiTable();
    try {
      const value: unknown = await grok.functions.call('Forge:trainModel', {table: iris, target: 'Species',
        method: 'XGBoost', options: {seed: 11, search: {metric: 'accuracy', budgetSeconds: 20, maxTrials: 10}},
        showProgress: false});
      const result = isRecord(value) ? value : {};
      const search = searchRecordOf(result.search);
      if (search === undefined)
        throw new Error(`No search in ${JSON.stringify(result)}`);
      expect(`${search.mode} ${search.metric} ${search.budgetSeconds} ${search.seed}`, 'tune accuracy 20 11');
      expect(search.completed <= 10 && search.completed >= 1, true, `${search.completed} trials completed`);
      expectArray(search.leaders.map((l) => `${l.method} ${l.chosen}`), ['XGBoost true']);
      expect(search.holdOut === undefined, true, 'A hold out on 150 rows');
      const metrics = isRecord(result.metrics) ? result.metrics : {};
      expect(metrics.test === undefined, true, 'A test without a hold out');
      expectFloat(isRecord(metrics.validation) ? Number(metrics.validation.accuracy) : NaN, search.leaders[0].score,
        1e-9);
      const [run] = await runsOf([iris]);
      expect(sortedJson(searchRecordOf(run.search)), sortedJson(search));
      expect(run.hyperparameters?.iterations, search.leaders[0].hyperparameters.iterations);
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: ENGINE_TIMEOUT});

  test('Forge:trainModel keeps Predict probability during a search', async () => {
    const iris = await twoSpeciesIris();
    try {
      const result = await runTrainModel({table: iris, target: 'Species', method: 'XGBoost',
        options: {predictProbability: true, threshold: 0.4, seed: 3, search: {maxTrials: 2}}, showProgress: false});
      expect(result.task, 'classification');
      expect(result.search?.metric, 'balancedAccuracy');
      const run = (await runsOf([iris])).find((r) => r.id === result.runId);
      expectArray(run?.options?.postprocessingInfo ?? [], ['class-probabilities']);
      expect(run?.options?.binaryClassificationThreshold, 0.4);
      expect(run?.options?.thresholdMetric === undefined, true, 'An explicit threshold stored a metric');

      const optimized = await runTrainModel({table: iris, target: 'Species', method: 'XGBoost',
        options: {predictProbability: true, thresholdMetric: 'mcc', seed: 3, search: {maxTrials: 2}},
        showProgress: false});
      const optimizedRun = (await runsOf([iris])).find((r) => r.id === optimized.runId);
      const threshold = optimizedRun?.options?.binaryClassificationThreshold;
      expect(optimizedRun?.options?.thresholdMetric, 'mcc');
      expect(typeof threshold === 'number' && threshold > 0 && threshold < 1, true, `Threshold ${threshold}`);
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: ENGINE_TIMEOUT});

  test('Forge:trainModel refuses a bad search option', async () => {
    const iris = await apiTable();
    const refusals: [{[name: string]: unknown}, string][] = [
      [{budget: 1}, 'The option \'search\' needs an object with metric, budgetSeconds, maxTrials.'],
      [{metric: 'speed'}, 'The search metric \'speed\' is not one of: r2, rmse, mae, balancedAccuracy, accuracy, f1, ' +
        'mcc, auc, logLoss, brier.'],
      [{maxTrials: 0}, 'The option \'maxTrials\' needs a whole number of at least 1.'],
      [{metric: 'r2'}, 'The search metric \'r2\' does not apply to a classification. Choose one of: ' +
        'balancedAccuracy, accuracy, f1, mcc, auc, logLoss, brier.'],
    ];
    try {
      const wrong: string[] = [];
      for (const [search, message] of refusals) {
        const error = await runTrainModel({table: iris, target: 'Species', options: {search}, showProgress: false})
          .then(() => null, (e: unknown) => e);
        if (!(error instanceof ForgeError) || error.message !== message)
          wrong.push(`${JSON.stringify(search)}: ${String(error)}`);
      }
      expectArray(wrong, []);
      expect((await runsOf([iris])).length, 0, 'A refused search was recorded');
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: TIMEOUT});

  test('Forge:trainModel names the first failure of a search where every trial fails', async () => {
    const iris = await apiTable();
    try {
      const error = await runTrainModel({table: iris, target: 'Species', method: 'SVM',
        options: {seed: 3, search: {metric: 'auc'}}, showProgress: false}).then(() => null, (e: unknown) => e);
      const message = error instanceof ForgeError ? error.message : String(error);
      expect(message.startsWith('No trial completed: ') && message.includes('needs class probabilities'), true,
        message);
      expect((await runsOf([iris])).length, 0, 'A failed search was recorded');
    } finally {
      await deleteRunsOf([iris]);
    }
  }, {timeout: ENGINE_TIMEOUT});

  test('A saved search model carries search and test', async () => {
    const demog = await apiTable('System:DemoFiles/demog.csv');
    const name = `forge-test-search-${Date.now()}`;
    try {
      const result = await runTrainModel({table: demog, target: 'SEX', features: ['AGE', 'HEIGHT', 'WEIGHT'],
        method: 'XGBoost', options: {seed: 5, search: {maxTrials: 3, budgetSeconds: 60}}, name, showProgress: false});
      expect(typeof result.metrics.test?.accuracy, 'number');
      const rows = result.search?.holdOut?.rows ?? 0;
      expect(rows >= 1100 && rows <= 1200, true, `${rows} rows held out`);
      expect(result.rowCount, demog.rowCount - result.skippedRows - rows);
      const details = await runGetModel(result.modelId ?? '', []);
      expect(sortedJson(details.search), sortedJson(result.search));
      expect(sortedJson(details.metrics?.test), sortedJson(result.metrics.test));
    } finally {
      await deleteModelsNamed(name);
      await deleteRunsOf([demog]);
    }
  }, {timeout: ENGINE_TIMEOUT});

  test('Forge:trainModel is registered with its inputs and output', async () => {
    const func = DG.Func.find({package: 'Forge', name: 'trainModel'})[0];
    expectArray(func.inputs.map((p) => `${p.name}: ${p.propertyType}`), ['table: dataframe', 'target: string',
      'features: list', 'method: string', 'hyperparameters: map', 'options: map', 'name: string',
      'description: string', 'tags: string', 'storage: string', 'showProgress: bool']);
    expectArray(func.inputs.find((p) => p.name === 'storage')?.choices ?? [], ['none', 'reference', 'copy']);
    expectArray(func.outputs.map((p) => `${p.name}: ${p.propertyType}`), ['result: object']);
  });
});
