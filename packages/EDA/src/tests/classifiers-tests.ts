// Tests for classifiers

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';

import {applyTwice, classificationDataset, accuracy, demogFrame, expectGapsKept, expectProbabilities, expectRejected,
  GAP_ROWS, IRIS_SPECIES, irisMeasurements, madError, mulberry32, numericFeatures, predictColumn, r2, sameValues,
  seededRegression, splitTarget} from './utils';
import {FITS} from './baseline';
import {isSvmInteractive, trainSoftmax, trainSvm, trainXgboost} from '../crux/methods';
import {cruxPredict, readCruxManifest} from '../crux/crux-model';

const ROWS_K = 50;
const MIN_COLS = 2;
const COLS = 100;
const TIMEOUT = 60000;
const MIN_ACCURACY = 0.9;
/** Fixed seed so the Correctness datasets are reproducible (no random flakiness) */
const CORRECTNESS_SEED = 42;

const {softmax, xgboost, svm} = FITS;
const svmKernel = (df: DG.DataFrame, target: DG.Column, kernel: string) =>
  trainSvm(df, target, kernel, 1, 0, 3, 0, 0.1);

/** Seeded classification set: features in the frame, `Labels` as the target */
function classification(rows: number, cols: number, shift: boolean, seed?: number):
  {df: DG.DataFrame, target: DG.Column} {
  const df = classificationDataset(rows, cols, shift, seed);
  return {df, target: splitTarget(df, 'Labels')};
}

/** The target's classes as crux orders string classes */
const sortedClasses = (target: DG.Column) => [...target.categories].sort();

category('Softmax', () => {
  test(`Performance: ${ROWS_K}K samples, ${COLS} features`, async () => {
    const {df, target} = classification(ROWS_K * 1000, COLS, false);
    await predictColumn(df, await softmax(df, target));
  }, {timeout: TIMEOUT, benchmark: true});

  test('Correctness', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    const prediction = await predictColumn(df, await trainSoftmax(df, target, 500, 0.1, 1e-5));
    const acc = accuracy(target, prediction);
    expect(acc > MIN_ACCURACY, true, `Softmax failed, too small accuracy: ${acc}; expected: > ${MIN_ACCURACY}`);
  }, {timeout: TIMEOUT});

  test('Bool target', async () => {
    const {df, target} = await demogFrame('CONTROL', ['AGE', 'HEIGHT', 'WEIGHT'], true);
    expect(target.type, DG.COLUMN_TYPE.BOOL, 'CONTROL is expected to be a boolean column');

    const prediction = await applyTwice(df, await softmax(df, target), 'Softmax bool');
    expect(prediction.type, DG.COLUMN_TYPE.BOOL, 'prediction of a bool target must be a BOOL column');
  }, {timeout: TIMEOUT});

  test('C from penalty', async () => {
    const binary = await demogFrame('CONTROL', ['AGE', 'HEIGHT', 'WEIGHT'], true);
    expect(readCruxManifest(await trainSoftmax(binary.df, binary.target, 100, 0.1, 0.001)).params.C, 20);

    const multi = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    expect(multi.target.categories.length >= 3, true, 'dataset is expected to have 3+ categories');
    expect(readCruxManifest(await trainSoftmax(multi.df, multi.target, 100, 0.1, 0.001)).params.C, 10);

    await expectRejected([['penalty 0', () => trainSoftmax(multi.df, multi.target, 100, 0, 0.001)]]);
  }, {timeout: TIMEOUT});

  test('Class probabilities', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    await expectProbabilities(df, await softmax(df, target), sortedClasses(target));
  }, {timeout: TIMEOUT});
}); // Softmax

category('XGBoost', () => {
  test(`Performance: ${ROWS_K}K samples, ${COLS} features`, async () => {
    const {df, target} = classification(ROWS_K * 1000, COLS, false);
    await predictColumn(df, await xgboost(df, target));
  }, {timeout: TIMEOUT, benchmark: true});

  test('Correctness', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    const acc = accuracy(target, await predictColumn(df, await xgboost(df, target)));
    expect(acc > MIN_ACCURACY, true, `XGBoost failed, too small accuracy: ${acc}; expected: > ${MIN_ACCURACY}`);
  }, {timeout: TIMEOUT});

  test('Feature column order', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    const bytes = await xgboost(df, target);
    const before = await predictColumn(df, bytes);
    const reversed = DG.DataFrame.fromColumns(df.columns.toList().slice().reverse());
    expect(accuracy(before, await predictColumn(reversed, bytes)), 1, 'XGBoost prediction depends on column order');
  }, {timeout: TIMEOUT});

  test('Bool target', async () => {
    const {df, target} = await demogFrame('CONTROL', ['AGE', 'HEIGHT', 'WEIGHT'], false);
    expect(target.type, DG.COLUMN_TYPE.BOOL, 'CONTROL is expected to be a boolean column');

    const prediction = await applyTwice(df, await xgboost(df, target), 'XGBoost bool');
    expect(prediction.type, DG.COLUMN_TYPE.BOOL, 'prediction of a bool target must be a BOOL column');
    expect(accuracy(target, prediction) > 0.5, true, 'XGBoost bool target: training accuracy is no better than chance');
  }, {timeout: TIMEOUT});

  test('Regression target (float)', async () => {
    const df = seededRegression(1000, 4, CORRECTNESS_SEED);
    const target = splitTarget(df, 'y');
    const prediction = await applyTwice(df, await xgboost(df, target), 'XGBoost regression');
    expect(prediction.type, DG.COLUMN_TYPE.FLOAT, 'regression prediction must be a float column');
    expect(madError(target, prediction) < 0.5, true, 'XGBoost regression: training error is too large');
  }, {timeout: TIMEOUT});

  test('Multiclass target (3+ unordered categories)', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    expect(target.categories.length >= 3, true, 'dataset is expected to have 3+ categories');

    const prediction = await predictColumn(df, await xgboost(df, target));
    expect(prediction.type, DG.COLUMN_TYPE.STRING, 'multiclass prediction must be a string column');
    for (const cat of prediction.categories) {
      expect(target.categories.includes(cat), true,
        `predicted category "${cat}" is not among the target categories`);
    }
    const acc = accuracy(target, prediction);
    expect(acc > MIN_ACCURACY, true, `multiclass accuracy too small: ${acc}`);
  }, {timeout: TIMEOUT});

  test('Int target and int features', async () => {
    const samples = 500;
    const x = Int32Array.from({length: samples}, (_, i) => i % 100);
    const y = Int32Array.from({length: samples}, (_, i) => 2 * (i % 100) + 3);
    const df = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('x', x)]);
    const target = DG.Column.fromInt32Array('y', y);

    const prediction = await predictColumn(df, await xgboost(df, target));
    expect(prediction.type, DG.COLUMN_TYPE.FLOAT, 'int target prediction must be a float column');
    const score = r2(target, prediction);
    expect(score > 0.9, true, `int regression: poor training fit, R2 ${score}`);
  }, {timeout: TIMEOUT});

  test('Raw buffer longer than column (capacity)', async () => {
    const rows = 100;
    const capacity = 256;

    const exactX = Float32Array.from({length: rows}, (_, i) => i / rows);
    const exactY = Float32Array.from(exactX, (x) => 3 * x + 1);
    const bigX = new Float32Array(capacity).fill(1e6);
    const bigY = new Float32Array(capacity).fill(1e6);
    bigX.set(exactX);
    bigY.set(exactY);

    const over = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x', bigX, rows)]);
    const overTarget = DG.Column.fromFloat32Array('y', bigY, rows);
    expect(overTarget.length, rows, 'column length must be the logical length');
    const exact = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x', exactX, rows)]);

    const overPred = await predictColumn(over, await xgboost(over, overTarget));
    expect(overPred.length, rows, 'prediction length must equal the column length');
    const exactPred = await predictColumn(exact, await xgboost(exact, DG.Column.fromFloat32Array('y', exactY, rows)));
    expect(sameValues(overPred, exactPred), true, 'capacity tail leaked into computation');
  }, {timeout: TIMEOUT});

  test('Validation rejects bad inputs', async () => {
    const {df, target} = classification(100, MIN_COLS, true, CORRECTNESS_SEED);
    await expectRejected([
      ['iterations = 0', () => trainXgboost(df, target, 0, 0.3, 6, 1, 0)],
      ['eta = 5', () => trainXgboost(df, target, 10, 5, 6, 1, 0)],
      ['maxDepth = 0', () => trainXgboost(df, target, 10, 0.3, 0, 1, 0)],
      ['lambda < 0', () => trainXgboost(df, target, 10, 0.3, 6, -1, 0)],
      ['single category', () => xgboost(df, DG.Column.fromStrings('y', new Array(target.length).fill('only')))],
    ]);
  }, {timeout: TIMEOUT});

  test('Hyperparameters in the manifest', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    const params = readCruxManifest(await xgboost(df, target)).params;
    expect(params.n_estimators, 20);
    expect(params.learning_rate, 0.3);
    expect(params.max_depth, 6);
    expect(params.reg_lambda, 1);
    expect(params.reg_alpha, 0);
  }, {timeout: TIMEOUT});

  test('Missing feature values are kept', async () => {
    await expectGapsKept(xgboost, await grok.dapi.files.readCsv('System:DemoFiles/iris.csv'));
  }, {timeout: TIMEOUT});

  test('Class probabilities', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    await expectProbabilities(df, await xgboost(df, target), sortedClasses(target));

    const iris = await grok.dapi.files.readCsv('System:DemoFiles/iris.csv');
    const gapped = irisMeasurements(iris, GAP_ROWS);
    const out = await expectProbabilities(gapped, await xgboost(gapped, iris.col('Species')!), IRIS_SPECIES);
    for (const col of out.columns.toList())
      expect(col.stats.missingValueCount, 0, `${col.name}: rows with gaps must be predicted`);
  }, {timeout: TIMEOUT});

  test('Benchmarks (demog)', async () => {
    const lines: string[] = [];
    const bench = async (name: string, df: DG.DataFrame, targetName: string, featureNames: string[]) => {
      const features = df.clone(null, featureNames);
      const target = df.col(targetName)!;

      let t0 = performance.now();
      const bytes = await xgboost(features, target);
      const fitMs = performance.now() - t0;

      t0 = performance.now();
      await predictColumn(features, bytes);
      const firstMs = performance.now() - t0;

      t0 = performance.now();
      for (let k = 0; k < 5; ++k)
        await predictColumn(features, bytes);
      const cachedMs = (performance.now() - t0) / 5;

      lines.push(`BENCH ${name} fit=${fitMs.toFixed(0)}ms ` +
        `first-predict=${firstMs.toFixed(1)}ms cached-predict=${cachedMs.toFixed(1)}ms`);
    };

    const NUM = ['AGE', 'HEIGHT', 'WEIGHT'];
    for (const rows of [10000, 50000, 100000]) {
      const df = grok.data.demo.demog(rows);
      await bench(`reg-AGE-${rows / 1000}k`, df, 'AGE', ['HEIGHT', 'WEIGHT']);
      await bench(`binary-SEX-${rows / 1000}k`, df, 'SEX', NUM);
      await bench(`multi-RACE-${rows / 1000}k`, df, 'RACE', NUM);
      await bench(`bool-CONTROL-${rows / 1000}k`, df, 'CONTROL', NUM);
    }

    for (const [rows, cols] of [[10000, 100], [50000, 20], [100000, 10]]) {
      const df = grok.data.demo.randomWalk(rows, cols + 1);
      const names = df.columns.names();
      await bench(`walk-reg-${rows / 1000}kx${cols}`, df, names[names.length - 1], names.slice(0, cols));
    }

    for (const line of lines)
      console.log(line);
    expect(lines.length, 15, 'all benchmark cases must complete');
  }, {timeout: 180000, benchmark: true});
}); // XGBoost

category('SVM', () => {
  test('Correctness (classification)', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    const acc = accuracy(target, await applyTwice(df, await svm(df, target), 'SVM'));
    expect(acc > MIN_ACCURACY, true, `SVM accuracy too small: ${acc}; expected > ${MIN_ACCURACY}`);
  }, {timeout: TIMEOUT});

  test('Multiclass (3+ categories)', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    expect(target.categories.length >= 3, true, 'dataset is expected to have 3+ categories');

    const prediction = await predictColumn(df, await svm(df, target));
    expect(prediction.type, DG.COLUMN_TYPE.STRING, 'multiclass prediction must be a string column');
    for (const cat of prediction.categories) {
      expect(target.categories.includes(cat), true,
        `predicted category "${cat}" is not among the target categories`);
    }
    expect(accuracy(target, prediction) > MIN_ACCURACY, true, 'multiclass accuracy too small');
  }, {timeout: TIMEOUT});

  test('Bool target', async () => {
    const rows = 300;
    const rng = mulberry32(CORRECTNESS_SEED);
    const x0 = new Float32Array(rows);
    const x1 = new Float32Array(rows);
    const y = new Array<boolean>(rows);
    for (let i = 0; i < rows; ++i) {
      x0[i] = rng();
      x1[i] = rng();
      y[i] = (x0[i] + x1[i]) > 1;
    }
    const df = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x0', x0), DG.Column.fromFloat32Array('x1', x1)]);
    const target = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'y', y);

    const prediction = await applyTwice(df, await svm(df, target), 'SVM bool');
    expect(prediction.type, DG.COLUMN_TYPE.BOOL, 'prediction of a bool target must be a BOOL column');
    expect(accuracy(target, prediction) > MIN_ACCURACY, true, 'bool training accuracy too low');
  }, {timeout: TIMEOUT});

  test('Regression (float)', async () => {
    const rows = 300;
    const rng = mulberry32(CORRECTNESS_SEED);
    const x0 = new Float32Array(rows);
    const x1 = new Float32Array(rows);
    const y = new Float32Array(rows);
    for (let i = 0; i < rows; ++i) {
      x0[i] = rng();
      x1[i] = rng();
      y[i] = 2 * x0[i] + 3 * x1[i];
    }
    const df = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x0', x0), DG.Column.fromFloat32Array('x1', x1)]);
    const target = DG.Column.fromFloat32Array('y', y);

    const prediction = await applyTwice(df, await svm(df, target), 'SVM regression');
    expect(prediction.type, DG.COLUMN_TYPE.FLOAT, 'regression prediction must be a float column');
    expect(r2(target, prediction) > 0.5, true, 'SVM regression training fit too poor');
  }, {timeout: TIMEOUT});

  test('Standardization is applied', async () => {
    const {df, target} = classification(200, MIN_COLS, true, CORRECTNESS_SEED);
    const raw0 = df.columns.byIndex(0).getRawData();
    const raw1 = df.columns.byIndex(1).getRawData();
    const n = target.length;
    const scaled = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('big', Float32Array.from({length: n}, (_, i) => raw0[i] * 1000)),
      DG.Column.fromFloat32Array('f1', Float32Array.from({length: n}, (_, i) => raw1[i]))]);

    const acc = accuracy(target, await predictColumn(scaled, await svm(scaled, target)));
    expect(acc > MIN_ACCURACY, true, `standardization not applied: accuracy ${acc}`);
  }, {timeout: TIMEOUT});

  test('Kernel round-trip (RBF & linear)', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    for (const kernel of ['RBF', 'Linear'])
      await applyTwice(df, await svmKernel(df, target, kernel), `kernel ${kernel}`);
  }, {timeout: TIMEOUT});

  test('Kernel and gamma in the manifest', async () => {
    const iris = await grok.dapi.files.readCsv('System:DemoFiles/iris.csv');
    const df = iris.clone(null, ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width']);
    const target = iris.col('Species')!;

    const rbf = readCruxManifest(await svm(df, target)).params;
    expect(rbf.kernel, 'rbf');
    expect(rbf.gamma, 1 / 4);
    expect(readCruxManifest(await svmKernel(df, target, 'Polynomial')).params.kernel, 'poly');
  }, {timeout: TIMEOUT});

  test('Validation rejects bad inputs', async () => {
    const {df, target} = classification(100, MIN_COLS, true, CORRECTNESS_SEED);
    const misFeat = DG.Column.fromFloat32Array('mf', new Float32Array(target.length).fill(1));
    misFeat.set(0, null);

    await expectRejected([
      ['C = 0', () => trainSvm(df, target, 'RBF', 0, 0, 3, 0, 0.1)],
      ['degree = 0', () => trainSvm(df, target, 'Polynomial', 1, 0, 0, 0, 0.1)],
      ['missing in feature', () => svm(DG.DataFrame.fromColumns([misFeat]), target)],
      ['single category', () => svm(df, DG.Column.fromStrings('y', new Array(target.length).fill('only')))],
    ]);

    const holed = DG.Column.fromFloat32Array('y', Float32Array.from(df.columns.byIndex(0).getRawData()));
    holed.set(0, null);
    expect(readCruxManifest(await svm(df, holed)).training.n_rows, target.length - 1, 'a null target row is skipped');
  }, {timeout: TIMEOUT});

  test('Raw buffer longer than column (capacity)', async () => {
    const rows = 100;
    const capacity = 256;
    const rng = mulberry32(CORRECTNESS_SEED);

    const exactX = Float32Array.from({length: rows}, () => rng());
    const target = DG.Column.fromStrings('y', Array.from(exactX, (x) => x > 0.5 ? 'A' : 'B'));
    const bigX = new Float32Array(capacity).fill(1e6);
    bigX.set(exactX);
    const over = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x', bigX, rows)]);
    const exact = DG.DataFrame.fromColumns([DG.Column.fromFloat32Array('x', exactX, rows)]);

    const overPred = await predictColumn(over, await svm(over, target));
    const exactPred = await predictColumn(exact, await svm(exact, target));
    expect(accuracy(overPred, exactPred), 1, 'capacity tail leaked into computation');
  }, {timeout: TIMEOUT});

  test('Interactive on winequality', async () => {
    const df = await grok.dapi.files.readCsv('System:DemoFiles/winequality.csv');
    const features = numericFeatures(df, 'type', ['quality']);
    const target = df.col('type')!;
    expect(isSvmInteractive(features, target), true,
      `winequality classification must be interactive (${df.rowCount} rows)`);
    expect(isSvmInteractive(numericFeatures(df, 'type'), df.col('quality')!), true,
      `winequality regression must be interactive (${df.rowCount} rows)`);

    const t0 = Date.now();
    await svm(DG.DataFrame.fromColumns(features.toList()), target);
    console.log(`winequality classification fit: ${Date.now() - t0} ms`);
  }, {timeout: TIMEOUT});

  test('No class probabilities', async () => {
    const {df, target} = classification(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    expectArray((await cruxPredict(df, await svm(df, target))).columns.names(), ['Prediction']);
  }, {timeout: TIMEOUT});
}); // SVM
