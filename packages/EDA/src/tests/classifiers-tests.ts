// Tests for classifiers

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {classificationDataset, accuracy, madError, mulberry32} from './utils';
import {SoftmaxClassifier} from '../softmax-classifier';
import {XGBooster} from '../xgbooster';

/** Expect that an async action throws (validation must reject before wasm). */
async function expectThrows(action: () => Promise<unknown>, label: string): Promise<void> {
  let thrown = false;
  try {
    await action();
  } catch {
    thrown = true;
  }
  expect(thrown, true, `${label}: expected an error, but none was thrown`);
}

const ROWS_K = 50;
const MIN_COLS = 2;
const COLS = 100;
const TIMEOUT = 8000;
const MIN_ACCURACY = 0.9;
/** Fixed seed so the Correctness datasets are reproducible (no random flakiness) */
const CORRECTNESS_SEED = 42;

category('Softmax', () => {
  test(`Performance: ${ROWS_K}K samples, ${COLS} features`, async () => {
    // Data
    const df = classificationDataset(ROWS_K * 1000, COLS, false);
    const features = df.columns;
    const target = features.byIndex(COLS);
    features.remove(target.name);

    // Fit & pack trained model
    const model = new SoftmaxClassifier({
      classesCount: target.categories.length,
      featuresCount: features.length,
    });
    await model.fit(features, target);
    const modelBytes = model.toBytes();

    // Unpack & apply model
    const unpackedModel = new SoftmaxClassifier(undefined, modelBytes);
    unpackedModel.predict(features);
  }, {timeout: TIMEOUT, benchmark: true});

  test('Correctness', async () => {
    // Prepare data — fixed seed for a reproducible, non-flaky dataset
    const df = classificationDataset(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    const features = df.columns;
    const target = features.byIndex(MIN_COLS);
    features.remove(target.name);

    // Fit & pack trained model
    const model = new SoftmaxClassifier({
      classesCount: target.categories.length,
      featuresCount: features.length,
    });

    // More iterations + tighter tolerance: the Rust GD kernel under-converges
    // on this tiny (50-row) separable set with the defaults, where the old C++
    // solver did not. The problem is linearly separable, so a converged model
    // clears the 0.9 bar deterministically.
    await model.fit(features, target, 1, 500, 0.1, 1e-5);
    const modelBytes = model.toBytes();

    // Unpack & apply model
    const unpackedModel = new SoftmaxClassifier(undefined, modelBytes);
    const prediction = unpackedModel.predict(features);

    // Evaluate accuracy
    const acc = accuracy(target, prediction);
    expect(
      acc > MIN_ACCURACY,
      true,
      `Softmax failed, too small accuracy: ${acc}; expected: <= ${MIN_ACCURACY}`,
    );
  }, {timeout: TIMEOUT});

  test('Bool target', async () => {
    // demog has a boolean column CONTROL — use it as a classification target
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    const target = df.col('CONTROL')!;
    expect(target.type, DG.COLUMN_TYPE.BOOL, 'CONTROL is expected to be a boolean column');

    const features = DG.DataFrame.fromColumns(
      ['AGE', 'HEIGHT', 'WEIGHT'].map((name) => df.col(name)!),
    ).columns;

    // Fit & pack trained model (bool target → 2 classes)
    const model = new SoftmaxClassifier({classesCount: 2, featuresCount: features.length});
    await model.fit(features, target);
    const before = model.predict(features);

    // The prediction must be a BOOL column, not strings 'true'/'false'
    expect(before.type, DG.COLUMN_TYPE.BOOL, 'prediction of a bool target must be a BOOL column');

    // Round-trip: the unpacked model reproduces the predictions exactly (covers the trailing flag byte)
    const unpackedModel = new SoftmaxClassifier(undefined, model.toBytes());
    const after = unpackedModel.predict(features);
    expect(after.type, DG.COLUMN_TYPE.BOOL, 'unpacked model must also predict a BOOL column');
    expect(accuracy(before, after), 1, 'Softmax bool pack/unpack changed predictions');
  }, {timeout: TIMEOUT});
}); // Softmax

category('XGBoost', () => {
  test(`Performance: ${ROWS_K}K samples, ${COLS} features`, async () => {
    // Data
    const df = classificationDataset(ROWS_K * 1000, COLS, false);
    const features = df.columns;
    const target = features.byIndex(COLS);
    features.remove(target.name);

    // Fit & pack trained model
    const model = new XGBooster();
    await model.fit(features, target);
    const modelBytes = model.toBytes();

    // Unpack & apply model
    const unpackedModel = new XGBooster(modelBytes);
    unpackedModel.predict(features);
  }, {timeout: TIMEOUT, benchmark: true});

  test('Correctness', async () => {
    // Prepare data — fixed seed for a reproducible, non-flaky dataset
    const df = classificationDataset(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    const features = df.columns;
    const target = features.byIndex(MIN_COLS);
    features.remove(target.name);

    // Fit & pack trained model
    const model = new XGBooster();

    await model.fit(features, target);
    const modelBytes = model.toBytes();

    // Unpack & apply model
    const unpackedModel = new XGBooster(modelBytes);
    const prediction = unpackedModel.predict(features);

    // Evaluate accuracy
    const acc = accuracy(target, prediction);
    expect(
      acc > MIN_ACCURACY,
      true,
      `XGBoost failed, too small accuracy: ${acc}; expected: <= ${MIN_ACCURACY}`,
    );
  }, {timeout: TIMEOUT});

  test('Bool target', async () => {
    // demog has a boolean column CONTROL — use it as a classification target
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    const target = df.col('CONTROL')!;
    expect(target.type, DG.COLUMN_TYPE.BOOL, 'CONTROL is expected to be a boolean column');

    const features = DG.DataFrame.fromColumns(
      ['AGE', 'HEIGHT', 'WEIGHT'].map((name) => df.col(name)!),
    ).columns;

    // Fit & pack trained model
    const model = new XGBooster();
    await model.fit(features, target);
    const before = model.predict(features);

    // The prediction must be a BOOL column, not strings 'true'/'false'
    expect(before.type, DG.COLUMN_TYPE.BOOL, 'prediction of a bool target must be a BOOL column');

    // Round-trip: the unpacked model reproduces the predictions exactly (covers WAS_BOOL flag)
    const unpackedModel = new XGBooster(model.toBytes());
    const after = unpackedModel.predict(features);
    expect(after.type, DG.COLUMN_TYPE.BOOL, 'unpacked model must also predict a BOOL column');
    expect(accuracy(before, after), 1, 'XGBoost bool pack/unpack changed predictions');

    // Sanity: the model fits the training data better than a constant predictor
    expect(
      accuracy(target, before) > 0.5,
      true,
      'XGBoost bool target: training accuracy is no better than chance',
    );
  }, {timeout: TIMEOUT});

  test('Regression target (float)', async () => {
    // Seeded data: regressionDataset() is unseeded (randomWalk + Math.random
    // coefficients), which made the madError threshold flaky.
    const rows = 1000;
    const cols = 4;
    const rng = mulberry32(CORRECTNESS_SEED);
    const raw = Array.from({length: cols}, () => {
      const arr = new Float32Array(rows);
      for (let i = 0; i < rows; ++i)
        arr[i] = rng();
      return arr;
    });
    const y = new Float32Array(rows);
    for (let i = 0; i < rows; ++i)
      y[i] = raw.reduce((sum, arr, j) => sum + arr[i] * (j + 1), 0);

    const features = DG.DataFrame.fromColumns(
      raw.map((arr, j) => DG.Column.fromFloat32Array(`x${j}`, arr))).columns;
    const target = DG.Column.fromFloat32Array('y', y);

    const model = new XGBooster();
    await model.fit(features, target);
    const before = model.predict(features);
    expect(before.type, DG.COLUMN_TYPE.FLOAT, 'regression prediction must be a float column');
    expect(
      madError(target, before) < 0.5,
      true,
      'XGBoost regression: training error is too large',
    );

    // Round-trip must reproduce predictions exactly
    const after = new XGBooster(model.toBytes()).predict(features);
    for (let i = 0; i < before.length; ++i) {
      if (before.get(i) !== after.get(i))
        throw new Error(`regression pack/unpack changed prediction at row ${i}`);
    }
  }, {timeout: TIMEOUT});

  test('Multiclass target (3+ unordered categories)', async () => {
    // classificationDataset yields up to 4 combined categories -> multi:softmax
    const df = classificationDataset(ROWS_K, MIN_COLS, true, CORRECTNESS_SEED);
    const features = df.columns;
    const target = features.byIndex(MIN_COLS);
    features.remove(target.name);
    expect(target.categories.length >= 3, true, 'dataset is expected to have 3+ categories');

    const model = new XGBooster();
    await model.fit(features, target);
    const prediction = model.predict(features);

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
    const x = new Int32Array(samples);
    const y = new Int32Array(samples);
    for (let i = 0; i < samples; ++i) {
      x[i] = i % 100;
      y[i] = 2 * (i % 100) + 3;
    }
    const features = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('x', x)]).columns;
    const target = DG.Column.fromInt32Array('y', y);

    const model = new XGBooster();
    await model.fit(features, target);
    const prediction = model.predict(features);
    expect(prediction.type, DG.COLUMN_TYPE.INT, 'int target prediction must be an int column');
    expect(accuracy(target, prediction) > 0.9, true, 'int regression: poor training fit');
  }, {timeout: TIMEOUT});

  test('Raw buffer longer than column (capacity)', async () => {
    // col.getRawData() may return an array LONGER than col.length; the glue
    // must slice to col.length - otherwise garbage rows / heap overrun.
    const rows = 100;
    const capacity = 256;

    const exactX = new Float32Array(rows);
    const exactY = new Float32Array(rows);
    for (let i = 0; i < rows; ++i) {
      exactX[i] = i / rows;
      exactY[i] = 3 * exactX[i] + 1;
    }

    // Oversized backing buffers with garbage in the tail
    const bigX = new Float32Array(capacity).fill(1e6);
    const bigY = new Float32Array(capacity).fill(1e6);
    bigX.set(exactX);
    bigY.set(exactY);

    const overFeatures = DG.DataFrame.fromColumns(
      [DG.Column.fromFloat32Array('x', bigX, rows)]).columns;
    const overTarget = DG.Column.fromFloat32Array('y', bigY, rows);
    expect(overTarget.length, rows, 'column length must be the logical length');

    const exactFeatures = DG.DataFrame.fromColumns(
      [DG.Column.fromFloat32Array('x', exactX, rows)]).columns;
    const exactTarget = DG.Column.fromFloat32Array('y', exactY, rows);

    const overModel = new XGBooster();
    await overModel.fit(overFeatures, overTarget);
    const overPred = overModel.predict(overFeatures);
    expect(overPred.length, rows, 'prediction length must equal the column length');

    const exactModel = new XGBooster();
    await exactModel.fit(exactFeatures, exactTarget);
    const exactPred = exactModel.predict(exactFeatures);

    // Identical data -> identical predictions; any difference means the
    // capacity tail leaked into training or prediction.
    for (let i = 0; i < rows; ++i) {
      if (overPred.get(i) !== exactPred.get(i))
        throw new Error(`capacity tail leaked into computation at row ${i}`);
    }
  }, {timeout: TIMEOUT});

  test('Validation rejects bad inputs before wasm', async () => {
    const df = classificationDataset(100, MIN_COLS, true, CORRECTNESS_SEED);
    const features = df.columns;
    const target = features.byIndex(MIN_COLS);
    features.remove(target.name);

    await expectThrows(() => new XGBooster().fit(features, target, 0), 'iterations = 0');
    await expectThrows(() => new XGBooster().fit(features, target, 10, 5), 'eta = 5');
    await expectThrows(() => new XGBooster().fit(features, target, 10, 0.3, 0), 'maxDepth = 0');
    await expectThrows(() => new XGBooster().fit(features, target, 10, 0.3, 6, -1), 'lambda < 0');

    // Missing values in the target must be rejected (the wasm build has no
    // exceptions - a bad label would abort the whole instance).
    const holedTarget = DG.Column.fromFloat32Array('y',
      new Float32Array(target.length).fill(1));
    holedTarget.set(0, null);
    await expectThrows(() => new XGBooster().fit(features, holedTarget), 'missing in target');

    // Single-category classification target is meaningless
    const oneCat = DG.Column.fromStrings('y', new Array(target.length).fill('only'));
    await expectThrows(() => new XGBooster().fit(features, oneCat), 'single category');
  }, {timeout: TIMEOUT});

  test('Pre-1.7.0 (v1) model container is rejected', async () => {
    // v1 support is removed: a container without the Version header field
    // must fail with a clear "retrain the model" error, not garbage output.
    const v1Header = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Type', [DG.COLUMN_TYPE.STRING as string]),
      DG.Column.fromInt32Array('Params count', new Int32Array([4])),
      DG.Column.fromInt32Array('Categories size', new Int32Array([0])),
    ]).toByteArray();

    const v1 = new Uint8Array(4 + v1Header.length + 64);
    new Uint32Array(v1.buffer, 0, 1)[0] = v1Header.length;
    v1.set(v1Header, 4);

    let message = '';
    try {
      new XGBooster(v1);
    } catch (err) {
      message = err instanceof Error ? err.message : String(err);
    }
    expect(message.includes('older EDA version'), true,
      `v1 container must be rejected with a clear error, got: "${message}"`);
  }, {timeout: TIMEOUT});

  test('Benchmarks (demog)', async () => {
    // End-to-end benchmarks through the full XGBooster stack (column packing,
    // training worker, handle cache) on realistic platform data.
    // Times are printed as BENCH lines to the console.
    const lines: string[] = [];
    const bench = async (name: string, df: DG.DataFrame,
      targetName: string, featureNames: string[]) => {
      const features = DG.DataFrame.fromColumns(
        featureNames.map((n) => df.col(n)!)).columns;
      const target = df.col(targetName)!;

      const model = new XGBooster();
      let t0 = performance.now();
      await model.fit(features, target);
      const fitMs = performance.now() - t0;

      t0 = performance.now();
      model.predict(features); // first predict: includes model load
      const firstMs = performance.now() - t0;

      t0 = performance.now();
      for (let k = 0; k < 5; ++k)
        model.predict(features); // cached-handle predicts
      const cachedMs = (performance.now() - t0) / 5;

      model.dispose();
      lines.push(`BENCH ${name} fit=${fitMs.toFixed(0)}ms ` +
        `first-predict=${firstMs.toFixed(1)}ms cached-predict=${cachedMs.toFixed(1)}ms`);
    };

    const NUM = ['AGE', 'HEIGHT', 'WEIGHT'];
    for (const rows of [10000, 50000, 100000]) {
      const df = grok.data.demo.demog(rows);
      // AGE as a regression target: demog has no missing values in it,
      // while HEIGHT/WEIGHT nulls in features are handled by XGBoost.
      await bench(`reg-AGE-${rows / 1000}k`, df, 'AGE', ['HEIGHT', 'WEIGHT']);
      await bench(`binary-SEX-${rows / 1000}k`, df, 'SEX', NUM);
      await bench(`multi-RACE-${rows / 1000}k`, df, 'RACE', NUM);
      await bench(`bool-CONTROL-${rows / 1000}k`, df, 'CONTROL', NUM);
    }

    // Wide numeric matrices via randomWalk - the shapes match the wasm-level
    // node benchmarks (smoke/bench.mjs), so the full-stack overhead (column
    // packing, worker transfer) is directly comparable to the bare module.
    for (const [rows, cols] of [[10000, 100], [50000, 20], [100000, 10]]) {
      const df = grok.data.demo.randomWalk(rows, cols + 1);
      const names = df.columns.names();
      await bench(`walk-reg-${rows / 1000}kx${cols}`, df,
        names[names.length - 1], names.slice(0, cols));
    }

    for (const line of lines)
      console.log(line);
    expect(lines.length, 15, 'all benchmark cases must complete');
  }, {timeout: 180000, benchmark: true});

  test('Handle cache and dispose', async () => {
    const df = classificationDataset(200, MIN_COLS, true, CORRECTNESS_SEED);
    const features = df.columns;
    const target = features.byIndex(MIN_COLS);
    features.remove(target.name);

    const model = new XGBooster();
    await model.fit(features, target);

    // Repeated predicts reuse the cached handle and must be identical
    const first = model.predict(features);
    for (let k = 0; k < 5; ++k) {
      const next = model.predict(features);
      expect(accuracy(first, next), 1, `repeated predict #${k} differs`);
    }

    // dispose() frees the handle; the next predict reloads from bytes
    model.dispose();
    const after = model.predict(features);
    expect(accuracy(first, after), 1, 'predict after dispose differs');
  }, {timeout: TIMEOUT});
}); // XGBoost
