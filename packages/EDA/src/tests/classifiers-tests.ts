// Tests for classifiers

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {classificationDataset, accuracy} from './utils';
import {SoftmaxClassifier} from '../softmax-classifier';
import {XGBooster} from '../xgbooster';

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
}); // XGBoost
