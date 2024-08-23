// Tests for classifiers

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';

import {classificationDataset, accuracy} from './utils';
import {SoftmaxClassifier} from '../softmax-classifier';
import {XGBooster} from '../xgbooster';

const ROWS_K = 50;
const MIN_COLS = 2;
const COLS = 100;
const TIMEOUT = 8000;
const MIN_ACCURACY = 0.9;

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
    // Prepare data
    const df = classificationDataset(ROWS_K, MIN_COLS, true);
    const features = df.columns;
    const target = features.byIndex(MIN_COLS);
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
    const prediction = unpackedModel.predict(features);

    // Evaluate accuracy
    const acc = accuracy(target, prediction);
    expect(
      acc > MIN_ACCURACY,
      true,
      `Softmax failed, too small accuracy: ${acc}; expected: <= ${MIN_ACCURACY}`,
    );
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
    // Prepare data
    const df = classificationDataset(ROWS_K, MIN_COLS, true);
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
}); // XGBoost
