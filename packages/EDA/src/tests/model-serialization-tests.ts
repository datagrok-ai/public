// Model pack/unpack (serialization) round-trip tests on real datasets.
//
// For each model we train, predict with the in-memory model, then pack to
// bytes, unpack into a fresh model and predict again. The strict check is
// that the unpacked model reproduces the original predictions exactly
// (Variant B = round-trip identity + a loose sanity gate on model quality,
// so a broken `fit` is also caught, not only a broken serialization).
//
// Datasets are loaded through the platform exactly like
// control-comparisons-tests.ts `grok.dapi.files.readCsv("System:DemoFiles/...")`.

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {SoftmaxClassifier} from '../softmax-classifier';
import {getLinearRegressionParams, getPredictionByLinearRegression,
  packLinearRegressionModel, unpackLinearRegressionModel} from '../regression';
import {PlsModel} from '../pls/pls-ml';
import {accuracy, madError} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';
const CARS = 'System:DemoFiles/cars.csv';
const WINE = 'System:DemoFiles/winequality.csv';

const TIMEOUT = 30000;
const PLS_COMPONENTS = 3;
// Predictions are Float32 from identical parameters; a lossless round-trip
// reproduces them bit-for-bit, so this is effectively zero.
const ROUNDTRIP_EPS = 1e-6;

/** Numeric feature columns of `df`, excluding the target and any `exclude` names. */
function numericFeatures(df: DG.DataFrame, target: string, exclude: string[] = []): DG.ColumnList {
  const cols = df.columns.toList().filter((c) =>
    (c.type === DG.COLUMN_TYPE.INT || c.type === DG.COLUMN_TYPE.FLOAT) &&
    c.name !== target && !exclude.includes(c.name));
  return DG.DataFrame.fromColumns(cols).columns;
}

/** Coefficient of determination R² of `pred` against `target`. */
function r2(target: DG.Column, pred: DG.Column): number {
  const n = target.length;
  const t = target.getRawData();
  const p = pred.getRawData();
  let mean = 0;
  for (let i = 0; i < n; ++i)
    mean += t[i];
  mean /= n;
  let ssRes = 0;
  let ssTot = 0;
  for (let i = 0; i < n; ++i) {
    ssRes += (t[i] - p[i]) ** 2;
    ssTot += (t[i] - mean) ** 2;
  }
  return ssTot === 0 ? 0 : 1 - ssRes / ssTot;
}

// ── softmax ──────────────────────────────────────────────────────────

async function softmaxRoundTrip(df: DG.DataFrame, target: DG.Column, features: DG.ColumnList,
  accuracyFloor: number): Promise<void> {
  const model = new SoftmaxClassifier({classesCount: target.categories.length, featuresCount: features.length});
  await model.fit(features, target);
  const before = model.predict(features);

  const unpacked = new SoftmaxClassifier(undefined, model.toBytes());
  const after = unpacked.predict(features);

  // Strict: unpacked model gives identical labels.
  expect(accuracy(before, after), 1, 'softmax pack/unpack changed predictions');
  // Sanity: the trained model actually classifies.
  const acc = accuracy(target, before);
  expect(acc > accuracyFloor, true, `softmax accuracy too low: ${acc} (expected > ${accuracyFloor})`);
}

// ── linear regression ────────────────────────────────────────────────

async function linregRoundTrip(target: DG.Column, features: DG.ColumnList, r2Floor: number): Promise<void> {
  const params = await getLinearRegressionParams(features, target);
  const before = getPredictionByLinearRegression(features, params);

  // Pack/unpack exactly as train/applyLinearRegression do.
  const restored = unpackLinearRegressionModel(packLinearRegressionModel(params, features.names()));
  const after = getPredictionByLinearRegression(features, restored.params, restored.names);

  expect(madError(before, after) < ROUNDTRIP_EPS, true, 'linreg pack/unpack changed predictions');
  expect(restored.names !== undefined, true, 'linreg pack/unpack lost the feature names');
  expect(restored.names!.join(), features.names().join(), 'linreg pack/unpack changed feature names');

  // Reordering the table must not disturb the prediction: features are matched
  // by name, not by position.
  const shuffled = DG.DataFrame.fromColumns(features.toList().slice().reverse()).columns;
  const afterShuffle = getPredictionByLinearRegression(shuffled, restored.params, restored.names);
  expect(madError(before, afterShuffle) < ROUNDTRIP_EPS, true, 'linreg prediction depends on column order');

  // A model saved before the names were stored is a bare Float32Array buffer.
  const legacy = unpackLinearRegressionModel(Uint8Array.from(new Uint8Array(params.buffer)));
  const afterLegacy = getPredictionByLinearRegression(features, legacy.params, legacy.names);
  expect(legacy.names === undefined, true, 'legacy linreg model reported feature names');
  expect(madError(before, afterLegacy) < ROUNDTRIP_EPS, true, 'legacy linreg model no longer loads');

  const score = r2(target, before);
  expect(score > r2Floor, true, `linreg R² too low: ${score} (expected > ${r2Floor})`);
}

// ── PLS1 ─────────────────────────────────────────────────────────────

async function plsRoundTrip(target: DG.Column, features: DG.ColumnList, r2Floor: number): Promise<void> {
  const model = new PlsModel();
  await model.fit(features, target, PLS_COMPONENTS);
  const before = model.predict(features);

  const unpacked = new PlsModel(model.toBytes());
  const after = unpacked.predict(features);

  expect(madError(before, after) < ROUNDTRIP_EPS, true, 'PLS pack/unpack changed predictions');
  const score = r2(target, before);
  expect(score > r2Floor, true, `PLS R² too low: ${score} (expected > ${r2Floor})`);
}

category('Model serialization', () => {
  test('Softmax — iris', async () => {
    const df = await grok.dapi.files.readCsv(IRIS);
    const target = df.col('Species')!;
    const features = numericFeatures(df, 'Species', ['col 1']);
    await softmaxRoundTrip(df, target, features, 0.8);
  }, {timeout: TIMEOUT});

  test('Softmax — winequality', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    const target = df.col('type')!;
    const features = numericFeatures(df, 'type', ['quality']);
    await softmaxRoundTrip(df, target, features, 0.8);
  }, {timeout: TIMEOUT});

  test('Linear regression — cars', async () => {
    const df = await grok.dapi.files.readCsv(CARS);
    const target = df.col('price')!;
    const features = numericFeatures(df, 'price');
    await linregRoundTrip(target, features, 0.4);
  }, {timeout: TIMEOUT});

  test('Linear regression — winequality', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    const target = df.col('quality')!;
    const features = numericFeatures(df, 'quality');
    await linregRoundTrip(target, features, 0.1);
  }, {timeout: TIMEOUT});

  test('PLS — cars', async () => {
    const df = await grok.dapi.files.readCsv(CARS);
    const target = df.col('price')!;
    const features = numericFeatures(df, 'price');
    await plsRoundTrip(target, features, 0.4);
  }, {timeout: TIMEOUT});

  test('PLS — winequality', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    const target = df.col('quality')!;
    const features = numericFeatures(df, 'quality');
    await plsRoundTrip(target, features, 0.1);
  }, {timeout: TIMEOUT});
}); // Model serialization
