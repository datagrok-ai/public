// Model pack/unpack (serialization) round-trip tests on real datasets.
//
// For each model we train, predict with the in-memory model, then pack to
// bytes, unpack into a fresh model and predict again. The strict check is
// that the unpacked model reproduces the original predictions exactly
// (Variant B = round-trip identity + a loose sanity gate on model quality,
// so a broken `fit` is also caught, not only a broken serialization).
//
// Datasets are loaded through the platform exactly like
// control-comparisons-tests.ts: `grok.dapi.files.readCsv("System:DemoFiles/...")`.

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {SoftmaxClassifier} from '../softmax-classifier';
import {getLinearRegressionParams, getPredictionByLinearRegression} from '../regression';
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

  // Pack/unpack exactly as train/applyLinearRegression do, through a copied
  // byte buffer so it is a genuine serialization round-trip.
  const bytes = Uint8Array.from(new Uint8Array(params.buffer));
  const restored = new Float32Array(bytes.buffer);
  const after = getPredictionByLinearRegression(features, restored);

  expect(madError(before, after) < ROUNDTRIP_EPS, true, 'linreg pack/unpack changed predictions');
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
