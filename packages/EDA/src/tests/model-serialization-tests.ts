// Serialization round-trip + apply-time feature matching, on real datasets.
// See expectFeatureMatching for the order/subset/missing checks.

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {SoftmaxClassifier} from '../softmax-classifier';
import {getLinearRegressionParams, getPredictionByLinearRegression,
  packLinearRegressionModel, unpackLinearRegressionModel} from '../regression';
import {PlsModel} from '../pls/pls-ml';
import {SVM, SvmFitOptions} from '../svm';
import {XGBooster} from '../xgbooster';
import {KernelType} from '../../wasm/svm';
import {accuracy, madError} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';
const CARS = 'System:DemoFiles/cars.csv';
const WINE = 'System:DemoFiles/winequality.csv';
const DEMOG = 'System:DemoFiles/demog.csv';

const TIMEOUT = 30000;
const PLS_COMPONENTS = 3;
// Effectively zero: identical params → identical Float32.
const ROUNDTRIP_EPS = 1e-6;

// ── data helpers ──────────────────────────────────────────────────────

/** Numeric feature columns of `df`, excluding the target and any `exclude` names. */
function numericFeatures(df: DG.DataFrame, target: string, exclude: string[] = []): DG.ColumnList {
  const cols = df.columns.toList().filter((c) =>
    (c.type === DG.COLUMN_TYPE.INT || c.type === DG.COLUMN_TYPE.FLOAT) &&
    c.name !== target && !exclude.includes(c.name));
  return DG.DataFrame.fromColumns(cols).columns;
}

/** Rows where every `keep` column is non-null. */
function nonNullRows(df: DG.DataFrame, keep: string[]): DG.DataFrame {
  const cols = keep.map((n) => df.col(n)!);
  return df.clone(DG.BitSet.create(df.rowCount, (i) => cols.every((c) => !c.isNone(i))));
}

/** Up to `n` evenly-strided rows (a prefix would drop a class in ordered data). */
function sampleRows(df: DG.DataFrame, n: number): DG.DataFrame {
  if (df.rowCount <= n)
    return df;
  const stride = Math.ceil(df.rowCount / n);
  return df.clone(DG.BitSet.create(df.rowCount, (i) => i % stride === 0));
}

/** Independent column copy — won't reparent the original. */
function copyCol(col: DG.Column): DG.Column {
  const n = col.length;
  switch (col.type) {
  case DG.COLUMN_TYPE.INT:
    return DG.Column.fromInt32Array(col.name, Int32Array.from(col.getRawData().subarray(0, n)), n);
  case DG.COLUMN_TYPE.FLOAT:
    return DG.Column.fromFloat32Array(col.name, Float32Array.from(col.getRawData().subarray(0, n)));
  default:
    return DG.Column.fromStrings(col.name, col.toList() as string[]);
  }
}

/** Deterministic decoy columns (unrelated columns the apply-time table may carry). */
function decoyColumns(rows: number, count: number): DG.Column[] {
  const out = new Array<DG.Column>(count);
  for (let k = 0; k < count; ++k) {
    const arr = new Float32Array(rows);
    for (let i = 0; i < rows; ++i)
      arr[i] = ((i * 7 + k * 13) % 100) * 0.1 - 5;
    out[k] = DG.Column.fromFloat32Array(`_decoy_${k}`, arr);
  }
  return out;
}

/** Wrap columns into a fresh table's column list (an apply-time table). */
function asTable(cols: DG.Column[]): DG.ColumnList {
  return DG.DataFrame.fromColumns(cols).columns;
}

/** Expect that a synchronous action throws (all model `predict`s are sync). */
function expectThrowsSync(action: () => unknown, label: string): void {
  let thrown = false;
  try {
    action();
  } catch {
    thrown = true;
  }
  expect(thrown, true, `${label}: expected an error, but none was thrown`);
}

// ── metrics ───────────────────────────────────────────────────────────

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

type Predict = (cols: DG.ColumnList) => DG.Column;
type Same = (a: DG.Column, b: DG.Column) => boolean;

/** Classification: identical labels. */
const sameLabels: Same = (a, b) => accuracy(a, b) === 1;
/** Regression: identical Float32 values (a lossless model reproduces them). */
const sameValues: Same = (a, b) => madError(a, b) < ROUNDTRIP_EPS;

/** Apply-time contract: order- and extra-column-invariant; missing feature throws. */
function expectFeatureMatching(predict: Predict, features: DG.ColumnList, before: DG.Column, same: Same): void {
  const cols = features.toList();
  const rows = cols[0].length;

  // (1) reversed order
  expect(same(before, predict(asTable(cols.map(copyCol).reverse()))), true,
    'prediction depends on column order (reversed)');

  // (2) rotated order — a permutation distinct from the reverse (>2 features)
  if (cols.length > 2) {
    const copies = cols.map(copyCol);
    expect(same(before, predict(asTable([...copies.slice(1), copies[0]]))), true,
      'prediction depends on column order (rotated)');
  }

  // (3) superset: pick features by name from a scrambled table with extra columns.
  const decoys = decoyColumns(rows, 2);
  const superset = [decoys[0], ...cols.map(copyCol).reverse()];
  superset.splice(1 + Math.floor(cols.length / 2), 0, decoys[1]);
  expect(same(before, predict(asTable(superset))), true, 'extra columns in the table disturb the prediction');

  // (4) a missing feature must throw, not silently mispredict
  expectThrowsSync(() => predict(asTable(cols.slice(1).map(copyCol))), 'a missing feature column must throw');
}

// ── softmax ──────────────────────────────────────────────────────────

async function softmaxRoundTrip(target: DG.Column, features: DG.ColumnList, accuracyFloor: number): Promise<void> {
  const model = new SoftmaxClassifier({classesCount: target.categories.length, featuresCount: features.length});
  await model.fit(features, target);
  const before = model.predict(features);

  const unpacked = new SoftmaxClassifier(undefined, model.toBytes());
  expect(accuracy(before, unpacked.predict(features)), 1, 'softmax pack/unpack changed predictions');
  expectFeatureMatching((c) => unpacked.predict(c), features, before, sameLabels);

  const acc = accuracy(target, before);
  expect(acc > accuracyFloor, true, `softmax accuracy too low: ${acc} (expected > ${accuracyFloor})`);
}

// ── linear regression ────────────────────────────────────────────────

async function linregRoundTrip(target: DG.Column, features: DG.ColumnList, r2Floor: number): Promise<void> {
  const params = await getLinearRegressionParams(features, target);
  const before = getPredictionByLinearRegression(features, params);

  // Pack/unpack exactly as train/applyLinearRegression do.
  const restored = unpackLinearRegressionModel(packLinearRegressionModel(params, features.names()));
  expect(madError(before, getPredictionByLinearRegression(features, restored.params, restored.names)) < ROUNDTRIP_EPS,
    true, 'linreg pack/unpack changed predictions');
  expect(restored.names !== undefined, true, 'linreg pack/unpack lost the feature names');
  expect(restored.names!.join(), features.names().join(), 'linreg pack/unpack changed feature names');

  expectFeatureMatching((c) => getPredictionByLinearRegression(c, restored.params, restored.names),
    features, before, sameValues);

  // Legacy model (bare Float32Array, no names) still loads — positional matching.
  const legacy = unpackLinearRegressionModel(Uint8Array.from(new Uint8Array(params.buffer)));
  expect(legacy.names === undefined, true, 'legacy linreg model reported feature names');
  expect(madError(before, getPredictionByLinearRegression(features, legacy.params, legacy.names)) < ROUNDTRIP_EPS,
    true, 'legacy linreg model no longer loads');

  const score = r2(target, before);
  expect(score > r2Floor, true, `linreg R² too low: ${score} (expected > ${r2Floor})`);
}

// ── PLS1 ─────────────────────────────────────────────────────────────

async function plsRoundTrip(target: DG.Column, features: DG.ColumnList, r2Floor: number,
  components: number = PLS_COMPONENTS): Promise<void> {
  const model = new PlsModel();
  await model.fit(features, target, Math.min(components, features.length));
  const before = model.predict(features);

  const unpacked = new PlsModel(model.toBytes());
  expect(madError(before, unpacked.predict(features)) < ROUNDTRIP_EPS, true, 'PLS pack/unpack changed predictions');
  expectFeatureMatching((c) => unpacked.predict(c), features, before, sameValues);

  const score = r2(target, before);
  expect(score > r2Floor, true, `PLS R² too low: ${score} (expected > ${r2Floor})`);
}

// ── SVM (libsvm) ─────────────────────────────────────────────────────

const SVM_OPTS: SvmFitOptions = {kernelType: KernelType.Rbf, C: 1, gamma: 0, degree: 3, coef0: 0, epsilon: 0.1};
// SMO is superlinear in the sample count; cap large datasets for the SVM tests.
const SVM_MAX_ROWS = 1500;

async function svmClassifyRoundTrip(target: DG.Column, features: DG.ColumnList,
  accuracyFloor: number): Promise<void> {
  const model = new SVM();
  await model.fit(features, target, SVM_OPTS);
  const before = model.predict(features);

  const unpacked = new SVM(model.toBytes());
  expect(accuracy(before, unpacked.predict(features)), 1, 'SVM pack/unpack changed predictions');
  expectFeatureMatching((c) => unpacked.predict(c), features, before, sameLabels);

  const acc = accuracy(target, before);
  expect(acc > accuracyFloor, true, `SVM accuracy too low: ${acc} (expected > ${accuracyFloor})`);
}

async function svmRegressRoundTrip(target: DG.Column, features: DG.ColumnList, r2Floor: number): Promise<void> {
  const model = new SVM();
  await model.fit(features, target, SVM_OPTS);
  const before = model.predict(features);

  const unpacked = new SVM(model.toBytes());
  expect(madError(before, unpacked.predict(features)) < ROUNDTRIP_EPS, true, 'SVM pack/unpack changed predictions');
  expectFeatureMatching((c) => unpacked.predict(c), features, before, sameValues);

  const score = r2(target, before);
  expect(score > r2Floor, true, `SVM R² too low: ${score} (expected > ${r2Floor})`);
}

// ── XGBoost ──────────────────────────────────────────────────────────

async function xgbClassifyRoundTrip(target: DG.Column, features: DG.ColumnList,
  accuracyFloor: number): Promise<void> {
  const model = new XGBooster();
  await model.fit(features, target);
  const before = model.predict(features);

  const unpacked = new XGBooster(model.toBytes());
  expect(accuracy(before, unpacked.predict(features)), 1, 'XGBoost pack/unpack changed predictions');
  expectFeatureMatching((c) => unpacked.predict(c), features, before, sameLabels);

  const acc = accuracy(target, before);
  expect(acc > accuracyFloor, true, `XGBoost accuracy too low: ${acc} (expected > ${accuracyFloor})`);
}

async function xgbRegressRoundTrip(target: DG.Column, features: DG.ColumnList, r2Floor: number): Promise<void> {
  const model = new XGBooster();
  await model.fit(features, target);
  const before = model.predict(features);

  const unpacked = new XGBooster(model.toBytes());
  expect(madError(before, unpacked.predict(features)) < ROUNDTRIP_EPS, true, 'XGBoost pack/unpack changed predictions');
  expectFeatureMatching((c) => unpacked.predict(c), features, before, sameValues);

  const score = r2(target, before);
  expect(score > r2Floor, true, `XGBoost R² too low: ${score} (expected > ${r2Floor})`);
}

// ── datasets ─────────────────────────────────────────────────────────

/** demog rows without nulls in the columns used, as {features, target}. */
async function demog(featureNames: string[], targetName: string, maxRows = Number.MAX_SAFE_INTEGER):
  Promise<{features: DG.ColumnList, target: DG.Column}> {
  const clean = nonNullRows(await grok.dapi.files.readCsv(DEMOG), [...featureNames, targetName]);
  const df = clean.rowCount > maxRows ? sampleRows(clean, maxRows) : clean;
  return {
    features: DG.DataFrame.fromColumns(featureNames.map((n) => df.col(n)!)).columns,
    target: df.col(targetName)!,
  };
}

const DEMOG_FEATS = ['AGE', 'HEIGHT', 'WEIGHT'];

category('Model serialization', () => {
  // ── softmax (classification) ──
  test('Softmax — iris', async () => {
    const df = await grok.dapi.files.readCsv(IRIS);
    await softmaxRoundTrip(df.col('Species')!, numericFeatures(df, 'Species', ['col 1']), 0.8);
  }, {timeout: TIMEOUT});

  test('Softmax — winequality', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await softmaxRoundTrip(df.col('type')!, numericFeatures(df, 'type', ['quality']), 0.8);
  }, {timeout: TIMEOUT});

  test('Softmax — demog (SEX)', async () => {
    const {features, target} = await demog(DEMOG_FEATS, 'SEX');
    await softmaxRoundTrip(target, features, 0.55);
  }, {timeout: TIMEOUT});

  // ── linear regression ──
  test('Linear regression — cars', async () => {
    const df = await grok.dapi.files.readCsv(CARS);
    await linregRoundTrip(df.col('price')!, numericFeatures(df, 'price'), 0.4);
  }, {timeout: TIMEOUT});

  test('Linear regression — winequality', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await linregRoundTrip(df.col('quality')!, numericFeatures(df, 'quality'), 0.1);
  }, {timeout: TIMEOUT});

  test('Linear regression — demog (WEIGHT)', async () => {
    const {features, target} = await demog(['AGE', 'HEIGHT'], 'WEIGHT');
    await linregRoundTrip(target, features, 0.2);
  }, {timeout: TIMEOUT});

  // ── PLS ──
  test('PLS — cars', async () => {
    const df = await grok.dapi.files.readCsv(CARS);
    await plsRoundTrip(df.col('price')!, numericFeatures(df, 'price'), 0.4);
  }, {timeout: TIMEOUT});

  test('PLS — winequality', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await plsRoundTrip(df.col('quality')!, numericFeatures(df, 'quality'), 0.1);
  }, {timeout: TIMEOUT});

  test('PLS — demog (WEIGHT)', async () => {
    const {features, target} = await demog(['AGE', 'HEIGHT'], 'WEIGHT');
    await plsRoundTrip(target, features, 0.2, 2);
  }, {timeout: TIMEOUT});

  // ── SVM (classification) ──
  test('SVM — iris', async () => {
    const df = await grok.dapi.files.readCsv(IRIS);
    await svmClassifyRoundTrip(df.col('Species')!, numericFeatures(df, 'Species', ['col 1']), 0.8);
  }, {timeout: TIMEOUT});

  test('SVM — winequality (type)', async () => {
    const df = sampleRows(await grok.dapi.files.readCsv(WINE), SVM_MAX_ROWS);
    await svmClassifyRoundTrip(df.col('type')!, numericFeatures(df, 'type', ['quality']), 0.8);
  }, {timeout: TIMEOUT});

  test('SVM — demog (SEX)', async () => {
    const {features, target} = await demog(DEMOG_FEATS, 'SEX', SVM_MAX_ROWS);
    await svmClassifyRoundTrip(target, features, 0.55);
  }, {timeout: TIMEOUT});

  // ── SVM (regression) ──
  test('SVM — winequality (quality)', async () => {
    const df = sampleRows(await grok.dapi.files.readCsv(WINE), SVM_MAX_ROWS);
    await svmRegressRoundTrip(df.col('quality')!, numericFeatures(df, 'quality'), 0.3);
  }, {timeout: TIMEOUT});

  // SVR defaults (C=1) need a small-magnitude target; cars price (~10^4) underfits.
  test('SVM — iris (Petal.Length)', async () => {
    const df = await grok.dapi.files.readCsv(IRIS);
    await svmRegressRoundTrip(df.col('Petal.Length')!, numericFeatures(df, 'Petal.Length', ['col 1']), 0.5);
  }, {timeout: TIMEOUT});

  // ── XGBoost (classification) ──
  test('XGBoost — iris', async () => {
    const df = await grok.dapi.files.readCsv(IRIS);
    await xgbClassifyRoundTrip(df.col('Species')!, numericFeatures(df, 'Species', ['col 1']), 0.8);
  }, {timeout: TIMEOUT});

  test('XGBoost — winequality (type)', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await xgbClassifyRoundTrip(df.col('type')!, numericFeatures(df, 'type', ['quality']), 0.8);
  }, {timeout: TIMEOUT});

  test('XGBoost — demog (SEX)', async () => {
    const {features, target} = await demog(DEMOG_FEATS, 'SEX');
    await xgbClassifyRoundTrip(target, features, 0.55);
  }, {timeout: TIMEOUT});

  // ── XGBoost (regression) ──
  test('XGBoost — cars (price)', async () => {
    const df = await grok.dapi.files.readCsv(CARS);
    await xgbRegressRoundTrip(df.col('price')!, numericFeatures(df, 'price'), 0.4);
  }, {timeout: TIMEOUT});

  test('XGBoost — winequality (quality)', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await xgbRegressRoundTrip(df.col('quality')!, numericFeatures(df, 'quality'), 0.3);
  }, {timeout: TIMEOUT});
}); // Model serialization
