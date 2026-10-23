// Serialization round-trip + apply-time feature matching, on real datasets.
// See expectFeatureMatching for the order/subset/missing checks.

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {isClassTarget} from '../crux/methods';
import {trainPlsRegression} from '../pls/pls-ml';
import {nonNullRows} from '../utils';
import {Fit, FITS} from './baseline';
import {accuracy, numericFeatures, copyCol, decoyColumns, asTable, r2, predictColumn, sameValues} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';
const CARS = 'System:DemoFiles/cars.csv';
const WINE = 'System:DemoFiles/winequality.csv';
const DEMOG = 'System:DemoFiles/demog.csv';

const TIMEOUT = 60000;
// Effectively zero: identical params → identical Float32.
const ROUNDTRIP_EPS = 1e-6;

// ── data helpers ──────────────────────────────────────────────────────

/** Up to `n` evenly-strided rows (a prefix would drop a class in ordered data). */
function sampleRows(df: DG.DataFrame, n: number): DG.DataFrame {
  if (df.rowCount <= n)
    return df;
  const stride = Math.ceil(df.rowCount / n);
  return df.clone(DG.BitSet.create(df.rowCount, (i) => i % stride === 0));
}

// ── metrics ───────────────────────────────────────────────────────────

type Predict = (cols: DG.ColumnList) => Promise<DG.Column>;
type Same = (a: DG.Column, b: DG.Column) => boolean;

/** Apply-time contract: order- and extra-column-invariant; missing feature throws. */
async function expectFeatureMatching(predict: Predict, features: DG.ColumnList, before: DG.Column,
  same: Same): Promise<void> {
  const cols = features.toList();
  const rows = cols[0].length;

  // (1) reversed order
  expect(same(before, await predict(asTable(cols.map(copyCol).reverse()))), true,
    'prediction depends on column order (reversed)');

  // (2) rotated order — a permutation distinct from the reverse (>2 features)
  if (cols.length > 2) {
    const copies = cols.map(copyCol);
    expect(same(before, await predict(asTable([...copies.slice(1), copies[0]]))), true,
      'prediction depends on column order (rotated)');
  }

  // (3) superset: pick features by name from a scrambled table with extra columns.
  const decoys = decoyColumns(rows, 2);
  const superset = [decoys[0], ...cols.map(copyCol).reverse()];
  superset.splice(1 + Math.floor(cols.length / 2), 0, decoys[1]);
  expect(same(before, await predict(asTable(superset))), true, 'extra columns in the table disturb the prediction');

  // (4) a missing feature must throw, not silently mispredict
  await expectExceptionAsync(async () => {
    await predict(asTable(cols.slice(1).map(copyCol)));
  }, (e) => e.message.includes('is missing'));
}

/** Trains on `features`, applies twice (identical), checks feature matching and the training score floor */
async function roundTrip(label: string, fit: Fit, target: DG.Column, features: DG.ColumnList,
  floor: number): Promise<void> {
  const classification = isClassTarget(target);
  // Classification: identical labels; regression: identical Float32 values (a lossless model reproduces them)
  const same: Same = (a, b) => sameValues(a, b, classification ? 0 : ROUNDTRIP_EPS);
  const bytes = await fit(DG.DataFrame.fromColumns(features.toList()), target);
  const predict: Predict = (cols) => predictColumn(DG.DataFrame.fromColumns(cols.toList()), bytes);

  const before = await predict(features);
  expect(same(before, await predict(features)), true, `${label}: a second apply changed predictions`);
  await expectFeatureMatching(predict, features, before, same);

  const score = classification ? accuracy(target, before) : r2(target, before);
  expect(score > floor, true, `${label}: ${classification ? 'accuracy' : 'R²'} too low: ${score} ` +
    `(expected > ${floor})`);
}

const {softmax, linreg, pls, svm, xgboost} = FITS;
const pls2: Fit = (df, target) => trainPlsRegression(df, target, 2);

// SMO is superlinear in the sample count; cap large datasets for the SVM tests.
const SVM_MAX_ROWS = 1500;

// ── datasets ─────────────────────────────────────────────────────────

/** demog rows without nulls in the columns used, as {features, target}. */
async function demog(featureNames: string[], targetName: string, maxRows = Number.MAX_SAFE_INTEGER):
  Promise<{features: DG.ColumnList, target: DG.Column}> {
  const all = await grok.dapi.files.readCsv(DEMOG);
  const clean = all.clone(nonNullRows(all.columns.byNames([...featureNames, targetName])));
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
    await roundTrip('softmax', softmax, df.col('Species')!, numericFeatures(df, 'Species', ['col 1']), 0.8);
  }, {timeout: TIMEOUT});

  test('Softmax — winequality', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await roundTrip('softmax', softmax, df.col('type')!, numericFeatures(df, 'type', ['quality']), 0.8);
  }, {timeout: TIMEOUT});

  test('Softmax — demog (SEX)', async () => {
    const {features, target} = await demog(DEMOG_FEATS, 'SEX');
    await roundTrip('softmax', softmax, target, features, 0.55);
  }, {timeout: TIMEOUT});

  // ── linear regression ──
  test('Linear regression — cars', async () => {
    const df = await grok.dapi.files.readCsv(CARS);
    await roundTrip('linreg', linreg, df.col('price')!, numericFeatures(df, 'price'), 0.4);
  }, {timeout: TIMEOUT});

  test('Linear regression — winequality', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await roundTrip('linreg', linreg, df.col('quality')!, numericFeatures(df, 'quality'), 0.1);
  }, {timeout: TIMEOUT});

  test('Linear regression — demog (WEIGHT)', async () => {
    const {features, target} = await demog(['AGE', 'HEIGHT'], 'WEIGHT');
    await roundTrip('linreg', linreg, target, features, 0.2);
  }, {timeout: TIMEOUT});

  // ── PLS ──
  test('PLS — cars', async () => {
    const df = await grok.dapi.files.readCsv(CARS);
    await roundTrip('PLS', pls, df.col('price')!, numericFeatures(df, 'price'), 0.4);
  }, {timeout: TIMEOUT});

  test('PLS — winequality', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await roundTrip('PLS', pls, df.col('quality')!, numericFeatures(df, 'quality'), 0.1);
  }, {timeout: TIMEOUT});

  test('PLS — demog (WEIGHT)', async () => {
    const {features, target} = await demog(['AGE', 'HEIGHT'], 'WEIGHT');
    await roundTrip('PLS', pls2, target, features, 0.2);
  }, {timeout: TIMEOUT});

  // ── SVM (classification) ──
  test('SVM — iris', async () => {
    const df = await grok.dapi.files.readCsv(IRIS);
    await roundTrip('SVM', svm, df.col('Species')!, numericFeatures(df, 'Species', ['col 1']), 0.8);
  }, {timeout: TIMEOUT});

  test('SVM — winequality (type)', async () => {
    const df = sampleRows(await grok.dapi.files.readCsv(WINE), SVM_MAX_ROWS);
    await roundTrip('SVM', svm, df.col('type')!, numericFeatures(df, 'type', ['quality']), 0.8);
  }, {timeout: TIMEOUT});

  test('SVM — demog (SEX)', async () => {
    const {features, target} = await demog(DEMOG_FEATS, 'SEX', SVM_MAX_ROWS);
    await roundTrip('SVM', svm, target, features, 0.55);
  }, {timeout: TIMEOUT});

  // ── SVM (regression) ──
  test('SVM — winequality (quality)', async () => {
    const df = sampleRows(await grok.dapi.files.readCsv(WINE), SVM_MAX_ROWS);
    await roundTrip('SVM', svm, df.col('quality')!, numericFeatures(df, 'quality'), 0.3);
  }, {timeout: TIMEOUT});

  // SVR defaults (C=1) need a small-magnitude target; cars price (~10^4) underfits.
  test('SVM — iris (Petal.Length)', async () => {
    const df = await grok.dapi.files.readCsv(IRIS);
    await roundTrip('SVM', svm, df.col('Petal.Length')!, numericFeatures(df, 'Petal.Length', ['col 1']), 0.5);
  }, {timeout: TIMEOUT});

  // ── XGBoost (classification) ──
  test('XGBoost — iris', async () => {
    const df = await grok.dapi.files.readCsv(IRIS);
    await roundTrip('XGBoost', xgboost, df.col('Species')!, numericFeatures(df, 'Species', ['col 1']), 0.8);
  }, {timeout: TIMEOUT});

  test('XGBoost — winequality (type)', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await roundTrip('XGBoost', xgboost, df.col('type')!, numericFeatures(df, 'type', ['quality']), 0.8);
  }, {timeout: TIMEOUT});

  test('XGBoost — demog (SEX)', async () => {
    const {features, target} = await demog(DEMOG_FEATS, 'SEX');
    await roundTrip('XGBoost', xgboost, target, features, 0.55);
  }, {timeout: TIMEOUT});

  // ── XGBoost (regression) ──
  test('XGBoost — cars (price)', async () => {
    const df = await grok.dapi.files.readCsv(CARS);
    await roundTrip('XGBoost', xgboost, df.col('price')!, numericFeatures(df, 'price'), 0.4);
  }, {timeout: TIMEOUT});

  test('XGBoost — winequality (quality)', async () => {
    const df = await grok.dapi.files.readCsv(WINE);
    await roundTrip('XGBoost', xgboost, df.col('quality')!, numericFeatures(df, 'quality'), 0.3);
  }, {timeout: TIMEOUT});
}); // Model serialization
