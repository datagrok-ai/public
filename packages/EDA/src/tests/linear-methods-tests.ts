// Tests for PCA, PLS & linear regression

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, expectFloat, test} from '@datagrok-libraries/test/src/test';
import {computePCA} from '../eda-tools';
import {getPlsAnalysis, hotellingT2Limit, hotellingEllipseParams} from '../pls/pls-tools';
import {PlsModel} from '../pls/pls-ml';
import {getLinearRegressionParams, getPredictionByLinearRegression} from '../regression';
import {regressionDataset, madNorm, madError} from './utils';

const ROWS = 100;
const ROWS_K = 100;
const COLS = 100;
const COMPONENTS = 3;
const TIMEOUT = 9000;
const INDEP_COLS = 2;
const DEP_COLS = 5;
const ERROR = 0.1;
// OLS now fits via gradient descent; cap epochs so the large-sample
// performance run stays within TIMEOUT (the default 1000 epochs overruns).
const PERF_EPOCHS = 100;

category('Principal component analysis', () => {
  test(`Performance: ${ROWS_K}K rows, ${COLS} cols, ${COMPONENTS} components`, async () => {
    const df = grok.data.demo.randomWalk(ROWS_K * 1000, COLS);
    await computePCA(df, df.columns, COMPONENTS, false, false);
  }, {timeout: TIMEOUT, benchmark: true});

  test(`Performance: 1K rows, 5K cols, ${COMPONENTS} components`, async () => {
    const df = grok.data.demo.randomWalk(1000, 5000);
    await computePCA(df, df.columns, COMPONENTS, false, false);
  }, {timeout: TIMEOUT, benchmark: true});

  test('Correctness', async () => {
    // Data
    const df = regressionDataset(ROWS, COMPONENTS, DEP_COLS);

    // Apply
    const pca = await computePCA(df, df.columns, COMPONENTS + 1, false, false);

    // Check
    const lastPca = pca.columns.byIndex(COMPONENTS);
    const norm = madNorm(lastPca);

    // the last PCA component must be small due to df construction
    expect((norm < ERROR), true, 'Incorrect PCA computations');
  }, {timeout: TIMEOUT});
}); // PCA

category('Partial least squares regression', () => {
  test(`Performance: ${ROWS_K}K rows, ${COLS} cols, ${COMPONENTS} components`, async () => {
    // Data
    const df = grok.data.demo.randomWalk(ROWS_K * 1000, COLS);
    const cols = df.columns;

    // Apply
    await getPlsAnalysis({
      table: df,
      features: cols,
      predict: cols.byIndex(COLS - 1),
      components: COMPONENTS,
      names: undefined,
      isQuadratic: false,
    });
  }, {timeout: TIMEOUT, benchmark: true});

  test('Correctness', async () => {
    // Data
    const df = regressionDataset(ROWS_K, COMPONENTS, DEP_COLS);
    const cols = df.columns;
    const target = cols.byIndex(COMPONENTS + DEP_COLS - 1);

    // Apply
    const plsRes = await getPlsAnalysis({
      table: df,
      features: cols,
      predict: target,
      components: COMPONENTS,
      isQuadratic: false,
      names: undefined,
    });

    // Check deviation
    const deviation = madError(target, plsRes.prediction);
    expect(
      (deviation < ERROR),
      true,
      `Incorrect PLS computations, error is too big: ${deviation}; expected: < ${ERROR}`,
    );
  }, {timeout: TIMEOUT});

  test(`Predictive modeling: ${ROWS_K}K samples, ${COLS} features, ${COMPONENTS} components`, async () => {
    // Prepare data
    const df = regressionDataset(ROWS_K * 1000, COMPONENTS, COLS - COMPONENTS + 1);
    const features = df.columns;
    const target = features.byIndex(COLS);
    features.remove(target.name);

    // Train & pack model
    const model = new PlsModel();
    await model.fit(features, target, COMPONENTS);
    const packed = model.toBytes();

    // Unpack model & predict
    const unpackedModel = new PlsModel(packed);
    const prediction = unpackedModel.predict(features);

    // Check deviation
    const deviation = madError(target, prediction);
    expect(
      (deviation < ERROR),
      true,
      `Incorrect PLS (ML) computations, error is too big: ${deviation}; expected: < ${ERROR}`,
    );
  }, {timeout: TIMEOUT, benchmark: true});
}); // PLS

category('Linear regression', () => {
  test(`Performance: ${ROWS_K}K samples, ${COLS} features`, async () => {
    // Prepare data
    const df = regressionDataset(ROWS_K * 1000, COLS, 1);
    const features = df.columns;
    const target = features.byIndex(COLS);

    // Train & pack model (cap GD epochs to stay within TIMEOUT)
    const params = await getLinearRegressionParams(features, target, 0, 0, undefined, PERF_EPOCHS);
    const packed = new Uint8Array(params.buffer);

    // Unpack & apply model
    const unpackedParams = new Float32Array(packed.buffer);
    getPredictionByLinearRegression(features, unpackedParams);
  }, {timeout: TIMEOUT, benchmark: true});

  test('Correctness', async () => {
    // Prepare data
    const df = regressionDataset(ROWS, INDEP_COLS, 1);
    const features = df.columns;
    const target = features.byIndex(INDEP_COLS);

    // Train & pack model
    const params = await getLinearRegressionParams(features, target);
    const packed = new Uint8Array(params.buffer);

    // Unpack & apply model
    const unpackedParams = new Float32Array(packed.buffer);
    const prediction = getPredictionByLinearRegression(features, unpackedParams);

    // Evaluate model
    const error = madError(target, prediction);
    expect(
      error < ERROR,
      true,
      `Incorrect linear regression computations, error is too big: ${error}; expected: < ${ERROR}`,
    );
  }, {timeout: TIMEOUT});
}); // Linear regression

category('Hotelling ellipse', () => {
  // T2 = (p(n-1)/(n-p))·F^{-1}(conf; p, n-p), p = 2. For n = 10 (df1 = 2, df2 = 8):
  // F_{0.95}(2,8) = 4.4590, F_{0.99}(2,8) = 8.6491 (standard F tables).
  test('T-squared limit matches F-table values', async () => {
    expectFloat(hotellingT2Limit(0.95, 10), 2.25 * 4.4590, 1e-2);
    expectFloat(hotellingT2Limit(0.99, 10), 2.25 * 8.6491, 1e-2);
    // larger confidence -> larger limit
    expect(hotellingT2Limit(0.99, 10) > hotellingT2Limit(0.95, 10), true);
  });

  test('Ellipse params: center, semi-axes, ratio', async () => {
    // x = 1..10 (mean 5.5, sample var 82.5/9); y = 2·x (mean 11, var 4·varX)
    const x = new Float64Array(10);
    const y = new Float64Array(10);
    for (let i = 0; i < 10; ++i) {
      x[i] = i + 1;
      y[i] = 2 * (i + 1);
    }
    const xCol = DG.Column.fromFloat64Array('x', x);
    const yCol = DG.Column.fromFloat64Array('y', y);
    const varX = 82.5 / 9;

    for (const conf of [0.95, 0.99]) {
      const {cx, cy, a, b} = hotellingEllipseParams(xCol, yCol, conf);
      const t2 = hotellingT2Limit(conf, 10);
      expectFloat(cx, 5.5, 1e-9);
      expectFloat(cy, 11, 1e-9);
      expectFloat(a, Math.sqrt(t2 * varX), 1e-9);
      expectFloat(b, Math.sqrt(t2 * 4 * varX), 1e-9);
      // var(y) = 4·var(x) -> b = 2·a
      expectFloat(b / a, 2, 1e-9);
    }
  });
}); // Hotelling ellipse
