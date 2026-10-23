// Tests for PCA, PLS & linear regression

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectExceptionAsync, expectFloat, test} from '@datagrok-libraries/test/src/test';
import {computePCA} from '../eda-tools';
import {getPlsAnalysis, hotellingT2Limit, hotellingEllipseParams} from '../pls/pls-tools';
import {trainPlsRegression} from '../pls/pls-ml';
import {trainLinearRegression} from '../regression';
import {readCruxManifest} from '../crux/crux-model';
import {nonNullRows} from '../utils';
import {regressionDataset, madNorm, madError, predictColumn, splitTarget} from './utils';

const ROWS = 100;
const ROWS_K = 100;
const COLS = 100;
const COMPONENTS = 3;
const TIMEOUT = 9000;
const CRUX_TIMEOUT = 60000;
const INDEP_COLS = 2;
const DEP_COLS = 5;
const ERROR = 0.1;

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
  }, {timeout: CRUX_TIMEOUT});

  test('Center and scale flags', async () => {
    const iris = await grok.dapi.files.readCsv('System:DemoFiles/iris.csv');
    const cols = DG.DataFrame.fromColumns(['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width']
      .map((name) => iris.col(name)!)).columns;

    for (const col of (await computePCA(iris, cols, 4, true, true)).columns) {
      const {avg, stdev} = col.stats;
      expect(Math.abs(avg) < 1e-5, true, `TT: component ${col.name} mean ${avg}`);
      expect(Math.abs(stdev - 1) < 1e-4, true, `TT: component ${col.name} stdev ${stdev}`);
    }

    for (const col of (await computePCA(iris, cols, 4, false, true)).columns) {
      const {stdev} = col.stats;
      expect(Math.abs(stdev - 1) < 1e-4, true, `FT: component ${col.name} stdev ${stdev}`);
    }

    const {stdev} = (await computePCA(iris, cols, 4, false, false)).columns.byIndex(0).stats;
    expect(Math.abs(stdev - 1) > 1e-2, true, `FF: component 1 stdev ${stdev}`);
  }, {timeout: CRUX_TIMEOUT});
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
  }, {timeout: CRUX_TIMEOUT});

  test(`Predictive modeling: ${ROWS_K}K samples, ${COLS} features, ${COMPONENTS} components`, async () => {
    const df = regressionDataset(ROWS_K * 1000, COMPONENTS, COLS - COMPONENTS + 1);
    const target = splitTarget(df, df.columns.byIndex(COLS).name);

    const bytes = await trainPlsRegression(df, target, COMPONENTS);
    const deviation = madError(target, await predictColumn(df, bytes));
    expect(
      (deviation < ERROR),
      true,
      `Incorrect PLS (ML) computations, error is too big: ${deviation}; expected: < ${ERROR}`,
    );
  }, {timeout: CRUX_TIMEOUT, benchmark: true});

  test('Rank error', async () => {
    const v = [1, 2, 3, 4, 5, 6];
    const df = DG.DataFrame.fromColumns([1, 2, 3].map((m) =>
      DG.Column.fromFloat32Array(`x${m}`, Float32Array.from(v, (x) => m * x))));
    const y = DG.Column.fromFloat32Array('y', Float32Array.from([1, 3, 2, 5, 4, 6]));
    await expectExceptionAsync(async () => {
      await trainPlsRegression(df, y, 2);
    }, (e) => e.message.includes('only 1 of 2 PLS components'));
  }, {timeout: CRUX_TIMEOUT});

  test('Missing target values', async () => {
    const df = regressionDataset(ROWS, COMPONENTS, 1);
    const target = splitTarget(df, df.columns.byIndex(COMPONENTS).name);
    for (const i of [0, 7, 42, ROWS - 1])
      target.set(i, null);

    const pred = await predictColumn(df, await trainPlsRegression(df, target, COMPONENTS));
    const mask = nonNullRows([target]);
    const ref = await predictColumn(df,
      await trainPlsRegression(df.clone(mask), target.clone(mask), COMPONENTS));

    for (let i = 0; i < ROWS; ++i) {
      expect(!pred.isNone(i) && Number.isFinite(pred.get(i)), true, `row ${i}: no prediction`);
      expectFloat(pred.get(i), ref.get(i), 1e-9, `row ${i}`);
    }
  }, {timeout: CRUX_TIMEOUT});
}); // PLS

category('Linear regression', () => {
  test(`Performance: ${ROWS_K}K samples, ${COLS} features`, async () => {
    const df = regressionDataset(ROWS_K * 1000, COLS, 1);
    const target = splitTarget(df, df.columns.byIndex(COLS).name);

    const bytes = await trainLinearRegression(df, target, 0, 0, 1000);
    await predictColumn(df, bytes);
  }, {timeout: CRUX_TIMEOUT, benchmark: true});

  test('Correctness', async () => {
    const df = regressionDataset(ROWS, INDEP_COLS, 1);
    const target = splitTarget(df, df.columns.byIndex(INDEP_COLS).name);

    const bytes = await trainLinearRegression(df, target, 0, 0, 1000);
    const error = madError(target, await predictColumn(df, bytes));
    expect(
      error < ERROR,
      true,
      `Incorrect linear regression computations, error is too big: ${error}; expected: < ${ERROR}`,
    );
  }, {timeout: CRUX_TIMEOUT});

  test('Model choice', async () => {
    const df = regressionDataset(ROWS, INDEP_COLS, 1);
    const target = splitTarget(df, df.columns.byIndex(INDEP_COLS).name);

    const linear = readCruxManifest(await trainLinearRegression(df, target, 0, 0, 1000));
    const ridge = readCruxManifest(await trainLinearRegression(df, target, 0, 0.5, 1000));
    const enet = readCruxManifest(await trainLinearRegression(df, target, 0.3, 0.1, 1000));

    expect(linear.model, 'linear');
    expect(ridge.model, 'ridge');
    expect(ridge.params.alpha, 0.5);
    expect(enet.model, 'elasticnet');
    expectFloat(enet.params.l1_ratio as number, 0.75, 1e-9);
    expectFloat(enet.params.alpha as number, 0.4, 1e-9);
    for (const manifest of [linear, ridge, enet])
      expect(manifest.pipeline?.[0].step, 'standard_scaler', `${manifest.model}: no standard scaler`);
  }, {timeout: CRUX_TIMEOUT});
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
