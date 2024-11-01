// Tests for PCA, PLS & linear regression

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {computePCA} from '../eda-tools';
import {getPlsAnalysis} from '../pls/pls-tools';
import {PlsModel} from '../pls/pls-ml';
import {getLinearRegressionParams, getPredictionByLinearRegression} from '../regression';
import {regressionDataset, madNorm, madError} from './utils';

const ROWS = 100;
const ROWS_K = 100;
const COLS = 100;
const COMPONENTS = 3;
const TIMEOUT = 8000;
const INDEP_COLS = 2;
const DEP_COLS = 5;
const ERROR = 0.1;

category('Principal component analysis', () => {
  test(`Performance: ${ROWS_K}K rows, ${COLS} cols, ${COMPONENTS} components`, async () => {
    const df = grok.data.demo.randomWalk(ROWS_K * 1000, COLS);
    await computePCA(df, df.columns, COMPONENTS, false, false, true);
  }, {timeout: TIMEOUT, benchmark: true});

  test('Correctness', async () => {
    // Data
    const df = regressionDataset(ROWS, COMPONENTS, DEP_COLS);

    // Apply
    const pca = await computePCA(df, df.columns, COMPONENTS + 1, false, false, true);

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

    // Train & pack model
    const params = await getLinearRegressionParams(features, target);
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
    const error = madError(prediction, prediction);
    expect(
      error < ERROR,
      true,
      `Incorrect linear regression computations, error is too big: ${error}; expected: < ${ERROR}`,
    );
  }, {timeout: TIMEOUT});
}); // Linear regression
