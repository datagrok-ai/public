import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {tableFromIPC} from 'apache-arrow';

import {cruxMl} from '../crux/crux-loader';
import {toArrowIpc} from '../crux/crux-bridge';
import {LEGACY_MODEL_ERROR, cruxFit, cruxPredict, cruxTerminate, readCruxManifest} from '../crux/crux-model';
import {accuracy, classificationDataset, copyCol, decoyColumns, numericFeatures, r2,
  regressionDataset} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';
const ROWS = 200;
const FIRST_TIMEOUT = 120000;
const TIMEOUT = 60000;

function regressionData(): {features: DG.Column[], y: DG.Column} {
  const df = regressionDataset(ROWS, 3, 1);
  return {features: [0, 1, 2].map((j) => copyCol(df.columns.byIndex(j))), y: copyCol(df.col('y0')!)};
}

function frame(cols: DG.Column[]): DG.DataFrame {
  return DG.DataFrame.fromColumns(cols.map(copyCol));
}

async function predictionColumn(df: DG.DataFrame, bytes: Uint8Array): Promise<DG.Column> {
  return (await cruxPredict(df, bytes)).columns.byIndex(0);
}

async function irisModel(): Promise<{features: DG.Column[], target: DG.Column, bytes: Uint8Array}> {
  const df = await grok.dapi.files.readCsv(IRIS);
  const features = numericFeatures(df, 'Species', ['col 1']).toList();
  const target = df.col('Species')!;
  const bytes = await cruxFit('classification', 'logistic', frame(features), target, {preprocess: 'standard'});
  return {features, target, bytes};
}

category('Crux bridge', () => {
  test('Module loads', async () => {
    const spec = JSON.parse(await cruxMl.fnDescribe('ml.fit'));
    expect(spec['meta.fn'], 'ml.fit');
    expect(spec.inputs.some((p: {name: string}) => p.name === 'x'), true, 'ml.fit has no x input');
  }, {timeout: FIRST_TIMEOUT});

  test('Column types', async () => {
    const {features, y} = regressionData();
    const raw = features[1].getRawData();
    features[1] = DG.Column.fromInt32Array(features[1].name,
      Int32Array.from({length: ROWS}, (_, i) => Math.round(raw[i] * 100)));

    const bytes = await cruxFit('regression', 'linear', frame(features), y);
    const pred = await predictionColumn(frame(features), bytes);
    expect(pred.type, DG.COLUMN_TYPE.FLOAT);
    expect(pred.length, ROWS);
    for (let i = 0; i < ROWS; ++i)
      expect(Number.isFinite(pred.get(i)), true, `row ${i}: ${pred.get(i)}`);
    const score = r2(y, pred);
    expect(score > 0.9, true, `R2 ${score}`);

    const text = DG.Column.fromStrings('text', Array.from({length: ROWS}, (_, i) => `t${i}`));
    await expectExceptionAsync(async () => {
      await cruxFit('regression', 'linear', frame([...features, text]), y);
    }, (e) => e.message.includes('"text"'));
  }, {timeout: TIMEOUT});

  test('Null and NaN', async () => {
    const {features, y} = regressionData();
    const dirty = features.map(copyCol);
    dirty[0].set(5, null);
    dirty[1].set(7, NaN);

    const x = tableFromIPC(toArrowIpc(dirty)).getChild('x')!;
    expect(Number.isNaN(x.get(5).get(0)), true, `null cell is ${x.get(5).get(0)}`);
    expect(Number.isNaN(x.get(7).get(1)), true, `NaN cell is ${x.get(7).get(1)}`);

    await expectExceptionAsync(async () => {
      await cruxFit('regression', 'linear', frame(dirty), y);
    }, (e) => e.message.includes('NaN'));
    await cruxFit('regression', 'xgboost', frame(dirty), y, {estimator: {n_estimators: 5}});

    const bytes = await cruxFit('regression', 'linear', frame(features), y);
    const pred = await predictionColumn(frame(dirty), bytes);
    expect(pred.isNone(5) && pred.isNone(7), true, 'rows with missing features must have no prediction');
    for (let i = 0; i < ROWS; ++i) {
      if (i !== 5 && i !== 7)
        expect(pred.isNone(i), false, `row ${i} has no prediction`);
    }

    const yNull = copyCol(y);
    yNull.set(3, null);
    const manifest = readCruxManifest(await cruxFit('regression', 'linear', frame(features), yNull));
    expect(manifest.training.n_rows, ROWS - 1);
  }, {timeout: TIMEOUT});

  test('Boolean target', async () => {
    const df = classificationDataset(300, 3, true, 42);
    const features = numericFeatures(df, 'Labels').toList();
    const target = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'flag',
      df.col('Labels')!.toList().map((l: string) => l[0] === 'A'));

    const bytes = await cruxFit('classification', 'logistic', frame(features), target);
    const pred = await predictionColumn(frame(features), bytes);
    expect(pred.type, DG.COLUMN_TYPE.BOOL);
    const acc = accuracy(target, pred);
    expect(acc > 0.95, true, `accuracy ${acc}`);
    expectArray(readCruxManifest(bytes).classes!, [false, true]);
  }, {timeout: TIMEOUT});

  test('String target', async () => {
    const {features, target, bytes} = await irisModel();
    const pred = await predictionColumn(frame(features), bytes);
    expect(pred.type, DG.COLUMN_TYPE.STRING);
    const acc = accuracy(target, pred);
    expect(acc > 0.9, true, `accuracy ${acc}`);

    const manifest = readCruxManifest(bytes);
    expect(manifest.pipeline![0].step, 'standard_scaler');
    expectArray(manifest.classes!, Array.from(new Set(target.toList())).sort());
  }, {timeout: TIMEOUT});

  test('Save and read back', async () => {
    const {features, bytes} = await irisModel();
    const manifest = readCruxManifest(bytes);
    expect(manifest.format, 'crux-model-v1');
    expectArray(manifest.eda!.features, features.map((c) => c.name));
    expect(manifest.n_features, 4);

    const before = await predictionColumn(frame(features), bytes);
    cruxTerminate();
    const reloaded = await predictionColumn(frame(features), bytes);
    expect(accuracy(before, reloaded), 1);

    const other = await cruxFit('regression', 'linear', frame(features.slice(1)), features[0]);
    const otherPred = await predictionColumn(frame(features), other);
    expect(otherPred.type, DG.COLUMN_TYPE.FLOAT, 'the other model');
    const score = r2(features[0], otherPred);
    expect(score > 0.8, true, `the other model: R2 ${score}`);
    const again = await predictionColumn(frame(features), bytes);
    expect(accuracy(before, again), 1, 'after switching models');

    const reversed = await predictionColumn(frame([...features].reverse()), bytes);
    expect(accuracy(before, reversed), 1, 'reversed column order');
    const superset = await predictionColumn(frame([...decoyColumns(before.length, 2), ...features]), bytes);
    expect(accuracy(before, superset), 1, 'extra columns');

    await expectExceptionAsync(async () => {
      await cruxPredict(frame(features.slice(1)), bytes);
    }, (e) => e.message.includes('is missing'));
  }, {timeout: TIMEOUT});

  test('Seed', async () => {
    const df = classificationDataset(ROWS, 3, true, 1);
    const features = numericFeatures(df, 'Labels').toList();
    const target = df.col('Labels')!;
    const params = {seed: 7, estimator: {n_estimators: 5}};

    const first = await cruxFit('classification', 'xgboost', frame(features), target, params);
    const second = await cruxFit('classification', 'xgboost', frame(features), target, params);
    const manifest = readCruxManifest(first);
    expect(manifest.params.random_state, 7, 'estimator.random_state');
    const state = (manifest as unknown as {state: {params: [string, string][]}}).state;
    expect(state.params.some(([k, v]) => k === 'random_state' && v === '7'), true,
      `XGBoost params ${JSON.stringify(state.params)}`);
    const a = await predictionColumn(frame(features), first);
    const b = await predictionColumn(frame(features), second);
    expect(accuracy(a, b), 1);
  }, {timeout: TIMEOUT});

  test('Legacy blob is rejected', async () => {
    const {features} = regressionData();
    await expectExceptionAsync(async () => {
      await cruxPredict(frame(features), Uint8Array.of(1, 2, 3));
    }, (e) => e.message === LEGACY_MODEL_ERROR);
  }, {timeout: TIMEOUT});

  after(async () => cruxTerminate());
});
