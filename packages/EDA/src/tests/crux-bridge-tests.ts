import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectArray, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {tableFromIPC} from 'apache-arrow';

import {cruxMl} from '../crux/crux-loader';
import {toArrowIpc} from '../crux/crux-bridge';
import {LEGACY_MODEL_ERROR, cruxFit, cruxPredict, cruxTerminate, readCruxManifest, runCrux} from '../crux/crux-model';
import {GAP_ROWS, IRIS_SPECIES, accuracy, classificationDataset, copyCol, decoyColumns, expectProbabilities,
  irisMeasurements, mulberry32, numericFeatures, predictColumn, r2, regressionDataset} from './utils';

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
    for (const name of ['ml.pls_analysis', 'ml.pca'])
      expect(JSON.parse(await cruxMl.fnDescribe(name))['meta.fn'], name);
  }, {timeout: FIRST_TIMEOUT});

  test('Column types', async () => {
    const {features, y} = regressionData();
    const raw = features[1].getRawData();
    features[1] = DG.Column.fromInt32Array(features[1].name,
      Int32Array.from({length: ROWS}, (_, i) => Math.round(raw[i] * 100)));

    const bytes = await cruxFit('regression', 'linear', frame(features), y);
    const pred = await predictColumn(frame(features), bytes);
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
    const pred = await predictColumn(frame(dirty), bytes);
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
    const pred = await predictColumn(frame(features), bytes);
    expect(pred.type, DG.COLUMN_TYPE.BOOL);
    const acc = accuracy(target, pred);
    expect(acc > 0.95, true, `accuracy ${acc}`);
    expectArray(readCruxManifest(bytes).classes!, [false, true]);
  }, {timeout: TIMEOUT});

  test('Missing string target', async () => {
    const df = classificationDataset(ROWS, 3, true, 42);
    const features = numericFeatures(df, 'Labels').toList();
    const target = df.col('Labels')!.clone();
    target.set(0, null);
    target.set(1, null);

    const y = tableFromIPC(toArrowIpc(features, target)).getChild('y')!;
    expect(y.nullCount, 2, 'missing labels must go as Arrow nulls');
    expect(y.isValid(0) || y.isValid(1), false, 'missing labels must go as Arrow nulls');

    const manifest = readCruxManifest(await cruxFit('classification', 'logistic', frame(features), target));
    expect(manifest.training.n_rows, ROWS - 2, 'rows with a missing target are skipped');
    expect((manifest.classes as string[]).includes(''), false, `classes ${JSON.stringify(manifest.classes)}`);
  }, {timeout: TIMEOUT});

  test('String target', async () => {
    const {features, target, bytes} = await irisModel();
    const pred = await predictColumn(frame(features), bytes);
    expect(pred.type, DG.COLUMN_TYPE.STRING);
    const acc = accuracy(target, pred);
    expect(acc > 0.9, true, `accuracy ${acc}`);

    const manifest = readCruxManifest(bytes);
    expect(manifest.pipeline![0].step, 'standard_scaler');
    expectArray(manifest.classes!, Array.from(new Set(target.toList())).sort());
  }, {timeout: TIMEOUT});

  test('Class probabilities', async () => {
    const {features, bytes} = await irisModel();
    await expectProbabilities(frame(features), bytes, IRIS_SPECIES);

    const df = classificationDataset(300, 3, true, 42);
    const flags = numericFeatures(df, 'Labels').toList();
    const target = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'flag',
      df.col('Labels')!.toList().map((l: string) => l[0] === 'A'));
    await expectProbabilities(frame(flags), await cruxFit('classification', 'logistic', frame(flags), target),
      ['false', 'true']);

    const regression = await cruxFit('regression', 'linear', frame(features.slice(1)), features[0]);
    expectArray((await cruxPredict(frame(features), regression)).columns.names(), ['Prediction']);

    const gapped = await expectProbabilities(irisMeasurements(await grok.dapi.files.readCsv(IRIS), GAP_ROWS), bytes,
      IRIS_SPECIES);
    for (const col of gapped.columns.toList()) {
      for (let i = 0; i < gapped.rowCount; ++i)
        expect(col.isNone(i), GAP_ROWS.includes(i), `${col.name}, row ${i}`);
    }
  }, {timeout: TIMEOUT});

  test('Class probabilities of rows the model cannot predict and of no rows', async () => {
    const {bytes} = await irisModel();
    const iris = await grok.dapi.files.readCsv(IRIS);
    const everyRow = Array.from({length: iris.rowCount}, (_, i) => i);
    const allGaps = await expectProbabilities(irisMeasurements(iris, everyRow), bytes, IRIS_SPECIES);
    for (const col of allGaps.columns.toList())
      expect(col.stats.missingValueCount, iris.rowCount, `${col.name} has a value`);

    const none = irisMeasurements(iris).clone(DG.BitSet.create(iris.rowCount, () => false));
    expect((await expectProbabilities(none, bytes, IRIS_SPECIES)).rowCount, 0);
  }, {timeout: TIMEOUT});

  test('Classes that differ only in case get no probability columns', async () => {
    const df = classificationDataset(200, 3, true, 42);
    const features = numericFeatures(df, 'Labels').toList();
    const target = DG.Column.fromStrings('case',
      df.col('Labels')!.toList().map((l: string) => l[0] === 'A' ? 'a' : 'A'));
    const bytes = await cruxFit('classification', 'logistic', frame(features), target);
    expectArray((await cruxPredict(frame(features), bytes)).columns.names(), ['Prediction']);
  }, {timeout: TIMEOUT});

  test('Save and read back', async () => {
    const {features, bytes} = await irisModel();
    const manifest = readCruxManifest(bytes);
    expect(manifest.format, 'crux-model-v1');
    expectArray(manifest.eda!.features, features.map((c) => c.name));
    expect(manifest.n_features, 4);

    const before = await predictColumn(frame(features), bytes);
    cruxTerminate();
    const reloaded = await predictColumn(frame(features), bytes);
    expect(accuracy(before, reloaded), 1);

    const other = await cruxFit('regression', 'linear', frame(features.slice(1)), features[0]);
    const otherPred = await predictColumn(frame(features), other);
    expect(otherPred.type, DG.COLUMN_TYPE.FLOAT, 'the other model');
    const score = r2(features[0], otherPred);
    expect(score > 0.8, true, `the other model: R2 ${score}`);
    const again = await predictColumn(frame(features), bytes);
    expect(accuracy(before, again), 1, 'after switching models');

    const reversed = await predictColumn(frame([...features].reverse()), bytes);
    expect(accuracy(before, reversed), 1, 'reversed column order');
    const superset = await predictColumn(frame([...decoyColumns(before.length, 2), ...features]), bytes);
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
    const a = await predictColumn(frame(features), first);
    const b = await predictColumn(frame(features), second);
    expect(accuracy(a, b), 1);
  }, {timeout: TIMEOUT});

  test('Legacy blob is rejected', async () => {
    const {features} = regressionData();
    await expectExceptionAsync(async () => {
      await cruxPredict(frame(features), Uint8Array.of(1, 2, 3));
    }, (e) => e.message === LEGACY_MODEL_ERROR);
  }, {timeout: TIMEOUT});

  test('No feature names', async () => {
    const {features, y} = regressionData();
    const out = await runCrux('ml.fit', toArrowIpc(features, y),
      {task: 'regression', model: 'linear', featurizer: 'vector', cv: 0, estimator: '{}'});
    const bytes = out.getChild('model')!.get(0) as Uint8Array;
    expect(readCruxManifest(bytes).eda === undefined, true, 'eda.features is still in the blob');
    await expectExceptionAsync(async () => {
      await cruxPredict(frame(features), bytes);
    }, (e) => e.message.includes('has no feature names'));
  }, {timeout: TIMEOUT});

  test('Failed model registration is retried', async () => {
    const {features, bytes} = await irisModel();
    const module = cruxMl as unknown as {post: (op: string, args: object, transfer?: Transferable[]) => Promise<any>};
    const post = module.post;
    let registrations = 0;
    module.post = function(op: string, args: object, transfer?: Transferable[]) {
      if (op === 'fnResource' && ++registrations === 1)
        return Promise.reject(new Error('registration failed'));
      return post.call(this, op, args, transfer);
    };
    try {
      cruxTerminate();
      await expectExceptionAsync(async () => {
        await cruxPredict(frame(features), bytes);
      }, (e) => e.message === 'registration failed');
      const pred = await predictColumn(frame(features), bytes);
      expect(registrations, 2, 'the model must be registered again after a failed registration');
      expect(pred.length, features[0].length);
    } finally {
      delete (module as {post?: unknown}).post;
    }
  }, {timeout: TIMEOUT});

  test('Live-model cap', async () => {
    const module = cruxMl as unknown as {terminate: () => void};
    const terminate = module.terminate;
    let terminations = 0;
    module.terminate = function() {
      ++terminations;
      terminate.call(this);
    };
    try {
      const rng = mulberry32(11);
      for (let k = 0; k < 33; ++k) {
        const x = DG.Column.fromFloat64Array('x', Float64Array.from({length: 10}, () => rng()));
        const y = DG.Column.fromFloat64Array('y', Float64Array.from(x.getRawData().subarray(0, 10), (v) => k + 2 * v));
        const pred = await predictColumn(frame([x]), await cruxFit('regression', 'linear', frame([x]), y));
        for (let i = 0; i < pred.length; ++i)
          expect(Number.isFinite(pred.get(i)), true, `model ${k}, row ${i}: ${pred.get(i)}`);
      }
      expect(terminations > 0, true, 'the worker was not restarted after 32 models');
    } finally {
      delete (module as {terminate?: unknown}).terminate;
    }
  }, {timeout: TIMEOUT});

  test('Arrow packing time', async () => {
    const lines: string[] = [];
    for (const [rows, cols] of [[100000, 10], [100000, 100], [20000, 500]]) {
      const df = grok.data.demo.randomWalk(rows, cols + 1);
      const all = df.columns.toList();
      const features = DG.DataFrame.fromColumns(all.slice(0, cols));
      const y = all[cols];

      let t0 = performance.now();
      toArrowIpc(features.columns.toList(), y);
      const packMs = performance.now() - t0;

      t0 = performance.now();
      await cruxFit('regression', 'xgboost', features, y, {estimator: {n_estimators: 20}});
      const xgbMs = performance.now() - t0;

      t0 = performance.now();
      await cruxFit('regression', 'linear', features, y);
      const linearMs = performance.now() - t0;

      lines.push(`BENCH pack ${rows / 1000}kx${cols} pack=${packMs.toFixed(0)}ms xgb-fit=${xgbMs.toFixed(0)}ms ` +
        `linear-fit=${linearMs.toFixed(0)}ms pack-share=${(100 * packMs / linearMs).toFixed(1)}%/` +
        `${(100 * packMs / xgbMs).toFixed(1)}%`);
    }
    for (const line of lines)
      console.log(line);
  }, {timeout: 600000, benchmark: true});
});
