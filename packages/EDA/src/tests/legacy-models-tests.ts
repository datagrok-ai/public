import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {cruxPredict} from '../crux/crux-model';
import {LEGACY_MODEL_ERROR, PREDICTION} from '../crux/model-common';
import {preparedLegacyModel, readLegacyModel} from '../legacy/legacy-models';
import {ubjsonToJson} from '../legacy/ubjson';
import {expectFeatureMatching, expectRejected, mulberry32, numericFeatures, predictColumn, sameValues} from './utils';

const FOLDER = 'System:AppData/Eda/legacy-models';
const TIMEOUT = 60000;

const FILES = ['linear-regression', 'pls-regression', 'softmax', 'xgboost', 'svm', 'xgboost-regression',
  'svm-regression', 'xgboost-bool', 'svm-bool', 'softmax-bool'];
/** Relative tolerance of FLOAT values where the prediction engine changed; every other fixture is exact */
const TOLERANCE: Record<string, {tol: number, reason: string}> = {
  'xgboost-regression': {tol: 1e-5, reason: 'crux XGBoost vs the old XGBoost wasm'},
  'svm-regression': {tol: 1e-5, reason: 'the TS LIBSVM port vs the old LIBSVM wasm'},
};

type Entry = {dataset: string, target: string, exclude: string[], type: string,
  prediction: (string | number | boolean | null)[]};
type Fixture = {bytes: Uint8Array, df: DG.DataFrame, features: DG.ColumnList, entry: Entry};

let reference: {fixtures: Record<string, Entry>};

async function fixture(file: string): Promise<Fixture> {
  const entry = reference.fixtures[file];
  const [df, bytes] = await Promise.all([grok.dapi.files.readCsv(entry.dataset),
    grok.dapi.files.readAsBytes(`${FOLDER}/${file}.bin`)]);
  return {bytes, df, entry, features: numericFeatures(df, entry.target, entry.exclude)};
}

function predict(features: DG.ColumnList, bytes: Uint8Array): Promise<DG.Column> {
  return predictColumn(DG.DataFrame.fromColumns(features.toList()), bytes);
}

/** `pred` equals the old prediction: exact, FLOAT within {@link TOLERANCE} where the engine changed */
function expectReference(file: string, entry: Entry, pred: DG.Column): void {
  expect(pred.name, PREDICTION, `${file}: column name`);
  expect(pred.type, entry.type, `${file}: column type`);
  expect(pred.length, entry.prediction.length, `${file}: rows`);
  const tolerance = TOLERANCE[file]?.tol ?? 0;
  const differ: number[] = [];
  let maxDiff = 0;
  for (let i = 0; i < pred.length; ++i) {
    const old = entry.prediction[i];
    const value = pred.isNone(i) ? null : pred.get(i);
    if (old === null || value === null || pred.type !== DG.COLUMN_TYPE.FLOAT) {
      if (old !== value)
        differ.push(i);
      continue;
    }
    const diff = Math.abs(value - Math.fround(old as number));
    maxDiff = Math.max(maxDiff, diff);
    if (diff > tolerance * Math.max(1, Math.abs(old as number)))
      differ.push(i);
  }
  const first = differ[0];
  expect(differ.length, 0, `${file}: ${differ.length} rows differ (max |diff| ${maxDiff}), first row ${first}: ` +
    `old ${entry.prediction[first]}, new ${first === undefined ? '' : pred.get(first)}`);
}

/** `[u32 size]...` followed by the parts */
function container(sizes: number[], ...parts: Uint8Array[]): Uint8Array {
  const head = 4 * sizes.length;
  const blob = new Uint8Array(head + parts.reduce((sum, part) => sum + part.length, 0));
  const view = new DataView(blob.buffer);
  for (const [i, size] of sizes.entries())
    view.setUint32(4 * i, size, true);
  let offset = head;
  for (const part of parts) {
    blob.set(part, offset);
    offset += part.length;
  }
  return blob;
}

class UbjsonWriter {
  bytes: number[] = [];

  private push(size: number, write: (view: DataView) => void): void {
    const view = new DataView(new ArrayBuffer(size));
    write(view);
    this.bytes.push(...new Uint8Array(view.buffer));
  }

  marker(...codes: number[]): this {
    this.bytes.push(...codes);
    return this;
  }

  int64(n: number): this {
    this.marker(0x4c);
    this.push(8, (v) => v.setBigInt64(0, BigInt(n)));
    return this;
  }

  key(text: string): this {
    this.int64(text.length);
    this.bytes.push(...new TextEncoder().encode(text));
    return this;
  }

  int32(n: number): this {
    this.push(4, (v) => v.setInt32(0, n));
    return this;
  }

  float32(x: number): this {
    this.push(4, (v) => v.setFloat32(0, x));
    return this;
  }

  float64(x: number): this {
    this.push(8, (v) => v.setFloat64(0, x));
    return this;
  }
}

category('Legacy models', () => {
  before(async () => {
    reference = JSON.parse(await grok.dapi.files.readAsText(`${FOLDER}/reference.json`));
  });

  for (const file of FILES) {
    test(file, async () => {
      const {bytes, features, entry} = await fixture(file);
      const out = await cruxPredict(DG.DataFrame.fromColumns(features.toList()), bytes);
      expect(out.columns.length, 1, `${file}: Prediction only`);
      expectReference(file, entry, out.columns.byIndex(0));
    }, {timeout: TIMEOUT});
  }

  test('PLS model without VIP', async () => {
    const {bytes, features, entry} = await fixture('pls-regression');
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    const modelSize = view.getUint32(0, true);
    const scoresSize = view.getUint32(4, true);
    const model = DG.DataFrame.fromByteArray(bytes.subarray(8, 8 + modelSize));
    expect(model.col('Variable Importance') !== null, true, 'the fixture has VIP');
    model.columns.remove('Variable Importance');
    const modelBytes = model.toByteArray();
    const blob = container([modelBytes.length, scoresSize], modelBytes,
      bytes.subarray(8 + modelSize, 8 + modelSize + scoresSize));

    const original = readLegacyModel(bytes);
    const stripped = readLegacyModel(blob);
    expect(stripped.kind, 'pls');
    expect(original.kind === 'pls' && stripped.kind === 'pls' && stripped.components === original.components, true,
      'components without VIP');
    expectReference('pls-regression', entry, await predict(features, blob));
  }, {timeout: TIMEOUT});

  test('Unknown blobs are rejected', async () => {
    const {features} = await fixture('svm');
    const rng = mulberry32(5);
    const softmax = (await fixture('softmax')).bytes;
    const softmaxSize = new DataView(softmax.buffer, softmax.byteOffset).getUint32(0, true);
    const v1Header = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Type', ['string']),
      DG.Column.fromInt32Array('Params count', Int32Array.of(3)),
      DG.Column.fromInt32Array('Categories size', Int32Array.of(0)),
      DG.Column.fromStrings('Categories', ['setosa']),
    ]).toByteArray();
    const blobs: [string, Uint8Array][] = [
      ['random', Uint8Array.from({length: 64}, () => Math.floor(rng() * 256))],
      ['bare Float32Array', new Uint8Array(Float32Array.from({length: 16}, (_, i) => i + 0.5).buffer)],
      ['frame with other columns', DG.DataFrame.fromColumns([DG.Column.fromStrings('a', ['x'])]).toByteArray()],
      ['XGBoost v1 header', container([v1Header.length], v1Header, new Uint8Array(64))],
      ['Softmax without feature names', softmax.slice(0, 4 + softmaxSize + 1)],
    ];
    await expectRejected(blobs.map(([label, blob]) => [label, () => predict(features, blob)]),
      (e) => e instanceof Error && e.message === LEGACY_MODEL_ERROR);
  }, {timeout: TIMEOUT});

  test('Prepared model is reused', async () => {
    const svm = await fixture('svm');
    const other = await fixture('svm-bool');
    await predict(svm.features, svm.bytes);
    const prepared = preparedLegacyModel(svm.bytes);
    await predict(svm.features, svm.bytes);
    expect(preparedLegacyModel(svm.bytes) === prepared, true, 'the second apply prepared the model again');
    await predict(other.features, other.bytes);
    await predict(svm.features, svm.bytes);
    expect(preparedLegacyModel(svm.bytes) !== prepared, true, 'another model must replace the cached one');
  }, {timeout: TIMEOUT});

  test('Features are matched by name', async () => {
    for (const file of ['xgboost', 'linear-regression']) {
      const {bytes, features} = await fixture(file);
      const baseline = await predict(features, bytes);
      await expectFeatureMatching((cols) => predict(cols, bytes), features, baseline, (a, b) => sameValues(a, b));
    }
  }, {timeout: TIMEOUT});

  test('Legacy PLS visualize', async () => {
    const {bytes, df} = await fixture('pls-regression');
    const price = df.col('price')!;
    const viewers = await grok.functions.call('Eda:visualizePLSRegression',
      {df, targetColumn: price, predictColumn: price, model: bytes});
    expect(viewers.length, 5, 'loadings, coefficients, VIP, explained variance, scores');
  }, {timeout: TIMEOUT});

  test('UBJSON to JSON', async () => {
    const w = new UbjsonWriter().marker(0x7b);
    w.key('a').marker(0x5b, 0x24, 0x6c, 0x23).int64(3).int32(1).int32(2).int32(3);
    w.key('b').marker(0x5b, 0x24, 0x64, 0x23).int64(2).float32(1.5).float32(2);
    w.key('c').marker(0x5b, 0x23).int64(2).marker(0x5a, 0x54);
    w.key('d').marker(0x53).int64(1).marker(0x78);
    w.key('e').marker(0x44).float64(0.25);
    w.key('f').marker(0x69, 0xfd);
    w.key('g').marker(0x55, 200);
    w.key('h').int64(2 ** 40);
    w.key('n').marker(0x5a);
    w.marker(0x7d);

    const json = ubjsonToJson(Uint8Array.from(w.bytes));
    expect(JSON.stringify(JSON.parse(json)), JSON.stringify({a: [1, 2, 3], b: [1.5, 2], c: [null, true], d: 'x',
      e: 0.25, f: -3, g: 200, h: 1099511627776, n: null}));
    expect(json.includes('2.0'), true, `a whole float keeps its fraction: ${json}`);
    expect(json.includes('"h":1099511627776'), true, `an integer has none: ${json}`);
    await expectExceptionAsync(async () => {
      ubjsonToJson(Uint8Array.of(0x48, 0x4c, 0, 0, 0, 0, 0, 0, 0, 1, 0x31));
    }, (e) => e.message.includes('unsupported marker'));
    await expectExceptionAsync(async () => {
      ubjsonToJson(Uint8Array.of(0x53, 0x4c, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff));
    }, (e) => e.message.includes('out of range'));
  });
});
