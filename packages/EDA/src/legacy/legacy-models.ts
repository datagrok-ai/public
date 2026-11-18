import * as DG from 'datagrok-api/dg';

import {LEGACY_MODEL_ERROR, PREDICTION, blobId, featureColumns, runCrux} from '../crux/model-common';
import {floatColumn, nullSentinel, toArrowIpc} from '../crux/crux-bridge';
import {LibsvmModel, parseLibsvmModel, svmPredict} from './libsvm';
import {ubjsonToJson} from './ubjson';

const SIZE_BYTES = 4;
const ALIGN = 4;
const MIN_BYTES = 16;
const BOOL_TRUE = 'true';
const BINARY_THRESHOLD = 0.5;
const AVGS = 'Avg-s';
const STDEVS = 'Stddev-s';
const OBJECTIVES = [0, 1, 2];

/** A pre-Crux EDA model (the containers with feature names); XGBoost and SVM keep their inner model bytes */
export type LegacyModel =
  {kind: 'linear', features: string[], params: Float32Array} |
  {kind: 'pls', features: string[], params: Float32Array, components: number} |
  {kind: 'softmax', features: string[], categories: string[], params: Float32Array[], avgs: Float32Array,
    stdevs: Float32Array, wasBool: boolean} |
  {kind: 'xgboost', features: string[], targetType: string, categories?: string[], wasBool: boolean,
    objective: number, model: Uint8Array} |
  {kind: 'svm', features: string[], targetType: string, categories?: string[], wasBool: boolean, avgs: Float32Array,
    stdevs: Float32Array, model: Uint8Array};

type SoftmaxModel = Extract<LegacyModel, {kind: 'softmax'}>;
type XgboostModel = Extract<LegacyModel, {kind: 'xgboost'}>;
type SvmModel = Extract<LegacyModel, {kind: 'svm'}>;

/** A legacy model ready to predict: the XGBoost model as JSON text, the SVM model parsed */
export type PreparedModel = Exclude<LegacyModel, {kind: 'xgboost' | 'svm'}> | (XgboostModel & {json: string}) |
  (SvmModel & {libsvm: LibsvmModel});

function legacyError(): Error {
  return new Error(LEGACY_MODEL_ERROR);
}

function frame(bytes: Uint8Array, offset: number, size: number): DG.DataFrame | null {
  if (size <= 0 || offset + size > bytes.length)
    return null;
  try {
    return DG.DataFrame.fromByteArray(bytes.subarray(offset, offset + size));
  } catch {
    return null;
  }
}

function u32(bytes: Uint8Array, offset: number): number {
  return offset + SIZE_BYTES > bytes.length ? 0 :
    new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength).getUint32(offset, true);
}

function names(df: DG.DataFrame | null, column: string): string[] {
  const col = df?.col(column) ?? null;
  if (col === null)
    throw legacyError();
  return col.toList() as string[];
}

function floats(df: DG.DataFrame, column: string): Float32Array {
  return new Float32Array(df.col(column)!.getRawData().subarray(0, df.rowCount));
}

function softmaxModel(bytes: Uint8Array, df: DG.DataFrame, size: number): SoftmaxModel {
  const classes = df.columns.toList().filter((col) => col.name !== AVGS && col.name !== STDEVS);
  const flag = SIZE_BYTES + size;
  const features = names(frame(bytes, flag + 1 + SIZE_BYTES, u32(bytes, flag + 1)), 'Features');
  if (classes.length < 1 || classes.some((col) => col.type !== DG.COLUMN_TYPE.FLOAT) ||
    df.rowCount !== features.length + 1)
    throw legacyError();
  return {kind: 'softmax', features, categories: classes.map((col) => col.name),
    params: classes.map((col) => floats(df, col.name)), avgs: floats(df, AVGS), stdevs: floats(df, STDEVS),
    wasBool: bytes[flag] === 1};
}

function containerModel(bytes: Uint8Array, header: DG.DataFrame, size: number): XgboostModel | SvmModel {
  const int = (name: string) => header.col(name) === null ? 0 : header.get(name, 0) as number;
  const targetType = header.get('Type', 0) as string;
  const wasBool = int('Was bool') === 1;
  let offset = SIZE_BYTES + size;

  const categoriesSize = int('Categories size');
  const categories = categoriesSize > 0 ? names(frame(bytes, offset, categoriesSize), 'Categories') : undefined;
  offset += categoriesSize;

  const featuresSize = int('Features size');
  const features = names(frame(bytes, offset, featuresSize), 'Features');
  offset += featuresSize;

  const statsSize = int('Stats size');
  const stats = frame(bytes, offset, statsSize);
  offset = Math.ceil((offset + statsSize) / ALIGN) * ALIGN;

  const modelSize = int('Model size');
  if (modelSize <= 0 || offset + modelSize > bytes.length)
    throw legacyError();
  const model = bytes.subarray(offset, offset + modelSize);

  if (header.columns.contains('Objective')) {
    const objective = int('Objective');
    if (!OBJECTIVES.includes(objective) || (objective !== 0 && categories === undefined))
      throw legacyError();
    return {kind: 'xgboost', features, targetType, categories, wasBool, objective, model};
  }
  if (stats === null || !stats.columns.contains(AVGS) || !stats.columns.contains(STDEVS) ||
    (targetType === DG.COLUMN_TYPE.STRING && categories === undefined))
    throw legacyError();
  const end = model.indexOf(0);
  return {kind: 'svm', features, targetType, categories, wasBool, avgs: floats(stats, AVGS),
    stdevs: floats(stats, STDEVS), model: end < 0 ? model : model.subarray(0, end)};
}

/** Recognises the containers of the pre-Crux Linear Regression, PLS Regression, Softmax, XGBoost (v2) and SVM; any
 * other blob, and every blob without feature names, throws {@link LEGACY_MODEL_ERROR} */
export function readLegacyModel(bytes: Uint8Array): LegacyModel {
  if (bytes.length < MIN_BYTES)
    throw legacyError();
  const size = u32(bytes, 0);

  const sized = frame(bytes, SIZE_BYTES, size);
  if (sized !== null) {
    const columns = sized.columns;
    if (columns.contains(AVGS) && columns.contains(STDEVS))
      return softmaxModel(bytes, sized, size);
    if (columns.contains('Type') && columns.contains('Model size') && columns.contains('Version'))
      return containerModel(bytes, sized, size);
    throw legacyError();
  }

  const pls = frame(bytes, 2 * SIZE_BYTES, size);
  if (pls !== null && pls.columns.contains('Feature names') && pls.columns.contains('Regression Coefficients')) {
    const features = names(pls, 'Feature names').slice(0, -1);
    const components = pls.columns.length - 2 - (pls.columns.contains('Variable Importance') ? 1 : 0);
    if (components < 1 || pls.rowCount !== features.length + 1)
      throw legacyError();
    return {kind: 'pls', features, params: floats(pls, 'Regression Coefficients'), components};
  }

  const linear = frame(bytes, 0, bytes.length);
  if (linear !== null && linear.columns.contains('features') && linear.columns.contains('values')) {
    const features = names(linear, 'features').slice(0, -1);
    if (linear.rowCount !== features.length + 1)
      throw legacyError();
    return {kind: 'linear', features, params: floats(linear, 'values')};
  }

  throw legacyError();
}

function prepare(bytes: Uint8Array): PreparedModel {
  const model = readLegacyModel(bytes);
  try {
    if (model.kind === 'xgboost')
      return {...model, json: ubjsonToJson(model.model)};
    if (model.kind === 'svm')
      return {...model, libsvm: parseLibsvmModel(new TextDecoder().decode(model.model))};
  } catch {
    throw legacyError();
  }
  return model;
}

let lastPrepared: {id: string, model: PreparedModel} | null = null;

/** The prepared model of a legacy blob, cached for the last blob */
export function preparedLegacyModel(bytes: Uint8Array): PreparedModel {
  const id = blobId(bytes);
  if (lastPrepared?.id !== id)
    lastPrepared = {id, model: prepare(bytes)};
  return lastPrepared.model;
}

function linearPrediction(cols: DG.Column[], params: Float32Array): DG.Column {
  const featuresCount = cols.length;
  const samplesCount = cols[0].length;
  const prediction = new Float32Array(samplesCount);

  let rawData = cols[0].getRawData();
  const bias = params[featuresCount];
  let weight = params[0];

  for (let i = 0; i < samplesCount; ++i)
    prediction[i] = bias + weight * rawData[i];

  for (let j = 1; j < featuresCount; ++j) {
    rawData = cols[j].getRawData();
    weight = params[j];

    for (let i = 0; i < samplesCount; ++i)
      prediction[i] += weight * rawData[i];
  }

  return DG.Column.fromFloat32Array(PREDICTION, prediction, samplesCount);
}

/** One f32 row per sample, each feature standardised with the stored stats (0 when the stdev is not positive) */
function standardizedRows(cols: DG.Column[], avgs: Float32Array, stdevs: Float32Array): Float32Array[] {
  const m = cols[0].length;
  const n = cols.length;

  const X = new Array<Float32Array>(m);
  for (let i = 0; i < m; ++i)
    X[i] = new Float32Array(n);

  for (let j = 0; j < n; ++j) {
    const raw = cols[j].getRawData();
    const avg = avgs[j];
    const stdev = stdevs[j];

    if (stdev > 0) {
      for (let i = 0; i < m; ++i)
        X[i][j] = (raw[i] - avg) / stdev;
    } else {
      for (let i = 0; i < m; ++i)
        X[i][j] = 0;
    }
  }

  return X;
}

function softmaxPrediction(cols: DG.Column[], model: SoftmaxModel): DG.Column {
  const X = standardizedRows(cols, model.avgs, model.stdevs);
  const m = X.length;
  const n = cols.length;
  const c = model.categories.length;
  const Z = new Float32Array(c);
  const argMaxes = new Int32Array(m);

  for (let j = 0; j < m; ++j) {
    const xBuf = X[j];

    for (let i = 0; i < c; ++i) {
      const wBuf = model.params[i];
      let sum = wBuf[n];

      for (let k = 0; k < n; ++k)
        sum += wBuf[k] * xBuf[k];

      Z[i] = Math.exp(sum);
    }

    let max = Z[0];
    let argMax = 0;

    for (let k = 1; k < c; ++k) {
      if (max < Z[k]) {
        max = Z[k];
        argMax = k;
      }
    }

    argMaxes[j] = argMax;
  }

  return classColumn(model.categories, model.wasBool, m, (i) => argMaxes[i]);
}

/** The class of each row by its category index (clamped as the old code did); null where `index` is null */
function classColumn(categories: string[], wasBool: boolean, rows: number, index: (i: number) => number | null):
  DG.Column {
  const last = categories.length - 1;
  const labels = Array.from({length: rows}, (_, i) => {
    const k = index(i);
    const label = k === null ? null : categories[Math.max(0, Math.min(k, last))];
    return wasBool && label !== null ? label === BOOL_TRUE : label;
  });
  return DG.Column.fromList(wasBool ? DG.COLUMN_TYPE.BOOL : DG.COLUMN_TYPE.STRING, PREDICTION, labels);
}

/** Regression values; INT / BIG_INT targets are rounded as the old XGBoost and SVM did */
function regressionColumn(values: DG.Column, targetType: string): DG.Column {
  const n = values.length;
  const raw = (values.getRawData() as Float32Array).subarray(0, n);
  if (targetType === DG.COLUMN_TYPE.INT) {
    return DG.Column.fromInt32Array(PREDICTION,
      Int32Array.from(raw, (v) => v === DG.FLOAT_NULL ? DG.INT_NULL : Math.round(v)), n);
  }
  if (targetType === DG.COLUMN_TYPE.BIG_INT) {
    return DG.Column.fromList(DG.COLUMN_TYPE.BIG_INT, PREDICTION,
      Array.from(raw, (v) => v === DG.FLOAT_NULL ? null : BigInt(Math.round(v))));
  }
  return values;
}

async function xgboostPrediction(cols: DG.Column[], model: XgboostModel & {json: string}): Promise<DG.Column> {
  const out = await runCrux('ml.xgboost_predict', toArrowIpc(cols), {model: model.json});
  const values = floatColumn(out, 'prediction', PREDICTION);
  if (model.objective === 0)
    return regressionColumn(values, model.targetType);
  const raw = values.getRawData() as Float32Array;
  const binary = model.objective === 1;
  return classColumn(model.categories!, model.wasBool, values.length, (i) => raw[i] === DG.FLOAT_NULL ? null :
    binary ? (raw[i] >= BINARY_THRESHOLD ? 1 : 0) : Math.round(raw[i]));
}

function svmPrediction(cols: DG.Column[], model: SvmModel & {libsvm: LibsvmModel}): DG.Column {
  const X = standardizedRows(cols, model.avgs, model.stdevs);
  const nRows = X.length;
  const prediction = new Float32Array(nRows);
  for (let i = 0; i < nRows; ++i)
    prediction[i] = svmPredict(model.libsvm, X[i]);

  return model.targetType === DG.COLUMN_TYPE.STRING ?
    classColumn(model.categories!, model.wasBool, nRows, (i) => Math.round(prediction[i])) :
    regressionColumn(DG.Column.fromFloat32Array(PREDICTION, prediction, nRows), model.targetType);
}

/** Applies a pre-Crux model blob; features are matched by the stored names; one `Prediction` column */
export async function legacyPredict(df: DG.DataFrame, bytes: Uint8Array): Promise<DG.DataFrame> {
  const model = preparedLegacyModel(bytes);
  const columns = featureColumns(df, model.features);
  for (const col of columns)
    nullSentinel(col);

  let prediction: DG.Column;
  switch (model.kind) {
  case 'linear':
  case 'pls':
    prediction = linearPrediction(columns, model.params);
    break;
  case 'softmax':
    prediction = softmaxPrediction(columns, model);
    break;
  case 'xgboost':
    prediction = await xgboostPrediction(columns, model);
    break;
  case 'svm':
    prediction = svmPrediction(columns, model);
    break;
  }
  return DG.DataFrame.fromColumns([prediction]);
}
