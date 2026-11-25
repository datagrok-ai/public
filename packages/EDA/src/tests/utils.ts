import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {expect, expectArray, expectExceptionAsync} from '@datagrok-libraries/test/src/test';

import {cruxPredict, readCruxManifest} from '../crux/crux-model';
import {isClassTarget} from '../crux/methods';
import {nonNullRows} from '../utils';
import type {Fit} from './baseline';

const TRESHOLD = 0.5;
const SHIFT = 1;
const LIMIT = 2;

const MAX_INT = 10;
const MAX_FLOAT = 10;
const CATEGORIES = ['Alpha', 'Beta', 'Gamma', 'Delta'];

/**
 * Deterministic PRNG (mulberry32). Returns a `Math.random`-compatible function
 * (0 <= x < 1). Pass a fixed seed to a dataset generator to make it reproducible
 * bit-for-bit and remove flakiness; omit it to keep the legacy random behaviour.
 */
export function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return function() {
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Check lengths of columns */
function checkLen(target: DG.Column, prediction: DG.Column): void {
  if (target.length !== prediction.length)
    throw new Error(`Non-equal elements count: ${target.length} vs. ${prediction.length}`);
}

/** Return dataframe for testing regression & linear methods */
export function regressionDataset(samples: number, features: number, dependent: number): DG.DataFrame {
  // create main features
  const df = grok.data.demo.randomWalk(samples, features);
  const cols = df.columns;

  const coefs = new Float32Array(features);

  // main features raw data
  const raw = new Array<Float32Array>(features);
  for (let j = 0; j < features; ++j)
    raw[j] = cols.byIndex(j).getRawData() as Float32Array;

  // dependent features raw data
  for (let j = 0; j < dependent; ++j) {
    const arr = new Float32Array(samples);

    // generate coefficients
    for (let k = 0; k < features; ++k)
      coefs[k] = Math.random();

    for (let i = 0; i < samples; ++i) {
      for (let k = 0; k < features; ++k)
        arr[i] += coefs[k] * raw[k][i];
    }

    cols.add(DG.Column.fromFloat32Array(`y${j}`, arr));
  }

  return df;
} // pcaTestDf

/** Seeded regression set: uniform features `x0`, `x1`, ... and the target `y = sum (j + 1) * xj` */
export function seededRegression(rows: number, features: number, seed: number): DG.DataFrame {
  const rng = mulberry32(seed);
  const raw = Array.from({length: features}, () => Float32Array.from({length: rows}, () => rng()));
  const y = Float32Array.from({length: rows}, (_, i) => raw.reduce((sum, arr, j) => sum + arr[i] * (j + 1), 0));
  return DG.DataFrame.fromColumns([...raw.map((arr, j) => DG.Column.fromFloat32Array(`x${j}`, arr)),
    DG.Column.fromFloat32Array('y', y)]);
}

/** Max absolute deviation of the column */
export function madNorm(col: DG.Column): number {
  let mad = 0;
  const rows = col.length;
  const raw = col.getRawData();

  for (let i = 0; i < rows; ++i)
    mad = Math.max(mad, Math.abs(raw[i]));

  return mad;
}

/** Max absolute deviation error */
export function madError(target: DG.Column, prediction: DG.Column): number {
  checkLen(target, prediction);

  let mad = 0;
  const rows = target.length;
  const raw1 = target.getRawData();
  const raw2 = prediction.getRawData();

  for (let i = 0; i < rows; ++i)
    mad = Math.max(mad, Math.abs(raw1[i] - raw2[i]));

  return mad;
}

/**
 * Return dataframe for testing classifiers.
 * Pass `seed` to generate a reproducible dataset (kills flakiness); omit it for
 * the legacy `Math.random` behaviour.
 */
export function classificationDataset(
  samples: number, features: number, useShift: boolean, seed?: number): DG.DataFrame {
  const rng = seed === undefined ? Math.random : mulberry32(seed);
  const labels = new Array<string>(samples);
  const raw = new Array<Float32Array>(features);

  for (let j = 0; j < features; ++j) {
    const arr = new Float32Array(samples);

    for (let i = 0; i < samples; ++i)
      arr[i] = rng();

    raw[j] = arr;
  }

  const df = DG.DataFrame.fromColumns(raw.map((arr, idx) => DG.Column.fromFloat32Array(`#${idx}`, arr)));

  for (let i = 0; i < samples; ++i)
    labels[i] = raw.slice(0, LIMIT).map((arr) => (arr[i] > TRESHOLD) ? 'A' : 'B').join('');

  df.columns.add(DG.Column.fromStrings('Labels', labels));

  if (useShift) {
    for (let j = 0; j < features; ++j) {
      for (let i = 0; i < samples; ++i)
        raw[j][i] += (raw[j][i] > 0 ? 1 : -1) * SHIFT;
    }
  }

  return df;
} // classificationDataset

/** Return accuracy */
export function accuracy(target: DG.Column, prediction: DG.Column): number {
  checkLen(target, prediction);

  let correctPredictions = 0;
  const rows = target.length;

  if (rows < 1)
    return 1;

  for (let i = 0; i < rows; ++i) {
    if (target.get(i) === prediction.get(i))
      ++correctPredictions;
  }

  return correctPredictions / rows;
}

/** Same type and values; numbers within `tolerance` when it is positive */
export function sameValues(a: DG.Column, b: DG.Column, tolerance: number = 0): boolean {
  if (a.type !== b.type || a.length !== b.length)
    return false;
  for (let i = 0; i < a.length; ++i) {
    if (tolerance > 0 ? !(Math.abs(a.get(i) - b.get(i)) < tolerance) : a.get(i) !== b.get(i))
      return false;
  }
  return true;
}

/** Coefficient of determination R² of `pred` against `target`. */
export function r2(target: DG.Column, pred: DG.Column): number {
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

/** Numeric feature columns of `df`, excluding the target and any `exclude` names. */
export function numericFeatures(df: DG.DataFrame, target: string, exclude: string[] = []): DG.ColumnList {
  const cols = df.columns.toList().filter((c) =>
    (c.type === DG.COLUMN_TYPE.INT || c.type === DG.COLUMN_TYPE.FLOAT) &&
    c.name !== target && !exclude.includes(c.name));
  return asTable(cols);
}

/** Independent column copy — won't reparent the original. */
export function copyCol(col: DG.Column): DG.Column {
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
export function decoyColumns(rows: number, count: number): DG.Column[] {
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
export function asTable(cols: DG.Column[]): DG.ColumnList {
  return DG.DataFrame.fromColumns(cols).columns;
}

/** Removes the column `name` from `df` (leaving the features) and returns it */
export function splitTarget(df: DG.DataFrame, name: string): DG.Column {
  const target = df.col(name)!;
  df.columns.remove(name);
  return target;
}

/** The `Prediction` column of a Crux model applied to `df` */
export async function predictColumn(df: DG.DataFrame, bytes: Uint8Array): Promise<DG.Column> {
  return (await cruxPredict(df, bytes)).columns.byIndex(0);
}

/** `cruxPredict` of a classifier: `Prediction` then one `proba_<class>` FLOAT column per `classes` entry; every
 * predicted row sums to 1 and its predicted class has the largest probability */
export async function expectProbabilities(df: DG.DataFrame, bytes: Uint8Array, classes: string[]):
  Promise<DG.DataFrame> {
  const out = await cruxPredict(df, bytes);
  expectArray(out.columns.names(), ['Prediction', ...classes.map((c) => `proba_${c}`)]);
  const prediction = out.columns.byIndex(0);
  const probabilities = out.columns.toList().slice(1);
  for (const col of probabilities) {
    expect(col.type, DG.COLUMN_TYPE.FLOAT, `${col.name} must be a FLOAT column`);
    expect(col.length, df.rowCount, `${col.name} length`);
  }

  for (let i = 0; i < df.rowCount; ++i) {
    if (prediction.isNone(i)) {
      for (const col of probabilities)
        expect(col.isNone(i), true, `row ${i}: ${col.name} must be null where the prediction is null`);
      continue;
    }
    const sum = probabilities.reduce((s, col) => s + col.get(i), 0);
    expect(Math.abs(sum - 1) < 1e-5, true, `row ${i}: the probabilities sum to ${sum}`);
    const predicted = out.col(`proba_${String(prediction.get(i))}`)!.get(i);
    for (const col of probabilities)
      expect(predicted >= col.get(i), true, `row ${i}: ${col.name} exceeds the predicted class's probability`);
  }
  return out;
}

export type Predict = (cols: DG.ColumnList) => Promise<DG.Column>;
export type Same = (a: DG.Column, b: DG.Column) => boolean;

// Effectively zero: identical params → identical Float32.
const ROUNDTRIP_EPS = 1e-6;

/** Apply-time contract: order- and extra-column-invariant; missing feature throws. */
export async function expectFeatureMatching(predict: Predict, features: DG.ColumnList, before: DG.Column,
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
export async function roundTrip(label: string, fit: Fit, target: DG.Column, features: DG.ColumnList,
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

const IRIS_MEASUREMENTS = ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width'];
/** The iris Species classes, in crux's (sorted) order */
export const IRIS_SPECIES = ['setosa', 'versicolor', 'virginica'];
export const GAP_ROWS = [3, 10];

/** Copies of the four iris measurements; `Sepal.Width` is null at `gaps` */
export function irisMeasurements(iris: DG.DataFrame, gaps: number[] = []): DG.DataFrame {
  const features = IRIS_MEASUREMENTS.map((name) => copyCol(iris.col(name)!));
  for (const row of gaps)
    features[1].set(row, null);
  return DG.DataFrame.fromColumns(features);
}

/** [fit] of Species by the iris measurements with gaps at `GAP_ROWS` trains on every row and predicts every row */
export async function expectGapsKept(fit: Fit, iris: DG.DataFrame): Promise<void> {
  const gapped = irisMeasurements(iris, GAP_ROWS);
  const bytes = await fit(gapped, iris.col('Species')!);
  expect(readCruxManifest(bytes).training.n_rows, 150, 'rows with gaps must be trained on');
  expect((await predictColumn(gapped, bytes)).stats.missingValueCount, 0, 'rows with gaps must be predicted');
}

/** Applies the model twice: the predictions must be identical */
export async function applyTwice(df: DG.DataFrame, bytes: Uint8Array, label: string): Promise<DG.Column> {
  const before = await predictColumn(df, bytes);
  expect(sameValues(before, await predictColumn(df, bytes)), true, `${label}: predictions differ`);
  return before;
}

export async function demogFrame(target: string, features: string[], nonNull: boolean):
  Promise<{df: DG.DataFrame, target: DG.Column}> {
  const demog = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
  const df = nonNull ? demog.clone(nonNullRows(demog.columns.byNames([...features, target]))) : demog;
  return {df: df.clone(null, features), target: df.col(target)!};
}

/** Every action must fail (with an error `check` accepts, when given); the message lists the ones that did not */
export async function expectRejected(cases: [string, () => Promise<unknown>][], check?: (e: any) => boolean):
  Promise<void> {
  const missed: string[] = [];
  for (const [label, action] of cases) {
    try {
      await action();
      missed.push(label);
    } catch (e) {
      if (check !== undefined && !check(e))
        missed.push(`${label} (${e instanceof Error ? e.message : e})`);
    }
  }
  expect(missed.length, 0, `not rejected: ${missed.join(', ')}`);
}

/** Return dataframe with missing values */
export function dataWithMissingVals(rows: number, intCols: number, floatCols: number,
  strCols: number, misValCount: number): {df: DG.DataFrame, misValsIds: Map<string, number[]>} {
  const catsCount = CATEGORIES.length;
  const cols = [];
  let idx = 0;

  const misValsIds = new Map<string, number[]>();

  for (let j = 0; j < intCols; ++j) {
    const arr = new Int32Array(rows);
    const name = `int #${j + 1}`;
    const indeces: number[] = [];

    for (let i = 0; i < rows; ++i)
      arr[i] = Math.floor(Math.random() * MAX_INT);

    for (let k = 0; k < misValCount; ++k) {
      idx = Math.floor(rows * Math.random());
      arr[idx] = DG.INT_NULL;
      indeces.push(idx);
    }

    cols.push(DG.Column.fromInt32Array(name, arr));
    misValsIds.set(name, indeces);
  }

  for (let j = 0; j < floatCols; ++j) {
    const arr = new Float32Array(rows);
    const name = `float #${j + 1}`;
    const indeces: number[] = [];

    for (let i = 0; i < rows; ++i)
      arr[i] = Math.random() * MAX_FLOAT;

    for (let k = 0; k < misValCount; ++k) {
      idx = Math.floor(rows * Math.random());
      arr[idx] = DG.FLOAT_NULL;
      indeces.push(idx);
    }

    cols.push(DG.Column.fromFloat32Array(name, arr));
    misValsIds.set(name, indeces);
  }

  for (let j = 0; j < strCols; ++j) {
    const arr = new Array<string>(rows);
    const name = `str #${j + 1}`;
    const indeces: number[] = [];

    for (let i = 0; i < rows; ++i)
      arr[i] = CATEGORIES[Math.floor(Math.random() * catsCount)];

    const col = DG.Column.fromStrings(name, arr);

    for (let k = 0; k < misValCount; ++k) {
      idx = Math.floor(rows * Math.random());
      col.set(idx, null);
      indeces.push(idx);
    }

    cols.push(col);
    misValsIds.set(name, indeces);
  }

  return {
    df: DG.DataFrame.fromColumns(cols),
    misValsIds: misValsIds,
  };
} // tableWithMissingVals
