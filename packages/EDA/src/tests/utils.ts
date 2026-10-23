import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {cruxPredict} from '../crux/crux-model';

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
