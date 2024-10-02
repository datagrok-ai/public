import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const TRESHOLD = 0.5;
const SHIFT = 1;
const LIMIT = 2;

const MAX_INT = 10;
const MAX_FLOAT = 10;
const CATEGORIES = ['Alpha', 'Beta', 'Gamma', 'Delta'];

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

/** Return dataframe for testing classifiers */
export function classificationDataset(samples: number, features: number, useShift: boolean): DG.DataFrame {
  const labels = new Array<string>(samples);
  const raw = new Array<Float32Array>(features);

  for (let j = 0; j < features; ++j) {
    const arr = new Float32Array(samples);

    for (let i = 0; i < samples; ++i)
      arr[i] = Math.random();

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
