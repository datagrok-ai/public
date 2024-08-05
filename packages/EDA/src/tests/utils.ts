import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Return test dataframe with the specified number of dependent columns */
export function testDf(samples: number, features: number, dependent: number): DG.DataFrame {
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
export function madError(col1: DG.Column, col2: DG.Column): number {
  let mad = 0;
  const rows = col1.length;

  if (rows !== col2.length)
    throw new Error(`Error compuation failed, non-equal elements count: col1 - ${col1.length}, col2 - ${col2.length}`);

  const raw1 = col1.getRawData();
  const raw2 = col2.getRawData();

  for (let i = 0; i < rows; ++i)
    mad = Math.max(mad, Math.abs(raw1[i] - raw2[i]));

  return mad;
}
