import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';
import {findColumn} from '../utils/column-detection';
import {GroupAssignment} from '../analysis/experiment-setup';

/** Returns log2-prefixed intensity column names. */
export function getIntensityColumns(df: DG.DataFrame): string[] {
  return df.columns.toList()
    .filter((c) => c.semType === SEMTYPE.INTENSITY && c.name.startsWith('log2('))
    .map((c) => c.name);
}

/** Computes the mean of non-null values across columns at a given row index.
 *  Returns NaN if all values are null. */
function groupMean(cols: DG.Column[], rowIdx: number): number {
  let sum = 0;
  let count = 0;
  for (const col of cols) {
    if (!col.isNone(rowIdx)) {
      sum += col.get(rowIdx) as number;
      count++;
    }
  }
  return count > 0 ? sum / count : NaN;
}

/** Ensures a column exists fresh -- removes it if present, then adds a new float column. */
function ensureFreshFloat(df: DG.DataFrame, name: string): DG.Column {
  if (df.columns.contains(name))
    df.columns.remove(name);
  return df.columns.addNewFloat(name);
}

/** Computes MA values (M = log2 ratio, A = average intensity) and adds M, A columns to df. */
export function computeMA(df: DG.DataFrame, groups: GroupAssignment): void {
  const g1Cols = groups.group1.columns.map((n) => df.col(n)!);
  const g2Cols = groups.group2.columns.map((n) => df.col(n)!);

  const nRows = df.rowCount;
  const mArr = new Float32Array(nRows);
  const aArr = new Float32Array(nRows);
  for (let i = 0; i < nRows; i++) {
    const g1Mean = groupMean(g1Cols, i);
    const g2Mean = groupMean(g2Cols, i);
    if (isNaN(g1Mean) || isNaN(g2Mean)) {
      mArr[i] = DG.FLOAT_NULL;
      aArr[i] = DG.FLOAT_NULL;
    } else {
      mArr[i] = g2Mean - g1Mean;
      aArr[i] = (g1Mean + g2Mean) / 2;
    }
  }
  ensureFreshFloat(df, 'M').init((i) => mArr[i]);
  ensureFreshFloat(df, 'A').init((i) => aArr[i]);
}

/** Computes a moving-average trend line for the MA plot.
 *  Approximates loess by sorting on A and averaging M within a sliding window.
 *  Adds an MA_trend column to df. */
export function computeLoessTrend(df: DG.DataFrame, windowFraction: number = 0.1): void {
  const aCol = df.col('A');
  const mCol = df.col('M');
  if (!aCol || !mCol) {
    ensureFreshFloat(df, 'MA_trend');
    return;
  }

  const aRaw = aCol.getRawData() as Float32Array | Float64Array;
  const mRaw = mCol.getRawData() as Float32Array | Float64Array;

  // Collect non-null (A, M) pairs with original row indices
  const points: {idx: number; a: number; m: number}[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    if (aRaw[i] !== DG.FLOAT_NULL && mRaw[i] !== DG.FLOAT_NULL)
      points.push({idx: i, a: aRaw[i], m: mRaw[i]});
  }

  const trendArr = new Float32Array(df.rowCount);
  trendArr.fill(DG.FLOAT_NULL);

  // Guard: need at least 3 non-null points
  if (points.length < 3) {
    ensureFreshFloat(df, 'MA_trend').init((i) => trendArr[i]);
    return;
  }

  // Sort by A value
  points.sort((a, b) => a.a - b.a);

  const halfWindow = Math.max(Math.round(points.length * windowFraction / 2), 5);

  // For each point, compute mean M in a window centered on it, writing to trendArr
  for (let i = 0; i < points.length; i++) {
    const lo = Math.max(0, i - halfWindow);
    const hi = Math.min(points.length - 1, i + halfWindow);
    let sum = 0;
    let count = 0;
    for (let j = lo; j <= hi; j++) {
      sum += points[j].m;
      count++;
    }
    trendArr[points[i].idx] = sum / count;
  }

  ensureFreshFloat(df, 'MA_trend').init((i) => trendArr[i]);
}

/** Computes coefficient of variation (CV = sd/mean) on raw (non-log) intensities
 *  within a group of columns. Adds cvColName and meanColName float columns to df. */
export function computeCV(
  df: DG.DataFrame,
  groupCols: string[],
  cvColName: string,
  meanColName: string,
): void {
  const rawArrays = groupCols
    .map((n) => df.col(n))
    .filter((c) => c != null && c.type === DG.COLUMN_TYPE.FLOAT)
    .map((c) => c!.getRawData() as Float32Array | Float64Array);

  const nRows = df.rowCount;
  const cvArr = new Float32Array(nRows);
  const meanArr = new Float32Array(nRows);

  // Welford's online algorithm for mean + sample variance.
  // The naive `(sumSq - sum*sum/count) / (count-1)` form suffers from
  // catastrophic cancellation when the values are large and the spread
  // is small — which is exactly the regime for raw intensities (often >1e6
  // with tight biological CVs).
  for (let i = 0; i < nRows; i++) {
    let count = 0;
    let mean = 0;
    let m2 = 0;
    for (const raw of rawArrays) {
      const v = raw[i];
      if (v !== DG.FLOAT_NULL) {
        // Exponentiate log2 values back to raw intensities for CV computation
        const val = Math.pow(2, v);
        count++;
        const delta = val - mean;
        mean += delta / count;
        m2 += delta * (val - mean);
      }
    }
    if (count < 1) {
      cvArr[i] = DG.FLOAT_NULL;
      meanArr[i] = DG.FLOAT_NULL;
    } else if (count < 2) {
      cvArr[i] = DG.FLOAT_NULL;
      meanArr[i] = mean;
    } else {
      const variance = m2 / (count - 1);
      const sd = Math.sqrt(Math.max(0, variance));
      cvArr[i] = mean > 0 ? sd / mean : DG.FLOAT_NULL;
      meanArr[i] = mean;
    }
  }

  ensureFreshFloat(df, cvColName).init((i) => cvArr[i]);
  ensureFreshFloat(df, meanColName).init((i) => meanArr[i]);
}

/** Creates a binary missingness matrix DataFrame (1 = present, 0 = missing).
 *  Includes protein ID column if available. */
export function createMissingnessMatrix(df: DG.DataFrame, intensityCols: string[]): DG.DataFrame {
  const columns: DG.Column[] = [];

  // Include protein ID column if available
  const labelCol = findColumn(df, SEMTYPE.PROTEIN_ID, ['protein id', 'accession']);
  if (labelCol)
    columns.push(labelCol.clone());

  // Binary columns: 1 = present, 0 = missing.
  // Write directly into a fresh Int32Array, then wrap as the column.
  for (const colName of intensityCols) {
    const srcCol = df.col(colName)!;
    const srcRaw = srcCol.type === DG.COLUMN_TYPE.FLOAT
      ? srcCol.getRawData() as Float32Array | Float64Array
      : null;
    const data = new Int32Array(df.rowCount);
    if (srcRaw) {
      for (let i = 0; i < df.rowCount; i++)
        data[i] = srcRaw[i] === DG.FLOAT_NULL ? 0 : 1;
    } else {
      for (let i = 0; i < df.rowCount; i++)
        data[i] = srcCol.isNone(i) ? 0 : 1;
    }
    columns.push(DG.Column.fromInt32Array(colName, data));
  }

  const result = DG.DataFrame.fromColumns(columns);
  result.name = 'Missing Values';
  return result;
}

/** Creates a long-format DataFrame with columns: ProteinId, Sample, Intensity.
 *  Skips null intensity values to keep the DataFrame compact. */
export function unpivotIntensities(df: DG.DataFrame, intensityCols: string[]): DG.DataFrame {
  const proteinIdCol = findColumn(df, SEMTYPE.PROTEIN_ID, ['protein id', 'accession']);

  const ids: string[] = [];
  const samples: string[] = [];
  const intensities: number[] = [];

  for (let i = 0; i < df.rowCount; i++) {
    const proteinId = proteinIdCol ? (proteinIdCol.get(i) as string ?? String(i)) : String(i);
    for (const colName of intensityCols) {
      const col = df.col(colName)!;
      if (!col.isNone(i)) {
        ids.push(proteinId);
        samples.push(colName);
        intensities.push(col.get(i) as number);
      }
    }
  }

  const result = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('ProteinId', ids),
    DG.Column.fromStrings('Sample', samples),
    DG.Column.fromFloat32Array('Intensity', new Float32Array(intensities)),
  ]);
  result.name = 'Intensity Distributions';
  return result;
}

/** Creates a sample-level DataFrame with columns: Sample, MissingPct, Group.
 *  Each row represents one intensity column with its missing value percentage. */
export function computeMissingBarData(
  df: DG.DataFrame,
  intensityCols: string[],
  groups: GroupAssignment,
): DG.DataFrame {
  const sampleNames: string[] = [];
  const missingPcts: number[] = [];
  const groupNames: string[] = [];

  const g1Set = new Set(groups.group1.columns);
  const g2Set = new Set(groups.group2.columns);

  for (const colName of intensityCols) {
    const col = df.col(colName)!;
    let nullCount = 0;
    for (let i = 0; i < df.rowCount; i++) {
      if (col.isNone(i))
        nullCount++;
    }

    sampleNames.push(colName);
    missingPcts.push((nullCount / df.rowCount) * 100);

    if (g1Set.has(colName))
      groupNames.push(groups.group1.name);
    else if (g2Set.has(colName))
      groupNames.push(groups.group2.name);
    else
      groupNames.push('Other');
  }

  const result = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Sample', sampleNames),
    DG.Column.fromFloat32Array('MissingPct', new Float32Array(missingPcts)),
    DG.Column.fromStrings('Group', groupNames),
  ]);
  result.name = 'Missing Values per Sample';
  return result;
}
