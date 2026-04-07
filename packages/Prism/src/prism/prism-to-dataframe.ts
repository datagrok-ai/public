import * as DG from 'datagrok-api/dg';

import {PrismSheet, PrismAnalysis} from './prism-types';


/** Number of Y columns per dataset for each data format. */
const FORMAT_COLS_PER_DATASET: Record<string, number> = {
  'y_single': 1,
  'y_replicates': -1, // special: uses replicatesCount
  'y_sd': 2,
  'y_se': 2,
  'y_cv': 2,
  'y_sd_n': 3,
  'y_se_n': 3,
  'y_cv_n': 3,
  'y_plus_minus': 3,
  'y_high_low': 3,
};

/** Column suffixes for each data format. */
const FORMAT_SUFFIXES: Record<string, string[]> = {
  'y_sd': ['Mean', 'SD'],
  'y_se': ['Mean', 'SEM'],
  'y_cv': ['Mean', '%CV'],
  'y_sd_n': ['Mean', 'SD', 'N'],
  'y_se_n': ['Mean', 'SEM', 'N'],
  'y_cv_n': ['Mean', '%CV', 'N'],
  'y_plus_minus': ['Mean', '+Error', '-Error'],
  'y_high_low': ['Mean', 'Upper Limit', 'Lower Limit'],
};

/** Parses a cell value as a float, returning NaN for empty/invalid values. */
function parseValue(s: string): number {
  if (s === '' || s === undefined || s === null)
    return NaN;
  const v = parseFloat(s);
  return v;
}

/** Returns the number of X columns in the CSV (0 or 1, before the Y data). */
function getXColumnCount(sheet: PrismSheet): number {
  return (sheet.xTitle != null || sheet.rowTitlesPresent) ? 1 : 0;
}

/** Computes the expected number of Y data columns in the CSV. */
function getExpectedYCols(sheet: PrismSheet): number {
  const format = sheet.dataFormat;
  const colsPerDs = FORMAT_COLS_PER_DATASET[format] ?? 1;
  const dsCount = Math.max(sheet.yColumnTitles.length, 1);

  if (format === 'y_replicates')
    return dsCount * sheet.replicatesCount;

  return dsCount * colsPerDs;
}

/** Generates column names for Y data based on format and dataset titles. */
function generateYColumnNames(sheet: PrismSheet): string[] {
  const format = sheet.dataFormat;
  const titles = sheet.yColumnTitles.length > 0
    ? sheet.yColumnTitles
    : ['Values'];
  const names: string[] = [];

  if (format === 'y_replicates') {
    for (const title of titles) {
      if (sheet.replicatesCount === 1)
        names.push(title || 'Values');
      else {
        for (let i = 1; i <= sheet.replicatesCount; i++)
          names.push(`${title || 'Values'}_${i}`);
      }
    }
  }
  else if (format === 'y_single') {
    for (const title of titles)
      names.push(title || 'Values');
  }
  else {
    const suffixes = FORMAT_SUFFIXES[format] ?? ['Values'];
    for (const title of titles) {
      for (const suffix of suffixes)
        names.push(`${title || 'Values'}_${suffix}`);
    }
  }

  return names;
}


/** Converts a PrismSheet to a Datagrok DataFrame. */
export function prismSheetToDataFrame(sheet: PrismSheet): DG.DataFrame {
  const rowCount = sheet.data.length;
  const df = DG.DataFrame.create(rowCount);
  df.name = sheet.title || 'Data';

  const xColCount = getXColumnCount(sheet);
  const yColNames = generateYColumnNames(sheet);

  // Add X / row titles column
  if (xColCount > 0) {
    const xName = sheet.xTitle || (sheet.rowTitlesPresent ? 'Row' : 'X');
    const xValues = sheet.data.map((row) => row[0] ?? '');

    // Try to parse as numeric
    const numericValues = xValues.map(parseValue);
    const allNumeric = numericValues.every((v, i) => !isNaN(v) || xValues[i] === '');

    if (allNumeric && !sheet.rowTitlesPresent) {
      const col = df.columns.addNewFloat(xName);
      for (let i = 0; i < rowCount; i++)
        col.set(i, isNaN(numericValues[i]) ? DG.FLOAT_NULL : numericValues[i]);
    }
    else {
      const col = df.columns.addNewString(xName);
      for (let i = 0; i < rowCount; i++)
        col.set(i, xValues[i]);
    }
  }

  // Add Y data columns
  for (let ci = 0; ci < yColNames.length; ci++) {
    const csvIdx = xColCount + ci;
    const col = df.columns.addNewFloat(yColNames[ci]);
    for (let ri = 0; ri < rowCount; ri++) {
      const row = sheet.data[ri];
      if (csvIdx < row.length) {
        const v = parseValue(row[csvIdx]);
        col.set(ri, isNaN(v) ? DG.FLOAT_NULL : v);
      }
      else
        col.set(ri, DG.FLOAT_NULL);
    }
  }

  // Apply transformations for special formats
  applyFormatTransformations(df, sheet);

  return df;
}

/** Applies post-parse transformations for %CV and high/low formats. */
function applyFormatTransformations(df: DG.DataFrame, sheet: PrismSheet): void {
  const format = sheet.dataFormat;

  if (format === 'y_cv' || format === 'y_cv_n') {
    // %CV columns store raw ratio; transform to percentage: %CV = raw / mean * 100
    for (let ci = 0; ci < df.columns.length; ci++) {
      const col = df.columns.byIndex(ci);
      if (col.name.endsWith('_%CV') || col.name.endsWith('_%CV')) {
        const meanColIdx = ci - 1;
        if (meanColIdx >= 0) {
          const meanCol = df.columns.byIndex(meanColIdx);
          for (let ri = 0; ri < df.rowCount; ri++) {
            const raw = col.get(ri);
            const mean = meanCol.get(ri);
            if (raw !== DG.FLOAT_NULL && mean !== DG.FLOAT_NULL && mean !== 0)
              col.set(ri, (raw / mean) * 100);
          }
        }
      }
    }
  }

  if (format === 'y_high_low') {
    // Upper Limit = Mean + raw, Lower Limit = Mean - raw
    for (let ci = 0; ci < df.columns.length; ci++) {
      const col = df.columns.byIndex(ci);
      if (col.name.endsWith('_Upper Limit')) {
        const meanColIdx = ci - 1;
        if (meanColIdx >= 0) {
          const meanCol = df.columns.byIndex(meanColIdx);
          for (let ri = 0; ri < df.rowCount; ri++) {
            const raw = col.get(ri);
            const mean = meanCol.get(ri);
            if (raw !== DG.FLOAT_NULL && mean !== DG.FLOAT_NULL)
              col.set(ri, mean + raw);
          }
        }
      }
      else if (col.name.endsWith('_Lower Limit')) {
        const meanColIdx = ci - 2;
        if (meanColIdx >= 0) {
          const meanCol = df.columns.byIndex(meanColIdx);
          for (let ri = 0; ri < df.rowCount; ri++) {
            const raw = col.get(ri);
            const mean = meanCol.get(ri);
            if (raw !== DG.FLOAT_NULL && mean !== DG.FLOAT_NULL)
              col.set(ri, mean - raw);
          }
        }
      }
    }
  }
}


/** Converts a PrismAnalysis to a Datagrok DataFrame. */
export function prismAnalysisToDataFrame(analysis: PrismAnalysis): DG.DataFrame {
  const rowCount = analysis.resultData.length;
  const df = DG.DataFrame.create(rowCount);
  df.name = analysis.title || 'Analysis';

  if (rowCount === 0)
    return df;

  const colCount = analysis.resultData[0].length;

  // First column is typically row labels
  if (colCount > 0) {
    const col = df.columns.addNewString('Parameter');
    for (let ri = 0; ri < rowCount; ri++)
      col.set(ri, analysis.resultData[ri][0] ?? '');
  }

  // Remaining columns use titles from analysis metadata
  for (let ci = 1; ci < colCount; ci++) {
    const title = (ci - 1 < analysis.columnTitles.length)
      ? (analysis.columnTitles[ci - 1] || `Col ${ci}`)
      : `Col ${ci}`;

    // Try numeric first
    const values = analysis.resultData.map((row) => row[ci] ?? '');
    const numValues = values.map(parseValue);
    const allNumeric = numValues.every((v, i) => !isNaN(v) || values[i] === '');

    if (allNumeric) {
      const col = df.columns.addNewFloat(title);
      for (let ri = 0; ri < rowCount; ri++)
        col.set(ri, isNaN(numValues[ri]) ? DG.FLOAT_NULL : numValues[ri]);
    }
    else {
      const col = df.columns.addNewString(title);
      for (let ri = 0; ri < rowCount; ri++)
        col.set(ri, values[ri]);
    }
  }

  return df;
}
