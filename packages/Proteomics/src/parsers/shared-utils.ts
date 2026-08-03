import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';

/** Creates log2-transformed copies of the specified columns.
 * Original columns get SEMTYPE.INTENSITY. New log2 columns get SEMTYPE.INTENSITY
 * and names prefixed with 'log2(...)'.
 * Zero, negative, and null values produce DG.FLOAT_NULL in the log2 column. */
export function log2TransformColumns(df: DG.DataFrame, colNames: string[]): void {
  let nonFiniteCount = 0;
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    col.semType = SEMTYPE.INTENSITY;
    const log2Name = `log2(${name})`;
    const log2Col = df.columns.addNewFloat(log2Name);
    log2Col.init((i) => {
      if (col.isNone(i)) return DG.FLOAT_NULL;
      const val = Number(col.get(i));
      if (isNaN(val)) { nonFiniteCount++; return DG.FLOAT_NULL; }
      return val > 0 ? Math.log2(val) : DG.FLOAT_NULL;
    });
    log2Col.semType = SEMTYPE.INTENSITY;
  }
  if (nonFiniteCount > 0)
    console.warn(`log2TransformColumns: ${nonFiniteCount} non-numeric value(s) dropped to null`);
}

/** Copies intensity columns with log2() prefix for already-transformed data.
 * No math transformation -- just copies values into log2-named columns
 * for downstream compatibility. Sets both original and log2 column
 * semType = SEMTYPE.INTENSITY. */
export function copyAsLog2Columns(df: DG.DataFrame, colNames: string[]): void {
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    col.semType = SEMTYPE.INTENSITY;
    const log2Name = `log2(${name})`;
    const log2Col = df.columns.addNewFloat(log2Name);
    log2Col.init((i) => col.isNone(i) ? DG.FLOAT_NULL : Number(col.get(i)));
    log2Col.semType = SEMTYPE.INTENSITY;
  }
}

/** Creates a primary (first semicolon-delimited entry) column if semicolons
 * are detected in source data. Skips creation entirely if no semicolons found. */
export function addPrimaryColumnIfNeeded(
  df: DG.DataFrame, sourceColName: string, newName: string, semType: string,
): void {
  const srcCol = df.col(sourceColName);
  if (!srcCol) return;

  // Scan for semicolons -- skip creation if none found
  let hasSemicolon = false;
  for (let i = 0; i < df.rowCount && !hasSemicolon; i++) {
    const val = srcCol.get(i);
    if (typeof val === 'string' && val.includes(';'))
      hasSemicolon = true;
  }
  if (!hasSemicolon) return;

  const newCol = df.columns.addNewString(newName);
  newCol.init((i) => {
    const val = srcCol.get(i);
    if (typeof val !== 'string' || val.length === 0) return null;
    const idx = val.indexOf(';');
    return idx >= 0 ? val.substring(0, idx) : val;
  });
  newCol.semType = semType;
}

/** Detects whether intensity columns appear to be log2-transformed.
 * Samples up to 200 non-null numeric values across the specified columns.
 * Heuristic: values >= 1000 suggest raw intensities; values in [0, 30] suggest log2. */
export function detectLog2Status(
  df: DG.DataFrame, colNames: string[],
): {isLog2: boolean; message: string} {
  let inLog2Range = 0;
  let inRawRange = 0;
  let total = 0;
  const maxSamples = 200;

  // Only count values that fall in either bucket — ambiguous "middle" values
  // (e.g. 31..999) carry no signal about scale, and counting them in `total`
  // makes the ratios at L101-105 unreliable on heterogeneous datasets.
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    if (col.type !== DG.COLUMN_TYPE.FLOAT && col.type !== DG.COLUMN_TYPE.INT && col.type !== DG.COLUMN_TYPE.BIG_INT) continue;
    for (let i = 0; i < df.rowCount && total < maxSamples; i++) {
      if (col.isNone(i)) continue;
      const val = Number(col.get(i));
      if (val >= 0 && val <= 30) { inLog2Range++; total++; }
      else if (val >= 1000) { inRawRange++; total++; }
    }
    if (total >= maxSamples) break;
  }

  if (total === 0)
    return {isLog2: false, message: 'No intensity values found'};
  if (inRawRange / total > 0.5)
    return {isLog2: false, message: 'Data appears to be raw intensities'};
  if (inLog2Range / total > 0.8)
    return {isLog2: true, message: 'Data appears to be already log2-transformed'};
  return {isLog2: false, message: 'Could not determine data scale'};
}

/** Detects delimiter by checking the first line for tab vs comma prevalence.
 * Returns '\t' or ','. */
export function detectDelimiter(text: string): string {
  const newlineIdx = text.indexOf('\n');
  const firstLine = newlineIdx >= 0 ? text.substring(0, newlineIdx) : text;
  const tabs = (firstLine.match(/\t/g) || []).length;
  const commas = (firstLine.match(/,/g) || []).length;
  return tabs > commas ? '\t' : ',';
}

/** Scans string columns for names containing protein-related keywords.
 * Returns the first matching column, or null. */
export function autoSuggestProteinIdColumn(df: DG.DataFrame): DG.Column | null {
  const keywords = ['protein', 'accession', 'uniprot'];
  for (const col of df.columns.toList()) {
    if (col.type !== DG.COLUMN_TYPE.STRING) continue;
    const lower = col.name.toLowerCase();
    for (const kw of keywords) {
      if (lower.includes(kw))
        return col;
    }
  }
  return null;
}

/** Scans string columns for names containing gene-related keywords.
 * Returns the first matching column, or null. */
export function autoSuggestGeneNameColumn(df: DG.DataFrame): DG.Column | null {
  const keywords = ['gene', 'symbol', 'gene_name', 'genename'];
  for (const col of df.columns.toList()) {
    if (col.type !== DG.COLUMN_TYPE.STRING) continue;
    const lower = col.name.toLowerCase();
    for (const kw of keywords) {
      if (lower.includes(kw))
        return col;
    }
  }
  return null;
}

/** Scans numeric column names for intensity-related keywords.
 * Returns all matching column names. */
export function autoSuggestIntensityColumns(df: DG.DataFrame): string[] {
  const keywords = ['intensity', 'lfq', 'ibaq', 'tmt', 'reporter', 'abundance'];
  const result: string[] = [];
  for (const col of df.columns.toList()) {
    if (col.type !== DG.COLUMN_TYPE.FLOAT && col.type !== DG.COLUMN_TYPE.INT && col.type !== DG.COLUMN_TYPE.BIG_INT) continue;
    const lower = col.name.toLowerCase();
    for (const kw of keywords) {
      if (lower.includes(kw)) {
        result.push(col.name);
        break;
      }
    }
  }
  return result;
}
