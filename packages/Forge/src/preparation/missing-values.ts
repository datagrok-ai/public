import * as DG from 'datagrok-api/dg';
import {defaultValuesOf} from '../engines/engine';
import {ForgeError} from '../forge-error';
import {columnRowCopy, compactCategories, rowCopy} from './row-copy';
import {releaseFrame, sharedFrame} from './shared-frame';

export const MISSING_VALUES_MODES = ['skip', 'impute'] as const;
export type MissingValuesMode = typeof MISSING_VALUES_MODES[number];
interface ImputeSettings { neighbors: number; distance: string }
export interface MissingValuesSettings { mode: MissingValuesMode; impute?: ImputeSettings }
export interface MissingColumn { name: string; count: number }

export interface PreparedData {
  features: DG.DataFrame;
  target?: DG.Column;
  /** Rows of the input that were kept, or null when every row was kept. */
  keptRows: DG.BitSet | null;
  skippedRows: number;
  imputedColumns: string[];
  /** Rows the imputation could not fill, skipped and counted in `skippedRows`. */
  failedRows: number;
}

const IMPUTE_FUNCTION = {package: 'Eda', name: 'knnImpute'};
const IMPUTE_DATA_INPUTS = ['table', 'columns', 'features'];
const IMPUTABLE_TYPES: string[] = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.STRING,
  DG.COLUMN_TYPE.DATE_TIME, DG.COLUMN_TYPE.QNUM];
const NO_IMPUTATION = 'Imputation is not available: the EDA package has no knnImpute function.';
const TWO_FEATURES = 'Impute needs at least two features. Choose Skip rows or check more features.';

export function missingColumnsOf(columns: DG.Column[]): MissingColumn[] {
  return columns
    .map((c) => ({name: c.name, count: c.stats.missingValueCount}))
    .filter((m) => m.count > 0);
}

export function imputeFunction(): DG.Func | undefined {
  return DG.Func.find(IMPUTE_FUNCTION)[0];
}

export function imputeSettingsOf(func: DG.Func): DG.Property[] {
  return func.inputs.filter((p) => !IMPUTE_DATA_INPUTS.includes(p.name));
}

/** Why the chosen handling cannot run on these features; empty for Skip rows and for features without gaps. */
export function missingValuesProblems(features: DG.Column[], settings: MissingValuesSettings): string[] {
  if (settings.mode !== 'impute')
    return [];
  return imputeProblems(features, missingColumnsOf(features).map((m) => m.name), imputeFunction());
}

/** {@link missingValuesProblems} of Impute for [features] whose gapped columns are [missing]. */
function imputeProblems(features: DG.Column[], missing: string[], func: DG.Func | undefined): string[] {
  if (missing.length === 0)
    return [];
  if (func === undefined)
    return [NO_IMPUTATION];
  const problems = features
    .filter((c) => c.type === DG.COLUMN_TYPE.BOOL && missing.includes(c.name))
    .map((c) => `Impute cannot fill the yes/no column '${c.name}'. Choose Skip rows.`);
  if (imputableColumns(features).length < 2)
    problems.push(TWO_FEATURES);
  return problems;
}

/** The result's `features` may share columns with [features]: release both with `releaseFrame` after use. */
export async function prepareMissingValues(features: DG.DataFrame, target: DG.Column | undefined,
  settings: MissingValuesSettings): Promise<PreparedData> {
  const rowCount = features.rowCount;
  if (settings.mode === 'skip') {
    const gapped = [...features.columns.toList(), ...(target === undefined ? [] : [target])].filter(hasGaps);
    return {...keep(features, target, keptRowsOf(gapped, rowCount)), imputedColumns: [], failedRows: 0};
  }

  const func = imputeFunction();
  const columns = features.columns.toList();
  const gapped = missingColumnsOf(columns).map((m) => m.name);
  const problems = imputeProblems(columns, gapped, func);
  if (problems.length > 0)
    throw new ForgeError(problems.join(' '));
  const withTarget = keep(features, target, keptRowsOf(target !== undefined && hasGaps(target) ? [target] : [],
    rowCount));
  // Skipping the rows without a target can take away a column's only gaps.
  const missing = withTarget.keptRows === null ? gapped :
    missingColumnsOf(withTarget.features.columns.toList()).map((m) => m.name);
  if (missing.length === 0 || func === undefined)
    return {...withTarget, imputedColumns: [], failedRows: 0};

  // The imputer writes in place: a frame of the user's columns gets clones of the columns it fills.
  const isShared = withTarget.keptRows === null;
  const imputeFrame = isShared ? sharedFrame(withTarget.features.columns.toList()
    .map((c) => missing.includes(c.name) ? c.clone() : c)) : withTarget.features;
  let isReturned = false;
  try {
    await func.apply({...defaultValuesOf(imputeSettingsOf(func)), ...settings.impute, table: imputeFrame,
      columns: missing, features: imputableColumns(imputeFrame.columns.toList()).map((c) => c.name)});
    // The filled columns are Forge's copies; a filled text column still lists the empty category.
    for (const name of missing)
      compactCategories(imputeFrame.getCol(name));
    const unfilled = keep(imputeFrame, withTarget.target,
      keptRowsOf(missing.map((name) => imputeFrame.getCol(name)), imputeFrame.rowCount));
    isReturned = unfilled.features === imputeFrame;
    const failedRows = unfilled.skippedRows;
    return {features: unfilled.features, target: unfilled.target,
      keptRows: combinedRows(withTarget.keptRows, unfilled.keptRows, rowCount),
      skippedRows: withTarget.skippedRows + failedRows, imputedColumns: missing, failedRows};
  } finally {
    if (isShared && !isReturned)
      releaseFrame(imputeFrame);
  }
}

function hasGaps(col: DG.Column): boolean {
  return col.stats.missingValueCount > 0;
}

function imputableColumns(columns: DG.Column[]): DG.Column[] {
  return columns.filter((c) => IMPUTABLE_TYPES.includes(c.type));
}

/** Rows without a missing value in any of [gapped], or null when no row has one. */
function keptRowsOf(gapped: DG.Column[], rowCount: number): DG.BitSet | null {
  if (gapped.length === 0)
    return null;
  const kept = DG.BitSet.create(rowCount, (r) => !gapped.some((c) => c.isNone(r)));
  return kept.trueCount === rowCount ? null : kept;
}

function keep(features: DG.DataFrame, target: DG.Column | undefined, keptRows: DG.BitSet | null):
  Pick<PreparedData, 'features' | 'target' | 'keptRows' | 'skippedRows'> {
  if (keptRows === null)
    return {features, target, keptRows, skippedRows: 0};
  return {features: rowCopy(features, keptRows), target: target === undefined ? undefined :
    columnRowCopy(target, keptRows), keptRows, skippedRows: features.rowCount - keptRows.trueCount};
}

/** [second] selects among the rows [first] kept; the result selects the same rows of the original [rowCount]. */
function combinedRows(first: DG.BitSet | null, second: DG.BitSet | null, rowCount: number): DG.BitSet | null {
  if (first === null || second === null)
    return first ?? second;
  const firstRows = first.getSelectedIndexes();
  const combined = DG.BitSet.create(rowCount);
  for (const i of second.getSelectedIndexes())
    combined.set(firstRows[i], true, false);
  return combined;
}
