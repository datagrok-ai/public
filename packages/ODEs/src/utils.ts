// Utitlities

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MISC, PATH, ERROR_MSG, INPUTS_DF} from './ui-constants';
import {CONTROL_EXPR} from './constants';

/** Return max absolute deviation between the corresponding float values of 2 dataframes */
export function error(df1: DG.DataFrame, df2: DG.DataFrame): number {
  let mad = 0;

  const cols1 = df1.columns.toList();
  const cols2 = df2.columns.toList();

  const rowCount = df1.rowCount;

  if (rowCount !== df2.rowCount)
    throw new Error('Dataframes have non-equal rows counts');

  const colCount = cols1.length;

  if (colCount !== cols2.length)
    throw new Error('Dataframes have non-equal columns counts');

  for (let i = 1; i < colCount; ++i) { // we skip the argument column
    const c1 = cols1[i];
    const c2 = cols2[i];

    if ((c1.type !== DG.COLUMN_TYPE.FLOAT) || (c2.type !== DG.COLUMN_TYPE.FLOAT))
      continue;

    const a1 = c1.getRawData();
    const a2 = c2.getRawData();

    for (let j = 0; j < rowCount; ++j)
      mad = Math.max(mad, Math.abs(a1[j] - a2[j]));
  }

  return mad;
} // error

/** Return unused IVP-file name */
export function unusedFileName(name: string, files: string[]): string {
  if (!files.includes(`${name}.${MISC.IVP_EXT}`))
    return name;

  let num = 1;

  while (files.includes(`${name}(${num}).${MISC.IVP_EXT}`))
    ++num;

  return `${name}(${num})`;
}

/** Return dataframe with the specified number of last rows */
export function getTableFromLastRows(df: DG.DataFrame, maxRows: number): DG.DataFrame {
  if (df.rowCount <= maxRows)
    return df;

  const cols: DG.Column[] = [];

  for (const col of df.columns)
    cols.push(DG.Column.fromList(col.type, col.name, col.toList().slice(-maxRows)));

  return DG.DataFrame.fromColumns(cols);
}

/** Return table with inputs */
export async function getInputsTable(inputsPath: string): Promise<DG.DataFrame | null> {
  const idx = inputsPath.indexOf('/');
  let path: string;

  if (idx < 0) {
    grok.shell.warning(`${ERROR_MSG.INPUTS_PATH}'${CONTROL_EXPR.INPUTS}'-line`);
    return null;
  }

  // Create platform path
  if (inputsPath.includes(PATH.MY_FILES))
    path = `${grok.shell.user.login}:${PATH.HOME}${inputsPath.slice(idx)}`;
  else
    path = `${PATH.SYSTEM}:${inputsPath}`;

  // Load table with inputs
  try {
    const dfs = await grok.dapi.files.readBinaryDataFrames(path);
    const inputsDf = dfs[0];
    const cols = inputsDf.columns;

    const msg = `${ERROR_MSG.INCORRECT_INPUT_DF} '${path}'.`;

    // check columns count
    if (cols.length < INPUTS_DF.MIN_COLS_COUNT) {
      grok.shell.warning(`${msg} At least two columns are needed.`);
      return null;
    }

    // check nulls & numerical cols
    let numColsCount = 0;

    for (const col of cols) {
      if (col.stats.missingValueCount > 0) {
        grok.shell.warning(`${msg} ${ERROR_MSG.INPUT_DF_MIS_VALS}`);
        return null;
      }

      if (col.isNumerical)
        ++numColsCount;
    }

    if (numColsCount === 0) {
      grok.shell.warning(`${msg} ${ERROR_MSG.NO_NUM_COLS}`);
      return null;
    }

    // check column with names of inputs
    if (cols.byIndex(INPUTS_DF.INP_NAMES_IDX).type !== DG.COLUMN_TYPE.STRING) {
      grok.shell.warning(`${msg} ${ERROR_MSG.INPUT_DF_1_COL}`);
      return null;
    }

    return inputsDf;
  } catch (err) {
    grok.shell.warning(`${ERROR_MSG.INPUT_TABLE}'${CONTROL_EXPR.INPUTS}'-line`);
  }

  return null;
} // getInputsTable
