import * as DG from 'datagrok-api/dg';
import {IssueDetail} from '../types/validation-result';
import {COLUMN_FROM_DM_TAG} from '../constants/constants';

export function fixISO8601Date(value: any): string | null {
  if (value == null || value === '' || value === 'null' || value === 'Not in dataset')
    return null;

  const strValue = String(value).trim();
  if (strValue === '')
    return null;

  const iso8601Regex = /^\d{4}-\d{2}-\d{2}(T\d{2}:\d{2}:\d{2}(\.\d{1,3})?(Z|[+-]\d{2}:\d{2})?)?$/;
  if (iso8601Regex.test(strValue))
    return strValue;

  const date = new Date(strValue);
  if (!isNaN(date.getTime())) {
    const year = date.getFullYear();
    const month = String(date.getMonth() + 1).padStart(2, '0');
    const day = String(date.getDate()).padStart(2, '0');

    const hasTime = strValue.includes(':') || strValue.includes('T');
    if (hasTime) {
      const hours = String(date.getHours()).padStart(2, '0');
      const minutes = String(date.getMinutes()).padStart(2, '0');
      const seconds = String(date.getSeconds()).padStart(2, '0');
      return `${year}-${month}-${day}T${hours}:${minutes}:${seconds}`;
    } else
      return `${year}-${month}-${day}`;
  }

  return null;
}

export function fixSTRESNFromSTRESC(strescValue: any): number | null {
  if (strescValue == null || strescValue === '' || strescValue === 'null' || strescValue === 'Not in dataset')
    return null;

  const strValue = String(strescValue).trim();
  if (strValue === '')
    return null;

  const cleaned = strValue.replace(/[^\d.\-+Ee]/g, '');
  const parsed = parseFloat(cleaned);
  if (!isNaN(parsed) && isFinite(parsed))
    return parsed;

  return null;
}

function fixISO8601Format(
  df: DG.DataFrame,
  issueDetails: IssueDetail[],
): {df: DG.DataFrame, colsToFix: string[], colsOrder: string[]} {
  const affectedRows = new Set<number>();
  const affectedVariables = new Set<string>();
  const colsToFix: string[] = [];

  for (const issueDetail of issueDetails) {
    const row = issueDetail.row;
    if (row !== '' && row !== 'Not in dataset' && typeof row !== 'string') {
      const rowIdx = typeof row === 'number' ? row : parseInt(String(row));
      if (!isNaN(rowIdx) && rowIdx >= 1 && rowIdx < df.rowCount) {
        affectedRows.add(rowIdx - 1);
        if (Array.isArray(issueDetail.variables)) {
          issueDetail.variables.forEach((v: string) => {
            if (df.columns.names().includes(v) && (!df.col(v)!.getTag(COLUMN_FROM_DM_TAG) || df.name === 'dm'))
              affectedVariables.add(v);
          });
        }
      }
    }
  }

  const previewDf = df.clone(df.filter, Array.from(affectedVariables));
  const colsOrder = [];
  for (const varName of affectedVariables) {
    if (varName.startsWith('$') || varName.startsWith('define_') || varName.startsWith('library_'))
      continue;

    const originalCol = previewDf.col(varName);
    if (!originalCol)
      continue;
    colsOrder.push(varName);
    colsToFix.push(varName);

    const fixColName = `${varName}_fix`;
    const fixCol = previewDf.columns.addNewString(fixColName);
    colsOrder.push(fixColName);

    for (let i = 0; i < previewDf.rowCount; i++) {
      const originalValue = originalCol.get(i);
      const fixedValue = fixISO8601Date(originalValue);
      if (fixedValue !== null && fixedValue !== originalValue)
        fixCol.set(i, fixedValue);
    }
  }

  return {df: previewDf, colsToFix, colsOrder};
}

function fixSTRESNConversion(
  df: DG.DataFrame,
  issueDetails: IssueDetail[],
): {df: DG.DataFrame, colsToFix: string[], colsOrder: string[]} {
  const strescColName = `${df.name.toUpperCase()}STRESC`;
  const stresnColName = `${df.name.toUpperCase()}STRESN`;

  const strescCol = df.col(strescColName);
  const stresnCol = df.col(stresnColName);
  if (!strescCol || !stresnCol)
    return {df: DG.DataFrame.create(), colsToFix: [], colsOrder: []};

  const colsToFix = [stresnColName];
  const fixColName = `${stresnColName}_fix`;
  const colsOrder = [stresnColName, fixColName];

  const previewDf = df.clone(df.filter, [stresnColName]);
  const fixCol = previewDf.columns.addNewFloat(fixColName);

  let counter = 0;
  for (let i = 0; i < df.filter.length; i++) {
    if (df.filter.get(i)) {
      const rescValue = df.get(strescColName, i);
      const numericValue = fixSTRESNFromSTRESC(rescValue);
      if (numericValue !== null && numericValue !== rescValue)
        fixCol.set(counter, numericValue);
      counter++;
    }
  }
  previewDf.col(stresnColName)?.setTag('format', 'full precision');
  fixCol.setTag('format', 'full precision');

  return {df: previewDf, colsToFix, colsOrder};
}

export const validationFixFunctions: {
  [coreId: string]: (df: DG.DataFrame, issueDetails: IssueDetail[]) => {df: DG.DataFrame, colsToFix: string[],
    colsOrder: string[]}
} = {
  'CORE-000547': fixISO8601Format,
  'CORE-000353': fixISO8601Format,
  'CORE-000542': fixSTRESNConversion,
  'CORE-000863': fixSTRESNConversion,
};
