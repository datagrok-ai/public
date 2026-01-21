/* eslint-disable valid-jsdoc */
import * as DG from 'datagrok-api/dg';
import {IssueDetail} from '../types/validation-result';
import { COLUMN_FROM_DM_TAG } from '../constants/constants';

/**
 * Fixes ISO 8601 date/datetime format issues
 * Attempts to parse various date formats and convert them to ISO 8601 format
 * @param value - The date value to fix
 * @return Fixed ISO 8601 formatted date string, or null if parsing fails
 */
export function fixISO8601Date(value: any): string | null {
  if (value == null || value === '' || value === 'null' || value === 'Not in dataset')
    return null;

  const strValue = String(value).trim();
  if (strValue === '')
    return null;

  // Check if already in ISO 8601 format
  const iso8601Regex = /^\d{4}-\d{2}-\d{2}(T\d{2}:\d{2}:\d{2}(\.\d{1,3})?(Z|[+-]\d{2}:\d{2})?)?$/;
  if (iso8601Regex.test(strValue))
    return strValue;

  // Try to parse common date formats
  const date = new Date(strValue);
  if (!isNaN(date.getTime())) {
    // Format as ISO 8601 date (YYYY-MM-DD)
    const year = date.getFullYear();
    const month = String(date.getMonth() + 1).padStart(2, '0');
    const day = String(date.getDate()).padStart(2, '0');

    // Check if original value had time component
    const hasTime = strValue.includes(':') || strValue.includes('T');
    if (hasTime) {
      const hours = String(date.getHours()).padStart(2, '0');
      const minutes = String(date.getMinutes()).padStart(2, '0');
      const seconds = String(date.getSeconds()).padStart(2, '0');
      return `${year}-${month}-${day}T${hours}:${minutes}:${seconds}`;
    } else
      return `${year}-${month}-${day}`;
  }

  // Try common date patterns
  const patterns = [
    // MM/DD/YYYY or MM-DD-YYYY
    /^(\d{1,2})[\/\-](\d{1,2})[\/\-](\d{4})(\s+(\d{1,2}):(\d{2})(:(\d{2}))?)?$/,
    // DD/MM/YYYY or DD-MM-YYYY
    /^(\d{1,2})[\/\-](\d{1,2})[\/\-](\d{4})(\s+(\d{1,2}):(\d{2})(:(\d{2}))?)?$/,
    // YYYY/MM/DD or YYYY-MM-DD
    /^(\d{4})[\/\-](\d{1,2})[\/\-](\d{1,2})(\s+(\d{1,2}):(\d{2})(:(\d{2}))?)?$/,
  ];

  for (const pattern of patterns) {
    const match = strValue.match(pattern);
    if (match) {
      let year: string; let month: string; let day: string;
      if (match[3] && match[3].length === 4) {
        // MM/DD/YYYY format
        year = match[3];
        month = match[1].padStart(2, '0');
        day = match[2].padStart(2, '0');
      } else if (match[1] && match[1].length === 4) {
        // YYYY/MM/DD format
        year = match[1];
        month = match[2].padStart(2, '0');
        day = match[3].padStart(2, '0');
      } else {
        // Try DD/MM/YYYY (assume first part is day if > 12)
        const first = parseInt(match[1]);
        if (first > 12) {
          // DD/MM/YYYY
          year = match[3];
          month = match[2].padStart(2, '0');
          day = match[1].padStart(2, '0');
        } else {
          // MM/DD/YYYY
          year = match[3];
          month = match[1].padStart(2, '0');
          day = match[2].padStart(2, '0');
        }
      }

      if (match[5] && match[6]) {
        // Has time component
        const hours = match[5].padStart(2, '0');
        const minutes = match[6].padStart(2, '0');
        const seconds = match[8] ? match[8].padStart(2, '0') : '00';
        return `${year}-${month}-${day}T${hours}:${minutes}:${seconds}`;
      } else
        return `${year}-${month}-${day}`;
    }
  }

  return null;
}

/**
 * Fixes STRESC to STRESN conversion issues
 * Parses numeric values from STRESC and returns the numeric value for STRESN
 * @param strescValue - The STRESC (character result) value
 * @return Parsed numeric value, or null if not numeric
 */
export function fixSTRESNFromSTRESC(strescValue: any): number | null {
  if (strescValue == null || strescValue === '' || strescValue === 'null' || strescValue === 'Not in dataset')
    return null;

  const strValue = String(strescValue).trim();
  if (strValue === '')
    return null;

  // Remove common non-numeric characters but keep decimal point and minus sign
  const cleaned = strValue.replace(/[^\d.\-+Ee]/g, '');

  // Try to parse as float
  const parsed = parseFloat(cleaned);
  if (!isNaN(parsed) && isFinite(parsed))
    return parsed;

  return null;
}


/**
 * Fixes ISO 8601 date/datetime format issues for multiple columns and rows
 * @param df - The dataframe containing the data
 * @param domain - The domain name (e.g., 'LB', 'VS', 'OM')
 * @param issueDetails - Array of issue details to fix
 * @return Map of original column names to their fix columns
 */
function fixISO8601Format(
  df: DG.DataFrame,
  issueDetails: IssueDetail[],
): {df: DG.DataFrame, colsToFix: string[], colsOrder: string[]} {
  const affectedRows = new Set<number>();
  const affectedVariables = new Set<string>();
  const colsToFix: string[] = [];

  // Collect affected rows and variables from issue details
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
  // Create fix columns for each affected variable
  const colsOrder = [];
  for (const varName of affectedVariables) {
    // Skip metadata variables
    if (varName.startsWith('$') || varName.startsWith('define_') || varName.startsWith('library_'))
      continue;

    const originalCol = previewDf.col(varName);
    if (!originalCol)
      continue;
    colsOrder.push(varName);
    colsToFix.push(varName);

    // Add fix column to preview
    const fixColName = `${varName}_fix`;
    const fixCol = previewDf.columns.addNewString(fixColName);
    colsOrder.push(fixColName);

    // Apply fix function to affected rows
    for (let i = 0; i < previewDf.rowCount; i++) {
      const originalValue = originalCol.get(i);
      const fixedValue = fixISO8601Date(originalValue);
      if (fixedValue !== null && fixedValue !== originalValue)
        fixCol.set(i, fixedValue);
    }
  }

  return {df: previewDf, colsToFix, colsOrder};
}

/**
 * Fixes STRESC/STRESN numeric conversion issues
 * @param df - The dataframe containing the data
 * @param domain - The domain name (e.g., 'LB', 'VS', 'OM')
 * @param issueDetails - Array of issue details to fix
 * @return Map of original column names to their fix columns
 */
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
  // Check if STRESN column exists, create fix column if needed
  const fixColName = `${stresnColName}_fix`;
  const colsOrder = [stresnColName, fixColName];

  const previewDf = df.clone(df.filter, [stresnColName]);
  const fixCol = previewDf.columns.addNewFloat(fixColName);

  // Apply fix function to affected rows
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
  previewDf.col(stresnColName).setTag('format', 'full precision');
  fixCol.setTag('format', 'full precision');

  return {df: previewDf, colsToFix, colsOrder};
}

/**
 * Dictionary mapping core rule IDs to their corresponding fix functions
 * Each fix function takes a dataframe, domain name, and array of issue details,
 * and returns a map of original column names to fix columns
 * Based on Issue_Summary analysis from validation results
 */
export const validationFixFunctions: {
  [coreId: string]: (df: DG.DataFrame, issueDetails: IssueDetail[]) => {df: DG.DataFrame, colsToFix: string[],
    colsOrder: string[]}
} = {
  // ISO 8601 date/datetime format errors
  'CORE-000547': fixISO8601Format, // Variable value is not in correct ISO 8601 date or datetime format
  // 'CORE-000548': fixISO8601Format, // PCDUR, PCEVLINT or PCELTM is not in ISO 8601 Duration format
  'CORE-000353': fixISO8601Format, // Invalid date/datetime or duration format

  // STRESC/STRESN numeric conversion errors
  // STRESC is numeric but STRESN is not populated or not equal to STRESC
  // (works for LB, OM, PP, etc.)
  'CORE-000542': fixSTRESNConversion,
  // STRESC is populated with a numeric value, but STRESN is not populated
  'CORE-000863': fixSTRESNConversion,
};


