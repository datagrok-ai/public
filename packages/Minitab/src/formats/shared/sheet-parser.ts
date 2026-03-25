import {MwxColumn, MwxColumnType} from '../mwx/mwx-types';


/** Parses the columns from a Minitab sheet.json `Data` object. Shared by MWX and MPX. */
export function parseSheetColumns(data: any): {columns: MwxColumn[], rowCount: number} {
  const rawColumns: any[] = data.Columns ?? [];
  const columns: MwxColumn[] = [];
  let maxRows = 0;

  for (const rawCol of rawColumns) {
    const col = parseColumn(rawCol);
    if (col) {
      maxRows = Math.max(maxRows, col.values.length);
      columns.push(col);
    }
  }

  for (const col of columns) {
    while (col.values.length < maxRows)
      col.values.push(null);
  }

  return {columns, rowCount: maxRows};
}


function parseColumn(raw: any): MwxColumn | null {
  const body = raw.WorksheetVarBody;
  if (!body)
    return null;

  const name = body.Name ?? '';
  const description = body.Desc ?? '';
  const varData = body.VarData;
  if (!varData)
    return null;

  const varDataBody = varData.VarDataBody;
  if (!varDataBody)
    return null;

  const type = detectColumnType(body, varData);
  const values = extractValues(type, varDataBody, body);

  let categories: {[key: string]: number} | undefined;
  const ordering = varData.Ordering;
  if (ordering?.TextOrder && ordering.TextOrder.length > 0) {
    categories = {};
    for (const entry of ordering.TextOrder)
      categories[entry.Key] = entry.Value;
  }

  return {name, description, type, values, categories};
}


function detectColumnType(body: any, varData: any): MwxColumnType {
  const formatKey = body.Format?.Key;
  if (formatKey === 4 || formatKey === 5)
    return 'datetime';

  const dataType = varData.Data;
  if (dataType === 1)
    return 'numeric';
  if (dataType === 2)
    return 'text';

  const vdb = varData.VarDataBody;
  if (vdb?.HasNumericData)
    return 'numeric';
  if (vdb?.HasTextData)
    return 'text';

  return 'text';
}


function extractValues(type: MwxColumnType, varDataBody: any, body: any): (number | string | null)[] {
  if (type === 'datetime' && varDataBody.HasNumericData && varDataBody.NumericData)
    return varDataBody.NumericData.map((v: number | null) => convertMinitabDate(v));

  if (type === 'numeric' && varDataBody.HasNumericData && varDataBody.NumericData)
    return varDataBody.NumericData as (number | null)[];

  if (varDataBody.HasTextData && varDataBody.TextData) {
    const missingLabel = body.Format?.Value?.MissingValueLabel || '*';
    return varDataBody.TextData.map((v: string) => {
      if (v === missingLabel || v === '*')
        return null;
      return v;
    });
  }

  return [];
}


/** Converts Minitab date value (days since 1899-12-30, like Excel) to ISO string. */
function convertMinitabDate(value: number | null): string | null {
  if (value == null || isNaN(value))
    return null;
  const msPerDay = 86400000;
  const epoch = new Date(1899, 11, 30).getTime();
  return new Date(epoch + value * msPerDay).toISOString();
}
