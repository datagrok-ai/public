import * as DG from 'datagrok-api/dg';
import {MwxWorksheet, MwxColumn} from './mwx-types';


/** Converts a parsed MWX worksheet to a Datagrok DataFrame with full metadata tags. */
export function mwxWorksheetToDataFrame(ws: MwxWorksheet, format: string = 'mwx'): DG.DataFrame {
  const df = DG.DataFrame.create(ws.rowCount);
  df.name = ws.name;

  // Table-level metadata
  df.setTag('source.format', format);
  df.setTag('minitab.format', format);
  if (ws.version)
    df.setTag('minitab.version', String(ws.version));
  if (ws.worksheetId)
    df.setTag('minitab.worksheetId', ws.worksheetId);

  for (const col of ws.columns) {
    const dgCol = createColumn(col, ws.rowCount);
    setColumnTags(dgCol, col);
    df.columns.add(dgCol);
  }

  return df;
}


function setColumnTags(dgCol: DG.Column, col: MwxColumn): void {
  if (col.description)
    dgCol.setTag('description', col.description);

  dgCol.setTag('minitab.type', col.type);

  if (col.format) {
    dgCol.setTag('minitab.formatKey', String(col.format.formatKey));
    if (col.format.autoFormat != null)
      dgCol.setTag('minitab.autoFormat', String(col.format.autoFormat));
    if (col.format.numDecPlaces != null)
      dgCol.setTag('format', `#.${'0'.repeat(col.format.numDecPlaces)}`);
    if (col.format.minValue != null)
      dgCol.setTag('minitab.minValue', String(col.format.minValue));
    if (col.format.maxValue != null)
      dgCol.setTag('minitab.maxValue', String(col.format.maxValue));
    if (col.format.charCount != null)
      dgCol.setTag('minitab.charCount', String(col.format.charCount));
  }

  if (col.categories) {
    const ordered = Object.entries(col.categories)
      .sort((a, b) => a[1] - b[1])
      .map(([key]) => key);
    dgCol.setTag('minitab.categoryOrder', JSON.stringify(ordered));
  }
}


function createColumn(col: MwxColumn, rowCount: number): DG.Column {
  switch (col.type) {
  case 'numeric':
    return DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, col.name,
      col.values.map((v) => v == null ? DG.FLOAT_NULL : v as number));
  case 'datetime':
    return createDateTimeColumn(col, rowCount);
  default:
    return DG.Column.fromList(DG.COLUMN_TYPE.STRING, col.name,
      col.values.map((v) => v == null ? '' : String(v)));
  }
}


function createDateTimeColumn(col: MwxColumn, rowCount: number): DG.Column {
  const dgCol = DG.Column.dateTime(col.name, rowCount);
  for (let i = 0; i < rowCount; i++) {
    const v = col.values[i];
    if (v == null)
      dgCol.set(i, null);
    else
      dgCol.set(i, new Date(v as string).getTime() * 1000); // Datagrok uses microseconds
  }
  return dgCol;
}
