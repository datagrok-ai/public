import * as DG from 'datagrok-api/dg';
import {MwxWorksheet, MwxColumn} from './mwx-types';


/** Converts a parsed MWX worksheet to a Datagrok DataFrame. */
export function mwxWorksheetToDataFrame(ws: MwxWorksheet): DG.DataFrame {
  const df = DG.DataFrame.create(ws.rowCount);
  df.name = ws.name;

  for (const col of ws.columns) {
    const dgCol = createColumn(col, ws.rowCount);
    if (col.description)
      dgCol.setTag('description', col.description);
    df.columns.add(dgCol);
  }

  return df;
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
