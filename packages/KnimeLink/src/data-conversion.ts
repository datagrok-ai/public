import * as DG from 'datagrok-api/dg';

/** Convert a Datagrok DataFrame to KNIME Container Input (Table) JSON format. */
export function dataFrameToKnimeTable(df: DG.DataFrame): {'table-spec': any[]; 'table-data': any[][]} {
  const cols = df.columns.toList();
  const tableSpec = cols.map((col) => ({[col.name]: dgTypeToKnimeType(col.type)}));
  const dateTimeCols = new Set(cols.filter((c) => c.type === DG.COLUMN_TYPE.DATE_TIME).map((c) => c.name));
  const tableData: any[][] = [];
  for (let r = 0; r < df.rowCount; r++) {
    const row: any[] = [];
    for (const col of cols) {
      let val = col.get(r);
      // KNIME localdatetime expects "yyyy-MM-ddTHH:mm:ss.SSS" without timezone
      // col.get() returns dayjs for datetime columns
      if (dateTimeCols.has(col.name) && val != null)
        val = val.format ? val.format('YYYY-MM-DDTHH:mm:ss.SSS') : String(val).replace(/Z$/, '').replace(/[+-]\d{2}:\d{2}$/, '');
      row.push(val);
    }
    tableData.push(row);
  }
  return {'table-spec': tableSpec, 'table-data': tableData};
}

function dgTypeToKnimeType(dgType: string): string {
  switch (dgType) {
    case DG.COLUMN_TYPE.INT: return 'int';
    case DG.COLUMN_TYPE.FLOAT:
    case DG.COLUMN_TYPE.QNUM: return 'double';
    case DG.COLUMN_TYPE.BOOL: return 'boolean';
    case DG.COLUMN_TYPE.DATE_TIME: return 'localdatetime';
    case DG.COLUMN_TYPE.BIG_INT: return 'long';
    default: return 'string';
  }
}

/** Convert KNIME output (array of flat objects) to a Datagrok DataFrame. */
export function knimeTableToDataFrame(data: any[]): DG.DataFrame {
  if (!data || data.length === 0)
    return DG.DataFrame.create();
  return DG.DataFrame.fromObjects(data)!;
}

/** Convert KNIME table-spec/table-data format to a Datagrok DataFrame. */
export function knimeSpecDataToDataFrame(tableSpec: any[], tableData: any[][]): DG.DataFrame {
  if (!tableSpec || tableSpec.length === 0)
    return DG.DataFrame.create();
  // Each spec entry is {columnName: type}. KNIME may emit duplicate column names
  // (e.g. via Joiner / Column Appender); dedupe via columns.getUnusedName so we
  // don't silently drop columns when building objects below.
  const dedupeDf = DG.DataFrame.create();
  const colNames: string[] = [];
  for (const s of tableSpec) {
    const original = Object.keys(s)[0];
    const unique = dedupeDf.columns.getUnusedName(original);
    dedupeDf.columns.addNew(unique, DG.COLUMN_TYPE.STRING);
    colNames.push(unique);
  }
  const objects: any[] = [];
  for (const row of tableData) {
    const obj: any = {};
    for (let c = 0; c < colNames.length; c++)
      obj[colNames[c]] = row[c];
    objects.push(obj);
  }
  return DG.DataFrame.fromObjects(objects)!;
}

/** Map a KNIME column type string to the compatible Datagrok column types. */
export function knimeTypeToDgColumnTypes(knimeType: string): string[] {
  switch (knimeType) {
    case 'int': return [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.BIG_INT];
    case 'long': return [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.BIG_INT];
    case 'double': return [DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.QNUM];
    case 'boolean': return [DG.COLUMN_TYPE.BOOL];
    case 'localdatetime':
    case 'localdate':
    case 'localtime':
    case 'zoneddatetime': return [DG.COLUMN_TYPE.DATE_TIME];
    default: return [DG.COLUMN_TYPE.STRING];
  }
}
