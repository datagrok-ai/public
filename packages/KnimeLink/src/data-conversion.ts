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
export function knimeTableToDataFrame(data: any[], name?: string): DG.DataFrame {
  if (!data || data.length === 0)
    return DG.DataFrame.create();
  return DG.DataFrame.fromObjects(data)!;
}

/** Convert KNIME table-spec/table-data format to a Datagrok DataFrame. */
export function knimeSpecDataToDataFrame(tableSpec: any[], tableData: any[][]): DG.DataFrame {
  if (!tableSpec || tableSpec.length === 0)
    return DG.DataFrame.create();
  // Each spec entry is {columnName: type} — extract column names
  const colNames = tableSpec.map((s) => Object.keys(s)[0]);
  const objects: any[] = [];
  for (const row of tableData) {
    const obj: any = {};
    for (let c = 0; c < colNames.length; c++)
      obj[colNames[c]] = row[c];
    objects.push(obj);
  }
  return DG.DataFrame.fromObjects(objects)!;
}

/** Convert a File to a base64 string for KNIME Container Input (File). */
export function fileToBase64(file: File): Promise<string> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => {
      const result = reader.result as string;
      resolve(result.split(',')[1]);
    };
    reader.onerror = () => reject(reader.error);
    reader.readAsDataURL(file);
  });
}

/** Decode a base64 string to a Blob. */
export function base64ToBlob(base64: string, mimeType: string = 'application/octet-stream'): Blob {
  const bytes = atob(base64);
  const buffer = new Uint8Array(bytes.length);
  for (let i = 0; i < bytes.length; i++)
    buffer[i] = bytes.charCodeAt(i);
  return new Blob([buffer], {type: mimeType});
}

/** Try to parse KNIME execution output into DataFrames. */
export function parseKnimeOutputs(outputs: {[key: string]: any}): {tables: DG.DataFrame[]; variables: {[key: string]: any}} {
  const tables: DG.DataFrame[] = [];
  const variables: {[key: string]: any} = {};

  for (const key of Object.keys(outputs)) {
    const val = outputs[key];
    if (Array.isArray(val) && val.length > 0 && typeof val[0] === 'object') {
      const df = knimeTableToDataFrame(val, key);
      df.name = key;
      tables.push(df);
    }
    else if (typeof val === 'object' && val !== null && val['table-spec'] && val['table-data']) {
      const df = knimeSpecDataToDataFrame(val['table-spec'], val['table-data']);
      df.name = key;
      tables.push(df);
    }
    else
      variables[key] = val;
  }

  return {tables, variables};
}
