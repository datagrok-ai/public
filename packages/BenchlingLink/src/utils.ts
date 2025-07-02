import * as DG from 'datagrok-api/dg';

function inferColumnType(value: any): DG.COLUMN_TYPE {
  if (typeof value === 'boolean') return DG.COLUMN_TYPE.BOOL;
  if (typeof value === 'bigint') return DG.COLUMN_TYPE.BIG_INT;
  if (typeof value === 'number') {
    if (Number.isInteger(value)) return DG.COLUMN_TYPE.INT;
    return DG.COLUMN_TYPE.FLOAT;
  }
  if (typeof value === 'string') {
    // Try to infer QNUM (quantitative number) if possible
    if (!isNaN(Number(value))) {
      if (value.includes('.')) return DG.COLUMN_TYPE.FLOAT;
      return DG.COLUMN_TYPE.INT;
    }
    return DG.COLUMN_TYPE.STRING;
  }
  return DG.COLUMN_TYPE.STRING;
}

export function dataFrameFromObjects(objects: Record<string, any>[]): DG.DataFrame {
  if (!objects || objects.length === 0)
    return DG.DataFrame.create();
  const keys = Object.keys(objects[0]);
  const columns: DG.Column[] = [];
  for (const key of keys) {
    let firstValue = objects.find(obj => obj[key] !== undefined && obj[key] !== null)?.[key];
    if (firstValue === undefined || firstValue === null) firstValue = '';
    if (Array.isArray(firstValue)) {
      columns.push(DG.Column.fromList(DG.COLUMN_TYPE.STRING, key, objects.map(obj => Array.isArray(obj[key]) ? obj[key].join(',') : '')));
    } else if (typeof firstValue === 'object' && firstValue !== null) {
      continue;
    } else {
      const colType = inferColumnType(firstValue);
      columns.push(DG.Column.fromList(colType, key, objects.map(obj => obj[key] !== undefined && obj[key] !== null ? obj[key] : null)));
    }
  }
  return DG.DataFrame.fromColumns(columns);
} 