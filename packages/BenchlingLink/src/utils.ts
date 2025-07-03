import * as DG from 'datagrok-api/dg';

function inferColumnType(value: any): DG.COLUMN_TYPE {
  if (typeof value === 'boolean') return DG.COLUMN_TYPE.BOOL;
  if (typeof value === 'bigint') return DG.COLUMN_TYPE.BIG_INT;
  if (typeof value === 'number') {
    if (Number.isInteger(value)) return DG.COLUMN_TYPE.INT;
    return DG.COLUMN_TYPE.FLOAT;
  }
  if (typeof value === 'string') {
    // Handle empty string first
    if (value.trim() === '') return DG.COLUMN_TYPE.STRING;
    
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

  // Helper function to flatten objects and convert arrays to strings
  const flattenObject = (obj: Record<string, any>, prefix = ''): Record<string, any> => {
    const result: Record<string, any> = {};
    
    for (const key in obj) {
      if (obj.hasOwnProperty(key)) {
        const value = obj[key];
        const newKey = prefix ? `${prefix}.${key}` : key;
        
        if (typeof value === 'object' && value !== null && !Array.isArray(value)) {
          // Recursively flatten nested objects
          Object.assign(result, flattenObject(value, newKey));
        } else {
          // Convert arrays to comma-separated strings
          result[newKey] = Array.isArray(value) ? value.join(',') : value;
        }
      }
    }
    
    return result;
  };

  // First flatten all objects in the array
  const flattenedObjects = objects.map(obj => flattenObject(obj));

  // Now get all unique keys from all flattened objects
  const allKeys = new Set<string>();
  for (const obj of flattenedObjects) {
    Object.keys(obj).forEach(key => allKeys.add(key));
  }

  // Create columns for each key
  for (const key of Array.from(allKeys)) {
    // Find first non-null, non-undefined value to infer type
    let firstValue = flattenedObjects.find(obj => obj[key] !== undefined && obj[key] !== null)?.[key];
    if (firstValue === undefined || firstValue === null) firstValue = '';

    const colType = inferColumnType(firstValue);
    const arr = flattenedObjects.map(obj => obj[key] !== undefined && obj[key] !== null ? obj[key] : null);
    columns.push(DG.Column.fromList(colType, key, arr));
  }

  return DG.DataFrame.fromColumns(columns);
}