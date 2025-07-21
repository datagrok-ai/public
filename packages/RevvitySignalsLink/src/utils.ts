import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { SignalsSearchParams, SignalsSearchQuery } from './signalsSearchQuery';

export const STORAGE_NAME = 'RevvitySignalsSearch';
export const QUERY_KEY = 'query';
export const PARAMS_KEY = 'params';

export function loadRevvitySearchQueryAndParams(): { loadedQuery: SignalsSearchQuery | null, loadedParams: SignalsSearchParams | null } {
  // Try to load saved query and params
  let loadedQuery: SignalsSearchQuery | null = null;
  let loadedParams: SignalsSearchParams | null = null;
  try {
    const savedQuery = grok.userSettings.getValue(STORAGE_NAME, QUERY_KEY);
    if (savedQuery)
      loadedQuery = JSON.parse(savedQuery);
  } catch (e) {
    loadedQuery = null;
  }
  try {
    const savedParams = grok.userSettings.getValue(STORAGE_NAME, PARAMS_KEY);
    if (savedParams)
      loadedParams = JSON.parse(savedParams);
  } catch (e) {
    loadedParams = null;
  }
  return { loadedQuery, loadedParams }
}


function inferColumnType(value: any): DG.COLUMN_TYPE {
  if (typeof value === 'boolean') return DG.COLUMN_TYPE.BOOL;
  if (typeof value === 'bigint') return DG.COLUMN_TYPE.BIG_INT;
  if (typeof value === 'number') {
    if (Number.isInteger(value)) return DG.COLUMN_TYPE.INT;
    return DG.COLUMN_TYPE.FLOAT;
  }
  // All strings become STRING type, no number inference
  if (typeof value === 'string') return DG.COLUMN_TYPE.STRING;
  return DG.COLUMN_TYPE.STRING;
}

export function dataFrameFromObjects(objects: Record<string, any>[]): DG.DataFrame {
  if (!objects || objects.length === 0)
    return DG.DataFrame.create();

  // Helper function to process array of objects - extract id or identifier filed and concatenate
  const processObjectArray = (arr: any[]): string => {
    if (!Array.isArray(arr)) return '';
    
    return arr.map(obj => {
      if (typeof obj !== 'object' || obj === null) return String(obj);
      
      // Try to find identifier field
      let identifier = obj.id;
      if (identifier === undefined) {
        const identifierKey = Object.keys(obj).find(key => 
          key.toLowerCase().includes('identifier')
        );
        identifier = identifierKey ? obj[identifierKey] : undefined;
      }
      
      return identifier !== undefined ? String(identifier) : '';
    }).filter(Boolean).join(',');
  };

  // Helper function to flatten objects recursively
  const flattenObject = (obj: Record<string, any>, prefix = ''): Record<string, any> => {
    const result: Record<string, any> = {};
    
    for (const key in obj) {
      if (obj.hasOwnProperty(key)) {
        const value = obj[key];
        const newKey = prefix ? `${prefix}.${key}` : key;
        
        if (Array.isArray(value) && value.length > 0 && typeof value[0] === 'object') {
          // Handle array of objects specially
          result[newKey] = processObjectArray(value);
        }
        else if (typeof value === 'object' && value !== null && !Array.isArray(value)) {
          // Recursively flatten nested objects
          Object.assign(result, flattenObject(value, newKey));
        } else {
          // Convert regular arrays to comma-separated strings
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
  const columns: DG.Column[] = [];
  for (const key of Array.from(allKeys)) {
    // Find first non-null, non-undefined value to infer type
    let firstValue = flattenedObjects.find(obj => obj[key] !== undefined && obj[key] !== null)?.[key];
    if (firstValue === undefined || firstValue === null) firstValue = '';

    const colType = inferColumnType(firstValue);
    columns.push(DG.Column.fromList(colType, key,
      flattenedObjects.map(obj => obj[key] !== undefined && obj[key] !== null ? obj[key] : null)));
  }

  return DG.DataFrame.fromColumns(columns);
}