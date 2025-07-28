import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { SignalsSearchParams, SignalsSearchQuery } from './signalsSearchQuery';
import * as ui from 'datagrok-api/ui';
import { getRevvityUsers } from './users';
import { RevvityUser } from './revvityApi';
import { MOL_COL_NAME } from './compounds';

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


export function widgetFromObject(obj: Record<string, any>): DG.Widget {
  const map: Record<string, string> = {};
  const accordions: DG.Accordion[] = [];

  for (const key of Object.keys(obj)) {
    const value = obj[key];
    if (value === null || typeof value === 'string' || typeof value === 'number' || typeof value === 'boolean') {
      map[key] = value === null ? '' : String(value);
    } else if (Array.isArray(value)) {
      if (value.every(v => v === null || typeof v === 'string' || typeof v === 'number' || typeof v === 'boolean')) {
        map[key] = value.join(',');
      } else if (value.every(v => typeof v === 'object' && v !== null && !Array.isArray(v))) {
        // Array of objects: show as widgets in accordion
        const acc = ui.accordion();
        if (Array.isArray(value) && value?.length > 0) {
          acc.addPane(key, () => ui.divV(value.map(obj => widgetFromObject(obj).root)), true);
        } else {
          acc.addPane(key, () => ui.label('No data'), true);
        }
        accordions.push(acc);
      } else {
        // Mixed/other: show as JSON
        map[key] = JSON.stringify(value);
      }
    } else if (typeof value === 'object') {
      // Nested object: create accordion pane
      const acc = ui.accordion();
      acc.addPane(key, () => widgetFromObject(value).root, true);
      accordions.push(acc);
    }
  }

  const table = ui.tableFromMap(map);
  const children = [table, ...accordions.map(acc => acc.root)];
  return new DG.Widget(ui.divV(children));
}


const ENTITY_FIELDS_TO_EXCLUDE = ['type', 'isTemplate', 'tags', 'eid'];
const TAGS_TO_EXCLUDE = ['system.Keywords'];
const OWNER_FIELDS = ['owner', 'createdBy', 'editedBy'];
const FIRST_COL_NAMES = [MOL_COL_NAME, 'id', 'name'];

export async function transformData(data: Record<string, any>[]): Promise<Record<string, any>[]> {
  const items = [];
  for (const item of data) {
        const result: any = { id: item.id };
    const attrs = item.attributes || {};

    for (const [key, value] of Object.entries(attrs)) {
      if (ENTITY_FIELDS_TO_EXCLUDE.includes(key))
        continue;
      //handle fileds related to users
      if (OWNER_FIELDS.includes(key)) {
        const user = (await getRevvityUsers())![value as string] as RevvityUser;
        result[key] = user.firstName && user.lastName ? `${user.firstName} ${user.lastName}` : user.userName ?? value;
        continue;
      }
      if (Array.isArray(value)) {
        result[key] = value.join('; ');
      } else {
        result[key] = value;
      }
    }

    // Handle tags
    if (attrs.tags) {
      for (const [tagKey, tagValue] of Object.entries(attrs.tags)) {
        if (TAGS_TO_EXCLUDE.includes(tagKey))
          continue;
        const fieldName = tagKey.includes('.') ? tagKey.split('.').slice(1).join('.') : tagKey;
        if (Array.isArray(tagValue)) {
          result[fieldName] = tagValue.join('; ');
        } else {
          result[fieldName] = tagValue;
        }
      }
    }
    items.push(result);
  }
  return items;
}

export async function reorderColummns(df: DG.DataFrame) {
  const colNames = df.columns.names();
  const newColOrder = [];
  for (const colName of FIRST_COL_NAMES) {
    const index = colNames.indexOf(colName);
    if (index > -1) {
      colNames.splice(index, 1);
      newColOrder.push(colName);
    }
  }
  df.columns.setOrder(newColOrder.concat(colNames));
}