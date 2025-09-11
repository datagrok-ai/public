import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { SignalsSearchParams, SignalsSearchQuery } from './signals-search-query';
import * as ui from 'datagrok-api/ui';
import { getRevvityUsersWithMapping } from './users';
import { queryEntityById, RevvityApiResponse, RevvityData, RevvityUser } from './revvity-api';

import { awaitCheck } from '@datagrok-libraries/utils/src/test';
import { compoundTypeAndViewNameMapping, ENTITY_FIELDS_TO_EXCLUDE, FIELDS_SECTION_NAME, FIELDS_TO_EXCLUDE_FROM_WIDGET, FIRST_COL_NAMES, MOL_COL_NAME, MOLECULAR_FORMULA_FIELD_NAME, PARAMS_KEY, QUERY_KEY, STORAGE_NAME, SUBMITTER_FIELD_NAME, TABS_TO_EXCLUDE_FROM_WIDGET, TAGS_TO_EXCLUDE, USER_FIELDS } from './constants';


function extractNameFromJavaObjectString(javaString: string): string {
  try {
    // Extract firstName
    const firstNameMatch = javaString.match(/firstName=([^,}]+)/);
    const firstName = firstNameMatch ? firstNameMatch[1].trim() : '';
    
    // Extract lastName
    const lastNameMatch = javaString.match(/lastName=([^,}]+)/);
    const lastName = lastNameMatch ? lastNameMatch[1].trim() : '';
    
    if (firstName && lastName) {
      return `${firstName} ${lastName}`;
    } else if (firstName) {
      return firstName;
    } else if (lastName) {
      return lastName;
    } else {
      // Fallback: try to extract userName
      const userNameMatch = javaString.match(/userName=([^,}]+)/);
      return userNameMatch ? userNameMatch[1].trim() : 'Unknown User';
    }
  } catch (e) {
    return 'Unknown User';
  }
}

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


export function processField(map: Record<string, any>, key: string, value: any): void {
  if (typeof value === 'object' || FIELDS_TO_EXCLUDE_FROM_WIDGET.includes(key)) {
    // Skip objects except for fields
    return;
  } else if (Array.isArray(value))
    map[key] = value.join(', ');
  else {
    //handle molecular formula fields since they contain <sub> tags
    if (key.includes(MOLECULAR_FORMULA_FIELD_NAME))
      value = (value as string).replaceAll('<sub>', '').replaceAll('</sub>', '');
    else if (key === SUBMITTER_FIELD_NAME) {
      if (typeof value === 'string' && value.includes('firstName=')) {
        value = extractNameFromJavaObjectString(value);
      } else {
        value = 'Unknown User';
      }
    }
    map[key] = value === null ? '' : value;
  }
}


export function createRevvityResponseWidget(response: RevvityApiResponse, idSemValue: DG.SemanticValue<string>): HTMLDivElement {
  const data = response.data as RevvityData;
  if (!data || !data.attributes) {
    return ui.divText('No data available');
  }

  const attributes = data.attributes;
  const relationships = data.relationships || {};
  const included = response.included || [];

  // Create the main table from attributes
  const mainMap: Record<string, any> = {};

  // Process regular attributes
  for (const [key, value] of Object.entries(attributes)) {
    if (key === FIELDS_SECTION_NAME) continue; // Handle fields separately
    processField(mainMap, key, value as string);
  }

  // Process fields
  if (attributes.fields) {
    for (const [fieldName, fieldObj] of Object.entries(attributes.fields)) {
      if ((fieldObj as any)['value']) {
        processField(mainMap, fieldName, (fieldObj as any)['value'] as string);
      }
    }
  }
  const widgetMapWithAddIcons: { [key: string]: HTMLElement } = {};
  Object.keys(mainMap).forEach((key) => {
      const propDiv = ui.divH([mainMap[key]], {style: {'position': 'relative'}});
        if (idSemValue.cell.dataFrame && idSemValue.cell.column) {
      const addColumnIcon = ui.iconFA('plus', async () => {
        calculatePropForColumn(key, idSemValue.cell.dataFrame, idSemValue.cell.column);
      }, `Calculate ${key} for the whole table`);
      addColumnIcon.classList.add('revvity-add-prop-as-column-icon');
      propDiv.prepend(addColumnIcon);
      ui.tools.setHoverVisibility(propDiv, [addColumnIcon]);
      widgetMapWithAddIcons[key] = propDiv;
    }
  });

  const mainTable = ui.tableFromMap(widgetMapWithAddIcons);
  const accordions: DG.Accordion[] = [];

  // Create accordion tabs for relationships
  for (const [relationshipName, relationshipData] of Object.entries(relationships)) {
    if (TABS_TO_EXCLUDE_FROM_WIDGET.includes(relationshipName))
      continue;
    const relationship = relationshipData as any;
    if (!relationship.data) continue;

    const relationshipItems: any[] = [];

    // Handle single relationship or array of relationships
    const items = Array.isArray(relationship.data) ? relationship.data : [relationship.data];

    for (const item of items) {
      // Find corresponding entity in included section
      const includedEntity = included.find((inc: any) => inc.id === item.id);
      if (includedEntity && includedEntity.attributes) {
        const entityMap: Record<string, string> = {};

        // Process attributes
        for (const [key, value] of Object.entries(includedEntity.attributes)) {
          if (key === FIELDS_SECTION_NAME) continue; // Handle fields separately
          processField(entityMap, key, value as string);
        }

        // Process fields if present
        if (includedEntity.attributes.fields) {
          for (const [fieldName, fieldObj] of Object.entries(includedEntity.attributes.fields)) {
            processField(entityMap, fieldName, (fieldObj as any)['value'] as string);
          }
        }
        relationshipItems.push(entityMap);
      }
    }

    if (relationshipItems.length > 0) {
      const acc = ui.accordion();
      acc.addPane(relationshipName, () => {
        if (relationshipItems.length === 1)
          return ui.tableFromMap(relationshipItems[0]);
        else {
          const innerAcc = ui.accordion();
          for (const relationshipItem of relationshipItems) {
            innerAcc.addPane(USER_FIELDS.includes(relationshipName) ? relationshipItem.userId ?? 'User' : relationshipItem.name ?? relationshipItem.id, () => {
              return ui.tableFromMap(relationshipItem);
            });
          }
          return innerAcc.root;
        }
      });
      accordions.push(acc);
    }
  }

  const children = [mainTable, ...accordions.map(acc => acc.root)];
  return ui.divV(children);
}

export async function transformData(data: Record<string, any>[]): Promise<Record<string, any>[]> {
  const items = [];
  for (const item of data) {
        const result: any = { id: item.id };
    const attrs = item.attributes || {};

    for (const [key, value] of Object.entries(attrs)) {
      if (ENTITY_FIELDS_TO_EXCLUDE.includes(key))
        continue;
      //handle fileds related to users
      if (USER_FIELDS.includes(key)) {
        const user = (await getRevvityUsersWithMapping())![value as string].revvityUser as RevvityUser;
        result[key] = user.userId;
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

export function getCompoundTypeByViewName(viewName: string): string {
  const res = compoundTypeAndViewNameMapping.filter((it) => it.viewName === viewName);
  return !res.length ? viewName : res[0].compoundType;
}

export function getViewNameByCompoundType(compoundType: string): string {
  const res = compoundTypeAndViewNameMapping.filter((it) => it.compoundType === compoundType);
  return !res.length ? compoundType : res[0].viewName;
}

export async function calculatePropForColumn(propName: string, df: DG.DataFrame, idsCol: DG.Column) {
  const newColName = df.columns.getUnusedName(propName);
  const newCol = df.columns.addNewString(newColName);
  const promises = idsCol.toList().map((id, index) =>
    getEntityPropByEntityId(propName, id)
      .then(res => newCol.set(index, res))
      .catch(() => { })
  );
  await Promise.allSettled(promises);

}

export async function getEntityPropByEntityId(propName: string, id: string): Promise<any> {
  let prop: any = null;
  try {
  const entity = await queryEntityById(id);
  if (entity.data) {
    const attributes = (entity.data as RevvityData).attributes;
    if (attributes) {
      //look for prop in attributes
      if(Object.keys(attributes).includes(propName))
        prop = attributes[propName];
      //if not in attributes - look in fields
      if(attributes['fields'] && Object.keys(attributes['fields']))
        prop = attributes['fields'][propName].value;
      const tempMap = {[propName]: prop};
      processField({[propName]: prop}, propName, prop);
      prop = tempMap[propName];
    }
  }
  } catch(e: any) {
    grok.shell.error(e.message ?? e);
  }
  return prop;
}