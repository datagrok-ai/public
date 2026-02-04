import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { SignalsSearchParams, SignalsSearchQuery } from './signals-search-query';
import * as ui from 'datagrok-api/ui';
import { getRevvityUsers, getUserStringIdById } from './users';
import { queryEntityById, RevvityApiResponse, RevvityData, RevvityUser, search } from './revvity-api';

import { awaitCheck } from '@datagrok-libraries/test/src/test';
import { compoundTypeAndViewNameMapping, ENTITY_FIELDS_TO_EXCLUDE, FIELDS_SECTION_NAME, FIELDS_TO_EXCLUDE_FROM_CORPORATE_ID_WIDGET, FIELDS_TO_EXCLUDE_FROM_WIDGET, FIRST_COL_NAMES, LAST_COL_NAMES, MOL_COL_NAME, MOLECULAR_FORMULA_FIELD_NAME, NAME, PARAMS_KEY, QUERY_KEY, STORAGE_NAME, SUBMITTER_FIELD_NAME, TABS_TO_EXCLUDE_FROM_WIDGET, TAGS_TO_EXCLUDE, USER_FIELDS } from './constants';
import { getRevvityLibraries } from './libraries';
import { u2 } from '@datagrok-libraries/utils/src/u2';
import { _package } from './package';
import { currentQueryBuilderConfig, filterProperties, runSearchQuery } from './search-utils';
import { funcs } from './package-api';
import { ComplexCondition, Operators } from '@datagrok-libraries/utils/src/query-builder/query-builder';
import { openedView, updateView } from './view-utils';
import dayjs from 'dayjs';


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

export async function processAttribute(attributeName: string, value: any) {
  //handle user fields
  if (USER_FIELDS.includes(attributeName)) {
    return await getUserStringIdById(value as string);
  }
  //handle arrays considering duplicates
  if (Array.isArray(value)) {
    return Array.from(new Set([1, 1, 2])).join(',');
  }
  return value;
}

export function createRevvityWidgetByEntityId(response: RevvityApiResponse, idSemValue: DG.SemanticValue<string>): HTMLDivElement {
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
    const propDiv = ui.divH([mainMap[key]], { style: { 'position': 'relative' } });
    if (idSemValue.cell?.dataFrame && idSemValue.cell?.column) {
      const addColumnIcon = ui.iconFA('plus', async () => {
       // calculatePropForColumn(key, idSemValue.cell.dataFrame, idSemValue.cell.column);
      }, `Calculate ${key} for the whole table`);
      addColumnIcon.classList.add('revvity-add-prop-as-column-icon');
      propDiv.prepend(addColumnIcon);
      ui.tools.setHoverVisibility(propDiv, [addColumnIcon]);
    }
    widgetMapWithAddIcons[key] = propDiv;
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

export function createRevvityWidgetByCorporateId(response: RevvityApiResponse, idSemValue: DG.SemanticValue<string>): HTMLElement {
  let data = response.data;

    const createAttributesPane = async (item: RevvityData<any>, div: HTMLDivElement, entityType: string) => {
    //move all tags to attributes
    if (item.attributes.tags) {
      for (const [key, value] of Object.entries(item.attributes.tags)) {
        if (FIELDS_TO_EXCLUDE_FROM_CORPORATE_ID_WIDGET.includes(key)) continue;
        const splittedKey = key.split('.');
        const keyWithoutTagName = splittedKey.length > 1 ? splittedKey[1] : splittedKey[0];
        item.attributes[keyWithoutTagName] = value;
      }
    }
    delete item.attributes['tags'];

    //create attributes object for widget
    const fieldsForWidget: { [key: string]: any } = {};
    const fieldsToTypesMap: { [key: string]: 'string' | 'double' | 'bool' } = {};
    for (const [key, value] of Object.entries(item.attributes)) {
      if (FIELDS_TO_EXCLUDE_FROM_CORPORATE_ID_WIDGET.includes(key)) continue;
      fieldsForWidget[key] = await processAttribute(key, value);
      fieldsToTypesMap[key] = typeof value === 'number' ? 'double' : typeof value === 'boolean' ? 'bool' : 'string';
    }

    const widgetMapWithAddIcons: { [key: string]: HTMLElement } = {};
    Object.keys(fieldsForWidget).forEach((key) => {
      const propDiv = ui.divH([fieldsForWidget[key]], { style: { 'position': 'relative' } });
      if (idSemValue.cell?.dataFrame && idSemValue.cell?.column) {
        const addColumnIcon = ui.iconFA('plus', async () => {
          calculatePropForColumn(key, fieldsToTypesMap[key], entityType, idSemValue.cell.dataFrame, idSemValue.cell.column);
        }, `Calculate ${key} for the whole table`);
        addColumnIcon.classList.add('revvity-add-prop-as-column-icon');
        propDiv.prepend(addColumnIcon);
        ui.tools.setHoverVisibility(propDiv, [addColumnIcon]);
      }
      widgetMapWithAddIcons[key] = propDiv;
    });

    const table = ui.tableFromMap(widgetMapWithAddIcons);
    ui.setUpdateIndicator(div, false);
    div.append(table);
  }

  if (!data || Array.isArray(data) && !data.length)
    return ui.divText('No data available');
  if (!Array.isArray(data))
    data = [data];
  if (data.length === 1) {
    const itemType = data[0].attributes.type ?? data[0].attributes.eid ? data[0].attributes.eid.split(':')[0] : null;
    if (!itemType)
        return ui.divText('Entity of unknown type');;
    const div = ui.div();
    createAttributesPane(data[0], div, itemType);
    return div;
  }

  const acc = ui.accordion();

  for (let item of data) {
    if (item.attributes) {
      const itemType = item.attributes.type ?? item.attributes.eid ? item.attributes.eid.split(':')[0] : null;
      if (!itemType)
        continue;
      acc.addPane(itemType, () => {
        const div = ui.div();
        ui.setUpdateIndicator(div, true);
        createAttributesPane(item, div, itemType);
        return div;
      })
    }
  }

  return acc.root;
}

export async function transformData(data: Record<string, any>[], libId: string, entityType: string): Promise<DG.DataFrame> {
  if (!data.length)
    return DG.DataFrame.create();
  const users = await getRevvityUsers();
  const columnsData: {[key: string]: {type: string, data: any[]}} = {};
  const createCol = (key: string, counter: number) => {
    //find column type
    const props = filterProperties[`${libId}|${entityType}`];
    const prop = props.filter((it) => it.name === key);
    if (!columnsData[key]) {
      columnsData[key] = { type: prop.length ? prop[0].type : DG.TYPE.STRING, data: [] };
      for (let j = 0; j < counter; j++)
        columnsData[key].data.push(undefined);
    }
  }
  for (let i = 0; i < data.length; i++) {
    if (!columnsData.id)
      columnsData.id = {type: DG.TYPE.STRING, data: []};
    columnsData.id.data.push(data[i].id);

    const attrs = data[i].attributes || {};

    for (const [key, value] of Object.entries(attrs)) {
      if (ENTITY_FIELDS_TO_EXCLUDE.includes(key))
        continue;
      createCol(key, i);
      //handle fileds related to users
      if (USER_FIELDS.includes(key)) {
        let userArr = users?.filter((it) => it.userId === value);
        let userString = 'unknown user';
        if (userArr?.length) {
          const user = userArr[0];
          userString = user.email ?? user.userName ??
          (user.firstName && user.lastName ? `${user.firstName} ${user.lastName}` : user.userId ?? 'unknown user');
        }
        columnsData[key].data.push(userString);
        continue;
      }
      if (Array.isArray(value)) {
        columnsData[key].data.push(value.join('; '));
      } else {
        columnsData[key].data.push(value);
      }
    }

    // Handle tags
    if (attrs.tags) {
      for (const [tagKey, tagValue] of Object.entries(attrs.tags)) {
        if (TAGS_TO_EXCLUDE.includes(tagKey))
          continue;
      createCol(tagKey, i);
        if (Array.isArray(tagValue)) {
          columnsData[tagKey].data.push(tagValue.join('; '));
        } else {
          columnsData[tagKey].data.push(tagValue);
        }
      }
    }
    
    //in case item lacked some of the keys, fill in corresponding columns with undefined values
    for (const key of Object.keys(columnsData)) {
      if (columnsData[key].data.length < i + 1)
        columnsData[key].data.push(undefined);
    }
  }

  const columns: DG.Column[] = [];
  for (let colName of Object.keys(columnsData)) {
    if (columnsData[colName].type === 'datetime')
      columnsData[colName].data = columnsData[colName].data.map((it) => dayjs(it));
    const col = DG.Column.fromList(columnsData[colName].type as 'string' | 'int' | 'double' | 'bool' | 'qnum' | 'datetime',
      colName, columnsData[colName].data);
    columns.push(col);
  }
  return DG.DataFrame.fromColumns(columns);
}

export async function reorderColumns(df: DG.DataFrame) {
  const colNames = df.columns.names().map((it) => it.toLowerCase());
  const createColNamesArr = (colsToReorder: string[]) => {
    const reorderedNames = [];
    for (const colName of colsToReorder) {
      const index = colNames.indexOf(colName.toLowerCase());
      if (index > -1) {
        colNames.splice(index, 1);
        reorderedNames.push(colName);
      }
    }
    return reorderedNames;
  } 
  const firstCols = createColNamesArr(FIRST_COL_NAMES);
  const lastCols = createColNamesArr(LAST_COL_NAMES);

  df.columns.setOrder(firstCols.concat(colNames).concat(lastCols));
}

export function getCompoundTypeByViewName(viewName: string): string {
  const res = compoundTypeAndViewNameMapping.filter((it) => it.viewName === viewName);
  return !res.length ? viewName : res[0].compoundType;
}

export function getViewNameByCompoundType(compoundType: string): string {
  const res = compoundTypeAndViewNameMapping.filter((it) => it.compoundType === compoundType);
  return !res.length ? compoundType : res[0].viewName;
}

export async function calculatePropForColumn(propName: string, colType: "string" | "double" | "bool", entityType: string, df: DG.DataFrame, idsCol: DG.Column) {
  const newColName = df.columns.getUnusedName(propName);
  const newCol = df.columns.addNew(newColName, colType);
  const promises = idsCol.toList().map((id, index) =>
    getEntityPropByCorporateEntityId(propName, id, entityType)
      .then(res => newCol.set(index, res.toString()))
      .catch(() => { })
  );
  await Promise.allSettled(promises);
}

export async function getEntityPropByCorporateEntityId(propName: string, name: string, type: string): Promise<any> {
  const query = {
    "query": {
      "$and": [
        {
          "$match": {
            "field": "type",
            "value": type,
            "mode": "keyword"
          }
        },
        {
          "$match": {
            "field": "name",
            "value": name,
            "mode": "keyword"
          }
        }
      ]
    }
  }
  let prop: any = null;
  try {
    const entity = await search(query);
    if (entity.data) {
      const entityData = Array.isArray(entity.data) ? entity.data[0] : entity.data;
      const attributes = (entityData as RevvityData).attributes;
      if (attributes) {
        //look for prop in attributes
        if (Object.keys(attributes).includes(propName))
          prop = await processAttribute(propName, attributes[propName]);
        else if (attributes.tags) { //look for prop in tags
          for (const tag of Object.keys(attributes.tags)) {
            const splittedTag = tag.split('.');
            const tagName = splittedTag.length > 1 ? splittedTag[1] : splittedTag[0];
            if (tagName === propName) {
              prop = await processAttribute(propName, attributes.tags[tag]);
              break;
            }
          }
        }
      }
    }
  } catch (e: any) {
    grok.shell.error(e.message ?? e);
  }
  return prop;
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

export function getAppHeader(): HTMLElement {
  const appHeader = u2.appHeader({
      iconPath: _package.webRoot + '/img/revvity.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/RevvitySignalsLink/README.md',
      description: '- Connect to your Revvity account\n' +
        '- Run, save, and share searches\n' +
        '- Visualize & analyze assay data: assets, batches, and more',
      appTitle: 'RevvityLink',
      appSubTitle: 'Access and analyze your assay data',
      bottomLine: true
    });
  return appHeader;
}

let latestLabel: string | undefined = undefined; //workaround since widget is re-created once search is performed and dataframe in the tv is chenged
export async function createWidgetByRevvityLabel(idSemValue: DG.SemanticValue<string>) {

  let div = ui.divText(`No data found for ${idSemValue.value}`);
  let searchConfig = currentQueryBuilderConfig;
  if (!searchConfig) {
    const libs = (await getRevvityLibraries()).filter((it) => it.name.toLowerCase() === 'compounds');
    if (libs.length)
      searchConfig = { libId: libs[0].id, libName: 'compounds', entityType: 'batch' } //use compounds|batches as default
  }
  if (idSemValue.cell?.column) { //now can search only through materials library - the last true parameter
    const terms: string[] = await funcs.getTermsForField(idSemValue.cell?.column.name,
      searchConfig!.entityType, searchConfig!.libId, true);
    if (terms.length) {
      const initValue = latestLabel ?? idSemValue.value;
      latestLabel = undefined;
      if (terms.includes(initValue)) {
        const inputName = idSemValue.cell?.column.getTag('friendlyName') ?? idSemValue.cell?.column.name;

        const labelInput = ui.input.choice(inputName, { value: initValue, items: terms, nullable: false });

        const searchButton = ui.button('Search', async () => {
          latestLabel = labelInput.value!;
          const viewToUpdate = openedView ?? grok.shell.tv;
          if (viewToUpdate) {
            ui.setUpdateIndicator(viewToUpdate.root, true, `Searching ${inputName} = ${labelInput.value}`);
            const queryBuilderCondition: ComplexCondition = {
              logicalOperator: Operators.Logical.and,
              conditions: [
                { field: idSemValue.cell?.column.name, operator: Operators.EQ, value: labelInput.value }]
            };
            const df = await runSearchQuery(searchConfig!.libId, searchConfig!.entityType, queryBuilderCondition);
            const filtersDiv = Array.from(document.getElementsByClassName('revvity-signals-filter-panel'));
            updateView(viewToUpdate as DG.TableView, df, searchConfig!.entityType, searchConfig!.libName,
              searchConfig!.libId, filtersDiv.length ? filtersDiv[0] as HTMLDivElement : undefined);
            ui.setUpdateIndicator(viewToUpdate.root, false);
            if (searchConfig?.qb)
              searchConfig?.qb.loadCondition(queryBuilderCondition);
          }
        });

        div = ui.divV([
          labelInput,
          searchButton
        ]);
      }
    }
  }

  return div;
}
