/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Plate} from '../plate/plate';
import {Matcher} from './numeric_matcher';
import {Subject} from 'rxjs';
import * as api from '../package-api';

export const events: Subject<CrudEvent> = new Subject();

/** Events emitted by the plates CRUD layer. */
export type CrudEvent = {
  on: 'before' | 'after';
  eventType: 'created' | 'read' | 'updated' | 'deleted';
  objectType: TYPE;
  object: Partial<PlateTemplate> | Partial<Plate> | Partial<PlateProperty>;
  cancel?: boolean; // when set by the listener at the 'before' event, the operation is cancelled
}

export enum TYPE {
  PLATE = 'Plate',
  WELL = 'PlateWell',
  TEMPLATE = 'PlateTemplate',
  PROPERTY = 'PlateProperty',
}

export const entityTypes: string[] = [TYPE.PLATE, TYPE.WELL, TYPE.TEMPLATE, TYPE.PROPERTY];

export type PlateType = {
  id: number;
  name: string;
  rows: number;
  cols: number;
  maxVolume?: number;
}

export interface PlateProperty extends DG.IProperty {
  id: number;
  name: string;
  type: string;
  template_id?: number;
  scope?: 'plate' | 'well';
  min?: number;
  max?: number;
}

export type PropertyCondition = {
  property: PlateProperty;
  matcher: Matcher;
}

export type PlateQuery = {
  plateMatchers: PropertyCondition[];
  wellMatchers: PropertyCondition[];
}

export type PlateTemplate = {
  id: number;
  name: string;
  description: string;
  plate_layout_id: number;
  plateProperties: Partial<PlateProperty>[];
  wellProperties: Partial<PlateProperty>[];
}

// export const wellProperties: PlateProperty[] = [
//   {id: 1000, name: 'Volume', type: DG.COLUMN_TYPE.FLOAT},
//   {id: 1001, name: 'Concentration', type: DG.COLUMN_TYPE.FLOAT},
//   {id: 1002, name: 'Sample', type: DG.COLUMN_TYPE.STRING},
//   {id: 1003, name: 'Well Role', type: DG.COLUMN_TYPE.STRING}
// ];

export let allProperties: PlateProperty[] = [];
export let plateTemplates: PlateTemplate[] = [];
export let plateTypes: PlateType[] = [
  {id: 1, name: 'Generic 96 wells', rows: 8, cols: 12},
  {id: 2, name: 'Generic 384 wells', rows: 16, cols: 24},
  {id: 3, name: 'Generic 1536 wells', rows: 32, cols: 48},
];

export const plateUniquePropertyValues: DG.DataFrame = DG.DataFrame.create();
export const wellUniquePropertyValues: DG.DataFrame = DG.DataFrame.create();

export const plateDbColumn: {[key: string]: string} = {
  [DG.COLUMN_TYPE.FLOAT]: 'value_num',
  [DG.COLUMN_TYPE.INT]: 'value_num',
  [DG.COLUMN_TYPE.BOOL]: 'value_bool',
  [DG.COLUMN_TYPE.STRING]: 'value_string',
};


let _initialized = false;
export async function initPlates(force: boolean = false) {
  if (_initialized && !force)
    return;

  if (!_initialized)
    events.subscribe((event) => grok.shell.info(`${event.on} ${event.eventType} ${event.objectType}`));

  // Fetch all templates and all properties
  plateTemplates = (await api.queries.getPlateTemplates()).toJson();
  allProperties = (await api.queries.getProperties()).toJson();
  plateTypes = (await api.queries.getPlateTypes()).toJson();

  // Populate each template with its properties based on template_id and scope
  for (const template of plateTemplates) {
    template.plateProperties = allProperties.filter(
      (p) => p.template_id === template.id && p.scope === 'plate'
    );
    template.wellProperties = allProperties.filter(
      (p) => p.template_id === template.id && p.scope === 'well'
    );
  }

  _initialized = true;
}


/** Creates a new plate for a user to edit. Does not add it to the database. */
export async function createNewPlateForTemplate(plateType: PlateType, plateTemplate: PlateTemplate): Promise<Plate> {
  // ADD THIS CHECK
  if (!plateType)
    throw new Error('Cannot create plate: plateType is undefined. Are plate_types seeded in the database?');

  if (!plateTemplate.plate_layout_id) {
    const plate = new Plate(plateType.rows, plateType.cols);
    for (const property of plateTemplate.wellProperties)
      plate.data.columns.addNew(property.name!, property.type! as DG.ColumnType);

    return plate;
  }

  const p = await getPlateById(plateTemplate.plate_layout_id);
  p.id = undefined;
  return p;
}


export function findProp(props: PlateProperty[], id: number | string): PlateProperty {
  if (typeof id === 'number')
    return props.find((p) => p.id == id)!;
  else
    return props.find((p) => p.name.toLowerCase() == id.toLowerCase())!;
}

export function getUniquePropertyValues(prop: PlateProperty, df: DG.DataFrame): string[] {
  const nameCol = df.col('name')!;
  const valueCol = df.col('value_string')!;

  return df.rows
    .where((i) => nameCol.get(i) == prop.name)
    .map((i) => valueCol.get(i))
    .toArray();
}


export function getPlateUniquePropertyValues(prop: PlateProperty): string[] {
  return getUniquePropertyValues(prop, plateUniquePropertyValues);
}

export function getWellUniquePropertyValues(prop: PlateProperty): string[] {
  return getUniquePropertyValues(prop, wellUniquePropertyValues);
}

function getValueType(x: any): string {
  if (typeof x === 'string')
    return DG.TYPE.STRING;
  if (typeof x === 'number')
    return DG.TYPE.FLOAT;
  if (typeof x === 'boolean')
    return DG.TYPE.BOOL;
  throw new Error('Not supported type');
}


function sqlStr(s?: string | number): string {
  if (!s)
    return 'null';
  if (typeof s == 'string')
    return `'${s}'`;
  else
    return `${s}`;
}


export async function getPlateById(id: number): Promise<Plate> {
  const df: DG.DataFrame = await grok.functions.call('Curves:getWellValuesById', {id: id});
  const plate = Plate.fromDbDataFrame(df);
  events.next({on: 'after', eventType: 'read', objectType: TYPE.PLATE, object: plate});
  return plate;
}


/**
 * Build SQL that
 *   1) keeps only plates whose plate-level properties match the user's filters
 *   2) returns every plate-level property (name + value) for those plates
 *      packed into a single JSONB column called "properties".
 *
 *   The JSON looks like:
 *   {
 *     "Volume":        "100.0",
 *     "Concentration": "10.5",
 *     "Sample":        "ABC-123",
 *     "Well Role":     "DMSO",
 *     ...
 *   }
 */
function getPlateSearchSql(query: PlateQuery): string {
  // ---------- 1. Build the boolean conditions for the filter ----------
  const existsClauses: string[] = [];

  // Handle plate-level property conditions
  for (const condition of query.plateMatchers) {
    const dbColumn = plateDbColumn[condition.property.type];
    existsClauses.push(`
      EXISTS (
        SELECT 1
        FROM plates.plate_details pd_f
        WHERE pd_f.plate_id = p.id
          AND pd_f.property_id = ${condition.property.id}
          AND (${condition.matcher.toSql(`pd_f.${dbColumn}`)})
      )`);
  }

  // Handle well-level property conditions
  for (const condition of query.wellMatchers) {
    const dbColumn = plateDbColumn[condition.property.type];
    existsClauses.push(`
      EXISTS (
        SELECT 1
        FROM plates.plate_well_values pwv
        WHERE pwv.plate_id = p.id
          AND pwv.property_id = ${condition.property.id}
          AND (${condition.matcher.toSql(`pwv.${dbColumn}`)})
      )`);
  }

  const whereFilter =
    existsClauses.length > 0 ? `WHERE ${existsClauses.join(' AND ')}` : '';

  // ---------- 2. Final query – filter first, then collect all properties ----------
  return `
WITH filtered_plates AS (
  SELECT p.id, p.barcode, p.description
  FROM plates.plates p
  ${whereFilter}
)
SELECT
  fp.id         AS plate_id,
  fp.barcode,
  fp.description,
  COALESCE(
    jsonb_object_agg(
      pr.name,
      CASE
        WHEN pd.value_string   IS NOT NULL THEN to_jsonb(pd.value_string)
        WHEN pd.value_num      IS NOT NULL THEN to_jsonb(pd.value_num)
        WHEN pd.value_bool     IS NOT NULL THEN to_jsonb(pd.value_bool)
        WHEN pd.value_uuid     IS NOT NULL THEN to_jsonb(pd.value_uuid)
        WHEN pd.value_datetime IS NOT NULL THEN to_jsonb(pd.value_datetime)
        ELSE to_jsonb(NULL::text)
      END
    ) FILTER (WHERE pr.name IS NOT NULL),
    '{}'::jsonb
  )::text AS properties
FROM filtered_plates fp
LEFT JOIN plates.plate_details pd ON pd.plate_id = fp.id
LEFT JOIN plates.properties     pr ON pr.id = pd.property_id
GROUP BY fp.id, fp.barcode, fp.description
ORDER BY fp.id;
`.trim();
}

/**
 * Build SQL that
 *   1) keeps only wells whose plate-level and well-level properties match the user's filters
 *   2) returns every well-level property (name + value) for those wells
 *      packed into a single JSONB column called "properties".
 *
 *   The JSON looks like:
 *   {
 *     "Volume":        "100.0",
 *     "Concentration": "10.5",
 *     "Sample":        "ABC-123",
 *     "Well Role":     "DMSO",
 *     ...
 *   }
 */
function getWellSearchSql(query: PlateQuery): string {
  // ---------- 1. Build the boolean conditions for the filter ----------
  const whereClauses: string[] = [];

  // Handle plate-level property conditions (using EXISTS for parent plate)
  for (const condition of query.plateMatchers) {
    const dbColumn = plateDbColumn[condition.property.type];
    whereClauses.push(`
      EXISTS (
        SELECT 1
        FROM plates.plate_details pd_f
        WHERE pd_f.plate_id = pwv.plate_id
          AND pd_f.property_id = ${condition.property.id}
          AND (${condition.matcher.toSql(`pd_f.${dbColumn}`)})
      )`);
  }

  // Handle well-level property conditions (direct filter on wells)
  for (const condition of query.wellMatchers) {
    const dbColumn = plateDbColumn[condition.property.type];
    whereClauses.push(`
      EXISTS (
        SELECT 1
        FROM plates.plate_well_values pwv_f
        WHERE pwv_f.plate_id = pwv.plate_id
          AND pwv_f.row = pwv.row
          AND pwv_f.col = pwv.col
          AND pwv_f.property_id = ${condition.property.id}
          AND (${condition.matcher.toSql(`pwv_f.${dbColumn}`)})
      )`);
  }

  const whereFilter =
    whereClauses.length > 0 ? `WHERE ${whereClauses.join(' AND ')}` : '';

  // ---------- 2. Final query – filter first, then collect all well properties ----------
  return `
WITH filtered_plates AS (
  SELECT p.id, p.barcode, p.description
  FROM plates.plates p
  ${whereFilter}
)
SELECT
  fp.id         AS plate_id,
  fp.barcode,
  fp.description,
  COALESCE(
    jsonb_object_agg(
      pr.name,
      CASE
        WHEN pd.value_string   IS NOT NULL THEN to_jsonb(pd.value_string)
        WHEN pd.value_num      IS NOT NULL THEN to_jsonb(pd.value_num)
        WHEN pd.value_bool     IS NOT NULL THEN to_jsonb(pd.value_bool)
        WHEN pd.value_uuid     IS NOT NULL THEN to_jsonb(pd.value_uuid)
        WHEN pd.value_datetime IS NOT NULL THEN to_jsonb(pd.value_datetime)
        ELSE to_jsonb(NULL::text)
      END
    ) FILTER (WHERE pr.name IS NOT NULL),
    '{}'::jsonb
  )::text AS properties
FROM filtered_plates fp
LEFT JOIN plates.plate_details pd ON pd.plate_id = fp.id
LEFT JOIN plates.properties     pr ON pr.id      = pd.property_id
GROUP BY fp.id, fp.barcode, fp.description
ORDER BY fp.id;
`.trim();
}

export async function queryWells(query: PlateQuery): Promise<DG.DataFrame> {
  const df = await grok.data.db.query('Curves:Plates', getWellSearchSql(query));
  DG.Utils.jsonToColumns(df.col('properties')!);
  df.col('properties')!.name = '~properties';
  return df;
}

export async function queryPlates(query: PlateQuery): Promise<DG.DataFrame> {
  const df = await grok.data.db.query('Curves:Plates', getPlateSearchSql(query));
  DG.Utils.jsonToColumns(df.col('properties')!);
  df.col('properties')!.name = '~properties';
  return df;
}

export async function createProperty(prop: Partial<PlateProperty>): Promise<PlateProperty> {
  // Need to pass ALL property fields to the backend
  prop.id = await grok.functions.call('Curves:createProperty', {
    propertyName: prop.name!,
    valueType: prop.type!,
    templateId: prop.template_id,
    scope: prop.scope,
    choices: prop.choices ? JSON.stringify(prop.choices) : null, // Add this
    min: prop.min, // Add this
    max: prop.max, // Add this
    // Add other fields as needed
  });

  events.next({on: 'after', eventType: 'created', objectType: TYPE.PROPERTY, object: prop});
  allProperties.push(prop as PlateProperty);
  return prop as PlateProperty;
}

/**
 * Saves the plate to the database.
 * When a property is not found in the plateProperties array, it is created automatically
 * if autoCreateProperties is true; otherwise, an error is thrown.
*/
export async function savePlate(plate: Plate, options?: { autoCreateProperties?: boolean }) {
  const autoCreateProperties = options?.autoCreateProperties ?? true;
  await initPlates();

  const plateSql =
    `insert into plates.plates(plate_type_id, barcode)
     values(${plate.plateTypeId}, ${sqlStr(plate.barcode)})
     returning id`;
  plate.id = (await grok.data.db.query('Curves:Plates', plateSql)).get('id', 0);

  // Get global properties (template_id is null) from the unified cache
  const globalPlateProperties = allProperties.filter((p) => p.template_id == null && p.scope === 'plate');
  const globalWellProperties = allProperties.filter((p) => p.template_id == null && p.scope === 'well');

  // FIRST: register new GLOBAL plate level properties
  for (const layer of Object.keys(plate.details)) {
    const prop = findProp(globalPlateProperties, layer);
    if (autoCreateProperties && !prop) {
      const valueType = getValueType(plate.details[layer]);
      // Create a global (template_id: null) property with a 'plate' scope
      await createProperty({name: layer, type: valueType, scope: 'plate'});
      grok.shell.info('Global plate property created: ' + layer);
    } else if (!prop) {
      throw new Error(`Global property ${layer} not found`);
    }
  }

  // SECOND: register new GLOBAL well level properties
  for (const layer of plate.getLayerNames()) {
    const prop = findProp(globalWellProperties, layer);
    if (autoCreateProperties && !prop) {
      const col = plate.data.col(layer);
      // Create a global (template_id: null) property with a 'well' scope
      await createProperty({name: layer, type: col!.type, scope: 'well'});
      grok.shell.info('Global well property created: ' + layer);
    } else if (!prop) {
      throw new Error(`Global property ${layer} not found`);
    }
  }

  // THIRD: NOW execute the SQL
  const sql = getPlateInsertSql(plate);
  await grok.data.db.query('Curves:Plates', sql);

  events.next({on: 'after', eventType: 'created', objectType: TYPE.PLATE, object: plate});
  grok.shell.info('Plate saved');
}

export async function savePlateAsTemplate(plate: Plate, template: PlateTemplate) {
  await savePlate(plate);
  const sql = `update plates.templates set plate_layout_id = ${plate.id} where id = ${template.id}`;
  await grok.data.db.query('Curves:Plates', sql);
}


function getPlateInsertSql(plate: Plate): string {
  let sql = 'insert into plates.plate_wells(plate_id, row, col) values ' +
    plate.wells.map((pw) => `  (${plate.id}, ${pw.row}, ${pw.col})`).toArray().join(',\n') + ';';

  // Find properties from the unified 'allProperties' cache, filtering by scope
  for (const layer of Object.keys(plate.details)) {
    const property = allProperties.find((p) => p.scope === 'plate' && p.name.toLowerCase() == layer.toLowerCase())!;
    const dbCol = plateDbColumn[property.type];

    sql += `\n insert into plates.plate_details(plate_id, property_id, ${dbCol}) values ` +
      `(${plate.id}, ${property.id}, ${sqlStr(plate.details[layer])});`;
  }

  // well data
  for (const layer of plate.getLayerNames()) {
    const property = allProperties.find((p) => p.scope === 'well' && p.name.toLowerCase() == layer.toLowerCase())!;
    const dbCol = plateDbColumn[property.type];
    sql += `\n\n insert into plates.plate_well_values(plate_id, row, col, property_id, ${dbCol}) values\n` +
      plate.wells
        .map((pw) => `  (${plate.id}, ${pw.row}, ${pw.col}, ${property.id}, ${sqlStr(pw[layer])})`)
        .toArray().join(',\n') + ';';
  }

  return sql;
}


export async function createPlateTemplate(template: Partial<PlateTemplate>): Promise<PlateTemplate> {
  template.id = await grok.functions.call('Curves:createTemplate', {
    name: template.name,
    description: template.description
  });

  const createdPlateProperties: PlateProperty[] = [];
  const createdWellProperties: PlateProperty[] = [];

  // Create plate properties with ALL their fields
  for (const property of template.plateProperties ?? []) {
    property.template_id = template.id;
    property.scope = 'plate';
    const newProp = await createProperty(property); // This should now pass all fields
    createdPlateProperties.push(newProp);
  }

  // Create well properties with ALL their fields
  for (const property of template.wellProperties ?? []) {
    property.template_id = template.id;
    property.scope = 'well';
    const newProp = await createProperty(property);
    createdWellProperties.push(newProp);
  }

  // Assign the complete property objects back
  template.plateProperties = createdPlateProperties;
  template.wellProperties = createdWellProperties;

  plateTemplates.push(template as PlateTemplate);
  events.next({on: 'after', eventType: 'created', objectType: TYPE.TEMPLATE, object: template});
  return template as PlateTemplate;
}
/** Checks if a plate template's properties are being used in any plates */
export async function plateTemplatePropertiesUsed(template: PlateTemplate): Promise<boolean> {
  // Check if any plates are using this template's properties by looking for
  // property IDs that belong to this template.
  const sql = `
    SELECT EXISTS (
      SELECT 1 FROM plates.plate_details pd
      WHERE pd.property_id IN (SELECT id FROM plates.properties WHERE template_id = ${template.id})
      UNION
      SELECT 1 FROM plates.plate_well_values pwv
      WHERE pwv.property_id IN (SELECT id FROM plates.properties WHERE template_id = ${template.id})
    ) as is_used`;

  const result = await grok.data.db.query('Curves:Plates', sql);
  return result.get('is_used', 0);
}

/** Deletes a plate template and its associated properties */
export async function deletePlateTemplate(template: PlateTemplate): Promise<void> {
  // The database cascade on the foreign key from properties to templates
  // should handle deleting the properties. If not, you'd delete them first.
  // A simple approach is to delete properties owned by the template, then the template itself.
  const deleteSql = `
    DELETE FROM plates.properties WHERE template_id = ${template.id};
    DELETE FROM plates.templates WHERE id = ${template.id};
  `;

  await grok.data.db.query('Curves:Plates', deleteSql);
  events.next({on: 'after', eventType: 'deleted', objectType: TYPE.TEMPLATE, object: template});
  await initPlates(true); // Refresh the cache
}
