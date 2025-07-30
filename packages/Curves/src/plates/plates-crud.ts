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
  type: string; //DG.COLUMN_TYPE;
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

export let wellProperties: PlateProperty[] = [
  {id: 1000, name: 'Volume', type: DG.COLUMN_TYPE.FLOAT},
  {id: 1001, name: 'Concentration', type: DG.COLUMN_TYPE.FLOAT},
  {id: 1002, name: 'Sample', type: DG.COLUMN_TYPE.STRING},
  {id: 1003, name: 'Well Role', type: DG.COLUMN_TYPE.STRING}
];

export let plateProperties: PlateProperty[] = [];
export let plateTemplates: PlateTemplate[] = [];

export let plateTypes: PlateType[] = [
  {id: 1, name: 'Generic 96 wells', rows: 8, cols: 12},
  {id: 2, name: 'Generic 384 wells', rows: 16, cols: 24},
  {id: 3, name: 'Generic 1536 wells', rows: 32, cols: 48},
];

export let plateUniquePropertyValues: DG.DataFrame = DG.DataFrame.create();
export let wellUniquePropertyValues: DG.DataFrame = DG.DataFrame.create();

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

  plateTemplates = (await api.queries.getPlateTemplates()).toJson();
  plateProperties = (await api.queries.getPlateLevelProperties()).toJson();
  wellProperties = (await api.queries.getWellLevelProperties()).toJson();
  plateTypes = (await api.queries.getPlateTypes()).toJson();

  plateUniquePropertyValues = await api.queries.getUniquePlatePropertyValues();
  wellUniquePropertyValues = await api.queries.getUniqueWellPropertyValues();

  const templateWellProperties = (await api.queries.getTemplateWellProperties()).toJson();
  const templatePlateProperties = (await api.queries.getTemplatePlateProperties()).toJson();
  for (const template of plateTemplates) {
    template.wellProperties = templateWellProperties
      .filter((p) => p.template_id == template.id).map((p) => findProp(wellProperties, p.property_id));
    template.plateProperties = templatePlateProperties
      .filter((p) => p.template_id == template.id).map((p) => findProp(plateProperties, p.property_id));
  }

  _initialized = true;
}


/** Creates a new plate for a user to edit. Does not add it to the database. */
export async function createNewPlateForTemplate(plateType: PlateType, plateTemplate: PlateTemplate): Promise<Plate> {
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
  jsonb_pretty(
    jsonb_object_agg(
      pr.name,
      CASE
        WHEN pd.value_string   IS NOT NULL THEN to_jsonb(pd.value_string)
        WHEN pd.value_num      IS NOT NULL THEN to_jsonb(pd.value_num)
        WHEN pd.value_bool     IS NOT NULL THEN to_jsonb(pd.value_bool)
        WHEN pd.value_uuid     IS NOT NULL THEN to_jsonb(pd.value_uuid)
        WHEN pd.value_datetime IS NOT NULL THEN to_jsonb(pd.value_datetime)
        ELSE to_jsonb(NULL::text)               -- should never happen
      END
    )
  ) AS properties
FROM filtered_plates fp
JOIN plates.plate_details pd ON pd.plate_id = fp.id
JOIN plates.properties     pr ON pr.id      = pd.property_id
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
WITH filtered_wells AS (
  SELECT DISTINCT pwv.plate_id, pwv.row, pwv.col
  FROM plates.plate_well_values pwv
  JOIN plates.plates p ON p.id = pwv.plate_id
  ${whereFilter}
)
SELECT
  fw.plate_id,
  p.barcode,
  p.description AS plate_description,
  fw.row,
  fw.col,
  jsonb_pretty(
    jsonb_object_agg(
      pr.name,
      CASE
        WHEN pwv.value_string IS NOT NULL THEN to_jsonb(pwv.value_string)
        WHEN pwv.value_num    IS NOT NULL THEN to_jsonb(pwv.value_num)
        WHEN pwv.value_bool   IS NOT NULL THEN to_jsonb(pwv.value_bool)
        ELSE to_jsonb(NULL::text)               -- should never happen
      END
    )
  ) AS properties
FROM filtered_wells fw
JOIN plates.plates p ON p.id = fw.plate_id
JOIN plates.plate_well_values pwv ON pwv.plate_id = fw.plate_id
                                   AND pwv.row = fw.row
                                   AND pwv.col = fw.col
JOIN plates.properties pr ON pr.id = pwv.property_id
GROUP BY fw.plate_id, p.barcode, p.description, fw.row, fw.col
ORDER BY fw.plate_id, fw.row, fw.col;
`.trim();
}

export async function queryWells(query: PlateQuery): Promise<DG.DataFrame> {
  const df = await grok.data.db.query('Admin:Plates', getWellSearchSql(query));
  DG.Utils.jsonToColumns(df.col('properties')!);
  df.col('properties')!.name = '~properties';
  return df;
}

export async function queryPlates(query: PlateQuery): Promise<DG.DataFrame> {
  const df = await grok.data.db.query('Admin:Plates', getPlateSearchSql(query));
  DG.Utils.jsonToColumns(df.col('properties')!);
  df.col('properties')!.name = '~properties';
  return df;
}

export async function createProperty(prop: Partial<PlateProperty>): Promise<PlateProperty> {
  // prop.id = await api.queries.createProperty(prop.name!, prop.type!);
  prop.id = await grok.functions.call('Curves:createProperty', {propertyName: prop.name!, valueType: prop.type});
  events.next({on: 'after', eventType: 'created', objectType: TYPE.PROPERTY, object: prop});
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
  plate.id = (await grok.data.db.query('Admin:Plates', plateSql)).get('id', 0);

  // register new plate level properties
  for (const layer of Object.keys(plate.details)) {
    const prop = findProp(plateProperties, layer);
    if (autoCreateProperties && !prop) {
      const valueType = getValueType(plate.details[layer]);
      plateProperties.push(await createProperty({name: layer, type: valueType}));
      grok.shell.info('Plate layer created: ' + layer);
    } else if (!prop) { throw new Error(`Property ${layer} not found in plateProperties`); }
  }

  // register new well level properties
  for (const layer of plate.getLayerNames()) {
    const prop = findProp(wellProperties, layer);
    if (autoCreateProperties && !prop) {
      const col = plate.data.col(layer);
      wellProperties.push(await createProperty({name: layer, type: col!.type}));
      grok.shell.info('Well layer created: ' + layer);
    } else if (!prop) { throw new Error(`Property ${layer} not found in wellProperties`); }
  }

  await grok.data.db.query('Admin:Plates', getPlateInsertSql(plate));
  events.next({on: 'after', eventType: 'created', objectType: TYPE.PLATE, object: plate});
  grok.shell.info('Plate saved');
}

export async function savePlateAsTemplate(plate: Plate, template: PlateTemplate) {
  await savePlate(plate);
  const sql = `update plates.templates set plate_layout_id = ${plate.id} where id = ${template.id}`;
  await grok.data.db.query('Admin:Plates', sql);
}


function getPlateInsertSql(plate: Plate): string {
  let sql = 'insert into plates.plate_wells(plate_id, row, col) values ' +
    plate.wells.map((pw) => `  (${plate.id}, ${pw.row}, ${pw.col})`).toArray().join(',\n') + ';';
  for (const layer of Object.keys(plate.details)) {
    const property = plateProperties.find((p) => p.name.toLowerCase() == layer.toLowerCase())!;
    const dbCol = plateDbColumn[property.type];

    sql += `\n insert into plates.plate_details(plate_id, property_id, ${dbCol}) values ` +
      `(${plate.id}, ${property.id}, ${sqlStr(plate.details[layer])});`;
  }

  // well data
  for (const layer of plate.getLayerNames()) {
    const property = wellProperties.find((p) => p.name.toLowerCase() == layer.toLowerCase())!;
    const dbCol = plateDbColumn[property.type];
    sql += `\n\n insert into plates.plate_well_values(plate_id, row, col, property_id, ${dbCol}) values\n` +
      plate.wells
        .map((pw) => `  (${plate.id}, ${pw.row}, ${pw.col}, ${property.id}, ${sqlStr(pw[layer])})`)
        .toArray().join(',\n') + ';';
  }

  return sql;
}

export async function createPlateTemplate(template: Partial<PlateTemplate>): Promise<PlateTemplate> {
  const results = await grok.functions.call('Curves:createTemplate', {name: template.name, description: template.description})!;

  template.id = results;

  let sql = '';

  for (let property of template.plateProperties ?? []) {
    property = await createProperty(property);
    plateProperties.push(property as PlateProperty);
    sql += `insert into plates.template_plate_properties(template_id, property_id) values (${template.id}, ${property.id});\n`;
  }

  for (let property of template.wellProperties ?? []) {
    property = await createProperty(property);
    wellProperties.push(property as PlateProperty);
    sql += `insert into plates.template_well_properties(template_id, property_id) values (${template.id}, ${property.id});\n`;
  }

  await grok.data.db.query('Admin:Plates', sql);
  plateTemplates.push(template as PlateTemplate);

  events.next({on: 'after', eventType: 'created', objectType: TYPE.TEMPLATE, object: template});
  return template as PlateTemplate;
}

/** Checks if a plate template's properties are being used in any plates */
export async function plateTemplatePropertiesUsed(template: PlateTemplate): Promise<boolean> {
  // Check if any plates are using this template's properties
  const sql = `
    SELECT EXISTS (
      SELECT 1 FROM plates.plates p
      JOIN plates.plate_details pd ON p.id = pd.plate_id
      JOIN plates.template_plate_properties tpp ON pd.property_id = tpp.property_id
      WHERE tpp.template_id = ${template.id}
      UNION
      SELECT 1 FROM plates.plates p
      JOIN plates.plate_well_values pwv ON p.id = pwv.plate_id
      JOIN plates.template_well_properties twp ON pwv.property_id = twp.property_id
      WHERE twp.template_id = ${template.id}
    ) as is_used`;

  const result = await grok.data.db.query('Admin:Plates', sql);
  return result.get('is_used', 0);
}

/** Deletes a plate template and its associated properties */
export async function deletePlateTemplate(template: PlateTemplate): Promise<void> {
  // First delete the template's property associations
  const deletePropsSql = `
    DELETE FROM plates.template_plate_properties WHERE template_id = ${template.id};
    DELETE FROM plates.template_well_properties WHERE template_id = ${template.id};
    DELETE FROM plates.templates WHERE id = ${template.id};`;

  await grok.data.db.query('Admin:Plates', deletePropsSql);
  events.next({on: 'after', eventType: 'deleted', objectType: TYPE.TEMPLATE, object: template});
  await initPlates(true); // Refresh the templates cache
}
