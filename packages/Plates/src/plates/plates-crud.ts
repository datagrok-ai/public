/* eslint-disable camelcase */
/* eslint-disable prefer-const */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Plate} from '../plate/plate';
import {Matcher, NumericMatcher} from './matchers';
import {Subject} from 'rxjs';
import * as api from '../package-api';
import {AnalysisManager} from '../plate/analyses/analysis-manager';

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
  // template_id?: number;
  scope?: 'plate' | 'well';
  min?: number;
  max?: number;
}

export type PropertyCondition = {
  property: PlateProperty;
  matcher: Matcher | NumericMatcher;
}

export type PlateQuery = {
  plateMatchers: PropertyCondition[];
  wellMatchers: PropertyCondition[];
  analysisMatchers: AnalysisCondition[];
}

export type PlateTemplate = {
  id: number;
  name: string;
  description: string;
  plateProperties: Partial<PlateProperty>[];
  wellProperties: Partial<PlateProperty>[];
  required_props: [number, string][];
}
export type PlateTemplateInput = Partial<Omit<PlateTemplate, 'plateProperties' | 'wellProperties' | 'required_props' | 'id'>> & {
  name: string;
  plateProperties?: TemplatePropertyInput[];
  wellProperties?: TemplatePropertyInput[];
}
export type TemplatePropertyInput = Partial<PlateProperty> & {
  is_required?: boolean;
  default_value?: string | null;
}
export type AnalysisProperty = {
  name: string; // User-friendly name, e.g., "IC50"
  type: DG.TYPE; // Data type, e.g., DG.TYPE.FLOAT
}
export type AnalysisQuery = {
  analysisName: string;
  propertyMatchers: PropertyCondition[];
  group?: string | string[];
}

export type AnalysisCondition = {
  property?: AnalysisProperty;
  matcher?: Matcher;
  analysisName: string;
  group?: string | string[];
}


export let allProperties: PlateProperty[] = [];
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
export const plateDbJsonColumn = 'value_jsonb';


let _initialized = false;
export async function initPlates(force: boolean = false) {
  if (_initialized && !force)
    return;

  if (!_initialized)
    events.subscribe((event) => grok.shell.info(`${event.on} ${event.eventType} ${event.objectType}`));

  plateTemplates = (await api.queries.getPlateTemplates()).toJson();
  allProperties = (await api.queries.getProperties()).toJson();
  plateTypes = (await api.queries.getPlateTypes()).toJson();

  for (const template of plateTemplates) {
    const platePropsDf = await api.queries.getTemplatePlateProperties(template.id);
    // await grok.functions.call('Curves:getTemplatePlateProperties', {templateId: template.id});
    const wellPropsDf = await api.queries.getTemplateWellProperties(template.id);
    // await grok.functions.call('Curves:getTemplateWellProperties', {templateId: template.id});

    const plateProps: any[] = platePropsDf.toJson();
    const wellProps: any[] = wellPropsDf.toJson();

    template.plateProperties = plateProps;
    template.wellProperties = wellProps;

    template.required_props = [];
    for (const prop of [...plateProps, ...wellProps]) {
      if (prop.is_required)
        template.required_props.push([prop.id, prop.name]);
    }
  }

  _initialized = true;
}


/** Creates a new plate for a user to edit. Does not add it to the database. */
export async function createNewPlateForTemplate(plateType: PlateType, plateTemplate: PlateTemplate): Promise<Plate> {
  if (!plateType)
    throw new Error('Cannot create plate: plateType is undefined. Are plate_types seeded in the database?');

  // Always create a new empty plate based on the plate type
  const plate = new Plate(plateType.rows, plateType.cols);
  // Add columns for well properties
  for (const property of plateTemplate.wellProperties)
    plate.data.columns.addNew(property.name!, property.type! as DG.ColumnType);

  return plate;
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

export async function getAnalysisRunGroups(analysisType: string): Promise<string[]> {
  await initPlates();
  const df = await api.queries.getAnalysisRunGroups(analysisType);
  // await grok.functions.call('Curves:getAnalysisRunGroups', {analysisType: analysisType});
  return df.col('group')?.toList() ?? [];
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
  const df: DG.DataFrame = await api.queries.getWellValuesById(id);
  // await grok.functions.call('Curves:getWellValuesById', {id: id});
  const plate = Plate.fromDbDataFrame(df);
  events.next({on: 'after', eventType: 'read', objectType: TYPE.PLATE, object: plate});
  return plate;
}


function getPlateSearchSql(query: PlateQuery): string {
  const existsClauses: string[] = [];

  for (const condition of query.plateMatchers) {
    const dbColumn = plateDbColumn[condition.property.type];
    existsClauses.push(`
      EXISTS (
        SELECT 1 FROM plts.plate_details pd_f
        WHERE pd_f.plate_id = p.id
        AND pd_f.property_id = ${condition.property.id}
        AND (${condition.matcher.toSql(`pd_f.${dbColumn}`)})
      )`);
  }

  for (const condition of query.wellMatchers) {
    const dbColumn = plateDbColumn[condition.property.type];
    existsClauses.push(`
      EXISTS (
        SELECT 1 FROM plts.plate_well_values pwv
        WHERE pwv.plate_id = p.id
        AND pwv.property_id = ${condition.property.id}
        AND (${condition.matcher.toSql(`pwv.${dbColumn}`)})
      )`);
  }

  const analysisConditionsByType = new Map<string, AnalysisCondition[]>();
  for (const condition of query.analysisMatchers) {
    if (!analysisConditionsByType.has(condition.analysisName))
      analysisConditionsByType.set(condition.analysisName, []);
    analysisConditionsByType.get(condition.analysisName)!.push(condition);
  }

  for (const [analysisName, conditions] of analysisConditionsByType.entries()) {
    const groupCondition = conditions.find((c) => c.group);
    const selectedGroup = groupCondition?.group;
    const propertyConditions = conditions.filter((c) => c.property);

    let groupsArray: string = '';
    let hasGroupFilter = false;

    if (selectedGroup) {
      const groups = Array.isArray(selectedGroup) ? selectedGroup : [selectedGroup];
      if (groups.length > 0) {
        hasGroupFilter = true;
        groupsArray = groups.map((g) => `'${g.replace(/'/g, '\'\'')}'`).join(',');
      }
    }
    let analysisSubClauses: string[] = [];

    for (const condition of propertyConditions) {
      const prop = allProperties.find((p) => p.name === condition.property!.name);
      if (!prop) continue;

      const dbColumn = plateDbColumn[prop.type] ?? plateDbJsonColumn;

      let subClause = `
        EXISTS (
          SELECT 1 FROM plts.analysis_results res
          WHERE res.analysis_run_id = ar.id
          AND res.property_id = ${prop.id}
          AND (${condition.matcher!.toSql(`res.${dbColumn}`)})
      `;

      if (hasGroupFilter)
        subClause += ` AND res.group_combination && ARRAY[${groupsArray}]`;


      subClause += `)`;
      analysisSubClauses.push(subClause);
    }

    let analysisClause = `
      EXISTS (
        SELECT 1 FROM plts.analysis_runs ar
        WHERE ar.plate_id = p.id AND ar.analysis_type = '${analysisName}'
    `;

    if (hasGroupFilter && propertyConditions.length === 0)
      analysisClause += ` AND ar.groups && ARRAY[${groupsArray}]`;


    if (analysisSubClauses.length > 0)
      analysisClause += ` AND ${analysisSubClauses.join(' AND ')}`;

    analysisClause += `)`;
    existsClauses.push(analysisClause);
  }

  const whereFilter = existsClauses.length > 0 ? `WHERE ${existsClauses.join(' AND ')}` : '';

  return `
WITH filtered_plates AS (
    SELECT p.id, p.barcode, p.description
    FROM plts.plates p
    ${whereFilter}
)
SELECT
    fp.id AS plate_id,
    fp.barcode,
    fp.description,
    COALESCE(
        jsonb_object_agg(
            pr.name,
            CASE
                WHEN pd.value_string IS NOT NULL THEN to_jsonb(pd.value_string)
                WHEN pd.value_num IS NOT NULL THEN to_jsonb(pd.value_num)
                WHEN pd.value_bool IS NOT NULL THEN to_jsonb(pd.value_bool)
                ELSE 'null'::jsonb
            END
        ) FILTER (WHERE pr.id IS NOT NULL),
        '{}'::jsonb
    )::text AS properties
FROM filtered_plates fp
LEFT JOIN plts.plate_details pd ON pd.plate_id = fp.id
LEFT JOIN plts.properties pr ON pr.id = pd.property_id
GROUP BY fp.id, fp.barcode, fp.description
ORDER BY fp.id;
    `.trim();
}


function getWellSearchSql(query: PlateQuery): string {
  const whereClauses: string[] = [];

  // Handle plate-level property conditions
  for (const condition of query.plateMatchers) {
    const dbColumn = plateDbColumn[condition.property.type];
    whereClauses.push(`
      EXISTS (
        SELECT 1
        FROM plts.plate_details pd_f
        WHERE pd_f.plate_id = pwv.plate_id
          AND pd_f.property_id = ${condition.property.id}
          AND (${condition.matcher.toSql(`pd_f.${dbColumn}`)})
      )`);
  }

  // Handle well-level property conditions
  for (const condition of query.wellMatchers) {
    const dbColumn = plateDbColumn[condition.property.type];
    whereClauses.push(`
      EXISTS (
        SELECT 1
        FROM plts.plate_well_values pwv_f
        WHERE pwv_f.plate_id = pwv.plate_id
          AND pwv_f.row = pwv.row
          AND pwv_f.col = pwv.col
          AND pwv_f.property_id = ${condition.property.id}
          AND (${condition.matcher.toSql(`pwv_f.${dbColumn}`)})
      )`);
  }

  const whereFilter = whereClauses.length > 0 ? `WHERE ${whereClauses.join(' AND ')}` : '';

  return `
SELECT
    pwv.plate_id,
    p.barcode,
    pwv.row,
    pwv.col,
    COALESCE(jsonb_object_agg(
        pr.name,
        CASE
            WHEN pwv.value_string IS NOT NULL THEN to_jsonb(pwv.value_string)
            WHEN pwv.value_num    IS NOT NULL THEN to_jsonb(pwv.value_num)
            WHEN pwv.value_bool   IS NOT NULL THEN to_jsonb(pwv.value_bool)
            ELSE 'null'::jsonb
        END
    ) FILTER (WHERE pr.id IS NOT NULL), '{}'::jsonb)::text AS properties
FROM
    plts.plate_well_values pwv
JOIN
    plts.plates p ON p.id = pwv.plate_id
LEFT JOIN
    plts.properties pr ON pr.id = pwv.property_id
${whereFilter}
GROUP BY
    pwv.plate_id, p.barcode, pwv.row, pwv.col
ORDER BY
    pwv.plate_id, pwv.row, pwv.col;
`.trim();
}

export async function queryWells(query: PlateQuery): Promise<DG.DataFrame> {
  const df = await grok.data.db.query('Plates:Plts', getWellSearchSql(query));
  DG.Utils.jsonToColumns(df.col('properties')!);
  df.col('properties')!.name = '~properties';
  return df;
}

export async function queryPlates(query: PlateQuery): Promise<DG.DataFrame> {
  const df = await grok.data.db.query('Plates:Plts', getPlateSearchSql(query));


  const propsCol = df.col('properties');
  if (propsCol && df.rowCount > 0) {
    try {
      DG.Utils.jsonToColumns(propsCol);
      propsCol.name = '~properties';
    } catch (e) {
      console.error('Failed to parse properties JSON:', e);
    }
  }

  return df;
}


export async function createProperty(prop: Partial<PlateProperty>): Promise<PlateProperty> {
  prop.id = await api.queries.createProperty(prop.name!, prop.type!, prop.scope!, prop.choices ? JSON.stringify(prop.choices) : null, prop.min ?? null, prop.max ?? null);
  // await grok.functions.call('Curves:createProperty', {
  //   propertyName: prop.name!,
  //   valueType: prop.type!,
  //   scope: prop.scope,
  //   choices: prop.choices ? JSON.stringify(prop.choices) : null,
  //   min: prop.min,
  //   max: prop.max,
  // });

  events.next({on: 'after', eventType: 'created', objectType: TYPE.PROPERTY, object: prop});
  allProperties.push(prop as PlateProperty);
  return prop as PlateProperty;
}

export async function savePlate(plate: Plate, options?: { autoCreateProperties?: boolean }) {
  const autoCreateProperties = options?.autoCreateProperties ?? true;
  await initPlates();

  const plateSql =
 `insert into plts.plates(plate_type_id, barcode)
  values(${plate.plateTypeId}, ${sqlStr(plate.barcode)})
 returning id`;
  plate.id = (await grok.data.db.query('Plates:Plts', plateSql)).get('id', 0);

  const globalPlateProperties = allProperties.filter((p) => p.scope === 'plate');
  const globalWellProperties = allProperties.filter((p) => p.scope === 'well');

  for (const layer of Object.keys(plate.details)) {
    const prop = findProp(globalPlateProperties, layer);
    if (autoCreateProperties && !prop) {
      const valueType = getValueType(plate.details[layer]);
      await createProperty({name: layer, type: valueType, scope: 'plate'});
      // grok.shell.info('Global plate property created: ' + layer);
    } else if (!prop) {
      throw new Error(`Global property ${layer} not found`);
    }
  }

  for (const col of plate.data.columns) {
    const prop = findProp(globalWellProperties, col.name);
    if (autoCreateProperties && !prop)
      await createProperty({name: col.name, type: col.type, scope: 'well'});
      // grok.shell.info('Global well property created: ' + col.name);
    else if (!prop)
      throw new Error(`Global property ${col.name} not found`);
  }

  await initPlates(true);
  const sql = getPlateInsertSql(plate);
  await grok.data.db.query('Plates:Plts', sql);

  events.next({on: 'after', eventType: 'created', objectType: TYPE.PLATE, object: plate});
  // grok.shell.info('Plate saved');
}

// This should be expanded on as a part of the "Template from data" feature, but currently uses obsolete layout field.
// export async function savePlateAsTemplate(plate: Plate, template: PlateTemplate) {
//   await savePlate(plate);
//   const sql = `update plts.templates set plate_layout_id = ${plate.id} where id = ${template.id}`;
//   await grok.data.db.query('Plates:Plts', sql);
// }


function getPlateInsertSql(plate: Plate): string {
  let sql = 'insert into plts.plate_wells(plate_id, row, col) values ' +
        plate.wells.map((pw) => `  (${plate.id}, ${pw.row}, ${pw.col})`).toArray().join(',\n') + ';';

  for (const layer of Object.keys(plate.details)) {
    const property = allProperties.find((p) => p.scope === 'plate' && p.name.toLowerCase() == layer.toLowerCase());

    if (!property) {
      console.warn(`Property '${layer}' not found in cache. Skipping save for this plate-level property.`);
      continue;
    }

    const dbCol = plateDbColumn[property.type];
    sql += `\n insert into plts.plate_details(plate_id, property_id, ${dbCol}) values ` +
           `(${plate.id}, ${property.id}, ${sqlStr(plate.details[layer])});`;
  }

  for (const layer of plate.getLayerNames()) {
    const property = allProperties.find((p) => p.scope === 'well' && p.name.toLowerCase() == layer.toLowerCase());

    if (!property) {
      console.warn(`Property '${layer}' not found in cache. Skipping save for this well-level property.`);
      continue;
    }

    const dbCol = plateDbColumn[property.type];
    sql += `\n\n insert into plts.plate_well_values(plate_id, row, col, property_id, ${dbCol}) values\n` +
           plate.wells
             .map((pw) => `  (${plate.id}, ${pw.row}, ${pw.col}, ${property.id}, ${sqlStr(pw[layer])})`)
             .toArray().join(',\n') + ';';
  }

  return sql;
}


export async function createPlateTemplate(template: PlateTemplateInput): Promise<PlateTemplate> {
  const templateId = await api.queries.createTemplate(template.name, template.description || '');
  // await grok.functions.call('Curves:createTemplate', {
  //   name: template.name,
  //   description: template.description
  // });

  const createdPlateProperties: Partial<PlateProperty>[] = [];
  const createdWellProperties: Partial<PlateProperty>[] = [];
  const required_props: [number, string][] = [];

  const getOrCreateGlobalProperty = async (prop: Partial<PlateProperty>): Promise<PlateProperty> => {
    await initPlates();
    let existingProp = allProperties.find((p) => p.name === prop.name && p.scope === prop.scope);
    if (existingProp)
      return existingProp;
    return await createProperty(prop);
  };

  for (const property of template.plateProperties ?? []) {
    property.scope = 'plate';
    const newProp = await getOrCreateGlobalProperty(property);

    const isRequired = property.is_required ?? false;
    await api.queries.addTemplatePlateProperty(templateId, newProp.id, isRequired, property.default_value ?? null);
    // await grok.functions.call('Curves:addTemplatePlateProperty', {
    //   templateId: templateId,
    //   propertyId: newProp.id,
    //   isRequired: isRequired,
    //   defaultValue: property.default_value ?? null
    // });

    createdPlateProperties.push(newProp);
    if (isRequired)
      required_props.push([newProp.id, newProp.name]);
  }

  for (const property of template.wellProperties ?? []) {
    property.scope = 'well';
    const newProp = await getOrCreateGlobalProperty(property);

    const isRequired = property.is_required ?? false;

    await api.queries.addTemplateWellProperty(templateId, newProp.id, isRequired, property.default_value ?? null);
    // await grok.functions.call('Curves:addTemplateWellProperty', {
    //   templateId: templateId,
    //   propertyId: newProp.id,
    //   isRequired: isRequired,
    //   defaultValue: property.default_value ?? null
    // });

    createdWellProperties.push(newProp);
    if (isRequired)
      required_props.push([newProp.id, newProp.name]);
  }

  const finalTemplate: PlateTemplate = {
    id: templateId,
    name: template.name!,
    description: template.description || '',
    plateProperties: createdPlateProperties,
    wellProperties: createdWellProperties,
    required_props: required_props,
  };

  plateTemplates.push(finalTemplate);
  events.next({on: 'after', eventType: 'created', objectType: TYPE.TEMPLATE, object: finalTemplate});
  return finalTemplate;
}
/** Checks if a plate template's properties are being used in any plates */
export async function plateTemplatePropertiesUsed(template: PlateTemplate): Promise<boolean> {
  // Check if any plates are using this template's properties by looking for
  // property IDs that belong to this template.
  const sql = `
    SELECT EXISTS (
      SELECT 1 FROM plts.plate_details pd
      WHERE pd.property_id IN (SELECT id FROM plts.properties WHERE template_id = ${template.id})
      UNION
      SELECT 1 FROM plts.plate_well_values pwv
      WHERE pwv.property_id IN (SELECT id FROM plts.properties WHERE template_id = ${template.id})
    ) as is_used`;

  const result = await grok.data.db.query('Plates:Plts', sql);
  return result.get('is_used', 0);
}

/** Deletes a plate template and its associated properties */
export async function deletePlateTemplate(template: PlateTemplate): Promise<void> {
  // The database cascade on the foreign key from properties to templates
  // should handle deleting the properties. If not, you'd delete them first.
  // A simple approach is to delete properties owned by the template, then the template itself.
  const deleteSql = `
    DELETE FROM plts.properties WHERE template_id = ${template.id};
    DELETE FROM plts.templates WHERE id = ${template.id};
  `;

  await grok.data.db.query('Plates:Plts', deleteSql);
  events.next({on: 'after', eventType: 'deleted', objectType: TYPE.TEMPLATE, object: template});
  await initPlates(true); // Refresh the cache
}

export async function createAnalysisRun(plateId: number, analysisType: string, groups: string[]): Promise<number> {
  if (!plateId)
    throw new Error('Cannot create analysis run: plateId is missing.');

  const runId: number = await api.queries.createAnalysisRun(plateId, analysisType, groups);
  // await grok.functions.call('Curves:createAnalysisRun', {
  //   plateId: plateId,
  //   analysisType: analysisType,
  //   groups: groups,
  // });

  if (typeof runId !== 'number' || runId <= 0)
    throw new Error(`Failed to create a valid analysis run. Received runId: ${runId}`);

  return runId;
}


export async function saveAnalysisRunParameter(params: {
    runId: number,
    propertyName: string,
    propertyType: DG.TYPE,
    value: any
}) {
  const prop = await getOrCreateProperty(params.propertyName, params.propertyType);
  const callParams: any = {
    analysisRunId: params.runId,
    propertyId: prop.id,
    valueString: null, valueNum: null, valueBool: null, valueJsonb: null
  };

  const dbColumnKey = Object.keys(plateDbColumn).find((key) => key === params.propertyType);
  if (!dbColumnKey)
    throw new Error(`Unsupported property type for parameter: ${params.propertyType}`);

  const dbColumn = plateDbColumn[dbColumnKey];
  if (dbColumn === 'value_string') callParams.valueString = String(params.value);
  else if (dbColumn === 'value_num') callParams.valueNum = params.value;
  else if (dbColumn === 'value_bool') callParams.valueBool = params.value;

  await api.queries.saveAnalysisRunParameter(callParams.analysisRunId, callParams.propertyId, callParams.valueString, callParams.valueNum, callParams.valueBool, callParams.valueJsonb);
  // await grok.functiosns.call('Curves:saveAnalysisRunParameter', callParams);
}

export async function saveAnalysisResult(params: {
    runId: number,
    propertyId: number,
    propertyName: string,
    propertyType: string,
    value: any,
    groupCombination: string[]
}): Promise<void> {
  const callParams: any = {
    analysisRunId: params.runId,
    propertyId: params.propertyId,
    groupCombination: params.groupCombination,
    valueString: null, valueNum: null, valueBool: null, valueJsonb: null
  };

  if (params.propertyName.toLowerCase().includes('curve') && typeof params.value === 'string') {
    callParams.valueJsonb = params.value;
  } else {
    const dbColumnKey = Object.keys(plateDbColumn).find((key) => key === params.propertyType);
    if (!dbColumnKey)
      throw new Error(`Unsupported property type for DB storage: ${params.propertyType}`);

    const dbColumn = plateDbColumn[dbColumnKey];
    if (dbColumn === 'value_string') callParams.valueString = String(params.value);
    else if (dbColumn === 'value_num') callParams.valueNum = params.value;
    else if (dbColumn === 'value_bool') callParams.valueBool = params.value;
  }

  await api.queries.saveAnalysisResult(callParams.analysisRunId, callParams.groupCombination, callParams.propertyId, callParams.valueString, callParams.valueNum, callParams.valueBool, callParams.valueJsonb);
  // await grok.functions.call('Curves:saveAnalysisResult', callParams);
}

export async function getOrCreateProperty(name: string, type: DG.TYPE, scope: 'plate' | 'well' = 'plate'): Promise<PlateProperty> {
  await initPlates();

  let prop = allProperties.find((p) => p.name === name && p.scope == scope);
  if (prop)
    return prop;

  const newProp = await createProperty({name: name, type: type, scope: scope});
  await initPlates(true);
  return allProperties.find((p) => p.id === newProp.id)!;
}

export async function queryAnalysesGeneric(query: AnalysisQuery): Promise<DG.DataFrame> {
  await initPlates();
  const selectedGroups = Array.isArray(query.group) ? query.group :
    (query.group ? [query.group] : []);

  let sqlQuery: string;

  if (selectedGroups.length > 0) {
    const groupsArray = selectedGroups.map((g) => `'${g.replace(/'/g, '\'\'')}'`).join(',');

    const crossJoinWithGroupFilter = `
      CROSS JOIN LATERAL (
        SELECT
          res.group_combination,
          jsonb_object_agg(
            prop.name,
            CASE
              WHEN res.value_jsonb IS NOT NULL THEN to_jsonb(res.value_jsonb::text)
              WHEN res.value_string IS NOT NULL THEN to_jsonb(res.value_string)
              WHEN res.value_num IS NOT NULL THEN to_jsonb(res.value_num)
              WHEN res.value_bool IS NOT NULL THEN to_jsonb(res.value_bool)
              ELSE 'null'::jsonb
            END
          )::text as properties
        FROM
          plts.analysis_results res
        JOIN
          plts.properties prop ON res.property_id = prop.id
        WHERE
          res.analysis_run_id = ar.id
          AND res.group_combination && ARRAY[${groupsArray}]
        GROUP BY
          res.group_combination
      ) as res_pivot`;

    const whereClauses: string[] = [];
    whereClauses.push(`ar.analysis_type = '${query.analysisName.replace(/'/g, '\'\'')}'`);

    for (const condition of query.propertyMatchers) {
      const prop = allProperties.find((p) => p.name === condition.property.name);
      if (!prop) continue;
      const dbColumn = plateDbColumn[prop.type] ?? plateDbJsonColumn;

      const existsClause = `
        EXISTS (
            SELECT 1
            FROM plts.analysis_results filter_res
            WHERE filter_res.analysis_run_id = ar.id
              AND filter_res.group_combination = res_pivot.group_combination
              AND filter_res.property_id = ${prop.id}
              AND (${condition.matcher.toSql(`filter_res.${dbColumn}`)})
        )`;

      whereClauses.push(existsClause);
    }

    const finalWhereClause = whereClauses.length > 0 ? `WHERE ${whereClauses.join(' AND ')}` : '';

    sqlQuery = `
      SELECT
        ar.id as run_id,
        p.id as plate_id,
        p.barcode,
        array_to_string(res_pivot.group_combination, ', ') as group_combination,
        res_pivot.properties
      FROM
        plts.analysis_runs ar
      JOIN
        plts.plates p ON ar.plate_id = p.id
      ${crossJoinWithGroupFilter}
      ${finalWhereClause}
      ORDER BY ar.id, res_pivot.group_combination;
    `;
  } else {
    const whereClauses: string[] = [];
    whereClauses.push(`ar.analysis_type = '${query.analysisName.replace(/'/g, '\'\'')}'`);

    for (const condition of query.propertyMatchers) {
      const prop = allProperties.find((p) => p.name === condition.property.name);
      if (!prop) continue;
      const dbColumn = plateDbColumn[prop.type] ?? plateDbJsonColumn;

      const existsClause = `
        EXISTS (
            SELECT 1
            FROM plts.analysis_results filter_res
            WHERE filter_res.analysis_run_id = ar.id
              AND filter_res.group_combination = res_pivot.group_combination
              AND filter_res.property_id = ${prop.id}
              AND (${condition.matcher.toSql(`filter_res.${dbColumn}`)})
        )`;

      whereClauses.push(existsClause);
    }

    const finalWhereClause = whereClauses.length > 0 ? `WHERE ${whereClauses.join(' AND ')}` : '';

    sqlQuery = `
      SELECT
        ar.id as run_id,
        p.id as plate_id,
        p.barcode,
        array_to_string(res_pivot.group_combination, ', ') as group_combination,
        res_pivot.properties
      FROM
        plts.analysis_runs ar
      JOIN
        plts.plates p ON ar.plate_id = p.id
      CROSS JOIN LATERAL (
        SELECT
          res.group_combination,
          jsonb_object_agg(
            prop.name,
            CASE
              WHEN res.value_jsonb IS NOT NULL THEN to_jsonb(res.value_jsonb::text)
              WHEN res.value_string IS NOT NULL THEN to_jsonb(res.value_string)
              WHEN res.value_num IS NOT NULL THEN to_jsonb(res.value_num)
              WHEN res.value_bool IS NOT NULL THEN to_jsonb(res.value_bool)
              ELSE 'null'::jsonb
            END
          )::text as properties
        FROM
          plts.analysis_results res
        JOIN
          plts.properties prop ON res.property_id = prop.id
        WHERE
          res.analysis_run_id = ar.id
        GROUP BY
          res.group_combination
      ) as res_pivot
      ${finalWhereClause}
      ORDER BY ar.id, res_pivot.group_combination;
    `;
  }

  console.log('Generated SQL Query:', sqlQuery);

  try {
    const df = await grok.data.db.query('Plates:Plts', sqlQuery);
    console.log('Query returned dataframe with', df.rowCount, 'rows');
    if (df.rowCount > 0)
      console.log('Group combinations in result:', df.getCol('group_combination')?.toList());

    return df;
  } catch (error) {
    console.error('Query failed:', error);
    throw error;
  }
}

// fallback
export async function queryAnalyses(query: AnalysisQuery): Promise<DG.DataFrame> {
  const analysis = AnalysisManager.instance.getAnalysis(query.analysisName);
  if (analysis && 'queryResults' in analysis)
    return await (analysis as any).queryResults(query);
  return queryAnalysesGeneric(query);
}
