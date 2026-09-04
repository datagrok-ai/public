/* eslint-disable camelcase */
/* eslint-disable prefer-const */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {Plate} from '../plate/plate';
import {Matcher, NumericMatcher} from './matchers';
import {Subject} from 'rxjs';
import {pltsDb, PlateDetailInsert, PlateWellValueInsert, AnalysisRunParameterInsert,
  AnalysisResultInsert, PltsTransactionOp} from '../generated/db';

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
  id: string;
  name: string;
  rows: number;
  cols: number;
  maxVolume?: number;
}

export interface PlateProperty extends DG.IProperty {
  id: string;
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
  id: string;
  name: string;
  description: string;
  plateProperties: Partial<PlateProperty>[];
  wellProperties: Partial<PlateProperty>[];
  required_props: [string, string][];
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
  {id: '', name: 'Generic 96 wells', rows: 8, cols: 12},
  {id: '', name: 'Generic 384 wells', rows: 16, cols: 24},
  {id: '', name: 'Generic 1536 wells', rows: 32, cols: 48},
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


const PAGE_SIZE = 5000;

async function queryAll<T>(page: (limit: number, offset: number) => PromiseLike<T[]>): Promise<T[]> {
  const res: T[] = [];
  for (let offset = 0; ; offset += PAGE_SIZE) {
    const rows = await page(PAGE_SIZE, offset);
    res.push(...rows);
    if (rows.length < PAGE_SIZE)
      return res;
  }
}

let _initialized = false;
export async function initPlates(force: boolean = false) {
  if (_initialized && !force)
    return;

  allProperties = await queryAll((limit, offset) => pltsDb.properties.query()
    .orderBy('created_on').orderBy('id').top(limit).skip(offset)) as unknown as PlateProperty[];
  plateTypes = (await queryAll((limit, offset) => pltsDb.plateTypes.query()
    .orderBy('created_on').orderBy('id').top(limit).skip(offset)))
    .map((pt) => ({id: pt.id, name: pt.name, rows: pt.rows, cols: pt.cols, maxVolume: pt.max_volume}));
  plateTemplates = (await queryAll((limit, offset) => pltsDb.templates.query()
    .orderBy('created_on').orderBy('id').top(limit).skip(offset)))
    .map((t) => ({id: t.id, name: t.name, description: t.description ?? '',
      plateProperties: [], wellProperties: [], required_props: []} as PlateTemplate));

  const templatePropRows = await queryAll((limit, offset) =>
    pltsDb.templateProperties.query().orderBy('id').top(limit).skip(offset));
  for (const template of plateTemplates) {
    const props = templatePropRows
      .filter((tp) => tp.template_id === template.id)
      .map((tp) => ({...allProperties.find((p) => p.id === tp.property_id),
        is_required: tp.is_required ?? false, default_value: tp.default_value ?? null}))
      .filter((p) => p.id != null)
      .sort((a, b) => a.name!.localeCompare(b.name!));

    template.plateProperties = props.filter((p) => p.scope === 'plate');
    template.wellProperties = props.filter((p) => p.scope === 'well');
    template.required_props = props.filter((p) => p.is_required).map((p) => [p.id!, p.name!]);
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


export function findProp(props: PlateProperty[], name: string): PlateProperty {
  return props.find((p) => p.name.toLowerCase() == name.toLowerCase())!;
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
  const runs = await queryAll((limit, offset) => pltsDb.analysisRuns.query()
    .where({analysis_type: analysisType}).select('groups').orderBy('id').top(limit).skip(offset));
  return [...new Set(runs.flatMap((r) => r.groups ?? []))].sort();
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


export async function getPlateById(id: string): Promise<Plate> {
  await initPlates();
  const df = await pltsDb.plateWellValues.query()
    .where({plate_id: id})
    .orderBy('property_id').orderBy('row').orderBy('col')
    .top(1000000).df();
  const plate = Plate.fromDbDataFrame(df);
  plate.id = id;
  events.next({on: 'after', eventType: 'read', objectType: TYPE.PLATE, object: plate});
  return plate;
}


function intersect(a: Set<string> | null, b: Set<string>): Set<string> {
  return a === null ? b : new Set([...a].filter((id) => b.has(id)));
}

function matcherCondition(condition: PropertyCondition) {
  const column = plateDbColumn[condition.property.type] ?? plateDbJsonColumn;
  const matcherCond = condition.matcher.toCondition(column);
  const propCond = DG.cond('property_id', '=', condition.property.id);
  return matcherCond == null ? propCond : DG.and(propCond, matcherCond);
}

const detailsFor = (condition: PropertyCondition) => queryAll((limit, offset) =>
  pltsDb.plateDetails.query().where(matcherCondition(condition))
    .select('plate_id').orderBy('id').top(limit).skip(offset));

const wellValuesFor = (condition: PropertyCondition) => queryAll((limit, offset) =>
  pltsDb.plateWellValues.query().where(matcherCondition(condition))
    .select('plate_id', 'row', 'col').orderBy('id').top(limit).skip(offset));

/** Intersects the key sets matched by each condition; null = no conditions (everything passes). */
async function matchingIds<T>(fetch: (c: PropertyCondition) => Promise<T[]>,
  conditions: PropertyCondition[], key: (r: T) => string): Promise<Set<string> | null> {
  let ids: Set<string> | null = null;
  for (const condition of conditions)
    ids = intersect(ids, new Set((await fetch(condition)).map(key)));
  return ids;
}

async function matchingAnalysisPlateIds(analysisMatchers: AnalysisCondition[]): Promise<Set<string> | null> {
  if (analysisMatchers.length === 0)
    return null;

  const analysisConditionsByType = new Map<string, AnalysisCondition[]>();
  for (const condition of analysisMatchers) {
    if (!analysisConditionsByType.has(condition.analysisName))
      analysisConditionsByType.set(condition.analysisName, []);
    analysisConditionsByType.get(condition.analysisName)!.push(condition);
  }

  let plateIds: Set<string> | null = null;
  for (const [analysisName, conditions] of analysisConditionsByType.entries()) {
    const selectedGroup = conditions.find((c) => c.group)?.group;
    const groups = selectedGroup ? (Array.isArray(selectedGroup) ? selectedGroup : [selectedGroup]) : [];
    const propertyConditions = conditions.filter((c) => c.property);

    let runs = await queryAll((limit, offset) => pltsDb.analysisRuns.query()
      .where({analysis_type: analysisName}).orderBy('id').top(limit).skip(offset));
    if (groups.length > 0 && propertyConditions.length === 0)
      runs = runs.filter((r) => (r.groups ?? []).some((g) => groups.includes(g)));

    let runIds = new Set(runs.map((r) => r.id));
    for (const condition of propertyConditions) {
      const prop = allProperties.find((p) => p.name === condition.property!.name);
      if (!prop) continue;
      const results = await queryAll((limit, offset) => pltsDb.analysisResults.query()
        .where(DG.and(DG.cond('analysis_run_id.analysis_type', '=', analysisName),
          matcherCondition({property: prop, matcher: condition.matcher!})))
        .select('analysis_run_id', 'group_combination').orderBy('id').top(limit).skip(offset));
      const passed = new Set(results
        .filter((r) => groups.length === 0 || (r.group_combination ?? []).some((g) => groups.includes(g)))
        .map((r) => r.analysis_run_id));
      runIds = intersect(runIds, passed);
    }

    const typePlateIds = new Set(runs.filter((r) => runIds.has(r.id)).map((r) => r.plate_id));
    plateIds = intersect(plateIds, typePlateIds);
  }
  return plateIds;
}

function propertyNamesById(): Map<string, string> {
  return new Map(allProperties.map((p) => [p.id, p.name]));
}

async function plateBarcodes(ids: Set<string> | null): Promise<Map<string, string>> {
  const rows = await queryAll((limit, offset) => {
    const q = pltsDb.plates.query().select('barcode');
    if (ids !== null)
      q.where(DG.cond('id', '=', [...ids]));
    return q.orderBy('id').top(limit).skip(offset);
  });
  return new Map(rows.map((p) => [p.id, p.barcode ?? '']));
}

function typedColumn(name: string, values: any[]): DG.Column {
  const sample = values.find((v) => v != null);
  const type = typeof sample === 'number' ?
    (values.every((v) => v == null || Number.isInteger(v)) ? DG.COLUMN_TYPE.INT : DG.COLUMN_TYPE.FLOAT) :
    typeof sample === 'boolean' ? DG.COLUMN_TYPE.BOOL : DG.COLUMN_TYPE.STRING;
  return DG.Column.fromList(type, name, values);
}

function dfFromRows(rows: {[key: string]: any}[]): DG.DataFrame {
  if (rows.length === 0)
    return DG.DataFrame.create();
  return DG.DataFrame.fromColumns(Object.keys(rows[0])
    .map((key) => typedColumn(key, rows.map((r) => r[key] ?? null))));
}

function pivotedDf(rows: ({[key: string]: any} & {props?: {[key: string]: any}})[]): DG.DataFrame {
  if (rows.length === 0)
    return DG.DataFrame.create();
  const df = dfFromRows(rows.map(({props, ...rest}) =>
    ({...rest, '~properties': JSON.stringify(props ?? {})})));
  for (const key of new Set(rows.flatMap((r) => Object.keys(r.props ?? {}))))
    df.columns.add(typedColumn(key, rows.map((r) => r.props?.[key] ?? null)));
  return df;
}

export async function queryWells(query: PlateQuery): Promise<DG.DataFrame> {
  await initPlates();
  const plateIds = await matchingIds(detailsFor, query.plateMatchers, (r) => r.plate_id);
  const wellKey = (r: {plate_id: string; row: number; col: number}) => `${r.plate_id}|${r.row}|${r.col}`;
  const wellIds = await matchingIds(wellValuesFor, query.wellMatchers, wellKey);

  let pivotPlateIds = plateIds;
  if (wellIds !== null)
    pivotPlateIds = intersect(pivotPlateIds, new Set([...wellIds].map((k) => k.split('|')[0])));
  if (pivotPlateIds !== null && pivotPlateIds.size === 0)
    return pivotedDf([]);

  const wellRows = await queryAll((limit, offset) => {
    const q = pltsDb.plateWellValues.query();
    if (pivotPlateIds !== null)
      q.where(DG.cond('plate_id', '=', [...pivotPlateIds]));
    return q.orderBy('id').top(limit).skip(offset);
  });
  const barcodes = await plateBarcodes(pivotPlateIds);
  const propNames = propertyNamesById();

  const wells = new Map<string, {plate_id: string; row: number; col: number; props: {[key: string]: any}}>();
  for (const r of wellRows) {
    const key = wellKey(r);
    if (wellIds !== null && !wellIds.has(key))
      continue;
    const well = wells.get(key) ?? {plate_id: r.plate_id, row: r.row, col: r.col, props: {}};
    const name = propNames.get(r.property_id);
    const value = r.value_string ?? r.value_num ?? r.value_bool;
    if (name && value != null)
      well.props[name] = value;
    wells.set(key, well);
  }

  const sorted = [...wells.values()]
    .sort((a, b) => a.plate_id.localeCompare(b.plate_id) || a.row - b.row || a.col - b.col);
  return pivotedDf(sorted.map((w) => ({plate_id: w.plate_id, barcode: barcodes.get(w.plate_id) ?? '',
    row: w.row, col: w.col, props: w.props})));
}

export async function queryPlates(query: PlateQuery): Promise<DG.DataFrame> {
  await initPlates();
  const idFilters = [
    await matchingIds(detailsFor, query.plateMatchers, (r) => r.plate_id),
    await matchingIds(wellValuesFor, query.wellMatchers, (r) => r.plate_id),
    await matchingAnalysisPlateIds(query.analysisMatchers),
  ];
  let matched: Set<string> | null = null;
  for (const ids of idFilters) {
    if (ids !== null)
      matched = intersect(matched, ids);
  }
  if (matched !== null && matched.size === 0)
    return pivotedDf([]);

  const plateRows = await queryAll((limit, offset) => {
    const q = pltsDb.plates.query();
    if (matched !== null)
      q.where(DG.cond('id', '=', [...matched]));
    return q.orderBy('created_on').orderBy('id').top(limit).skip(offset);
  });
  const detailRows = await queryAll((limit, offset) => {
    const q = pltsDb.plateDetails.query();
    if (matched !== null)
      q.where(DG.cond('plate_id', '=', [...matched]));
    return q.orderBy('id').top(limit).skip(offset);
  });

  const propNames = propertyNamesById();
  const detailsByPlate = new Map<string, {[key: string]: any}>();
  for (const d of detailRows) {
    const name = propNames.get(d.property_id);
    const value = d.value_string ?? d.value_num ?? d.value_bool;
    if (!name || value == null) continue;
    const props = detailsByPlate.get(d.plate_id) ?? {};
    props[name] = value;
    detailsByPlate.set(d.plate_id, props);
  }

  return pivotedDf(plateRows.map((p) => ({plate_id: p.id, barcode: p.barcode ?? '',
    description: p.description ?? '', props: detailsByPlate.get(p.id) ?? {}})));
}


export async function createProperty(prop: Partial<PlateProperty>): Promise<PlateProperty> {
  const report = (await pltsDb.properties.insert({
    name: prop.name!, type: prop.type! as any, scope: prop.scope!,
    choices: prop.choices ? JSON.stringify(prop.choices) : undefined,
    min: prop.min ?? undefined, max: prop.max ?? undefined,
  }))[0];
  prop.id = report.existingId ?? report.id;

  events.next({on: 'after', eventType: 'created', objectType: TYPE.PROPERTY, object: prop});
  allProperties.push(prop as PlateProperty);
  return prop as PlateProperty;
}

export async function savePlate(plate: Plate, options?: { autoCreateProperties?: boolean }) {
  const autoCreateProperties = options?.autoCreateProperties ?? true;
  await initPlates();

  if (!plate.plateTypeId)
    plate.plateTypeId = plateTypes.find((pt) => pt.rows === plate.rows && pt.cols === plate.cols)?.id;
  if (!plate.plateTypeId)
    throw new Error(`No plate type for ${plate.rows}x${plate.cols} plates. Are plate_types seeded in the database?`);

  const globalPlateProperties = allProperties.filter((p) => p.scope === 'plate');
  const globalWellProperties = allProperties.filter((p) => p.scope === 'well');

  for (const layer of Object.keys(plate.details)) {
    const prop = findProp(globalPlateProperties, layer);
    if (autoCreateProperties && !prop) {
      const valueType = getValueType(plate.details[layer]);
      await createProperty({name: layer, type: valueType, scope: 'plate'});
    } else if (!prop) {
      throw new Error(`Global property ${layer} not found`);
    }
  }

  for (const col of plate.data.columns) {
    const prop = findProp(globalWellProperties, col.name);
    if (autoCreateProperties && !prop)
      await createProperty({name: col.name, type: col.type, scope: 'well'});
    else if (!prop)
      throw new Error(`Global property ${col.name} not found`);
  }

  await initPlates(true);

  const detailOps: PltsTransactionOp[] = [];
  for (const layer of Object.keys(plate.details)) {
    const property = allProperties.find((p) => p.scope === 'plate' && p.name.toLowerCase() == layer.toLowerCase());
    if (!property) {
      console.warn(`Property '${layer}' not found in cache. Skipping save for this plate-level property.`);
      continue;
    }
    detailOps.push({op: 'insert', table: 'plate_detail',
      values: {plate_id: '$plate', property_id: property.id,
        ...({[plateDbColumn[property.type]]: plate.details[layer]} as Partial<PlateDetailInsert>)}});
  }

  const results = await pltsDb.transaction([
    {op: 'insert', table: 'plate', ref: 'plate',
      values: {plate_type_id: plate.plateTypeId, barcode: plate.barcode, template_id: plate.plateTemplateId}},
    ...detailOps,
  ]);
  plate.id = results[0].id;

  const wellRows: PlateWellValueInsert[] = [];
  for (const layer of plate.getLayerNames()) {
    const property = allProperties.find((p) => p.scope === 'well' && p.name.toLowerCase() == layer.toLowerCase());
    if (!property) {
      console.warn(`Property '${layer}' not found in cache. Skipping save for this well-level property.`);
      continue;
    }
    const valueColumn = plateDbColumn[property.type];
    for (const pw of plate.wells) {
      wellRows.push({plate_id: plate.id, row: pw.row, col: pw.col, property_id: property.id,
        ...({[valueColumn]: pw[layer] ?? undefined} as Partial<PlateWellValueInsert>)});
    }
  }
  if (wellRows.length > 0) {
    try {
      await pltsDb.plateWellValues.batch(wellRows);
    } catch (e) {
      await pltsDb.plates.delete(plate.id);
      plate.id = undefined;
      throw e;
    }
  }

  events.next({on: 'after', eventType: 'created', objectType: TYPE.PLATE, object: plate});
}


export async function createPlateTemplate(template: PlateTemplateInput): Promise<PlateTemplate> {
  const templateId = (await pltsDb.templates.insert({
    name: template.name, description: template.description || ''}))[0].id;

  const createdPlateProperties: Partial<PlateProperty>[] = [];
  const createdWellProperties: Partial<PlateProperty>[] = [];
  const required_props: [string, string][] = [];

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
    await pltsDb.templateProperties.upsert({template_id: templateId, property_id: newProp.id,
      is_required: isRequired, default_value: property.default_value ?? undefined});

    createdPlateProperties.push(newProp);
    if (isRequired)
      required_props.push([newProp.id, newProp.name]);
  }

  for (const property of template.wellProperties ?? []) {
    property.scope = 'well';
    const newProp = await getOrCreateGlobalProperty(property);

    const isRequired = property.is_required ?? false;
    await pltsDb.templateProperties.upsert({template_id: templateId, property_id: newProp.id,
      is_required: isRequired, default_value: property.default_value ?? undefined});

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
export async function deletePlateTemplate(template: PlateTemplate): Promise<void> {
  await pltsDb.templates.delete(template.id);
  events.next({on: 'after', eventType: 'deleted', objectType: TYPE.TEMPLATE, object: template});
  await initPlates(true); // Refresh the cache
}

export async function createAnalysisRun(plateId: string, analysisType: string, groups: string[]): Promise<string> {
  if (!plateId)
    throw new Error('Cannot create analysis run: plateId is missing.');

  const [report] = await pltsDb.analysisRuns.insert({plate_id: plateId, analysis_type: analysisType, groups: groups});
  return report.id;
}


export async function saveAnalysisRunParameter(params: {
    runId: string,
    propertyName: string,
    propertyType: DG.TYPE,
    value: any
}) {
  const prop = await getOrCreateProperty(params.propertyName, params.propertyType);
  const dbColumn = plateDbColumn[params.propertyType];
  if (!dbColumn)
    throw new Error(`Unsupported property type for parameter: ${params.propertyType}`);

  const values: AnalysisRunParameterInsert = {analysis_run_id: params.runId, property_id: prop.id};
  if (dbColumn === 'value_string') values.value_string = String(params.value);
  else if (dbColumn === 'value_num') values.value_num = params.value;
  else if (dbColumn === 'value_bool') values.value_bool = params.value;

  await pltsDb.analysisRunParameters.insert(values);
}

export async function saveAnalysisResult(params: {
    runId: string,
    propertyId: string,
    propertyName: string,
    propertyType: string,
    value: any,
    groupCombination: string[]
}): Promise<void> {
  const values: AnalysisResultInsert = {analysis_run_id: params.runId, property_id: params.propertyId,
    group_combination: params.groupCombination};

  if (params.propertyName.toLowerCase().includes('curve') && typeof params.value === 'string') {
    values.value_jsonb = params.value;
  } else {
    const dbColumn = plateDbColumn[params.propertyType];
    if (!dbColumn)
      throw new Error(`Unsupported property type for DB storage: ${params.propertyType}`);

    if (dbColumn === 'value_string') values.value_string = String(params.value);
    else if (dbColumn === 'value_num') values.value_num = params.value;
    else if (dbColumn === 'value_bool') values.value_bool = params.value;
  }

  await pltsDb.analysisResults.insert(values);
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

const SEP = String.fromCharCode(1);

export async function queryAnalysesGeneric(query: AnalysisQuery): Promise<DG.DataFrame> {
  await initPlates();
  const selectedGroups = Array.isArray(query.group) ? query.group :
    (query.group ? [query.group] : []);

  const runs = await queryAll((limit, offset) => pltsDb.analysisRuns.query()
    .where({analysis_type: query.analysisName}).orderBy('created_on').orderBy('id').top(limit).skip(offset));
  if (runs.length === 0)
    return dfFromRows([]);
  const results = await queryAll((limit, offset) => pltsDb.analysisResults.query()
    .where(DG.cond<any>('analysis_run_id.analysis_type', '=', query.analysisName))
    .orderBy('id').top(limit).skip(offset));

  const propNames = propertyNamesById();
  type Combo = {runId: string; combo: string[]; props: {[key: string]: any}};
  const comboKey = (runId: string, combo: string[]) => `${runId}${SEP}${combo.join(SEP)}`;
  const combos = new Map<string, Combo>();
  for (const r of results) {
    if (selectedGroups.length > 0 && !(r.group_combination ?? []).some((g) => selectedGroups.includes(g)))
      continue;
    const key = comboKey(r.analysis_run_id, r.group_combination ?? []);
    const combo = combos.get(key) ?? {runId: r.analysis_run_id, combo: r.group_combination ?? [], props: {}};
    const name = propNames.get(r.property_id);
    const value = r.value_jsonb ?? r.value_string ?? r.value_num ?? r.value_bool;
    if (name && value != null)
      combo.props[name] = value;
    combos.set(key, combo);
  }

  let filtered = [...combos.values()];
  for (const condition of query.propertyMatchers) {
    const prop = allProperties.find((p) => p.name === condition.property.name);
    if (!prop) continue;
    const matches = await queryAll((limit, offset) => pltsDb.analysisResults.query()
      .where(DG.and(DG.cond('analysis_run_id.analysis_type', '=', query.analysisName),
        matcherCondition({property: prop, matcher: condition.matcher})))
      .select('analysis_run_id', 'group_combination').orderBy('id').top(limit).skip(offset));
    const passed = new Set(matches.map((r) => comboKey(r.analysis_run_id, r.group_combination ?? [])));
    filtered = filtered.filter((c) => passed.has(comboKey(c.runId, c.combo)));
  }

  const runOrder = new Map(runs.map((r, i) => [r.id, i]));
  const plateByRun = new Map(runs.map((r) => [r.id, r.plate_id]));
  const barcodes = await plateBarcodes(new Set(runs.map((r) => r.plate_id)));
  filtered.sort((a, b) => (runOrder.get(a.runId)! - runOrder.get(b.runId)!) ||
    a.combo.join(', ').localeCompare(b.combo.join(', ')));

  return dfFromRows(filtered.map((c) => ({run_id: c.runId, plate_id: plateByRun.get(c.runId) ?? '',
    barcode: barcodes.get(plateByRun.get(c.runId) ?? '') ?? '',
    group_combination: c.combo.join(', '), properties: JSON.stringify(c.props)})));
}

export async function loadUniquePropertyValues(): Promise<void> {
  await initPlates();
  const stringPropNames = new Map(allProperties
    .filter((p) => p.type === DG.COLUMN_TYPE.STRING).map((p) => [p.id, p.name]));
  const stringPropIds = [...stringPropNames.keys()];
  const toDf = (rows: {property_id: string; value_string?: string}[]) => {
    const pairs = new Set<string>();
    for (const r of rows) {
      const name = stringPropNames.get(r.property_id);
      if (name && r.value_string != null)
        pairs.add(`${name}${SEP}${r.value_string}`);
    }
    const names: string[] = [];
    const values: string[] = [];
    for (const pair of pairs) {
      const [name, value] = pair.split(SEP);
      names.push(name);
      values.push(value);
    }
    return DG.DataFrame.fromColumns([
      DG.Column.fromStrings('name', names), DG.Column.fromStrings('value_string', values)]);
  };
  plateUniquePropertyValues = toDf(stringPropIds.length === 0 ? [] :
    await queryAll((limit, offset) => pltsDb.plateDetails.query()
      .where(DG.cond('property_id', '=', stringPropIds))
      .select('property_id', 'value_string').orderBy('id').top(limit).skip(offset)));
  wellUniquePropertyValues = toDf(stringPropIds.length === 0 ? [] :
    await queryAll((limit, offset) => pltsDb.plateWellValues.query()
      .where(DG.cond('property_id', '=', stringPropIds))
      .select('property_id', 'value_string').orderBy('id').top(limit).skip(offset)));
}
