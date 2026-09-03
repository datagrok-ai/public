/* eslint-disable camelcase */
/* eslint-disable prefer-const */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Plate} from '../plate/plate';
import {Matcher, NumericMatcher} from './matchers';
import {Subject} from 'rxjs';
import {pltsDb, PlateDetailsInsert, PlateWellValuesInsert, AnalysisRunParametersInsert,
  AnalysisResultsInsert, PltsTransactionOp} from '../generated/db';

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


let _initialized = false;
export async function initPlates(force: boolean = false) {
  if (_initialized && !force)
    return;

  if (!_initialized)
    events.subscribe((event) => grok.shell.info(`${event.on} ${event.eventType} ${event.objectType}`));

  allProperties = await pltsDb.propertieses.query().orderBy('created_on').top(10000) as unknown as PlateProperty[];
  plateTypes = (await pltsDb.plateTypeses.query().orderBy('created_on').top(1000))
    .map((pt) => ({id: pt.id, name: pt.name, rows: pt.rows, cols: pt.cols, maxVolume: pt.max_volume}));
  plateTemplates = (await pltsDb.templateses.query().orderBy('created_on').top(10000))
    .map((t) => ({id: t.id, name: t.name, description: t.description ?? '',
      plateProperties: [], wellProperties: [], required_props: []} as PlateTemplate));

  const templatePropRows = await pltsDb.templatePropertieses.query().top(100000);
  for (const template of plateTemplates) {
    const props = templatePropRows
      .filter((tp) => tp.template_id === template.id)
      .map((tp) => ({...allProperties.find((p) => p.id === tp.property_id),
        is_required: tp.is_required ?? false, default_value: tp.default_value ?? null}))
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
  const runs = await pltsDb.analysisRunses.query().where({analysis_type: analysisType}).top(100000);
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
  const df = await pltsDb.plateWellValueses.query()
    .where({plate_id: id})
    .orderBy('property_id').orderBy('row').orderBy('col')
    .top(1000000).df();
  const plate = Plate.fromDbDataFrame(df);
  plate.id = id;
  events.next({on: 'after', eventType: 'read', objectType: TYPE.PLATE, object: plate});
  return plate;
}


type ValueRow = {plate_id?: string; property_id: string; row?: number; col?: number;
  value_num?: number; value_string?: string; value_bool?: boolean; value_jsonb?: string};

function rowValue(row: ValueRow, prop: Partial<PlateProperty>): any {
  return (row as any)[plateDbColumn[prop.type!] ?? plateDbJsonColumn];
}

function intersect(a: Set<string> | null, b: Set<string>): Set<string> {
  return a === null ? b : new Set([...a].filter((id) => b.has(id)));
}

/** Intersects the row sets matched by each condition; null = no conditions (everything passes). */
function matchingIds<T extends ValueRow>(rows: T[], conditions: PropertyCondition[],
  key: (r: T) => string): Set<string> | null {
  let ids: Set<string> | null = null;
  for (const condition of conditions) {
    const passed = new Set<string>(rows
      .filter((r) => r.property_id === condition.property.id && condition.matcher.match(rowValue(r, condition.property)))
      .map(key));
    ids = intersect(ids, passed);
  }
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

    let runs = await pltsDb.analysisRunses.query().where({analysis_type: analysisName}).top(100000);
    if (groups.length > 0 && propertyConditions.length === 0)
      runs = runs.filter((r) => (r.groups ?? []).some((g) => groups.includes(g)));

    let runIds = new Set(runs.map((r) => r.id));
    if (propertyConditions.length > 0) {
      const results = (await pltsDb.analysisResultses.query().top(1000000))
        .filter((r) => runIds.has(r.analysis_run_id) &&
          (groups.length === 0 || (r.group_combination ?? []).some((g) => groups.includes(g))));
      for (const condition of propertyConditions) {
        const prop = allProperties.find((p) => p.name === condition.property!.name);
        if (!prop) continue;
        const passed = new Set<string>(results
          .filter((r) => r.property_id === prop.id && condition.matcher!.match(rowValue(r, prop)))
          .map((r) => r.analysis_run_id));
        runIds = intersect(runIds, passed);
      }
    }

    const typePlateIds = new Set(runs.filter((r) => runIds.has(r.id)).map((r) => r.plate_id));
    plateIds = intersect(plateIds, typePlateIds);
  }
  return plateIds;
}

function propertyNamesById(): Map<string, string> {
  return new Map(allProperties.map((p) => [p.id, p.name]));
}

function pivotedDf(rows: {[key: string]: any}[], options?: {explode?: boolean}): DG.DataFrame {
  const df = rows.length > 0 ? DG.DataFrame.fromObjects(rows)! : DG.DataFrame.create();
  const propsCol = df.col('properties');
  if ((options?.explode ?? true) && propsCol && df.rowCount > 0) {
    try {
      DG.Utils.jsonToColumns(propsCol);
      propsCol.name = '~properties';
    } catch (e) {
      console.error('Failed to parse properties JSON:', e);
    }
  }
  return df;
}

export async function queryWells(query: PlateQuery): Promise<DG.DataFrame> {
  await initPlates();
  const wellRows = await pltsDb.plateWellValueses.query().top(1000000);
  const detailRows = query.plateMatchers.length > 0 ? await pltsDb.plateDetailses.query().top(1000000) : [];
  const plateIds = matchingIds(detailRows, query.plateMatchers, (r) => r.plate_id);
  const wellKey = (r: ValueRow) => `${r.plate_id}|${r.row}|${r.col}`;
  const wellIds = matchingIds(wellRows, query.wellMatchers, wellKey);

  const barcodes = new Map((await pltsDb.plateses.query().top(1000000)).map((p) => [p.id, p.barcode ?? '']));
  const propNames = propertyNamesById();

  const wells = new Map<string, {plate_id: string; row: number; col: number; props: {[key: string]: any}}>();
  for (const r of wellRows) {
    if (plateIds !== null && !plateIds.has(r.plate_id))
      continue;
    const key = wellKey(r);
    if (wellIds !== null && !wellIds.has(key))
      continue;
    const well = wells.get(key) ?? {plate_id: r.plate_id, row: r.row, col: r.col, props: {}};
    const name = propNames.get(r.property_id);
    if (name)
      well.props[name] = r.value_string ?? r.value_num ?? r.value_bool ?? null;
    wells.set(key, well);
  }

  const sorted = [...wells.values()]
    .sort((a, b) => a.plate_id.localeCompare(b.plate_id) || a.row - b.row || a.col - b.col);
  return pivotedDf(sorted.map((w) => ({plate_id: w.plate_id, barcode: barcodes.get(w.plate_id) ?? '',
    row: w.row, col: w.col, properties: JSON.stringify(w.props)})));
}

export async function queryPlates(query: PlateQuery): Promise<DG.DataFrame> {
  await initPlates();
  const plateRows = await pltsDb.plateses.query().orderBy('created_on').top(1000000);
  const detailRows = await pltsDb.plateDetailses.query().top(1000000);

  const idFilters = [
    matchingIds(detailRows, query.plateMatchers, (r) => r.plate_id),
    query.wellMatchers.length > 0 ?
      matchingIds(await pltsDb.plateWellValueses.query().top(1000000), query.wellMatchers, (r) => r.plate_id) : null,
    await matchingAnalysisPlateIds(query.analysisMatchers),
  ];
  let filtered = plateRows;
  for (const ids of idFilters) {
    if (ids !== null)
      filtered = filtered.filter((p) => ids.has(p.id));
  }

  const propNames = propertyNamesById();
  const detailsByPlate = new Map<string, {[key: string]: any}>();
  for (const d of detailRows) {
    const name = propNames.get(d.property_id);
    if (!name) continue;
    const props = detailsByPlate.get(d.plate_id) ?? {};
    props[name] = d.value_string ?? d.value_num ?? d.value_bool ?? null;
    detailsByPlate.set(d.plate_id, props);
  }

  return pivotedDf(filtered.map((p) => ({plate_id: p.id, barcode: p.barcode ?? '',
    description: p.description ?? '', properties: JSON.stringify(detailsByPlate.get(p.id) ?? {})})));
}


export async function createProperty(prop: Partial<PlateProperty>): Promise<PlateProperty> {
  const report = (await pltsDb.propertieses.insert({
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
    detailOps.push({op: 'insert', table: 'plate_details',
      values: {plate_id: '$plate', property_id: property.id,
        ...({[plateDbColumn[property.type]]: plate.details[layer]} as Partial<PlateDetailsInsert>)}});
  }

  const results = await pltsDb.transaction([
    {op: 'insert', table: 'plates', ref: 'plate',
      values: {plate_type_id: plate.plateTypeId, barcode: plate.barcode, template_id: plate.plateTemplateId}},
    ...detailOps,
  ]);
  plate.id = results[0].id;

  const wellRows: PlateWellValuesInsert[] = [];
  for (const layer of plate.getLayerNames()) {
    const property = allProperties.find((p) => p.scope === 'well' && p.name.toLowerCase() == layer.toLowerCase());
    if (!property) {
      console.warn(`Property '${layer}' not found in cache. Skipping save for this well-level property.`);
      continue;
    }
    const valueColumn = plateDbColumn[property.type];
    for (const pw of plate.wells) {
      wellRows.push({plate_id: plate.id, row: pw.row, col: pw.col, property_id: property.id,
        ...({[valueColumn]: pw[layer] ?? undefined} as Partial<PlateWellValuesInsert>)});
    }
  }
  if (wellRows.length > 0) {
    try {
      await pltsDb.plateWellValueses.batch(wellRows);
    } catch (e) {
      await pltsDb.plateses.delete(plate.id);
      plate.id = undefined;
      throw e;
    }
  }

  events.next({on: 'after', eventType: 'created', objectType: TYPE.PLATE, object: plate});
}


export async function createPlateTemplate(template: PlateTemplateInput): Promise<PlateTemplate> {
  const templateId = (await pltsDb.templateses.insert({
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
    await pltsDb.templatePropertieses.upsert({template_id: templateId, property_id: newProp.id,
      is_required: isRequired, default_value: property.default_value ?? undefined});

    createdPlateProperties.push(newProp);
    if (isRequired)
      required_props.push([newProp.id, newProp.name]);
  }

  for (const property of template.wellProperties ?? []) {
    property.scope = 'well';
    const newProp = await getOrCreateGlobalProperty(property);

    const isRequired = property.is_required ?? false;
    await pltsDb.templatePropertieses.upsert({template_id: templateId, property_id: newProp.id,
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
/** Deletes a plate template; its template_property links go with it (cascade). */
export async function deletePlateTemplate(template: PlateTemplate): Promise<void> {
  await pltsDb.templateses.delete(template.id);
  events.next({on: 'after', eventType: 'deleted', objectType: TYPE.TEMPLATE, object: template});
  await initPlates(true); // Refresh the cache
}

export async function createAnalysisRun(plateId: string, analysisType: string, groups: string[]): Promise<string> {
  if (!plateId)
    throw new Error('Cannot create analysis run: plateId is missing.');

  const [report] = await pltsDb.analysisRunses.insert({plate_id: plateId, analysis_type: analysisType, groups: groups});
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

  const values: AnalysisRunParametersInsert = {analysis_run_id: params.runId, property_id: prop.id};
  if (dbColumn === 'value_string') values.value_string = String(params.value);
  else if (dbColumn === 'value_num') values.value_num = params.value;
  else if (dbColumn === 'value_bool') values.value_bool = params.value;

  await pltsDb.analysisRunParameterses.insert(values);
}

export async function saveAnalysisResult(params: {
    runId: string,
    propertyId: string,
    propertyName: string,
    propertyType: string,
    value: any,
    groupCombination: string[]
}): Promise<void> {
  const values: AnalysisResultsInsert = {analysis_run_id: params.runId, property_id: params.propertyId,
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

  await pltsDb.analysisResultses.insert(values);
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

  const runs = await pltsDb.analysisRunses.query()
    .where({analysis_type: query.analysisName}).orderBy('created_on').top(100000);
  const runIds = new Set(runs.map((r) => r.id));
  const results = (await pltsDb.analysisResultses.query().top(1000000))
    .filter((r) => runIds.has(r.analysis_run_id));

  const propNames = propertyNamesById();
  type Combo = {runId: string; combo: string[]; props: {[key: string]: any}; rows: ValueRow[]};
  const combos = new Map<string, Combo>();
  for (const r of results) {
    if (selectedGroups.length > 0 && !(r.group_combination ?? []).some((g) => selectedGroups.includes(g)))
      continue;
    const key = `${r.analysis_run_id}${SEP}${(r.group_combination ?? []).join(SEP)}`;
    const combo = combos.get(key) ?? {runId: r.analysis_run_id, combo: r.group_combination ?? [], props: {}, rows: []};
    const name = propNames.get(r.property_id);
    if (name)
      combo.props[name] = r.value_jsonb ?? r.value_string ?? r.value_num ?? r.value_bool ?? null;
    combo.rows.push({plate_id: '', property_id: r.property_id, value_num: r.value_num,
      value_string: r.value_string, value_bool: r.value_bool, value_jsonb: r.value_jsonb});
    combos.set(key, combo);
  }

  let filtered = [...combos.values()];
  for (const condition of query.propertyMatchers) {
    const prop = allProperties.find((p) => p.name === condition.property.name);
    if (!prop) continue;
    filtered = filtered.filter((c) =>
      c.rows.some((r) => r.property_id === prop.id && condition.matcher.match(rowValue(r, prop))));
  }

  const runOrder = new Map(runs.map((r, i) => [r.id, i]));
  const plateByRun = new Map(runs.map((r) => [r.id, r.plate_id]));
  const barcodes = new Map((await pltsDb.plateses.query().top(1000000)).map((p) => [p.id, p.barcode ?? '']));
  filtered.sort((a, b) => (runOrder.get(a.runId)! - runOrder.get(b.runId)!) ||
    a.combo.join(', ').localeCompare(b.combo.join(', ')));

  return pivotedDf(filtered.map((c) => ({run_id: c.runId, plate_id: plateByRun.get(c.runId) ?? '',
    barcode: barcodes.get(plateByRun.get(c.runId) ?? '') ?? '',
    group_combination: c.combo.join(', '), properties: JSON.stringify(c.props)})), {explode: false});
}

/** Loads distinct string values of plate- and well-level properties for search-form suggestions. */
export async function loadUniquePropertyValues(): Promise<void> {
  await initPlates();
  const stringPropNames = new Map(allProperties
    .filter((p) => p.type === DG.COLUMN_TYPE.STRING).map((p) => [p.id, p.name]));
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
  plateUniquePropertyValues = toDf(await pltsDb.plateDetailses.query().top(100000));
  wellUniquePropertyValues = toDf(await pltsDb.plateWellValueses.query().top(100000));
}
