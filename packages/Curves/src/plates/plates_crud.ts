import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Plate} from "../plate/plate";
import {Utils} from "datagrok-api/dg";
import {NumericMatcher, Matcher} from "./numeric_matcher";


export type PlateType = {
  id: number;
  name: string;
  rows: number;
  cols: number;
  maxVolume?: number;
}

export type PlateProperty = {
  id: number;
  name: string;
  unit?: string;
  value_type: string; //DG.COLUMN_TYPE;
}

export type PropertyCondition = {
  property: PlateProperty;
  matcher: Matcher;
}

export type PlateQuery = {
  plateMatchers: PropertyCondition[];
  wellMatchers: PropertyCondition[];
}


export let wellProperties: PlateProperty[] = [
  {id: 1000, name: 'Volume', value_type: DG.COLUMN_TYPE.FLOAT},
  {id: 1001, name: 'Concentration', value_type: DG.COLUMN_TYPE.FLOAT},
  {id: 1002, name: 'Sample', value_type: DG.COLUMN_TYPE.STRING},
  {id: 1003, name: 'Well Role', value_type: DG.COLUMN_TYPE.STRING}
]

export let plateProperties: PlateProperty[] = [
]

export let plateTypes: PlateType[] = [
  {id: 1, name: 'Generic 96 wells', rows: 8, cols: 12},
  {id: 2, name: 'Generic 384 wells', rows: 16, cols: 24},
  {id: 3, name: 'Generic 1536 wells', rows: 32, cols: 48},
]

export let plateUniquePropertyValues: DG.DataFrame = DG.DataFrame.create();
export let wellUniquePropertyValues: DG.DataFrame = DG.DataFrame.create();

export const plateDbColumn: {[key: string]: string} = {
  [DG.COLUMN_TYPE.FLOAT]: 'value_num',
  [DG.COLUMN_TYPE.INT]: 'value_num',
  [DG.COLUMN_TYPE.BOOL]: 'value_bool',
  [DG.COLUMN_TYPE.STRING]: 'value_string',
}

let _initialized = false;
export async function initPlates() {
  if (_initialized)
    return;
  plateProperties = (await grok.functions.call('Curves:getPlateLevelProperties') as DG.DataFrame).toJson();
  wellProperties = (await grok.functions.call('Curves:getWellLevelProperties') as DG.DataFrame).toJson();
  plateTypes = (await grok.functions.call('Curves:getPlateTypes') as DG.DataFrame).toJson();
  plateUniquePropertyValues = await grok.functions.call('Curves:getUniquePlatePropertyValues');
  wellUniquePropertyValues = await grok.functions.call('Curves:getUniqueWellPropertyValues');

  _initialized = true;
}

export function getPlateProperty(name: string): PlateProperty {
  return plateProperties.find(p => p.name === name)!;
}

export function getUniquePropertyValues(prop: PlateProperty, df: DG.DataFrame): string[] {
  const nameCol = df.col('name')!;
  const valueCol = df.col('value_string')!;

  return df.rows
    .where(i => nameCol.get(i) == prop.name)
    .map(i => valueCol.get(i))
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
  throw 'Not supported type';
}


function sqlStr(s?: string | number) {
  if (!s)
    return 'null';
  if (typeof s == 'string')
    return `'${s}'`;
  else
    return `${s}`;
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
    const dbColumn = plateDbColumn[condition.property.value_type];
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
    const dbColumn = plateDbColumn[condition.property.value_type];
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
    existsClauses.length > 0 ? `WHERE ${existsClauses.join(" AND ")}` : "";

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


export async function queryPlates(query: PlateQuery): Promise<DG.DataFrame> {
  const df = await grok.data.db.query('Admin:Plates', getPlateSearchSql(query));
  Utils.jsonToColumns(df.col('properties')!);
  df.col('properties')!.name = '~properties';
  return df;
}


/** Saves the plate to the database. */
export async function savePlate(plate: Plate){
  await initPlates();

  const plateSql =
    `insert into plates.plates(plate_type_id, barcode)
     values(${plate.plateTypeId}, ${sqlStr(plate.barcode)})
     returning id`;
  plate.id = (await grok.data.db.query('Admin:Plates', plateSql)).get('id', 0);   

  // register new plate level properties
  for (const layer of Object.keys(plate.details))
    if (!plateProperties.find(p => p.name.toLowerCase() == layer.toLowerCase())) {
      const valueType = getValueType(plate.details[layer]);
      const propertyId = await grok.functions.call('Curves:createProperty', {
        propertyName: layer,
        valueType: valueType
      });
      plateProperties.push({id: propertyId, name: layer, value_type: valueType});
      grok.shell.info('Plate layer created: ' + layer);
    }

  // register new well level properties
  for (const layer of plate.getLayerNames())
    if (!wellProperties.find(p => p.name.toLowerCase() == layer.toLowerCase())) {
      const col = plate.data.col(layer)
      const propertyId = await grok.functions.call('Curves:createProperty', {
        propertyName: layer,
        valueType: col?.type
      });
      wellProperties.push({id: propertyId, name: layer, value_type: col!.type});
      grok.shell.info('Well layer created: ' + layer);
    }
    
  await grok.data.db.query('Admin:Plates', getPlateInsertSql(plate));
  grok.shell.info('Plate saved');
}


function getPlateInsertSql(plate: Plate): string {
  let sql = 'insert into plates.plate_wells(plate_id, row, col) values '
  + plate.wells.map(pw => `  (${plate.id}, ${pw.row}, ${pw.col})`).toArray().join(',\n') + ';'

  // plate data
  for (const layer of Object.keys(plate.details)) {
    const property = plateProperties.find(p => p.name.toLowerCase() == layer.toLowerCase())!;
    const dbCol = plateDbColumn[property.value_type];

    sql += `\n insert into plates.plate_details(plate_id, property_id, ${dbCol}) values `
      + `(${plate.id}, ${property.id}, ${sqlStr(plate.details[layer])});`;
  }

  // well data
  for (const layer of plate.getLayerNames()) {
    const property = wellProperties.find(p => p.name.toLowerCase() == layer.toLowerCase())!;
    const dbCol = plateDbColumn[property.value_type];
    sql += `\n\n insert into plates.plate_well_values(plate_id, row, col, property_id, ${dbCol}) values\n`
      + plate.wells.map(pw => `  (${plate.id}, ${pw.row}, ${pw.col}, ${property.id}, ${sqlStr(pw[layer])})`).toArray().join(',\n') + ';'
  }

  return sql;
}


export async function __createDummyPlateData() {
  await initPlates();

  for (let i = 0; i < 50; i++) {
    const plate = Plate.demo();
    plate.details = {
      'Project': DG.Utils.random([
        'Allosteric Modulators of mGluR5',
        'Agonists for Orphan GPCR GPR139',
        'Glutaminase Inhibitors for TNBC']),
      'Stage': DG.Utils.random(['Lead generation', 'Lead optimization']),
      'Chemist': DG.Utils.random(['John Marlowski', 'Mary Hopton']),
      'Passed QC': DG.Utils.random([true, false]),
      'Z-score': Math.random() * 3,
    }
    await savePlate(plate);
  }

  for (let i = 0; i < 50; i++) {
    const plate = Plate.demo();
    plate.details = {
      'Project': DG.Utils.random([
        'NaV1.7 Blockers for Pain Relief',
        'TLR4 Antagonists for Sepsis',]),
      'Stage': DG.Utils.random(['Stage 1', 'Stage 2', 'Stage 3']),
      'Chemist': DG.Utils.random(['John Marlowski', 'Andrew Smith']),
      'Biologist': DG.Utils.random(['Anna Fei', 'Joan Dvorak']),
      'Cells': Math.ceil(Math.random() * 100)
    }
    await savePlate(plate);
  }

  grok.shell.info('100 plates saved');
}
