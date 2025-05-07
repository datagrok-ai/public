import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Plate} from "../plate/plate";


type PlateType = {
  id: number;
  name: string;
  rows: number;
  cols: number;
  maxVolume?: number;
}

type PlateProperty = {
  id: number;
  name: string;
  unit?: string;
  value_type: string; //DG.COLUMN_TYPE;
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
  _initialized = true;
}

function getValueType(x: any): string {
  if (typeof x === 'string')
    return DG.TYPE.STRING;
  if (typeof x === 'number')
    return DG.TYPE.FLOAT;
  if (typeof x === 'boolean')
    return DG.TYPE.BOOL;
  console.log(x);
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
  grok.shell.info('100 plates saved');
}
