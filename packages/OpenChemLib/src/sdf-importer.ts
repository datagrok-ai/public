/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';

export function _importSdf(bytes: Uint8Array): DG.DataFrame[] {
  let str = new TextDecoder().decode(bytes);
  // @ts-ignore
  let parser = new OCL.SDFileParser(str, null);
  let fieldNames = parser.getFieldNames(100);
  parser = new OCL.SDFileParser(str, fieldNames);
  let rows = 0;
  let data: any = {'molecule': []};

  while (parser.next()) {
    for (let field of fieldNames) {
      if (!data[field])
        data[field] = [];
      data[field][rows] = parser.getField(field);
      data['molecule'][rows] = parser.getNextMolFile();
    }
    rows++;
  }

  let df = DG.DataFrame.create(rows);
  for (let field of ['molecule'].concat(fieldNames))
    df.columns.add(DG.Column.fromStrings(field, data[field]));
  df.col('molecule')!.semType = DG.SEMTYPE.MOLECULE;
  df.col('molecule')!.tags[DG.TAGS.UNITS] = 'molblock';  DG.UNITS.Molecule.MOLBLOCK;

  return [df];
}