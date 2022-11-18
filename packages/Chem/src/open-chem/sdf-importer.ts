/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';

export function _importSdf(bytes: Uint8Array): DG.DataFrame[] {
  const str = new TextDecoder().decode(bytes);
  // @ts-ignore
  let parser = new OCL.SDFileParser(str, null);
  const fieldNames = parser.getFieldNames(100);
  parser = new OCL.SDFileParser(str, fieldNames);
  let rows = 0;
  const data: any = {'molecule': []};

  while (parser.next()) {
    for (const field of fieldNames) {
      if (!data[field])
        data[field] = [];
      data[field][rows] = parser.getField(field);
    }
    data['molecule'][rows] = parser.getNextMolFile();
    rows++;
  }

  const df = DG.DataFrame.create(rows);
  for (const field of ['molecule'].concat(fieldNames))
    df.columns.add(DG.Column.fromStrings(field, data[field]));
  df.col('molecule')!.semType = DG.SEMTYPE.MOLECULE;
  df.col('molecule')!.tags[DG.TAGS.UNITS] = 'molblock'; // DG.UNITS.Molecule.MOLBLOCK;

  return [df];
}
