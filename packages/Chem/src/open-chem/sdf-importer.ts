/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';

export function _importSdf(bytes: Uint8Array): DG.DataFrame[] {
  const str = new TextDecoder().decode(bytes);
  return _importSdfString(str);
}

export function _importSdfString(str: string): DG.DataFrame[] {
  // @ts-ignore
  let parser = new OCL.SDFileParser(str, null);
  const fieldNames = parser.getFieldNames(100);
  parser = new OCL.SDFileParser(str, fieldNames);
  let rows = 0;
  const data: any = {'molecule': []};
  const addName = !fieldNames.includes('Name');

  while (parser.next()) {
    for (const field of fieldNames) {
      if (!data[field])
        data[field] = [];
      data[field][rows] = parser.getField(field);
    }
    const molfile = parser.getNextMolFile();
    data['molecule'][rows] = molfile;
    if (addName)
      (data['Name'] ??= [])[rows] = molfile.split('\n', 1)[0].trim();
    rows++;
  }

  const df = DG.DataFrame.create(rows);
  for (const field of ['molecule'].concat(addName ? ['Name'] : [], fieldNames))
    df.columns.add(DG.Column.fromStrings(field, data[field]));
  df.col('molecule')!.semType = DG.SEMTYPE.MOLECULE;
  df.col('molecule')!.meta.units = 'molblock'; // DG.UNITS.Molecule.MOLBLOCK;

  return [df];
}
