/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full';

export function _importSmi(bytes: Uint8Array): DG.DataFrame[] {
  let str = new TextDecoder().decode(bytes);
  str = str.replace('#SMILES', 'SMILES');
  if(!str.includes('SMILES'))
    str = 'SMILES\n' + str;

  let df = DG.DataFrame.fromCsv(str);

  return [df];
}
