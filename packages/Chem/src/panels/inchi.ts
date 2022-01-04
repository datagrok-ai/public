import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from '../chem_common_rdkit';

/**
 * Adds InchI identifiers to dataframe.
 * @export
 * @param {DG.Column} col column with smiles.
 */
export function addInchis(col: DG.Column): void {
  const pi = DG.TaskBarProgressIndicator.create('Getting Inchi');
  const rdKitModule = getRdKitModule();
  const inchis = new Array(col.length);
  for (let i = 0; i < inchis.length; i++) {
    inchis[i] = rdKitModule.get_mol(col.get(i)).get_inchi();
  }

  const name = getName('inchi', col.dataFrame.columns.names());

  col.dataFrame.columns.add(DG.Column.fromList('string', name, inchis));
  pi.close();
}

/**
 * Adds InchI Keys to dataframe.
 * @export
 * @param {DG.Column} col column with smiles.
 */
export function addInchiKeys(col: DG.Column): void {
  const pi = DG.TaskBarProgressIndicator.create('Getting Inchi Keys');
  const rdKitModule = getRdKitModule();
  const inchiKeys = new Array(col.length);
  for (let i = 0; i < inchiKeys.length; i++) {
    inchiKeys[i] = rdKitModule.get_inchikey_for_inchi(rdKitModule.get_mol(col.get(i)).get_inchi());
  }

  const name = getName('inchi_key', col.dataFrame.columns.names());

  col.dataFrame.columns.add(DG.Column.fromList('string', name, inchiKeys));
  pi.close();
}

function getName(initialName: string, existingNames: string[]) {
  if (!existingNames.includes(initialName)) {
    return initialName;
  } else {
    let counter: number = 1;
    let newName: string = (' ' + initialName + '_' + counter).slice(1);
    while (existingNames.includes(newName)) {
      counter++;
      newName = (' ' + initialName + '_' + counter).slice(1);
    }

    return newName;
  }
}
