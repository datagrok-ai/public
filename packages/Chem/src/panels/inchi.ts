import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {_rdKitModule, getRdKitModule} from '../utils/chem-common-rdkit';
import {StringUtils} from "@datagrok-libraries/utils/src/string-utils";
import {RDMol} from "../rdkit-api";

function addDerived(col: DG.Column, description: string, extract: (mol: RDMol) => string) {
  const pi = DG.TaskBarProgressIndicator.create(description);
  const result = new Array(col.length);
  for (let i = 0; i < result.length; i++) {
    const mol = _rdKitModule.get_mol(col.get(i));
    result[i] = extract(mol);
    mol.delete();
  }

  const name = col.dataFrame.columns.getUnusedName('inchi');
  col.dataFrame.columns.add(DG.Column.fromList('string', name, result));
  pi.close();
}

/** Adds InchI identifiers for the specified molecular column. */
export function addInchis(col: DG.Column): void {
  addDerived(col, 'Getting Inchi', (m) => m.get_inchi());
}

/** Adds InchI keys identifiers for the specified molecular column. */
export function addInchiKeys(col: DG.Column): void {
  addDerived(col, 'Getting Inchi', (m) => _rdKitModule.get_inchikey_for_inchi(m.get_inchi()));
}
