import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {_rdKitModule} from '../utils/chem-common-rdkit';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

/** Adds a derived column, given a source column `col` and extraction function `extract`.
 * Handles progress indication, and molecule disposal. */
function addDerived(table: DG.DataFrame, col: DG.Column, description: string, extract: (mol: RDMol) => string, colName: string): DG.Column {
  const pi = DG.TaskBarProgressIndicator.create(description);
  const result = new Array(col.length);
  for (let i = 0; i < result.length; i++) {
    try{
      const mol = _rdKitModule.get_mol(col.get(i));
      result[i] = extract(mol);
      mol.delete();
    } catch (e: any) {
      result[i] = '';
    }
  }
  const name = table.columns.getUnusedName(colName);
  const resultColumn = table.columns.add(DG.Column.fromList('string', name, result));
  pi.close();
  return resultColumn;
}

/** Adds InchI identifiers for the specified molecular column. */
export function addInchis(table: DG.DataFrame, col: DG.Column): DG.Column {
  return addDerived(table, col, 'Getting Inchi', (m) => m.get_inchi(), 'inchi');
}

/** Adds InchI keys identifiers for the specified molecular column. */
export function addInchiKeys(table: DG.DataFrame, col: DG.Column): DG.Column {
  return addDerived(table, col, 'Getting Inchi', (m) => _rdKitModule.get_inchikey_for_inchi(m.get_inchi()), 'inchi_key');
}
