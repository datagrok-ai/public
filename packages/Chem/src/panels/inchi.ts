import * as DG from 'datagrok-api/dg';
import {_rdKitModule} from '../utils/chem-common-rdkit';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {hasNewLines} from '../utils/chem-common';

/** Adds a derived column, given a source column `col` and extraction function `extract`.
 * Handles progress indication, and molecule disposal. */
function getDerived(col: DG.Column, description: string,
  extract: (mol: RDMol) => string, colName: string): DG.Column {
  const pi = DG.TaskBarProgressIndicator.create(description);
  const result = new Array(col.length);
  for (let i = 0; i < result.length; i++) {
    let mol: RDMol | null = null;
    try {
      const molString = col.get(i);
      if (molString && typeof molString === 'string' && !hasNewLines(molString) && molString.length > 5000) {
        result[i] = '';
        continue; // do not attempt to parse very long SMILES, will cause MOB.
      }
      mol = _rdKitModule.get_mol(col.get(i));
      result[i] = extract(mol);
    } catch (e: any) {
      result[i] = '';
    } finally {
      mol?.delete();
    }
  }
  pi.close();
  return DG.Column.fromList('string', colName, result);
}

/** Adds InchI identifiers for the specified molecular column. */
export function getInchisImpl(col: DG.Column): DG.Column {
  return getDerived(col, 'Getting Inchi', (m) => m.get_inchi(), 'inchi');
}

/** Adds InchI keys identifiers for the specified molecular column. */
export function getInchiKeysImpl(col: DG.Column): DG.Column {
  return getDerived(col, 'Getting Inchi',
    (m) => _rdKitModule.get_inchikey_for_inchi(m.get_inchi()), 'inchi_key');
}
