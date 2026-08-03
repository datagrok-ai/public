import * as DG from 'datagrok-api/dg';
// type-only: rdkit-service imports package.ts, importing it as a value here would close an import cycle
import type {RdKitService} from '../rdkit-service/rdkit-service';
import {getRdKitService} from '../utils/chem-common-rdkit';

/** Builds a derived column out of a molecular column `col`, using the parallel RDKit service.
 * Handles progress indication. Malformed molecules result in an empty string. */
async function getDerived(col: DG.Column, description: string,
  extract: (service: RdKitService, molecules: string[]) => Promise<string[]>,
  colName: string): Promise<DG.Column> {
  const service = await getRdKitService();
  const pi = DG.TaskBarProgressIndicator.create(description);
  try {
    return DG.Column.fromList('string', colName, await extract(service, col.toList()));
  } finally {
    pi.close();
  }
}

/** Adds InchI identifiers for the specified molecular column. */
export async function getInchisImpl(col: DG.Column): Promise<DG.Column> {
  return getDerived(col, 'Getting Inchi', (service, molecules) => service.getInchis(molecules), 'inchi');
}

/** Adds InchI keys identifiers for the specified molecular column. */
export async function getInchiKeysImpl(col: DG.Column): Promise<DG.Column> {
  return getDerived(col, 'Getting Inchi', (service, molecules) => service.getInchiKeys(molecules), 'inchi_key');
}
