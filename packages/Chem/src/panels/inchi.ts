import * as DG from 'datagrok-api/dg';
import {inchisByEngine} from '../engine/worker-jobs';

/** Builds a derived column out of a molecular column `col` (the parallel workers of the Chemical Engine in effect).
 * Handles progress indication. Malformed molecules result in an empty string. */
async function getDerived(col: DG.Column, description: string,
  extract: (molecules: string[]) => Promise<string[]>, colName: string): Promise<DG.Column> {
  const pi = DG.TaskBarProgressIndicator.create(description);
  try {
    return DG.Column.fromList('string', colName, await extract(col.toList()));
  } finally {
    pi.close();
  }
}

/** Adds InchI identifiers for the specified molecular column (by the Chemical Engine in effect). */
export async function getInchisImpl(col: DG.Column): Promise<DG.Column> {
  return getDerived(col, 'Getting Inchi', (molecules) => inchisByEngine(molecules, false), 'inchi');
}

/** Adds InchI keys identifiers for the specified molecular column (by the Chemical Engine in effect). */
export async function getInchiKeysImpl(col: DG.Column): Promise<DG.Column> {
  return getDerived(col, 'Getting Inchi', (molecules) => inchisByEngine(molecules, true), 'inchi_key');
}
