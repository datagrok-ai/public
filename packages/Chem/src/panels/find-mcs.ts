import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {findMCS} from '../scripts-api';

/**
 * Finds most common substructure in molecules column and returns a variable.
 * @export
 * @param {DG.Column} col column with smiles.
 */
export async function getMcs(col: DG.Column): Promise<string> {
  if (col.length >= 100000) {
    grok.shell.error('Number of sructures exceeeds 100000');
  }

  const pi = DG.TaskBarProgressIndicator.create('Estimating MCS');
  const mcs: string = await findMCS(col.name, DG.DataFrame.fromColumns([col]));
  pi.close();
  return mcs;
}
