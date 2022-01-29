import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {findMCS} from '../scripts-api';
import {StringUtils} from "@datagrok-libraries/utils/src/string-utils";

/**
 * Finds most common substructure in molecules column and adds to dataframe.
 * @export
 * @param {DG.Column} col column with smiles.
 */
export async function addMcs(col: DG.Column): Promise<void> {
  if (col.length >= 100000) {
    grok.shell.error('Number of sructures exceeeds 100000');
    return;
  }

  const pi = DG.TaskBarProgressIndicator.create('Estimating MCS');
  const mcs: string = await findMCS(col.name, DG.DataFrame.fromColumns([col]));
  const name = col.dataFrame.columns.getUnusedName('MCS');
  const mcsCol = DG.Column.fromList('string', name, new Array(col.length).fill(mcs));
  mcsCol.semType = 'Molecule';
  mcsCol.setTag('cell.renderer', 'Molecule');
  col.dataFrame.columns.add(mcsCol);

  pi.close();
}
