import {delay} from '@datagrok-libraries/utils/src/test';
import {checkInputColumnUI} from './check-input-column';
import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';
import * as C from './constants';
import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export async function splitToMonomersUI(table: DG.DataFrame, seqCol: DG.Column<string>): Promise<DG.DataFrame> {
  // Delay is required for initial function dialog to close before starting invalidating of molfiles.
  // Otherwise, dialog is freezing
  await delay(10);
  if (!checkInputColumnUI(seqCol, 'Sequence space')) return table;

  const tempDf = splitAlignedSequences(seqCol);
  const originalDf = seqCol.dataFrame;
  for (const tempCol of tempDf.columns) {
    // TODO: GROK-11212
    // tempCol.setTag(DG.TAGS.CELL_RENDERER, C.SEM_TYPES.MONOMER);
    tempCol.semType = C.SEM_TYPES.MONOMER;
    tempCol.setTag(bioTAGS.alphabet, seqCol.getTag(bioTAGS.alphabet));
  }
  // Create the new data frame to enable platform to setup cell renderers
  const newDf = originalDf.join(tempDf, [], [], undefined, undefined, DG.JOIN_TYPE.LEFT, false);
  grok.shell.addTableView(newDf);

  return newDf;
}
