import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {ISeqHelper} from './seq-helper';


export function splitAlignedSequences(sequenceColumn: DG.Column<string>, seqHelper: ISeqHelper): DG.DataFrame {
  const getCol = (index: number): DG.Column<string> | null => columnList[index] ?? null;
  const createCol = (index: number): DG.Column<string> => {
    const positionCol = resultDf.columns.addNewString((index + 1).toString());
    columnList.push(positionCol);
    return positionCol;
  };

  const columnList: DG.Column<string>[] = [];
  const rowCount = sequenceColumn.length;
  const resultDf = DG.DataFrame.create(rowCount);

  const uh = seqHelper.getSeqHandler(sequenceColumn);
  for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
    const currentMonomerList = uh.getSplitted(rowIdx);
    for (let posIdx: number = 0; posIdx < currentMonomerList.length; ++posIdx) {
      const cm: string = currentMonomerList.getCanonical(posIdx);
      const col = getCol(posIdx) || createCol(posIdx);
      col.set(rowIdx, cm, false);
    }
  }

  return resultDf;
}
