import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {getSplitterForColumn} from './macromolecule/utils';


export function splitAlignedSequences(sequenceColumn: DG.Column<string>): DG.DataFrame {
  const splitter = getSplitterForColumn(sequenceColumn);
  const getCol = (index: number): DG.Column<string> | null => columnList[index] ?? null;
  const createCol = (index: number): DG.Column<string> => {
    const positionCol = resultDf.columns.addNewString((index + 1).toString());
    columnList.push(positionCol);
    return positionCol;
  };

  const columnList: DG.Column<string>[] = [];
  const rowCount = sequenceColumn.length;
  const resultDf = DG.DataFrame.create(rowCount);

  for (let rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
    const sequence = sequenceColumn.get(rowIndex);
    if (sequence == null)
      continue;

    const currentMonomerList = splitter(sequence);
    for (const [monomer, positionIndex] of wu.enumerate(currentMonomerList)) {
      const col = getCol(positionIndex) || createCol(positionIndex);
      col.set(rowIndex, monomer || '-', false);
    }
  }

  return resultDf;
}
