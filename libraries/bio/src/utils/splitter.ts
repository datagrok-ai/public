import * as DG from 'datagrok-api/dg';

import {getSplitterForColumn} from './macromolecule';

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
    currentMonomerList.forEach((monomer, positionIndex) => {
      const col = getCol(positionIndex) || createCol(positionIndex);
      col.set(rowIndex, monomer || '-', false);
    });
  }

  return resultDf;
}
