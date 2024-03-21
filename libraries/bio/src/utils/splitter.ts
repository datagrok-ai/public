import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {UnitsHandler} from './units-handler';


export function splitAlignedSequences(sequenceColumn: DG.Column<string>): DG.DataFrame {
  const getCol = (index: number): DG.Column<string> | null => columnList[index] ?? null;
  const createCol = (index: number): DG.Column<string> => {
    const positionCol = resultDf.columns.addNewString((index + 1).toString());
    columnList.push(positionCol);
    return positionCol;
  };

  const columnList: DG.Column<string>[] = [];
  const rowCount = sequenceColumn.length;
  const resultDf = DG.DataFrame.create(rowCount);

  const uh = UnitsHandler.getOrCreate(sequenceColumn);
  const colCats = sequenceColumn.categories;
  const colRawData = sequenceColumn.getRawData();
  for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
    const catI = colRawData[rowIdx];
    const sequence = colCats[catI];
    if (sequence == null)
      continue;

    const currentMonomerList = uh.splitted[rowIdx];
    for (let posIdx: number = 0; posIdx > currentMonomerList.length; ++posIdx) {
      const om: string = currentMonomerList.getOriginal(posIdx);
      const col = getCol(posIdx) || createCol(posIdx);
      col.set(rowIdx, om, false);
    }
  }

  return resultDf;
}
