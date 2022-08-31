import * as DG from 'datagrok-api/dg';

import {WebLogo} from '../viewers/web-logo';

export function splitAlignedSequences(peptideColumn: DG.Column<string>): DG.DataFrame {
  const splitter = WebLogo.getSplitterForColumn(peptideColumn);
  const colLen = peptideColumn.length;
  const resultDf = DG.DataFrame.create(colLen);
  let monomerList = splitter(peptideColumn.get(0)!);
  const columnList: DG.Column<string>[] = [];

  // create columns and fill the first row for faster values filling in the next loop
  for (let i = 0; i < monomerList.length; i++) {
    const col = resultDf.columns.addNewString((i + 1).toString());
    col.set(0, monomerList[i] || '-', false);
    columnList.push(col);
  }

  for (let rowIndex = 1; rowIndex < colLen; rowIndex++) {
    monomerList = splitter(peptideColumn.get(rowIndex)!);
    monomerList.forEach((monomer, colIndex) => columnList[colIndex].set(rowIndex, monomer || '-', false));
  }

  return resultDf;
}
  