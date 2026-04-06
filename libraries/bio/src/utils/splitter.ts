import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {ISeqHelper} from './seq-helper';
import {INotationProvider} from './macromolecule/types';
import {SeqTemps} from './macromolecule/seq-handler';
import {MONOMER_CANONICALIZER_FUNC_TAG} from './macromolecule/consts';


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
      const om: string = currentMonomerList.getOriginal(posIdx);
      const col = getCol(posIdx) || createCol(posIdx);
      col.set(rowIdx, om, false);
    }
  }

  // If the notation provider specifies a monomer canonicalizer, tag all monomer columns
  const notationProvider: INotationProvider | null =
    sequenceColumn.temp[SeqTemps.notationProvider] ?? null;
  const canonFuncName = notationProvider?.monomerCanonicalizerFuncName;
  if (canonFuncName)
    for (const col of columnList)
      col.setTag(MONOMER_CANONICALIZER_FUNC_TAG, canonFuncName);

  return resultDf;
}
