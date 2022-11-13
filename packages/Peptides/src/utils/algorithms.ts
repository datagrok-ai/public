import * as DG from 'datagrok-api/dg';

import * as C from './constants';
import * as type from './types';
import {getTypedArrayConstructor} from './misc';

//TODO: move out
export function findMutations(activityCol: DG.Column<number>, monomerColumns: DG.Column<string>[],
  settings: type.PeptidesSettings = {}): type.SubstitutionsInfo {
  const nCols = monomerColumns.length;
  if (nCols == 0)
    throw new Error(`PepAlgorithmError: Couldn't find any column of semType '${C.SEM_TYPES.MONOMER}'`);

  const substitutionsInfo: type.SubstitutionsInfo = new Map();
  const nRows = activityCol.length;
  for (let seq1Idx = 0; seq1Idx < nRows - 1; seq1Idx++) {
    for (let seq2Idx = seq1Idx + 1; seq2Idx < nRows; seq2Idx++) {
      let substCounter = 0;
      const activityValSeq1 = activityCol.get(seq1Idx)!;
      const activityValSeq2 = activityCol.get(seq2Idx)!;
      const delta = activityValSeq1 - activityValSeq2;
      if (Math.abs(delta) < (settings.minActivityDelta ?? 0))
        continue;

      let substCounterFlag = false;
      const tempData: { pos: string, seq1monomer: string, seq2monomer: string, seq1Idx: number, seq2Idx: number }[] =
        [];
      for (const currentPosCol of monomerColumns) {
        const seq1monomer = currentPosCol.get(seq1Idx)!;
        const seq2monomer = currentPosCol.get(seq2Idx)!;
        if (seq1monomer == seq2monomer)
          continue;

        substCounter++;
        substCounterFlag = substCounter > (settings.maxMutations ?? 1);
        if (substCounterFlag)
          break;

        tempData.push({
          pos: currentPosCol.name,
          seq1monomer: seq1monomer,
          seq2monomer: seq2monomer,
          seq1Idx: seq1Idx,
          seq2Idx: seq2Idx,
        });
      }

      if (substCounterFlag || substCounter == 0)
        continue;

      for (const tempDataElement of tempData) {
        const position = tempDataElement.pos;

        //Working with seq1monomer
        const seq1monomer = tempDataElement.seq1monomer;
        if (!substitutionsInfo.has(seq1monomer))
          substitutionsInfo.set(seq1monomer, new Map());

        let positionsMap = substitutionsInfo.get(seq1monomer)!;
        if (!positionsMap.has(position))
          positionsMap.set(position, new Map());

        let indexes = positionsMap.get(position)!;

        !indexes.has(seq1Idx) ? indexes.set(seq1Idx, [seq2Idx]) : (indexes.get(seq1Idx)! as number[]).push(seq2Idx);

        //Working with seq2monomer
        const seq2monomer = tempDataElement.seq2monomer;
        if (!substitutionsInfo.has(seq2monomer))
          substitutionsInfo.set(seq2monomer, new Map());

        positionsMap = substitutionsInfo.get(seq2monomer)!;
        if (!positionsMap.has(position))
          positionsMap.set(position, new Map());

        indexes = positionsMap.get(position)!;
        !indexes.has(seq2Idx) ? indexes.set(seq2Idx, [seq1Idx]) : (indexes.get(seq2Idx)! as number[]).push(seq1Idx);
      }
    }
  }

  const TypedArray = getTypedArrayConstructor(nRows);
  for (const positionMap of substitutionsInfo.values()) {
    for (const indexMap of positionMap.values()) {
      for (const [index, indexArray] of indexMap.entries())
        indexMap.set(index, new TypedArray(indexArray));
    }
  }

  return substitutionsInfo;
}
