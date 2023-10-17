import * as C from './constants';
import * as type from './types';
import {getTypedArrayConstructor} from './misc';
import {ParallelMutationCliffs} from './parallel-mutation-cliffs';

type MutationCliffInfo = {pos: string, seq1monomer: string, seq2monomer: string, seq1Idx: number, seq2Idx: number};

export async function findMutations(activityArray: type.RawData, monomerInfoArray: type.RawColumn[],
  settings: type.PeptidesSettings = {},
  targetOptions: {targetCol?: type.RawColumn | null, currentTarget?: string | null} = {},
): Promise<type.MutationCliffs> {
  const nCols = monomerInfoArray.length;
  if (nCols === 0)
    throw new Error(`PepAlgorithmError: Couldn't find any column of semType '${C.SEM_TYPES.MONOMER}'`);

  settings.minActivityDelta ??= 0;
  settings.maxMutations ??= 1;
  //const currentTargetIdx = targetOptions.targetCol?.cat!.indexOf(targetOptions.currentTarget!) ?? -1;

  //const substitutionsInfo: type.MutationCliffs = new Map();
  //const nRows = activityArray.length;

  const substitutionsInfo =
    await new ParallelMutationCliffs().calc(activityArray, monomerInfoArray, settings, targetOptions);
  return substitutionsInfo;
  /*
  for (let seq1Idx = 0; seq1Idx < nRows - 1; seq1Idx++) {
    if (currentTargetIdx !== -1 && targetOptions.targetCol?.rawData[seq1Idx] !== currentTargetIdx)
      continue;

    for (let seq2Idx = seq1Idx + 1; seq2Idx < nRows; seq2Idx++) {
      if (currentTargetIdx !== -1 && targetOptions.targetCol?.rawData[seq2Idx] !== currentTargetIdx)
        continue;

      let substCounter = 0;
      const activityValSeq1 = activityArray[seq1Idx];
      const activityValSeq2 = activityArray[seq2Idx];
      const delta = activityValSeq1 - activityValSeq2;
      if (Math.abs(delta) < settings.minActivityDelta)
        continue;

      let substCounterFlag = false;
      const tempData: MutationCliffInfo[] = new Array(monomerInfoArray.length);
      let tempDataIdx = 0;
      for (const monomerInfo of monomerInfoArray) {
        const seq1categoryIdx = monomerInfo.rawData[seq1Idx];
        const seq2categoryIdx = monomerInfo.rawData[seq2Idx];
        if (seq1categoryIdx === seq2categoryIdx)
          continue;

        substCounter++;
        substCounterFlag = substCounter > settings.maxMutations;
        if (substCounterFlag)
          break;

        tempData[tempDataIdx++] ={
          pos: monomerInfo.name,
          seq1monomer: monomerInfo.cat![seq1categoryIdx],
          seq2monomer: monomerInfo.cat![seq2categoryIdx],
          seq1Idx: seq1Idx,
          seq2Idx: seq2Idx,
        };
      }

      if (substCounterFlag || substCounter === 0)
        continue;

      // Separate processing loop in case substCOunter is 0 or out of restricted range to
      // prevent unnecessary computations
      for (let i = 0; i < tempDataIdx; i++) {
        const tempDataElement = tempData[i];
        //Working with seq1monomer
        const seq1monomer = tempDataElement.seq1monomer;
        if (!substitutionsInfo.has(seq1monomer))
          substitutionsInfo.set(seq1monomer, new Map());

        const position = tempDataElement.pos;

        let positionsMap = substitutionsInfo.get(seq1monomer)!;
        if (!positionsMap.has(position))
          positionsMap.set(position, new Map());

        let indexes = positionsMap.get(position)!;
        if (indexes.has(seq1Idx))
          (indexes.get(seq1Idx)! as number[]).push(seq2Idx);
        else
          indexes.set(seq1Idx, [seq2Idx]);

        //Working with seq2monomer
        const seq2monomer = tempDataElement.seq2monomer;
        if (!substitutionsInfo.has(seq2monomer))
          substitutionsInfo.set(seq2monomer, new Map());

        positionsMap = substitutionsInfo.get(seq2monomer)!;
        if (!positionsMap.has(position))
          positionsMap.set(position, new Map());

        indexes = positionsMap.get(position)!;
        if (indexes.has(seq2Idx))
          (indexes.get(seq2Idx)! as number[]).push(seq1Idx);
        else
          indexes.set(seq2Idx, [seq1Idx]);
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
  */
}
