import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {mmDistanceFunctionArgs} from '@datagrok-libraries/ml/src/macromolecule-distance-functions/types';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {getMonomerSubstitutionMatrix} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

export interface ISequenceSpaceResult {
  distance?: Float32Array;
  coordinates: DG.ColumnList;
}

export async function getEncodedSeqSpaceCol(
  seqCol: DG.Column, similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames, fingerprintType: string = 'Morgan'
): Promise<{seqList:string[], options: {[_:string]: any}}> {
// encodes sequences using utf charachters to also support multichar and non fasta sequences
  const ncUH = UnitsHandler.getOrCreate(seqCol);
  const seqList = seqCol.toList();
  const splitter = ncUH.getSplitter();
  const seqColLength = seqList.length;
  let charCodeCounter = 36;
  const charCodeMap = new Map<string, string>();
  for (let i = 0; i < seqColLength; i++) {
    const seq = seqList[i];
    if (seqList[i] === null || seqCol.isNone(i)) {
      seqList[i] = null;
      continue;
    }
    seqList[i] = '';
    const splittedSeq = splitter(seq);
    for (let j = 0; j < splittedSeq.length; j++) {
      const char = splittedSeq[j];
      if (!charCodeMap.has(char)) {
        charCodeMap.set(char, String.fromCharCode(charCodeCounter));
        charCodeCounter++;
      }
      seqList[i] += charCodeMap.get(char)!;
    }
  }
  let options = {};
  if (similarityMetric === MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE) {
    const monomers = Array.from(charCodeMap.keys());
    const monomerRes = await getMonomerSubstitutionMatrix(monomers, fingerprintType);
    // the susbstitution matrix contains similarity, but we need distances
    monomerRes.scoringMatrix.forEach((row, i) => {
      row.forEach((val, j) => {
        monomerRes.scoringMatrix[i][j] = 1 - val;
      });
    });
    const monomerHashToMatrixMap: {[_: string]: number} = {};
    Object.entries(monomerRes.alphabetIndexes).forEach(([key, value]) => {
      monomerHashToMatrixMap[charCodeMap.get(key)!] = value;
    });
    // sets distance function args in place.
    options = {scoringMatrix: monomerRes.scoringMatrix,
      alphabetIndexes: monomerHashToMatrixMap} satisfies mmDistanceFunctionArgs;
  } else if (similarityMetric === MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH) {
    const monomers = Array.from(charCodeMap.keys());
    const monomerRes = await getMonomerSubstitutionMatrix(monomers, fingerprintType);
    // the susbstitution matrix contains similarity, but we need distances
    // monomerRes.scoringMatrix.forEach((row, i) => {
    //   row.forEach((val, j) => {
    //     monomerRes.scoringMatrix[i][j] = 1 - val;
    //   });
    // });
    const monomerHashToMatrixMap: {[_: string]: number} = {};
    Object.entries(monomerRes.alphabetIndexes).forEach(([key, value]) => {
      monomerHashToMatrixMap[charCodeMap.get(key)!] = value;
    });
    // sets distance function args in place.
    options = {scoringMatrix: monomerRes.scoringMatrix,
      alphabetIndexes: monomerHashToMatrixMap} satisfies mmDistanceFunctionArgs;
  }
  return {seqList, options};
}
