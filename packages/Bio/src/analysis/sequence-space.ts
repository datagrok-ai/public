import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {mmDistanceFunctionArgs} from '@datagrok-libraries/ml/src/macromolecule-distance-functions/types';
import {getMonomerSubstitutionMatrix} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

import {_package} from '../package';

export interface ISequenceSpaceResult {
  distance?: Float32Array;
  coordinates: DG.ColumnList;
}

export interface DistanceFunctionParams {
  fingerprintType?: string;
  gapOpen?: number;
  gapExtend?: number;
}

export async function getEncodedSeqSpaceCol(
  seqCol: DG.Column,
  similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames,
  params: DistanceFunctionParams = {}
): Promise<{ seqList: string[], options: { [_: string]: any } }> {
  // encodes sequences using utf characters to also support multichar and non fasta sequences
  const rowCount = seqCol.length;
  const sh = _package.seqHelper.getSeqHandler(seqCol);
  const encList = Array<string>(rowCount);
  let charCodeCounter = 1; // start at 1, 0 is reserved for null.
  const charCodeMap = new Map<string, string>();
  const seqColCats = seqCol.categories;
  const seqColRawData = seqCol.getRawData();

  for (let rowIdx = 0; rowIdx < rowCount; rowIdx++) {
    const catI = seqColRawData[rowIdx];
    const seq = seqColCats[catI];
    if (seq === null || seqCol.isNone(rowIdx)) {
      //@ts-ignore
      encList[rowIdx] = null;
      continue;
    }
    encList[rowIdx] = '';
    const splittedSeq = sh.getSplitted(rowIdx);
    for (let j = 0; j < splittedSeq.length; j++) {
      const char = splittedSeq.getCanonical(j);
      if (!charCodeMap.has(char)) {
        charCodeMap.set(char, String.fromCharCode(charCodeCounter));
        charCodeCounter++;
      }
      encList[rowIdx] += charCodeMap.get(char)!;
    }
  }

  let options = {};//as mmDistanceFunctionArgs;

  // Handle fingerprint-based distance functions
  if (
    similarityMetric === MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE ||
    similarityMetric === MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH
  ) {
    const fingerprintType = params.fingerprintType || 'Morgan';
    const monomers = Array.from(charCodeMap.keys());
    const monomerRes = await getMonomerSubstitutionMatrix(monomers, fingerprintType);

    const monomerHashToMatrixMap: { [_: string]: number } = {};
    Object.entries(monomerRes.alphabetIndexes).forEach(([key, value]) => {
      monomerHashToMatrixMap[charCodeMap.get(key)!] = value;
    });

    // sets distance function args in place.
    const maxLength = encList.reduce((acc, val) => Math.max(acc, val?.length || 0), 0);
    options = {
      scoringMatrix: monomerRes.scoringMatrix,
      alphabetIndexes: monomerHashToMatrixMap,
      maxLength
    };
  }

  // Add gap penalties for Needleman-Wunsch
  if (similarityMetric === MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH) {
    options = {
      ...options,
      gapOpen: params.gapOpen ?? 1,
      gapExtend: params.gapExtend ?? 1,
    };
  }

  return {seqList: encList, options};
}
