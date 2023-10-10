import * as DG from 'datagrok-api/dg';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {BitArrayMetrics, StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {ISequenceSpaceParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {invalidateMols, MONOMERIC_COL_TAGS} from '../substructure-search/substructure-search';
import {mmDistanceFunctionArgs} from '@datagrok-libraries/ml/src/macromolecule-distance-functions/types';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {calculateMonomerSimilarity} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import * as grok from 'datagrok-api/grok';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

export interface ISequenceSpaceResult {
  distance?: Float32Array;
  coordinates: DG.ColumnList;
}

export async function sequenceSpace(spaceParams: ISequenceSpaceParams): Promise<ISequenceSpaceResult> {
  // code deprecated since seqCol is encoded
  /*    let preparedData: any;
  if (!(spaceParams.seqCol!.tags[DG.TAGS.UNITS] === 'HELM')) {
    const sep = spaceParams.seqCol.getTag(UnitsHandler.TAGS.separator);
    const sepFinal = sep ? sep === '.' ? '\\\.' : sep : '-';
    const regex = new RegExp(sepFinal, 'g');
    if (Object.keys(AvailableMetrics['String']).includes(spaceParams.similarityMetric))
      preparedData = spaceParams.seqCol.toList().map((v: string) => v.replace(regex, '')) as string[];
    else
      preparedData = spaceParams.seqCol.toList().map((v: string) => v.replace(regex, '')) as string[];
  } else {
    preparedData = spaceParams.seqCol.toList();
  }  */

  const sequenceSpaceResult = await reduceDimensinalityWithNormalization(
    spaceParams.seqCol.toList(),
    spaceParams.methodName,
    spaceParams.similarityMetric as StringMetrics | BitArrayMetrics,
    spaceParams.options);
  const cols: DG.Column[] = spaceParams.embedAxesNames.map(
    (name: string, index: number) => DG.Column.fromFloat32Array(name, sequenceSpaceResult.embedding[index]));
  return {distance: sequenceSpaceResult.distance, coordinates: new DG.ColumnList(cols)};
}

export async function sequenceSpaceByFingerprints(spaceParams: ISequenceSpaceParams): Promise<ISequenceSpaceResult> {
  if (spaceParams.seqCol.version !== spaceParams.seqCol.temp[MONOMERIC_COL_TAGS.LAST_INVALIDATED_VERSION])
    //we expect only string columns here
    await invalidateMols(spaceParams.seqCol as unknown as DG.Column<string>, false);

  const result = await grok.functions.call('Chem:getChemSpaceEmbeddings', {
    col: spaceParams.seqCol.temp[MONOMERIC_COL_TAGS.MONOMERIC_MOLS],
    methodName: spaceParams.methodName,
    similarityMetric: spaceParams.similarityMetric,
    xAxis: spaceParams.embedAxesNames[0],
    yAxis: spaceParams.embedAxesNames[1],
    options: spaceParams.options,
  });
  return result;
}

export async function getSequenceSpace(spaceParams: ISequenceSpaceParams,
  progressFunc?: (epochNum: number, epochsLength: number, embedding: number[][]) => void
): Promise<ISequenceSpaceResult> {
  const ncUH = UnitsHandler.getOrCreate(spaceParams.seqCol);

  //const distanceFName = ncUH.isMsa() ? MmDistanceFunctionsNames.HAMMING : MmDistanceFunctionsNames.LEVENSHTEIN;
  const seqList = spaceParams.seqCol.toList();

  const splitter = ncUH.getSplitter();
  const seqColLength = seqList.length;
  let charCodeCounter = 36;
  const charCodeMap = new Map<string, string>();
  for (let i = 0; i < seqColLength; i++) {
    const seq = seqList[i];
    if (seqList[i] === null || spaceParams.seqCol.isNone(i)) {
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

  if (spaceParams.similarityMetric === MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE) {
    const monomers = Array.from(charCodeMap.keys());
    const monomerRes = await calculateMonomerSimilarity(monomers);
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
    spaceParams.options.distanceFnArgs = {scoringMatrix: monomerRes.scoringMatrix,
      alphabetIndexes: monomerHashToMatrixMap} satisfies mmDistanceFunctionArgs;
  }

  const sequenceSpaceResult = await reduceDimensinalityWithNormalization(
    seqList,
    spaceParams.methodName,
    spaceParams.similarityMetric,
    spaceParams.options,
    true, progressFunc);
  const cols: DG.Column[] = spaceParams.embedAxesNames.map(
    (name: string, index: number) => DG.Column.fromFloat32Array(name, sequenceSpaceResult.embedding[index]));
  return {distance: sequenceSpaceResult.distance, coordinates: new DG.ColumnList(cols)};
}

export function getEmbeddingColsNames(df: DG.DataFrame) {
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it: string) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}
