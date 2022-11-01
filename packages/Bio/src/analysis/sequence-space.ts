import * as DG from 'datagrok-api/dg';
import {AvailableMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {BitArrayMetrics, StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {ISequenceSpaceParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {invalidateMols, MONOMERIC_COL_TAGS} from '../substructure-search/substructure-search';
import * as grok from 'datagrok-api/grok';

export interface ISequenceSpaceResult {
  distance: Matrix;
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
    await invalidateMols(spaceParams.seqCol, false);

  const result = await grok.functions.call('Chem:getChemSpaceEmbeddings', {
    col: spaceParams.seqCol.temp[MONOMERIC_COL_TAGS.MONOMERIC_MOLS],
    methodName: spaceParams.methodName,
    similarityMetric: spaceParams.similarityMetric,
    xAxis: spaceParams.embedAxesNames[0],
    yAxis: spaceParams.embedAxesNames[1]
  });
  return result;
}


export function getEmbeddingColsNames(df: DG.DataFrame) {
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it: string) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}
