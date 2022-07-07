import * as DG from 'datagrok-api/dg';
import {AvailableMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {BitArrayMetrics, StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {ISequenceSpaceParams} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';

export interface ISequenceSpaceResult {
  distance: Matrix;
  coordinates: DG.ColumnList;
}

export async function sequenceSpace(spaceParams: ISequenceSpaceParams): Promise<ISequenceSpaceResult> {
  let preparedData: any;
  if (!(spaceParams.seqCol!.tags[DG.TAGS.UNITS] === 'HELM')) {
    const sep = spaceParams.seqCol.getTag('separator');
    const sepFinal = sep ? sep === '.' ? '\\\.' : sep : '-';
    const regex = new RegExp(sepFinal, 'g');
    if (Object.keys(AvailableMetrics['String']).includes(spaceParams.similarityMetric))
      preparedData = spaceParams.seqCol.toList().map((v: string) => v.replace(regex, '')) as string[];
    else
      preparedData = spaceParams.seqCol.toList().map((v: string) => v.replace(regex, '')) as string[];
  } else {
    preparedData = spaceParams.seqCol.toList();
  }

  const sequenceSpaceResult = await reduceDimensinalityWithNormalization(
    preparedData,
    spaceParams.methodName,
    spaceParams.similarityMetric as StringMetrics | BitArrayMetrics,
    spaceParams.options);
  const cols: DG.Column[] = spaceParams.embedAxesNames.map(
    (name: string, index: number) => DG.Column.fromFloat32Array(name, sequenceSpaceResult.embedding[index]));
  return {distance: sequenceSpaceResult.distance, coordinates: new DG.ColumnList(cols)};
}


export function getEmbeddingColsNames(df: DG.DataFrame) {
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it: string) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}
