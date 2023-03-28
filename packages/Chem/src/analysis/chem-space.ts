import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
// TODO: clean up this module
import {chemGetFingerprints} from '../chem-searches';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {Fingerprint} from '../utils/chem-common';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {IReduceDimensionalityResult} from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {ISequenceSpaceParams, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';


export async function chemSpace(spaceParams: ISequenceSpaceParams): Promise<ISequenceSpaceResult> {
  const fpColumn = await chemGetFingerprints(spaceParams.seqCol, Fingerprint.Morgan);
  const chemSpaceResult: IReduceDimensionalityResult = await reduceDimensinalityWithNormalization(
    fpColumn,
    spaceParams.methodName,
    spaceParams.similarityMetric as BitArrayMetrics,
    spaceParams.options);
  const cols: DG.Column[] = spaceParams.embedAxesNames.map((name: string, index: number) => DG.Column.fromFloat32Array(name, chemSpaceResult.embedding[index]));
  return {distance: chemSpaceResult.distance, coordinates: new DG.ColumnList(cols)};
}

export function getEmbeddingColsNames(df: DG.DataFrame) {
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it: string) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}
