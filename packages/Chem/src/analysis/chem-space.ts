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
import { malformedDataWarning } from '../utils/malformed-data-utils';


export async function chemSpace(spaceParams: ISequenceSpaceParams): Promise<ISequenceSpaceResult> {
  const fpColumn = await chemGetFingerprints(spaceParams.seqCol, Fingerprint.Morgan);
  const malformedIdxs = malformedDataWarning(fpColumn);
  const chemSpaceResult: IReduceDimensionalityResult = await reduceDimensinalityWithNormalization(
    fpColumn,
    spaceParams.methodName,
    spaceParams.similarityMetric as BitArrayMetrics,
    spaceParams.options);
  malformedIdxs.forEach((idx: number | null) => {
    setNullForMalformedData(chemSpaceResult.embedding, idx!);
    if (chemSpaceResult.distance) {
      setNullForMalformedData(chemSpaceResult.distance, idx!);
    }
  })
  const cols: DG.Column[] = spaceParams.embedAxesNames.map((name: string, index: number) => DG.Column.fromFloat32Array(name, chemSpaceResult.embedding[index]));
  return {distance: chemSpaceResult.distance, coordinates: new DG.ColumnList(cols)};
}

function setNullForMalformedData(matrix: Matrix, idx: number) {
  for (const col of matrix) {
    col[idx] = DG.FLOAT_NULL; 
  };
}

export function getEmbeddingColsNames(df: DG.DataFrame) {
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it: string) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}
