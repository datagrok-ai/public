import * as DG from 'datagrok-api/dg';
// TODO: clean up this module
import {chemGetFingerprints} from '../chem-searches';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {Fingerprint} from '../utils/chem-common';
import { Matrix } from '@datagrok-libraries/utils/src/type-declarations';
import { IReduceDimensionalityResult } from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';

export interface IChemSpaceResult {
  distance: Matrix;
  coordinates: DG.ColumnList;
}

export async function chemSpace(molColumn: DG.Column, methodName: string, similarityMetric: string,
  axes: string[], options?: any): Promise<IChemSpaceResult> {

  const fpColumn = await chemGetFingerprints(molColumn, Fingerprint.Morgan);
  const chemSpaceResult: IReduceDimensionalityResult= await reduceDimensinalityWithNormalization(
    fpColumn,
    methodName,
    similarityMetric as BitArrayMetrics,
    options);
  const cols: DG.Column[] = axes.map((name, index) => DG.Column.fromFloat32Array(name, chemSpaceResult.embedding[index]))
  return {distance: chemSpaceResult.distance, coordinates: new DG.ColumnList(cols)};
}

export function getEmbeddingColsNames(df: DG.DataFrame){
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}
