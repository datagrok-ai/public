import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {chemGetFingerprints} from '../chem-searches';
import {Fingerprint} from '../utils/chem-common';
import {Matrix, Options} from '@datagrok-libraries/utils/src/type-declarations';
import {ISequenceSpaceParams, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {malformedDataWarning, setEmptyBitArraysForMalformed} from '../utils/malformed-data-utils';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {BitArrayMetrics, BitArrayMetricsNames, KnownMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {getNormalizedEmbeddings} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/embeddings-space';
import {multiColReduceDimensionality}
  from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {ITSNEOptions, IUMAPOptions} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/multi-column-dim-reducer';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';


export async function chemSpace(spaceParams: ISequenceSpaceParams,
  progressFunc?: (epochNum: number, epochsLength: number, embedding: number[][]) => void,
): Promise<ISequenceSpaceResult> {
  const fpColumn = await chemGetFingerprints(spaceParams.seqCol, Fingerprint.Morgan, false);
  const emptyAndMalformedIdxs = fpColumn.map((el: BitArray | null, idx: number) =>
    !el ? idx : null).filter((it) => it !== null);
  malformedDataWarning(fpColumn, spaceParams.seqCol);
  /* need to replace nulls with empty BitArrays since dimensionality reducing algorithmns
  fail in case fpColumn contains nulls. TODO: fix on dim reduction side */
  setEmptyBitArraysForMalformed(fpColumn);
  const chemSpaceResult = await getNormalizedEmbeddings([fpColumn], spaceParams.methodName,
    [spaceParams.similarityMetric], [1], 'MANHATTAN', {...spaceParams.options, distanceFnArgs: [{}]}, progressFunc);
  emptyAndMalformedIdxs.forEach((idx: number | null) => {
    setNullForEmptyAndMalformedData(chemSpaceResult, idx!);
  });
  const cols: DG.Column[] = spaceParams.embedAxesNames.map((name: string, index: number) =>
    DG.Column.fromFloat32Array(name, chemSpaceResult[index]));
  return {coordinates: new DG.ColumnList(cols)};
}

function setNullForEmptyAndMalformedData(matrix: Matrix, idx: number) {
  for (const col of matrix)
    col[idx] = DG.FLOAT_NULL;
}

export async function runChemSpace(table: DG.DataFrame, molecules: DG.Column, methodName: DimReductionMethods,
  similarityMetric: BitArrayMetrics = BitArrayMetricsNames.Tanimoto, plotEmbeddings: boolean,
  options?: (IUMAPOptions | ITSNEOptions) & Options, preprocessingFunction?: DG.Func,
  clusterEmbeddings?: boolean, tableView?: DG.TableView): Promise<DG.Viewer | undefined> {
  if (!preprocessingFunction)
    preprocessingFunction = DG.Func.find({name: 'getFingerprints', package: 'Chem'})[0];
  options ??= {};
  const res = await multiColReduceDimensionality(table, [molecules], methodName,
    [similarityMetric as KnownMetrics], [1], [preprocessingFunction], 'MANHATTAN',
    plotEmbeddings, clusterEmbeddings ?? false,
    /* dimRedOptions */ {...options, preprocessingFuncArgs: [options.preprocessingFuncArgs ?? {}]},
    /* uiOptions */{
      fastRowCount: options.fastRowCount ?? 10000,
      scatterPlotName: 'Chemical space',
      bypassLargeDataWarning: options?.[BYPASS_LARGE_DATA_WARNING],
      tableView: tableView,
    });
  return res;
}
