import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {chemGetFingerprints} from '../chem-searches';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {Fingerprint} from '../utils/chem-common';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {ISequenceSpaceParams, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {malformedDataWarning, setEmptyBitArraysForMalformed} from '../utils/malformed-data-utils';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {DimReductionMethods, IReduceDimensionalityResult, ITSNEOptions, IUMAPOptions}
  from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {dmLinearIndex} from '@datagrok-libraries/ml/src/distance-matrix';


export async function chemSpace(spaceParams: ISequenceSpaceParams): Promise<ISequenceSpaceResult> {
  const fpColumn = await chemGetFingerprints(spaceParams.seqCol, Fingerprint.Morgan, false);
  const emptyAndMalformedIdxs = fpColumn.map((el: BitArray | null, idx: number) =>
    !el ? idx : null).filter((it) => it !== null);
  malformedDataWarning(fpColumn, spaceParams.seqCol);
  /* need to replace nulls with empty BitArrays since dimensionality reducing algorithmns
  fail in case fpColumn contains nulls. TODO: fix on dim reduction side */
  setEmptyBitArraysForMalformed(fpColumn);
  const chemSpaceResult: IReduceDimensionalityResult = await reduceDimensinalityWithNormalization(
    fpColumn as BitArray[],
    spaceParams.methodName,
    spaceParams.similarityMetric,
    spaceParams.options);
  emptyAndMalformedIdxs.forEach((idx: number | null) => {
    setNullForEmptyAndMalformedData(chemSpaceResult.embedding, idx!);
    if (chemSpaceResult.distance)
      setNullForEmptyAndMalformedDistanceData(chemSpaceResult.distance, idx!, spaceParams.seqCol.length);
  });
  const cols: DG.Column[] = spaceParams.embedAxesNames.map((name: string, index: number) =>
    DG.Column.fromFloat32Array(name, chemSpaceResult.embedding[index]));
  return {distance: chemSpaceResult.distance, coordinates: new DG.ColumnList(cols)};
}

function setNullForEmptyAndMalformedData(matrix: Matrix, idx: number) {
  for (const col of matrix)
    col[idx] = DG.FLOAT_NULL;
}

function setNullForEmptyAndMalformedDistanceData(matrix: Float32Array, idx: number, length: number) {
  const linearIdx = dmLinearIndex(length);
  for (let i = 0; i < length; i++)
    matrix[linearIdx(idx, i)] = DG.FLOAT_NULL;
}

export function getEmbeddingColsNames(df: DG.DataFrame) {
  const axes = ['Embed_X', 'Embed_Y'];
  const colNameInd = df.columns.names().filter((it: string) => it.includes(axes[0])).length + 1;
  return axes.map((it) => `${it}_${colNameInd}`);
}

export async function runChemSpace(table: DG.DataFrame, molecules: DG.Column, methodName: DimReductionMethods,
  similarityMetric: BitArrayMetrics = BitArrayMetricsNames.Tanimoto, plotEmbeddings: boolean,
  options?: IUMAPOptions | ITSNEOptions):  Promise<DG.Viewer | undefined>{
  const embedColsNames = getEmbeddingColsNames(table);

  const chemSpaceParams = {
    seqCol: molecules,
    methodName: methodName,
    similarityMetric: similarityMetric as BitArrayMetrics,
    embedAxesNames: [embedColsNames[0], embedColsNames[1]],
    options: options,
  };
  const chemSpaceRes = await chemSpace(chemSpaceParams);
  const embeddings = chemSpaceRes.coordinates;

  for (const col of embeddings)
    table.columns.add(col);

  if (plotEmbeddings) {
    return grok.shell
      .tableView(table.name)
      .scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Chem space'});
  }
}
