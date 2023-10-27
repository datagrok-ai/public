import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {chemGetFingerprints} from '../chem-searches';
import {reduceDimensinalityWithNormalization} from '@datagrok-libraries/ml/src/sequence-space';
import {Fingerprint} from '../utils/chem-common';
import {Matrix, Options} from '@datagrok-libraries/utils/src/type-declarations';
import {ISequenceSpaceParams, ISequenceSpaceResult} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {malformedDataWarning, setEmptyBitArraysForMalformed} from '../utils/malformed-data-utils';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {DimReductionMethods, IReduceDimensionalityResult, ITSNEOptions, IUMAPOptions}
  from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {BitArrayMetrics, BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {dmLinearIndex} from '@datagrok-libraries/ml/src/distance-matrix';
import {DIMENSIONALITY_REDUCER_TERMINATE_EVENT}
  from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
//import {SHOW_SCATTERPLOT_PROGRESS} from '@datagrok-libraries/ml/src/functionEditors/seq-space-base-editor';


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
  const chemSpaceResult: IReduceDimensionalityResult = await reduceDimensinalityWithNormalization(
    fpColumn as BitArray[],
    spaceParams.methodName,
    spaceParams.similarityMetric,
    spaceParams.options, false, progressFunc); //to be changed to true
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
  options?: (IUMAPOptions | ITSNEOptions) & Options,
  progressF?: (percent: number) => void): Promise<DG.Viewer | undefined> {
  const embedColsNames = getEmbeddingColsNames(table);
  let scatterPlot: DG.ScatterPlotViewer | undefined = undefined;
  try {
    function progressFunc(_nEpoch: number, epochsLength: number, embeddings: number[][]) {
      let embedXCol: DG.Column | null = null;
      let embedYCol: DG.Column | null = null;
      if (!table.columns.names().includes(embedColsNames[0])) {
        embedXCol = table.columns.add(DG.Column.float(embedColsNames[0], table.rowCount));
        embedYCol = table.columns.add(DG.Column.float(embedColsNames[1], table.rowCount));
        if (plotEmbeddings) {
          scatterPlot = grok.shell
            .tableView(table.name)
            .scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Chem space'});
        }
      } else {
        embedXCol = table.columns.byName(embedColsNames[0]);
        embedYCol = table.columns.byName(embedColsNames[1]);
      }
      //if (options?.[SHOW_SCATTERPLOT_PROGRESS]) {
        scatterPlot?.root && ui.setUpdateIndicator(scatterPlot!.root, false);
        embedXCol.init((i) => embeddings[i][0]);
        embedYCol.init((i) => embeddings[i][1]);
      //}
      const progress = (_nEpoch / epochsLength * 100);
      progressF && progressF(progress);
    }
    const chemSpaceParams = {
      seqCol: molecules,
      methodName: methodName,
      similarityMetric: similarityMetric as BitArrayMetrics,
      embedAxesNames: [embedColsNames[0], embedColsNames[1]],
      options: options,
    };
    // const chemSpaceRes = await chemSpace(chemSpaceParams, progressFunc);
    // const embeddings = chemSpaceRes.coordinates;


    table.columns.add(DG.Column.float(embedColsNames[0], table.rowCount));
    table.columns.add(DG.Column.float(embedColsNames[1], table.rowCount));
    if (plotEmbeddings) {
      scatterPlot = grok.shell
        .tableView(table.name)
        .scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Chem space'});
      ui.setUpdateIndicator(scatterPlot.root, true);
    }
    let resolveF: Function | null = null;

    const sub = grok.events.onViewerClosed.subscribe((args) => {
      const v = args.args.viewer as unknown as DG.Viewer<any>;
      if (v?.getOptions()?.look?.title && scatterPlot?.getOptions()?.look?.title &&
          v?.getOptions()?.look?.title === scatterPlot?.getOptions()?.look?.title) {
        grok.events.fireCustomEvent(DIMENSIONALITY_REDUCER_TERMINATE_EVENT, {});
        sub.unsubscribe();
        resolveF?.();
      }
    });

    const chemSpaceResPromise = new Promise<ISequenceSpaceResult | undefined>(async (resolve) =>{
      resolveF = resolve;
      const r = await chemSpace(chemSpaceParams, progressFunc);
      resolve(r);
    });

    const chemSpaceRes = await chemSpaceResPromise;
    if (!chemSpaceRes)
      return undefined;
    const embeddings = chemSpaceRes.coordinates;
    if (!table.columns.names().includes(embedColsNames[0])) {
      for (const col of embeddings)
        table.columns.add(col);
    } else {
      for (const col of embeddings)
        table.columns.byName(col.name).init((i) => col.get(i));
    }

    if (plotEmbeddings) {
      if (!scatterPlot) {
        scatterPlot = grok.shell
          .tableView(table.name)
          .scatterPlot({x: embedColsNames[0], y: embedColsNames[1], title: 'Chem space'});
      }
      ui.setUpdateIndicator(scatterPlot.root, false);
      return scatterPlot;
    }
  } catch (e) {
    console.error(e);
  }
}
