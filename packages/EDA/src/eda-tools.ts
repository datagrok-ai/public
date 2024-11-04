// Exploratory data analysis (EDA) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_principalComponentAnalysisInWebWorker, _principalComponentAnalysis,
  _partialLeastSquareRegressionInWebWorker,
  _principalComponentAnalysisNipals, _principalComponentAnalysisNipalsInWebWorker,
} from '../wasm/EDAAPI';

import {checkWasmDimensionReducerInputs, checkUMAPinputs, checkTSNEinputs, NIPALS_PREFER_COLS_COUNT,
  getRowsOfNumericalColumnns, centerScaleDataFrame, extractNonConstantColsDf} from './utils';

// Principal components analysis (PCA)
export async function computePCA(table: DG.DataFrame, features: DG.ColumnList, components: number,
  toCenter: boolean, toScale: boolean): Promise<DG.DataFrame> {
  checkWasmDimensionReducerInputs(features, components);

  const rowCount = table.rowCount;

  // Extract non-const cols dataframe
  const nonConstData = extractNonConstantColsDf(features);
  const nonConstColsCount = nonConstData.columns.length;

  // Return zero columns if data is constant
  if (nonConstColsCount === 0) {
    const cols: DG.Column[] = [];

    for (let i = 0; i < components; ++i)
      cols.push(DG.Column.fromFloat32Array(`${i + 1}`, new Float32Array(rowCount).fill(0)));

    return DG.DataFrame.fromColumns(cols);
  }

  const zeroColsToAdd = (nonConstColsCount < components) ? (components - nonConstColsCount) : 0;
  const componentsToCompute = Math.min(components, nonConstColsCount);

  let output: DG.DataFrame | undefined = undefined;

  // PCA
  if (nonConstColsCount > NIPALS_PREFER_COLS_COUNT)
    output = await _principalComponentAnalysisNipalsInWebWorker(table, features, componentsToCompute);
  else {
    //try to apply the classic algorithm
    const res = await _principalComponentAnalysisInWebWorker(table, features, componentsToCompute);

    if (res !== -1) // the classic succeed
      output = centerScaleDataFrame(res, toCenter, toScale);
    else // the classic failed
      output = await _principalComponentAnalysisNipalsInWebWorker(table, features, componentsToCompute);
  }

  if (output === undefined)
    throw new Error('Failed to compute PCA');

  output = centerScaleDataFrame(output, toCenter, toScale);

  const cols = output.columns;
  const count = cols.length;

  // Add zero columns (with respect to the const cols count)
  for (let i = 0; i < zeroColsToAdd; ++i)
    cols.add(DG.Column.fromFloat32Array(`${count + i + 1}`, new Float32Array(rowCount).fill(0)));

  return output;
} // computePCA

// Uniform Manifold Approximation and Projection (UMAP)
export async function computeUMAP(features: DG.ColumnList, components: number, epochs: number,
  neighbors: number, minDist: number, spread: number): Promise<DG.DataFrame> {
  // check inputs
  checkUMAPinputs(features, components, epochs, neighbors, minDist, spread);

  // get row-by-row data
  const data = getRowsOfNumericalColumnns(features);

  let workerOutput: any;

  // UMAP in webworker
  const promise = new Promise((resolve, _reject) => {
    const worker = new Worker(new URL('workers/umap-worker.ts', import.meta.url));

    worker.postMessage({
      data: data,
      options: {
        nComponents: components,
        nEpochs: epochs,
        nNeighbors: neighbors,
        minDist: minDist,
        spread: spread,
      }});

    worker.onmessage = function(e) {
      worker.terminate();
      resolve(e.data.embeddings);
    };
  });

  await promise.then(
    (result) => {workerOutput = result;},
    (_error) => {throw new Error('applying UMAP fails.');},
  );

  const embeddings = workerOutput as number[][];
  const rowCount = embeddings.length;
  const range = [...Array(components).keys()];

  // Create output

  // columns data
  const umapColumnsData = range.map((_) => new Float32Array(rowCount));

  // perform transponation
  for (let i = 0; i < rowCount; ++i) {
    for (let j = 0; j < components; ++j)
      umapColumnsData[j][i] = embeddings[i][j];
  }

  return DG.DataFrame.fromColumns(range.map((i) =>
    DG.Column.fromFloat32Array('UMAP' + i.toString(), umapColumnsData[i]),
  ));
} // computeUMAP

// t-distributed stochastic neighbor embedding (t-SNE)
export async function computeTSNE(features: DG.ColumnList, components: number,
  learningRate: number, perplexity: number, iterations: number): Promise<DG.DataFrame> {
  // check inputs
  checkTSNEinputs(features, components, learningRate, perplexity, iterations);

  // get row-by-row data
  const data = getRowsOfNumericalColumnns(features);

  let workerOutput: any;

  // t-SNE in webworker
  const promise = new Promise((resolve, _reject) => {
    const worker = new Worker(new URL('workers/tsne-worker.ts', import.meta.url));

    worker.postMessage({
      data: data,
      options: {
        learningRate: learningRate,
        perplexity: perplexity,
        components: components,
        iterations: iterations,
      }});

    worker.onmessage = function(e) {
      worker.terminate();
      resolve(e.data.embeddings);
    };
  });

  await promise.then(
    (result) => {workerOutput = result;},
    (_error) => {throw new Error('applying t-SNE fails.');},
  );

  const embeddings = workerOutput as any[];

  const rowCount = embeddings.length;
  const range = [...Array(components).keys()];

  // Create output

  // columns data
  const umapColumnsData = range.map((_) => new Float32Array(rowCount));

  // perform transponation
  for (let i = 0; i < rowCount; ++i) {
    for (let j = 0; j < components; ++j)
      umapColumnsData[j][i] = embeddings[i][j];
  }

  return DG.DataFrame.fromColumns(range.map((i) =>
    DG.Column.fromFloat32Array('tSNE' + i.toString(), umapColumnsData[i]),
  ));
} // computeTSNE
