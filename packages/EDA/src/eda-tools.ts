// Exploratory data analysis (EDA) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_principalComponentAnalysisInWebWorker, _principalComponentAnalysis,
  _partialLeastSquareRegressionInWebWorker,
  _principalComponentAnalysisNipals,
} from '../wasm/EDAAPI';

import {checkWasmDimensionReducerInputs, checkUMAPinputs, checkTSNEinputs,
  getRowsOfNumericalColumnns} from './utils';

// Principal components analysis (PCA)
export async function computePCA(table: DG.DataFrame, features: DG.ColumnList, components: number,
  center: boolean, scale: boolean, nipals: boolean): Promise<DG.DataFrame> {
  checkWasmDimensionReducerInputs(features, components);

  const centerNum = center ? 1 : 0;
  const scaleNum = scale ? 1 : 0;

  try {
    const df = nipals ? _principalComponentAnalysisNipals(table, features, components) :
      _principalComponentAnalysis(table, features, components, centerNum, scaleNum);
    return df;
  } catch (e) {
    return DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('1', new Float32Array(table.rowCount)),
      DG.Column.fromFloat32Array('2', new Float32Array(table.rowCount)),
    ]);
  }

  return _principalComponentAnalysis(table, features, components, centerNum, scaleNum);
}

// Partial least square regression (PLS): TO REMOVE
export async function computePLS(
  table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, components: number,
): Promise<any> {
  // Inputs are checked in the same manner as in PCA, since the same computations are applied.
  checkWasmDimensionReducerInputs(features, components);

  return await _partialLeastSquareRegressionInWebWorker(table, features, predict, components);
}

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
