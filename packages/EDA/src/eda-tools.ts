// Exploratory data analysis (EDA) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {flatBatch, matrixIpc} from './crux/crux-bridge';
import {runCrux} from './crux/crux-model';

import {checkDimensionReducerInputs, checkUMAPinputs, checkTSNEinputs, getRowsOfNumericalColumnns,
  centerScaleDataFrame, extractNonConstantColsDf, standardizedRows} from './utils';

// Principal components analysis (PCA)
export async function computePCA(table: DG.DataFrame, features: DG.ColumnList, components: number,
  toCenter: boolean, toScale: boolean): Promise<DG.DataFrame> {
  checkDimensionReducerInputs(features, components);

  const rowCount = table.rowCount;
  const cols = extractNonConstantColsDf(features).columns.toList();
  const m = cols.length;
  const k = Math.min(components, m, rowCount);

  const scores: DG.Column[] = [];
  if (k > 0) {
    const out = await runCrux('ml.pca', matrixIpc(standardizedRows(cols), rowCount, m), {n_components: k,
      svd_solver: 'auto', random_state: 0});
    const flat = flatBatch(out, 'x');
    for (let a = 0; a < k; ++a) {
      const values = new Float32Array(rowCount);
      for (let i = 0; i < rowCount; ++i)
        values[i] = flat[i * k + a];
      scores.push(DG.Column.fromFloat32Array(`${a + 1}`, values));
    }
  }

  for (let a = k; a < components; ++a)
    scores.push(DG.Column.fromFloat32Array(`${a + 1}`, new Float32Array(rowCount)));

  return centerScaleDataFrame(DG.DataFrame.fromColumns(scores), toCenter, toScale);
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
