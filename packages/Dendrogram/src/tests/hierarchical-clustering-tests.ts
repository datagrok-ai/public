import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {category, test, expect, expectObject, expectArray, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {DistanceMetric} from '@datagrok-libraries/bio/src/trees';
import {DistanceMatrix} from '@datagrok-libraries/ml/src/distance-matrix';
import {ClusterMatrix} from '@datagrok-libraries/bio/src/trees';
import {getClusterMatrixWorker} from '@datagrok-libraries/math';

import {hierarchicalClusteringUI} from '../utils/hierarchical-clustering';

import {_package} from '../package-test';

/*
https://onecompiler.com/python

import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist


data = pd.DataFrame({'x': [8,6,5,1], 'y': [0,0,0,0]})
# data = pd.DataFrame({'x': [0,4,0], 'y': [0,0,3]})

distance_name = 'euclidean'

column_array = data[data.columns].to_numpy()
sys.stdout.write('column_array\n')
sys.stdout.write(str(column_array))
sys.stdout.write('\n\n')

dist1 = pdist(column_array, distance_name)
sys.stdout.write('dist1\n')
sys.stdout.write(str(dist1))
sys.stdout.write('\n\n')

dist_matrix = cdist(column_array, column_array)
sys.stdout.write('dist_matrix\n')
sys.stdout.write(str(dist_matrix))
sys.stdout.write('\n\n')

dist_list = pdist(dist_matrix, distance_name)

result = pd.DataFrame.from_dict({'distance': dist_list})
sys.stdout.write('distance\n')
sys.stdout.write(str(result))
*/


category('hierarchicalClustering', () => {
  // Single dimension for integer distances

  const tgt1Dist: number[] = [2, 3, 7, 1, 5, 4];
  // const tgt1NewickAverage = '(((2:1.00,1:1.00):1.50,0:2.50):2.83,3:5.33);';

  const tgt2Dist: number[] = [4, 3, 5];
  // const tgt2NewickAverage = '((2:3.00,0:3.00):1.50,1:4.50);';

  const tgt1ClusterMat: ClusterMatrix =
    {
      mergeRow1: new Int32Array([-2, -1, -4]),
      mergeRow2: new Int32Array([-3, 1, 2]),
      heightsResult: new Float32Array([1, 2.5, 5.3333])
    };

  const tgt2ClusterMat: ClusterMatrix =
    {
      mergeRow1: new Int32Array([-1, -2]),
      mergeRow2: new Int32Array([-3, 1]),
      heightsResult: new Float32Array([3, 4.5])
    };

  const AVERAGE_METHOD_CODE = 2;

  test('UI', async () => {
    const csv: string = await _package.files.readAsText('data/demog-short.csv');
    const dataDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    dataDf.name = 'testDemogShort';

    const tv: DG.TableView = grok.shell.addTableView(dataDf);
    await awaitCheck(() => {
      return $(tv.root).find('.d4-grid canvas').length > 0;
    }, 'The view grid canvas not found', 100);

    await hierarchicalClusteringUI(dataDf, ['HEIGHT'], DistanceMetric.Euclidean, 'average');
  });

  // test('hierarchicalClustering1', async () => {
  //   await _testHierarchicalClustering(data1, 'euclidean', 'average', tgt1NewickAverage);
  // });

  // test('hierarchicalClustering2', async () => {
  //   await _testHierarchicalClustering(data2, 'euclidean', 'average', tgt2NewickAverage);
  // });

  // async function _testHierarchicalClustering(
  //   csv: string, distance: string, linkage: string, tgtNewick: string,
  // ): Promise<void> {
  //   const th: ITreeHelper = new TreeHelper();
  //   const dataDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  //   const resTreeRoot: NodeType = await th.hierarchicalClustering(dataDf, distance, linkage);

  //   const tgtTreeRoot = parseNewick(tgtNewick);
  //   tgtTreeRoot.branch_length = 0;

  //   expectObject(resTreeRoot, tgtTreeRoot);
  // }


  // test('hierarchicalClusteringScript', async () => {
  //   const df: DG.DataFrame = DG.DataFrame.fromCsv(data1);
  //   const newick: string = await grok.functions.call('Dendrogram:hierarchicalClusteringScript',
  //     {data: df, distance_name: 'euclidean', linkage_name: 'average'});
  //   let k = 11;
  // });

  // test('hierarchicalClusterinfScript1', async () => {
  //   await _testHierarchicalClusteringScript(data1, 'euclidean', 'average', tgt1NewickAverage);
  // });

  // test('hierarchicalClusterinfScript2', async () => {
  //   await _testHierarchicalClusteringScript(data2, 'euclidean', 'average', tgt2NewickAverage);
  // });

  // async function _testHierarchicalClusteringScript(
  //   csv: string, distance: string, linkage: string, tgtNewick: string,
  // ): Promise<void> {
  //   const dataDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  //   const resNewick: string = await grok.functions.call('Dendrogram:hierarchicalClusteringScript',
  //     {data: dataDf, distance_name: distance, linkage_name: linkage});

  //   expect(resNewick, tgtNewick);
  // }

  async function _testDistanceScript(csv: string, tgtDist: number[], distM: number[][]): Promise<void> {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const t1: number = window.performance.now();
    const distDf: DG.DataFrame = await grok.functions.call('Dendrogram:distanceScript',
      {data: df, distance_name: 'euclidean'});
    const t2: number = window.performance.now();
    _package.logger.debug(`BsV: Tests: _testDistanceScript(), call Dendrogram:distanceScript ET: ${(t2 - t1)} ms`);
    const distCol: DG.Column = distDf.getCol('distance');
    const distData: Float32Array = distCol.getRawData() as Float32Array;
    const dist = new DistanceMatrix(distData, df.rowCount);

    expectArray(dist.data, tgtDist);
    for (let i = 0; i < distM.length; i++) {
      for (let j = 0; j < distM[i].length; j++)
        expect(dist.get(i, j), distM[i][j]);
    }
  }

  test('hierarchicalClusteringWasm1', async () => {
    await _testHierarchicalClusteringWasm(tgt1Dist, AVERAGE_METHOD_CODE, tgt1ClusterMat);
  });

  test('hierarchicalClusteringWasm2', async () => {
    await _testHierarchicalClusteringWasm(tgt2Dist, AVERAGE_METHOD_CODE, tgt2ClusterMat);
  });

  // test('hierarchicalClusteringWasmNoWorker1', async () => {
  //   await _testHierarchicalClusteringWasm(tgt1Dist, AVERAGE_METHOD_CODE, tgt1ClusterMat, false);
  // });

  // test('hierarchicalClusteringWasmNoWorker2', async () => {
  //   await _testHierarchicalClusteringWasm(tgt2Dist, AVERAGE_METHOD_CODE, tgt2ClusterMat, false);
  // });

  async function _testHierarchicalClusteringWasm(
    distA: number[], linkage: number, tgtClusterMatrix: ClusterMatrix,
  ) {
    //calculate number of observations from distance matrix length
    const n = (1 + Math.sqrt(1 + 4 * 2 * distA.length)) / 2;
    const distanceMatrix: DistanceMatrix = new DistanceMatrix(new Float32Array(distA));
    const clusterMatrix: ClusterMatrix = await getClusterMatrixWorker(
      distanceMatrix.data, n, linkage,
    );
    expectObject(clusterMatrix, tgtClusterMatrix);
  }
});
