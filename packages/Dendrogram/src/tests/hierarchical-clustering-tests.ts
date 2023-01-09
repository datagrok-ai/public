import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject, expectArray} from '@datagrok-libraries/utils/src/test';
import {hierarchicalClusteringUI} from '../utils/hierarchical-clustering';
import {_package} from '../package-test';
import {viewsTests} from './utils/views-tests';
import {DistanceMetric} from '@datagrok-libraries/utils/src/type-declarations';
import {DistanceMatrix, ITreeHelper, NodeType, parseNewick} from '@datagrok-libraries/bio';
import {TreeHelper} from '../utils/tree-helper';

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


category('hierarchicalClustering', viewsTests((ctx: { dfList: DG.DataFrame[], vList: DG.ViewBase[] }) => {
  // Single dimension for integer distances
  const data1: string = `x
8
6
5
1`;
  const tgt1Dist: number[] = [2, 3, 7, 1, 5, 4];
  const tgt1DistM: number[][] = [
    [0, 2, 3, 7],
    [2, 0, 1, 5],
    [3, 1, 0, 4],
    [7, 5, 4, 0]];
  const tgt1NewickAverage = '(((2:1.00,1:1.00):1.50,0:2.50):2.83,3:5.33);';

  const data2 = `x,y
0,0
4,0
0,3`;
  const tgt2Dist: number[] = [4, 3, 5];
  const tgt2DistM: number[][] = [
    [0, 4, 3],
    [4, 0, 5],
    [3, 5, 0]];
  const tgt2NewickAverage = '((2:3.00,0:3.00):1.50,1:4.50);';

  test('UI', async () => {
    const csv: string = await _package.files.readAsText('data/demog-short.csv');
    const dataDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    dataDf.name = 'testDemogShort';

    const tv: DG.TableView = grok.shell.addTableView(dataDf);

    ctx.vList.push(tv);
    ctx.dfList.push(dataDf);

    hierarchicalClusteringUI(dataDf, ['HEIGHT'], 'euclidean', 'average');
  });

  test('hierarchicalClustering1', async () => {
    await _testHierarchicalClustering(data1, 'euclidean', 'average', tgt1NewickAverage);
  });

  test('hierarchicalClustering2', async () => {
    await _testHierarchicalClustering(data2, 'euclidean', 'average', tgt2NewickAverage);
  });

  async function _testHierarchicalClustering(
    csv: string, distance: string, linkage: string, tgtNewick: string
  ): Promise<void> {
    const th: ITreeHelper = new TreeHelper();
    const dataDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const resTreeRoot: NodeType = await th.hierarchicalClustering(dataDf, distance, linkage);

    const tgtTreeRoot = parseNewick(tgtNewick);
    tgtTreeRoot.branch_length = 0;

    expectObject(resTreeRoot, tgtTreeRoot);
  }


  // test('hierarchicalClusteringScript', async () => {
  //   const df: DG.DataFrame = DG.DataFrame.fromCsv(data1);
  //   const newick: string = await grok.functions.call('Dendrogram:hierarchicalClusteringScript',
  //     {data: df, distance_name: 'euclidean', linkage_name: 'average'});
  //   let k = 11;
  // });

  test('hierarchicalClusterinfScript1', async () => {
    await _testHierarchicalClusteringScript(data1, 'euclidean', 'average', tgt1NewickAverage);
  });

  test('hierarchicalClusterinfScript2', async () => {
    await _testHierarchicalClusteringScript(data2, 'euclidean', 'average', tgt2NewickAverage);
  });

  async function _testHierarchicalClusteringScript(
    csv: string, distance: string, linkage: string, tgtNewick: string
  ): Promise<void> {
    const dataDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const resNewick: string = await grok.functions.call('Dendrogram:hierarchicalClusteringScript',
      {data: dataDf, distance_name: distance, linkage_name: linkage});

    expect(resNewick, tgtNewick);
  }


  test('distanceScript1', async () => {
    await _testDistanceScript(data1, tgt1Dist, tgt1DistM);
  });

  test('distanceScript2', async () => {
    await _testDistanceScript(data2, tgt2Dist, tgt2DistM);
  });

  async function _testDistanceScript(csv: string, tgtDist: number[], distM: number[][]): Promise<void> {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const distDf: DG.DataFrame = await grok.functions.call('Dendrogram:distanceScript',
      {data: df, distance_name: 'euclidean'});
    const distCol: DG.Column = distDf.getCol('distance');
    const distData: Float32Array = distCol.getRawData() as Float32Array;
    const dist = new DistanceMatrix(distData, df.rowCount);

    expectArray(dist.data, tgtDist);
    for (let i = 0; i < distM.length; i++) {
      for (let j = 0; j < distM[i].length; j++)
        expect(dist.get(i, j), distM[i][j]);
    }
  }


  test('hierarchicalClusteringByDistanceScript1', async () => {
    await _testHierarchicalClusteringByDistanceScript(tgt1Dist, 'average', tgt1NewickAverage);
  });

  test('hierarchicalClusteringByDistanceScript2', async () => {
    await _testHierarchicalClusteringByDistanceScript(tgt2Dist, 'average', tgt2NewickAverage);
  });

  async function _testHierarchicalClusteringByDistanceScript(
    distA: number[], linkage: string, tgtNewick: string
  ): Promise<void> {
    const distance: DistanceMatrix = new DistanceMatrix(new Float32Array(distA));
    const distanceCol: DG.Column = DG.Column.fromFloat32Array('distance', distance.data);
    const dataDf: DG.DataFrame = DG.DataFrame.fromColumns([distanceCol]);
    const resNewick: string = await grok.functions.call(
      'Dendrogram:hierarchicalClusteringByDistanceScript',
      {data: dataDf, size: distance.size, linkage_name: linkage});

    expect(resNewick, tgtNewick);
  }


  test('hierarchicalClusteringByDistance1', async () => {
    await _testHierarchicalClusteringByDistance(tgt1Dist, 'average', tgt1NewickAverage);
  });

  test('hierarchicalClusteringByDistance2', async () => {
    await _testHierarchicalClusteringByDistance(tgt2Dist, 'average', tgt2NewickAverage);
  });

  async function _testHierarchicalClusteringByDistance(
    distA: number[], linkage: string, tgtNewick: string
  ): Promise<void> {
    const th: ITreeHelper = new TreeHelper();
    const distance: DistanceMatrix = new DistanceMatrix(new Float32Array(distA));
    const resTreeRoot: NodeType = await th.hierarchicalClusteringByDistance(distance, linkage);

    const tgtTreeRoot: NodeType = parseNewick(tgtNewick);
    tgtTreeRoot.branch_length = 0;

    expectObject(resTreeRoot, tgtTreeRoot);
  }
}));
