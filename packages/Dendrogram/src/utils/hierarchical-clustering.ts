import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {injectTreeForGridUI2} from '../viewers/inject-tree-for-grid2';
import {DistanceMetric, isLeaf, LinkageMethod, NodeType} from '@datagrok-libraries/bio/src/trees';
import {parseNewick} from '@datagrok-libraries/bio/src/trees/phylocanvas';
import { DistanceMatrix } from '@datagrok-libraries/bio/src/trees/distance-matrix';
import { TreeHelper } from './tree-helper';
import { ITreeHelper } from '@datagrok-libraries/bio/src/trees/tree-helper';
import { ClusterMatrix } from '../wasm/clusterizerWasm';
import { getClusterMatrixWorker } from '../wasm/clustering-worker-creator';

/** Custom UI form for hierarchical clustering */
export async function hierarchicalClusteringUI2(df: DG.DataFrame): Promise<void> {
  async function exec(dlg: DG.Dialog) {
    const pi = DG.TaskBarProgressIndicator.create('Opening BioStructure* Viewer');
    try {
      // TODO: Get all params from input fields
      // hierarchicalClusteringUI(df, features, distance, linkage);
      throw new Error('Not implemented');
    } catch (err) {
      const errStr = errorToConsole(err);
      grok.shell.error(errStr);
    } finally {
      pi.close();
      dlg.close();
    }
  }

  // TODO: UI to select columns for distance
  const dlg: DG.Dialog = ui.dialog(
    {title: 'Hierarchical clustering'})
    .add(ui.div([])) // TODO: UI
    .addButton('Cancel', () => { dlg.close(); })
    .addButton('Ok', async () => { await exec(dlg); })
    .show();
}

/** Creates table view with injected tree of newick result */
export async function hierarchicalClusteringUI(
  df: DG.DataFrame,
  colNameList: string[],
  distance: DistanceMetric = DistanceMetric.Euclidean,
  linkage: string
): Promise<void> {

  const linkageCode = Object.values(LinkageMethod).findIndex(method => method === linkage);

  const colNameSet: Set<string> = new Set(colNameList);
  const [filteredDf, filteredIndexList]: [DG.DataFrame, Int32Array] =
    hierarchicalClusteringFilterDfForNulls(df, colNameSet);
  const th: ITreeHelper = new TreeHelper();

  let tv: DG.TableView = grok.shell.getTableView(df.name);
  if (filteredDf.rowCount != df.rowCount) {
    grok.shell.warning('Hierarchical clustering analysis on data filtered out for nulls.');
    tv = grok.shell.addTableView(filteredDf);
  }
  // TODO: Filter rows with nulls in selected columns
  const preparedDf = DG.DataFrame.fromColumns(
    filteredDf.columns.toList()
      .filter((col) => colNameSet.has(col.name))
      .map((col) => {
        let res: DG.Column;
        switch (col.type) {
        case DG.COLUMN_TYPE.DATE_TIME:
          // column of type 'datetime' getRawData() returns Float64Array
          const colData: Float64Array = col.getRawData() as Float64Array;
          res = DG.Column.float(col.name, col.length).init((rowI) => {
            return !col.isNone(rowI) ? colData[rowI] : null;
          });
          break;
        default:
          res = col;
        }
        return res;
      }));

  const distanceMatrix = await th.calcDistanceMatrix(preparedDf,
    preparedDf.columns.toList().map(col => col.name),
    distance);

  const clusterMatrixWorker = getClusterMatrixWorker({distMatArray: distanceMatrix!.data, n: preparedDf.rowCount, methodCode: linkageCode});
  const clusterMatrix = await clusterMatrixWorker;

  // const hcPromise = hierarchicalClusteringByDistanceExec(distanceMatrix!, linkage);
  // Replace rows indexes with filtered
  // newickStr returned with row indexes after filtering, so we need reversed dict { [fltIdx: number]: number}
  const fltRowIndexes: { [fltIdx: number]: number } = {};
  const fltRowCount: number = filteredDf.rowCount;
  for (let fltRowIdx: number = 0; fltRowIdx < fltRowCount; fltRowIdx++)
    fltRowIndexes[fltRowIdx] = filteredIndexList[fltRowIdx];

  const newickRoot: NodeType = parseClusterMatrix(clusterMatrix);
  // Fix branch_length for root node as required for hierarchical clustering result
  newickRoot.branch_length = 0;
  (function replaceNodeName(node: NodeType, fltRowIndexes: { [fltIdx: number]: number }) {
    const nodeFilteredIdx: number = parseInt(node.name);
    const nodeIdx: number = fltRowIndexes[nodeFilteredIdx];
    if (!isLeaf(node)) {
      for (const childNode of node.children!)
        replaceNodeName(childNode, fltRowIndexes);
    }
  })(newickRoot, fltRowIndexes);

  // empty clusterDf to stub injectTreeForGridUI2
  const clusterDf = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'cluster', [])]);
  injectTreeForGridUI2(tv.grid, newickRoot, undefined, 300);
}

export function hierarchicalClusteringFilterDfForNulls(
  df: DG.DataFrame, colNameSet: Set<string>
): [DG.DataFrame, Int32Array] {
  // filteredNullsDf to open new table view
  const colList: DG.Column[] = df.columns.toList().filter((col) => colNameSet.has(col.name));
  const filter: DG.BitSet = DG.BitSet.create(df.rowCount, (rowI: number) => {
    // TODO: Check nulls in columns of colNameList
    return colList.every((col) => !col.isNone(rowI));
  });
  const filteredDf: DG.DataFrame = df.clone(filter);
  const filteredIndexList: Int32Array = filter.getSelectedIndexes();
  return [filteredDf, filteredIndexList];
}

/** Runs script and returns newick result
 * @param {DG.DataFrame} preparedDf Data frame with features columns only
 */
export async function hierarchicalClusteringExec(
  preparedDf: DG.DataFrame, distance: string, linkage: string
): Promise<string> {
  // const newick: string = await hierarchicalClustering(selectedColumnsDf);
  // getNewick - is script scripts/hierarchicalClustering.py
  let newick = '';
  try {
    // Dendrogram:hierarchicalClusteringScript is the script at scripts/hierarchicalClustering.py
    newick = await grok.functions.call('Dendrogram:hierarchicalClusteringScript', {
      data: preparedDf,
      distance_name: distance,
      linkage_name: linkage,
    });
  } catch (err) {
    const errStr = errorToConsole(err);
    console.error('Dendrogram: hierarchicalClusteringExec() ' + `${errStr}`);
    throw err;
  }

  return newick;

  // // TODO: stub
  // return new Promise<string>((resolve, reject) => {
  //   const res: string = df.getTag(thTAGS.DF_NEWICK)!;
  //   resolve(res);
  // });
}

async function hierarchicalClusteringByDistanceExec(distance: DistanceMatrix, linkage: string): Promise<string> {
  const distanceCol: DG.Column = DG.Column.fromFloat32Array('distance', distance.data);
  const dataDf: DG.DataFrame = DG.DataFrame.fromColumns([distanceCol]);

  const newickStr: string = await grok.functions.call(
      'Dendrogram:hierarchicalClusteringByDistanceScript',
      {data: dataDf, size: distance.size, linkage_name: linkage});

  return newickStr;
}

function parseClusterMatrix(clusterMatrix:ClusterMatrix){
  /*
  clusert matrix is in R format, I.E. the indexings are 1-based.
  one of the reasons is that values in merge arrays are not always positive. if the value is negative
  it means that we are referencing a leaf node. otherwise we are referencing a cluster node.
  for example :
  1, -2, 0.1 would mean that the merge happened between cluster 0 and leaf node 1 with a distance of 0.1
  */

  function getSubTreeLength(node: NodeType) : number{
    if(isLeaf(node)){
      return (node.branch_length ?? 0);
    } else {
      return (node.branch_length ?? 0) + getSubTreeLength(node.children![0]);
    }
  }

  const clusters: NodeType[] = [];
  const {mergeRow1, mergeRow2, heightsResult} = clusterMatrix;
  for(let i = 0; i<heightsResult.length; i++){
    let left: NodeType, right: NodeType;
    if(mergeRow1[i] < 0){
      left = {name: (mergeRow1[i] * -1 - 1).toString(), branch_length: heightsResult[i]};
    } else {
      left = clusters[mergeRow1[i] - 1];
    }
    if(mergeRow2[i] < 0){
      right = {name: (mergeRow2[i] * -1 - 1).toString(), branch_length: heightsResult[i]};
    } else {
      right = clusters[mergeRow2[i] - 1];
    }

    const leftLength = getSubTreeLength(left);
    const rightLength = getSubTreeLength(right);

    if (leftLength > rightLength) {
      right.branch_length = leftLength - rightLength;
    }
    if(rightLength > leftLength){
      left.branch_length = rightLength - leftLength;
    }
    clusters.push({name: '', children: [left, right], branch_length:0});
  }
  return clusters[clusters.length - 1];
}
