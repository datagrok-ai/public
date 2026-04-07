import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import {DistanceMetric, isLeaf, NodeCuttedType, NodeType} from '@datagrok-libraries/bio/src/trees';
import {NO_NAME_ROOT, parseNewick} from '@datagrok-libraries/bio/src/trees/phylocanvas';
import {NEWICK_EMPTY} from '@datagrok-libraries/bio/src/trees/consts';
import {DistanceMatrix, DistanceMatrixService} from '@datagrok-libraries/ml/src/distance-matrix';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {ClusterMatrix} from '@datagrok-libraries/bio/src/trees';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {NumberMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {IntArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics/consts';
import { _package } from '../package';
import { getSeqHelper } from '@datagrok-libraries/bio/src/utils/seq-helper';

type DataNodeDict = { [nodeName: string]: number };

export const enum TAGS {
  DF_NEWICK = '.newick',
  DF_NEWICK_LEAF_COL_NAME = '.newickLeafColumn',
}

/** Implements ITreeHelper from bio lib
 * Common expected behavior: parseNewick and newickToDf returns name '[root]' for no-root-name newick string.
 * Because root node can have branch_length and root node name and parent node name for these row should be different.
 */
export class TreeHelper implements ITreeHelper {
  newickToDf(newick: string, dfName?: string, nodePrefix?: string): DG.DataFrame {
    const nodePrefixV: string = nodePrefix ?? '';

    let i = 0;

    const obj = parseNewick(newick);
    if (obj.name == NO_NAME_ROOT) obj.name = `${nodePrefixV}${NO_NAME_ROOT}`;

    const nodes: string[] = [];
    const parents: (string | null)[] = [];
    const leafs: boolean[] = [];
    const distances: (number | null)[] = [];
    const annotations: any[] = [];

    function traverse(obj: NodeType, parent: NodeType | null) {
      if (obj === null || typeof obj != 'object') return;

      const isRoot: boolean = parent == null;

      let name: string = obj.name;
      if (!name) {
        name = obj.name = `${nodePrefixV}node-${i}`;
        ++i;
      } else if (isRoot && obj.name == NO_NAME_ROOT) {
        name = `${nodePrefixV}${NO_NAME_ROOT}`;
      }

      if (!!obj.children) {
        const childrenCount = obj.children.length;
        for (let i = 0; i < childrenCount; i++)
          traverse(obj.children[i], obj);
      }

      nodes.push(name);
      distances.push(obj.branch_length ? obj.branch_length : null);
      // annotations.push(obj.annotation);
      parents.push(parent ? parent!.name : null);
      leafs.push(!obj.children || obj.children.length == 0);
    }

    traverse(obj, null);

    const nodeCol: DG.Column = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'node', nodes);
    const parentCol: DG.Column = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'parent', parents);
    const leafCol: DG.Column = DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'leaf', leafs);
    // Preventing semType detectors to interpret data
    nodeCol.semType = 'id';
    parentCol.semType = 'id';
    const columns = [nodeCol, parentCol, leafCol];

    if (distances.some((d) => d !== null))
      columns.push(DG.Column.fromList('double', 'distance', distances));

    if (annotations.some((a) => !!a))
      columns.push(DG.Column.fromList('string', 'annotation', annotations));

    const df = DG.DataFrame.fromColumns(columns);

    if (dfName)
      df.name = `df-${dfName}`;

    df.setTag('.newick', newick);
    df.setTag('.newickJson', JSON.stringify(obj));

    return df;
  }

  toNewick(node: NodeType | null): string {
    function toNewickInt(node: NodeType): string {
      const isLeaf = !node.children || node.children.length == 0;

      if (isLeaf) {
        return ([] as string[]).concat(
          node.name,
          node.hasOwnProperty('branch_length') ? `:${node.branch_length}` : []
        ).join('');
      } else {
        const childrenText = node.children!.map((childNode) => toNewickInt(childNode)).join(',');
        return ([] as string[]).concat(
          `(${childrenText})`,
          node.name,
          node.hasOwnProperty('branch_length') ? `:${node.branch_length}` : []
        ).join('');
      }
    }

    return !node ? NEWICK_EMPTY : `${toNewickInt(node)};`;
  }

  getLeafList<TNode extends NodeType>(node: TNode | null, list?: TNode[]): TNode[] {
    if (!node) return [];
    if (!list) list = [];

    if (isLeaf(node)) {
      list.push(node);
      return [node]; // node is a leaf
    } else {
      for (const child of node.children ?? [])
        this.getLeafList(child as TNode, list);
      //(node.children ?? []).forEach((child) => this.getLeafList(child as TNode, list));
      return list;
    }
  }

  getNodeList<TNode extends NodeType>(node: TNode | null, list?: TNode[]): TNode[] {
    if (!node) return [];
    if (!list) list = [];
    if (isLeaf(node)) {
      list.push(node);
      return [node]; // node is a leaf
    } else {
      // const childNodeListList: TNode[][] = node.children!
      //   .map((child) => { return this.getNodeList(child as TNode); });
      // const childNodeListList: TNode[][] = [];
      for (const child of node.children!)
        this.getNodeList(child as TNode, list);
      list.push(node);
      return list;
    }
  }

  includes<TNode extends NodeType>(node: TNode, sub: TNode): boolean {
    if (node == sub) return true;

    let res: boolean = false;
    for (const childNode of (node.children ?? [])) {
      if (this.includes(childNode, sub)) {
        res = true;
        break;
      }
    }
    return res;
  }

  filterTreeByLeaves(node: NodeType | null, leaves: { [name: string]: any }): NodeType | null {
    // copy node because phylocanvas.gl changes data structure completely
    if (!node) return null;
    const resNode = Object.assign({}, node); // shallow copy

    if (isLeaf(resNode)) {
      return resNode.name in leaves ? resNode : null;
    } else {
      resNode.children = [];
      for (const child of node.children!) {
        const resChild = this.filterTreeByLeaves(child, leaves);
        if (resChild) resNode.children.push(resChild);
      }
      // .map((child) => this.filterTreeByLeaves(child, leaves))
      // .filter((child) => child != null) as NodeType[];

      return resNode.children.length > 0 ? resNode : null;
    }
  }

  getNodesByLeaves<TNode extends NodeType>(node: TNode | null, leaves: { [name: string]: any }): TNode[] {
    if (!node) return [];

    if (isLeaf(node)) {
      return node.name in leaves ? [node] : [];
    } else {
      const children: TNode[] = node.children as TNode[] ?? [];
      const childrenRes: TNode[] = [];
      for (const child of children) {
        const childRes = this.getNodesByLeaves(child, leaves);
        childrenRes.push(...childRes);
      }
      //  children.map((child) => {
      //   return this.getNodesByLeaves(child, leaves);
      // }).flat();

      return (childrenRes.length == children.length &&
        children.every((child, i) => child == childrenRes[i])) ? [node] : childrenRes;
    }
  }


  /** Cuts tree, gets list of clusters as lists of leafs
   * @param {NodeType}node - root node of tree
   * @param {number}cutHeight - height of cut
   * @param {number}currentHeight - current height of node
   * @return {NodeType[]} - list of clusters as lists of leafs
   */
  treeCutAsLeaves(node: NodeType, cutHeight: number, currentHeight: number = 0): NodeType[] {
    let res: NodeType[];

    const nodeBranchLength = node.branch_length ?? 0;
    if ((currentHeight + nodeBranchLength) < cutHeight) {
      /* if node has no children, then an empty list will be returned */
      res = ([] as NodeType[]).concat(
        ...(node.children ?? [])
          .map((child) => this.treeCutAsLeaves(child, cutHeight, currentHeight + nodeBranchLength)));
    } else {
      res = [node];
    }

    return res;
  }

  /** Cuts tree, gets cutted tree of clusters as lists of leafs
   * @param {NodeType}node - root node of tree
   * @param {number}cutHeight - height of cut
   * @param {boolean}keepShorts - keep short branches
   * @param {number}currentHeight - current height of node
   * @return {NodeType} - cutted tree of clusters as lists of leafs*/
  treeCutAsTree(
    node: NodeType, cutHeight: number, keepShorts?: boolean, currentHeight?: number
  ): NodeType | null {
    const nodeBranchHeight = node.branch_length ?? 0;
    const currentHeightV: number = currentHeight ?? 0;
    const keepShortsV: boolean = keepShorts ?? false;
    if ((currentHeightV + nodeBranchHeight) < cutHeight) {
      if (!isLeaf(node)) {
        const res: NodeCuttedType = Object.assign({}, node) as NodeCuttedType;
        res.children = node.children!
          .map((child) => {
            return this.treeCutAsTree(child, cutHeight, keepShortsV, currentHeightV + nodeBranchHeight);
          })
          .filter((n) => n != null) as NodeType[];
        return res;
      } else {
        return keepShortsV ? node as NodeCuttedType : null;
      }
    } else {
      const res: NodeCuttedType = Object.assign({}, node) as NodeCuttedType;
      res.branch_length = cutHeight - currentHeightV; // shorten branch_length of the node to cut_length remains
      res.children = [];
      res.cuttedChildren = [{
        name: `${node.name}.cutted`,
        branch_length: nodeBranchHeight + currentHeightV - cutHeight,
        children: node.children ? [...node.children] : [],
      }];
      res.cuttedLeafNameList = this.getLeafList(node).map((l) => l.name);
      return res;
    }
  }

  /** Sets grid's row order and returns tree (root node) of nodes presented in data
   * @param {NodeType}tree - root node of tree
   * @param {DG.Grid}grid - grid to set order
   * @param {string}leafColName - name of column with leaf names
   * @param {boolean}removeMissingDataRows - remove rows with missing data
   * @return {[NodeType, string[]]} - tree (root node) of nodes presented in data and list of missed data nodes
   */
  setGridOrder(
    tree: NodeType | null, grid: DG.Grid, leafColName?: string,
    removeMissingDataRows: boolean = false
  ): [NodeType, string[]] {
    console.debug('Dendrogram.setGridOrder() start');

    const dataDf: DG.DataFrame = grid.dataFrame;
    if (!dataDf)
      throw new Error('DataFrame wasn\'t found');

    const missedDataNodeList: string[] = [];
    const missedTreeLeafList: string[] = [];

    // Build TREE node name dictionary
    const treeLeafList = this.getLeafList(tree); // ordered by tree
    const treeLeafDict: { [nodeName: string]: NodeType } = {};
    for (const treeLeaf of treeLeafList) {
      if (treeLeaf.name in treeLeafDict)
        throw new Error('Non unique key tree leaf name');

      treeLeafDict[treeLeaf.name] = treeLeaf;
    }

    const dataNodeCol: DG.Column = leafColName ? dataDf.getCol(leafColName) :
      DG.Column.fromList(DG.COLUMN_TYPE.INT, '<index>', wu.count(0).take(dataDf.rowCount).toArray());

    if (removeMissingDataRows) {
      for (let dataRowI = dataNodeCol.length - 1; dataRowI >= 0; dataRowI--) {
        const dataNodeName: string = dataNodeCol.get(dataRowI);
        if (!(dataNodeName in treeLeafDict)) {
          missedDataNodeList.push(dataNodeName);
          dataDf.rows.removeAt(dataRowI, 1, false);
        }
      }
    }
    // Build DATA node name dictionary
    const dataNodeDict: { [nodeName: string]: number } = {};
    for (let dataRowI = 0; dataRowI < dataNodeCol.length; dataRowI++) {
      const dataNodeName: string = dataNodeCol.get(dataRowI);

      const dataNodeOldFilter: boolean = dataDf.filter.get(dataRowI);
      const dataNodeFilter: boolean = (dataNodeName in treeLeafDict) && dataNodeOldFilter;
      if (dataNodeFilter) {
        if (dataNodeName in dataNodeDict)
          throw new Error('Non unique key data node name');

        dataNodeDict[dataNodeName] = dataRowI;
      } else { // !dataNodeFilter
        missedTreeLeafList.push(dataNodeName);

        // Hide grid data rows with node name not in tree leaf list
        if (dataNodeOldFilter)
          dataDf.filter.set(dataRowI, false, false);
      }
    }

    // TODO: Fire filter event (?)

    const resTree: NodeType = this.filterTreeByLeaves(tree, dataNodeDict)!;
    const resTreeLeafList = this.getLeafList(resTree);

    const order: number[] = new Array<number>(resTreeLeafList.length); // rowCount filtered for leaves
    for (let leafI = 0; leafI < resTreeLeafList.length; ++leafI) {
      const leafNodeName = resTreeLeafList[leafI].name;
      if (leafNodeName in dataNodeDict)
        order[leafI] = dataNodeDict[leafNodeName];
      else
        missedDataNodeList.push(leafNodeName);
    }

    grid.setRowOrder(order);

    // // missedWarnMessage
    // if (missedTreeLeafList.length > 0 || missedDataNodeList.length > 0) {
    //   const missedWarnMsg = Array<string>().concat(
    //     missedDataNodeList.length > 0 ?
    //       `Missed data nodes in tree: ${missedDataNodeList.map((n) => `'${n}'`).join(', ')}.` : [],
    //     missedTreeLeafList.length > 0 ?
    //       `Missed tree leaves in data: ${missedTreeLeafList.map((n) => `'${n}'`).join(', ')}. ` : [],
    //   ).join(' ');
    //   console.warn(missedWarnMsg);
    // }

    console.debug('Dendrogram.setGridOrder() ' + `resTreeLeafList.length = ${resTreeLeafList.length}`);

    console.debug('Dendrogram.setGridOrder() end');
    return [
      resTree,
      ([] as string[]).concat(missedDataNodeList, missedTreeLeafList),
    ];
  }

  markClusters(
    tree: NodeCuttedType, dataDf: DG.DataFrame, leafColName: string | null, clusterColName: string, na?: any
  ): void {
    const naValue = na ?? null;
    const clusterCol: DG.Column = dataDf.getCol(clusterColName);
    clusterCol.init((_rowI) => naValue);

    const dataNodeCol: DG.Column = leafColName ? dataDf.getCol(leafColName) :
      DG.Column.fromList(DG.COLUMN_TYPE.INT, '<index>', wu.count(0).take(dataDf.rowCount).toArray());

    const dataNodeDict: DataNodeDict = {};
    for (let dataRowI: number = 0; dataRowI < dataDf.rowCount; dataRowI++) {
      const dataNodeName = dataNodeCol.get(dataRowI);
      dataNodeDict[dataNodeName] = dataRowI;
    }

    const clusterList = this.getLeafList(tree);
    for (let clusterI: number = 1; clusterI < clusterList.length + 1; clusterI++) {
      const cluster: NodeCuttedType = clusterList[clusterI - 1] as NodeCuttedType;
      for (const leafName of cluster.cuttedLeafNameList ?? []) {
        const dataRowI = dataNodeDict[leafName];
        clusterCol.set(dataRowI, clusterI, false);
      }
    }
  }

  buildClusters(tree: NodeCuttedType, clusterDf: DG.DataFrame, _clusterColName: string, _leafColName?: string): void {
    for (let clusterRowI = clusterDf.rowCount - 1; clusterRowI >= 0; clusterRowI--)
      clusterDf.rows.removeAt(clusterRowI);

    const clusterList = this.getLeafList(tree);
    for (let clusterI: number = 1; clusterI < clusterList.length + 1; clusterI++) {
      const clusterNode = clusterList[clusterI - 1] as NodeCuttedType;
      const leafNameList: string[] = clusterNode.cuttedLeafNameList ?? [];
      const clusterLeavesStr = leafNameList.join(', ');
      clusterDf.rows.addNew([clusterI, clusterLeavesStr, leafNameList.length], false);
    }
  }

  /** Cuts tree at threshold (from root), returns array of subtrees, marks leafs for cluster in data
   * @param {NodeType}node - root node of tree
   * @param {number}cutHeight - threshold for cutting
   * @param {DG.DataFrame}dataDf - Dataframe
   * @param {number}leafColName - name of column with lef names
   * @param {string}clusterColName - name of column with cluster names
   * @param {any}na - Value of nulls*/
  cutTreeToGrid(
    node: NodeType, cutHeight: number, dataDf: DG.DataFrame,
    leafColName: string, clusterColName: string, na?: any
  ): void {
    const clusterList: NodeType[] = this.treeCutAsLeaves(node, cutHeight, 0);

    const naValue: any = na ?? null;
    const clusterCol: DG.Column = dataDf.getCol(clusterColName);
    // for (let rowI = 0; rowI < clusterCol.length; rowI++)
    //   clusterCol.set(rowI, na_value);
    clusterCol.init((_rowI) => naValue);

    /* A leaf with cumulative height less than threshold
       will not be included in nor marked as cluster */

    const dataNodeDict: DataNodeDict = {};
    for (let dataRowI: number = 0; dataRowI < dataDf.rowCount; dataRowI++) {
      const dataNodeName = dataDf.get(leafColName, dataRowI);
      dataNodeDict[dataNodeName] = dataRowI;
    }

    function markCluster(node: NodeType, cluster: number) {
      if (isLeaf(node)) {
        const nodeName = node.name;
        const dataRowI: number = dataNodeDict[nodeName];
        dataDf.set(clusterColName, dataRowI, cluster);
      } else {
        for (const childNode of node.children!)
          markCluster(childNode, cluster);
      }
    }

    // const leafCol: DG.Column = dataDf.getCol(leafColName);

    for (let clusterI = 1; clusterI <= clusterList.length; clusterI++)
      markCluster(clusterList[clusterI - 1], clusterI);
  }

  generateTree(size: number): NodeType {
    function placeNode(currentNode: NodeType, newNode: NodeType): void {
      if (currentNode.children!.length < 2) {
        currentNode.children!.push(newNode);
      } else {
        const rnd: number = Math.random();
        const tgtNodeI = Math.floor(rnd / (1 / currentNode.children!.length));
        const tgtNode = currentNode.children![tgtNodeI];
        placeNode(tgtNode, newNode);
      }
    }

    const root = {
      name: `node-0`,
      branch_length: Math.random(),
      children: [],
    };

    for (let nodeI = 1; nodeI < size; nodeI++) {
      const newNode = {
        name: `node-${nodeI}`,
        branch_length: Math.random(),
        children: [],
      };
      placeNode(root, newNode);
    }

    return root;
  }

  async hierarchicalClustering(data: DG.DataFrame, distance: string, linkage: string): Promise<NodeType> {
    const newickStr: string = await grok.functions.call(
      'Dendrogram:hierarchicalClusteringScript',
      {data: data, distance_name: distance, linkage_name: linkage});
    const treeRoot: NodeType = parseNewick(newickStr);
    if (!treeRoot.branch_length) treeRoot.branch_length = 0;
    return treeRoot;
  }

  async hierarchicalClusteringByDistance(distance: DistanceMatrix, linkage: string): Promise<NodeType> {
    const distanceCol: DG.Column = DG.Column.fromFloat32Array('distance', distance.data);
    const dataDf: DG.DataFrame = DG.DataFrame.fromColumns([distanceCol]);
    const newickStr: string = await grok.functions.call(
      'Dendrogram:hierarchicalClusteringByDistanceScript',
      {data: dataDf, size: distance.size, linkage_name: linkage});
    const treeRoot: NodeType = parseNewick(newickStr);
    if (!treeRoot.branch_length) treeRoot.branch_length = 0;
    return treeRoot;
  }

  async encodeSequences(seqs: DG.Column): Promise<string[]> {
    const seqHelper = await getSeqHelper();
    const ncSh = seqHelper.getSeqHandler(seqs);
    const seqList = seqs.toList();
    const seqColLength = seqList.length;
    let charCodeCounter = 36;
    const charCodeMap = new Map<string, string>();
    for (let rowIdx = 0; rowIdx < seqColLength; rowIdx++) {
      const seq = seqList[rowIdx];
      if (seq == null || seqs.isNone(rowIdx)) {
        seqList[rowIdx] = null;
        continue;
      }
      seqList[rowIdx] = '';
      const ss = ncSh.getSplitted(rowIdx);
      for (let j = 0; j < ss.length; j++) {
        const char = ss.getCanonical(j);
        if (!charCodeMap.has(char)) {
          charCodeMap.set(char, String.fromCharCode(charCodeCounter));
          charCodeCounter++;
        }
        seqList[rowIdx] += charCodeMap.get(char)!;
      }
    }
    return seqList;
  }

  async calcDistanceMatrix(
    df: DG.DataFrame, colNames: string[], method: DistanceMetric = DistanceMetric.Euclidean
  ) {
    // Output distance matrix. reusing it saves a lot of memory
    let out: DistanceMatrix | null = null;

    const distanceMatrixService = new DistanceMatrixService(true, false);
    const columns = colNames.map((name) => df.getCol(name));
    for (const col of columns) {
      let values: Float32Array;
      if (col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) {
        values = await distanceMatrixService.calc(col.getRawData(), NumberMetricsNames.Difference, false);
      } else if (col.semType === DG.SEMTYPE.MACROMOLECULE) {
        // Use Hamming distance when sequences are aligned
        const seqDistanceFunction: MmDistanceFunctionsNames = MmDistanceFunctionsNames.LEVENSHTEIN;
        const encodedSeqs = await this.encodeSequences(col);
        values = await distanceMatrixService.calc(encodedSeqs, seqDistanceFunction, false);
      } else if (col.semType === DG.SEMTYPE.MOLECULE) {
        const fingerPrintCol: DG.Column<DG.BitSet | null> =
          await grok.functions.call('Chem:getMorganFingerprints', {molColumn: col});
        const fingerPrintBitArrayCol = fingerPrintCol.toList().map((bs: DG.BitSet | null) =>
          bs ? bs.getBuffer() : null);
        values = await distanceMatrixService.calc(fingerPrintBitArrayCol, IntArrayMetricsNames.TanimotoIntArray, false);
      } else { throw new TypeError('Unsupported column type'); }

      if (!out) {
        out = new DistanceMatrix(values);
        if (columns.length > 1) {
          out.normalize();
          if (method === DistanceMetric.Euclidean)
            out.square();
        }
      } else {
        let newMat: DistanceMatrix | null = new DistanceMatrix(values);
        newMat.normalize();
        switch (method) {
          case DistanceMetric.Manhattan: {
            out.add(newMat);
            break;
          }
          default:
            newMat.square();
            out.add(newMat);
        }
        // remove reference
        newMat = null;
      }
    }
    distanceMatrixService.terminate();
    if (method === DistanceMetric.Euclidean && columns.length > 1)
      out?.sqrt();

    return out;
  }

  parseClusterMatrix(clusterMatrix: ClusterMatrix): NodeType {
    /*
    clusert matrix is in R format, I.E. the indexings are 1-based.
    one of the reasons is that values in merge arrays are not always positive. if the value is negative
    it means that we are referencing a leaf node. otherwise we are referencing a cluster node.
    for example :
    1, -2, 0.1 would mean that the merge happened between cluster 0 and leaf node 1 with a distance of 0.1
    */

    function getSubTreeLength(node: NodeType): number {
      //Tradeoff between performance and memory usage - in this case better to use less memory.
      //as wasm already takes up quite a bit
      function subTreeLength(children?: NodeType[]): number {
        return children && children.length ? (children[0].branch_length ?? 0) + subTreeLength(children[0].children) : 0;
      }

      if (isLeaf(node))
        return 0;
      else
        return subTreeLength(node.children);
    }


    const {mergeRow1, mergeRow2, heightsResult} = clusterMatrix;
    const clusters: NodeType[] = new Array<NodeType>(heightsResult.length);

    for (let i = 0; i < heightsResult.length; i++) {
      const left: NodeType = mergeRow1[i] < 0 ?
        {name: (mergeRow1[i] * -1 - 1).toString(), branch_length: heightsResult[i]} :
        clusters[mergeRow1[i] - 1];

      const right = mergeRow2[i] < 0 ?
        {name: (mergeRow2[i] * -1 - 1).toString(), branch_length: heightsResult[i]} :
        clusters[mergeRow2[i] - 1];

      const leftLength = getSubTreeLength(left);
      const rightLength = getSubTreeLength(right);

      left.branch_length = heightsResult[i] - leftLength;
      right.branch_length = heightsResult[i] - rightLength;

      clusters[i] = {name: '', children: [left, right], branch_length: 0};
    }
    return clusters[clusters.length - 1];
  }
}
