import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import {
  DistanceMatrix,
  isLeaf,
  ITreeHelper,
  NodeCuttedType,
  NodeType,
  parseNewick,
} from '@datagrok-libraries/bio';

type TreeLeafDict = { [nodeName: string]: NodeType };
type DataNodeDict = { [nodeName: string]: number };
type NodeNameCallback = (nodeName: string) => void;

export const enum TAGS {
  DF_NEWICK = '.newick',
  DF_NEWICK_LEAF_COL_NAME = '.newickLeafColumn',
}

/** Implements ITreeHelper from bio lib */
export class TreeHelper implements ITreeHelper {
  newickToDf(
    newick: string, dfName?: string, nodePrefix?: string, skipEmptyParentRoot?: boolean
  ): DG.DataFrame {
    const nodePrefixV: string = nodePrefix ?? '';
    const skipEmptyParentRootV: boolean = skipEmptyParentRoot ?? false;

    let parent: string | null = null;
    let i = 0;

    const obj = parseNewick(newick);

    const nodes: string[] = [];
    const parents: string[] = [];
    const leafs: boolean[] = [];
    const distances: (number | null)[] = [];
    const annotations: any[] = [];

    function traverse(obj: NodeType) {
      if (obj === null || typeof obj != 'object') return;

      const isRoot: boolean = obj.name == 'root';

      let name: string = obj.name;
      if (!name) {
        name = obj.name = `${nodePrefixV}node-${i}`;
        ++i;
      } else if (isRoot) {
        name = `${nodePrefixV}root`;
      }

      if (!isRoot || !skipEmptyParentRootV) {
        nodes.push(name);
        distances.push(obj.branch_length ? obj.branch_length : null);
        // annotations.push(obj.annotation);
        parents.push(parent!);
        leafs.push(!obj.children || obj.children.length == 0);
      }

      if (!obj.children) return;
      const childrenNum = obj.children.length;
      const prevParent = parent;
      parent = name;

      for (let i = 0; i < childrenNum; i++) {
        traverse(obj.children[i]);
        if (i === childrenNum - 1) parent = prevParent;
      }
    }

    traverse(obj);

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
          node.branch_length ? `:${node.branch_length}` : []
        ).join('');
      } else {
        const childrenText = node.children!.map((childNode) => toNewickInt(childNode)).join(',');
        return ([] as string[]).concat(
          `(${childrenText})`,
          node.name,
          node.branch_length ? `:${node.branch_length}` : []
        ).join('');
      }
    }

    return !node ? ';' : `${toNewickInt(node)};`;
  }

  getLeafList<TNode extends NodeType>(node: TNode): TNode[] {
    if (node == null) return [];

    if (isLeaf(node)) {
      return [node]; // node is a leaf
    } else {
      return ([] as TNode[]).concat(
        ...(node.children ?? []).map((child) => this.getLeafList(child as TNode)));
    }
  }

  getNodeList<TNode extends NodeType>(node: TNode): TNode[] {
    if (isLeaf(node)) {
      return [node]; // node is a leaf
    } else {
      const childNodeListList: TNode[][] = node.children!
        .map((child) => { return this.getNodeList(child as TNode); });
      return ([] as TNode[]).concat(
        [node],
        ...childNodeListList);
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

  filterTreeByLeaves(node: NodeType, leaves: { [name: string]: any }): NodeType | null {
    // copy node because phylocanvas.gl changes data structure completely
    const resNode = Object.assign({}, node); // shallow copy

    if (isLeaf(resNode)) {
      return resNode.name in leaves ? resNode : null;
    } else {
      resNode.children = node.children!
        .map((child) => this.filterTreeByLeaves(child, leaves))
        .filter((child) => child != null) as NodeType[];

      return resNode.children.length > 0 ? resNode : null;
    }
  }

  /***/
  getNodesByLeaves<TNode extends NodeType>(node: TNode, leaves: { [name: string]: any }): TNode[] {
    if (isLeaf(node)) {
      return node.name in leaves ? [node] : [];
    } else {
      const children: TNode[] = node.children as TNode[] ?? [];
      const childrenRes: TNode[] = children.map((child) => {
        return this.getNodesByLeaves(child, leaves);
      }).flat();

      return (childrenRes.length == children.length &&
        children.every((child, i) => child == childrenRes[i])) ? [node] : childrenRes;
    }
  }


  /** Cuts tree, gets list of clusters as lists of leafs */
  treeCutAsLeaves(
    node: NodeType, cutHeight: number, currentHeight: number = 0
  ): NodeType[] {
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

  /** Cuts tree, gets cutted tree of clusters as lists of leafs */
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

  /** Sets grid's row order and returns tree (root node) of nodes presented in data */
  setGridOrder(
    tree: NodeType, grid: DG.Grid, leafColName?: string,
    removeMissingDataRows: boolean = false
  ): [NodeType, string[]] {
    console.debug('Dendrogram.setGridOrder() start');

    const dataDf: DG.DataFrame = grid.dataFrame;

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
    let dataDfFilterModified: boolean = false;
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
        if (dataNodeOldFilter) {
          dataDf.filter.set(dataRowI, false, false);
          dataDfFilterModified = true;
        }
      }
    }

    // TODO: Fire filter event (?)

    const resTree: NodeType = this.filterTreeByLeaves(tree, dataNodeDict)!;
    const resTreeLeafList = this.getLeafList(resTree);

    const order: number[] = new Array<number>(resTreeLeafList.length); // rowCount filtered for leaves
    for (let leafI = 0; leafI < resTreeLeafList.length; ++leafI) {
      const leafNodeName = resTreeLeafList[leafI].name;
      if (leafNodeName in dataNodeDict) {
        order[leafI] = dataNodeDict[leafNodeName];
      } else {
        missedDataNodeList.push(leafNodeName);
        // TODO: consider to add empty row with leaf name
      }
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
    clusterCol.init((rowI) => { return naValue; });

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
      let k = 9;
    }
    let k = 11;
  }

  buildClusters(tree: NodeCuttedType, clusterDf: DG.DataFrame, clusterColName: string, leafColName?: string): void {
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

  /** Cuts tree at threshold (from root), returns array of subtrees, marks leafs for cluster in data */
  cutTreeToGrid(
    node: NodeType, cutHeight: number, dataDf: DG.DataFrame,
    leafColName: string, clusterColName: string, na?: any
  ): void {
    const clusterList: NodeType[] = this.treeCutAsLeaves(node, cutHeight, 0);

    const naValue: any = na ?? null;
    const clusterCol: DG.Column = dataDf.getCol(clusterColName);
    // for (let rowI = 0; rowI < clusterCol.length; rowI++)
    //   clusterCol.set(rowI, na_value);
    clusterCol.init((rowI) => { return naValue; });

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
}


