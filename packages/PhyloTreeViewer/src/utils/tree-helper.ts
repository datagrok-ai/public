import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import wu from 'wu';

import {PhylocanvasGL, Newick} from '@phylocanvas/phylocanvas.gl';


// function filterTree(node: NodeType, df: DG.DataFrame, leafColName?: string): NodeType {
//   // const grid :DG.Grid;
//   // grid.setRowOrder();
//   //
//   // df.filter.get
//   return null;
// }

type TreeLeafDict = { [nodeName: string]: bio.NodeType };
type DataNodeDict = { [nodeName: string]: number };
type NodeNameCallback = (nodeName: string) => void;

export class TreeHelper implements bio.ITreeHelper {
  toNewick(node: bio.NodeType | null): string {

    function toNewickInt(node: bio.NodeType): string {
      const isLeaf = !node.children || node.children.length == 0;

      if (isLeaf) {
        return Array<string>().concat(
          node.name,
          node.branch_length ? `:${node.branch_length}` : [],
        ).join('');
      } else {
        const childrenText = node.children!.map((childNode) => toNewickInt(childNode)).join(',');
        return Array<string>().concat(
          `(${childrenText})`,
          node.name,
          node.branch_length ? `:${node.branch_length}` : [],
        ).join('');
      }
    }

    return !node ? ';' : `${toNewickInt(node)};`;
  }

  getLeafList(node: bio.NodeType): bio.NodeType[] {
    if (node == null) return [];

    if (bio.isLeaf(node)) {
      return [node]; // node is a leaf
    } else {
      return Array<bio.NodeType>().concat(
        ...(node.children ?? []).map((child) => this.getLeafList(child)));
    }
  }

  getNodeList(node: bio.NodeType): bio.NodeType[] {
    if (bio.isLeaf(node)) {
      return [node]; // node is a leaf
    } else {
      const childNodeListList = node.children!.map((child) => this.getNodeList(child));
      return Array<bio.NodeType>().concat(
        [node],
        ...childNodeListList);
    }
  }

  treeFilterByLeaves(node: bio.NodeType, leaves: { [name: string]: any }): bio.NodeType | null {
    // copy node because phylocanvas.gl changes data structure completely
    const resNode = Object.assign({}, node); // shallow copy

    if (bio.isLeaf(resNode)) {
      return resNode.name in leaves ? resNode : null;
    } else {
      resNode.children = node.children!
        .map((child) => this.treeFilterByLeaves(child, leaves))
        .filter((child) => child != null) as bio.NodeType[];

      return resNode.children.length > 0 ? resNode : null;
    }
  }

  /** Cuts tree, gets list of clusters as lists of leafs */
  treeCutAsLeaves(
    node: bio.NodeType, cutHeight: number, currentHeight: number = 0
  ): bio.NodeType[] {
    let res: bio.NodeType[];

    const nodeBranchLength = node.branch_length ?? 0;
    if ((currentHeight + nodeBranchLength) < cutHeight) {
      /* if node has no children, then an empty list will be returned */
      res = new Array<bio.NodeType>().concat(
        ...(node.children ?? [])
          .map((child) => this.treeCutAsLeaves(child, cutHeight, currentHeight + nodeBranchLength)));
    } else {
      res = [node];
    }

    return res;
  }

  /** Cuts tree, gets cutted tree of clusters as lists of leafs */
  treeCutAsTree(node: bio.NodeType, cutHeight: number, currentHeight?: number): bio.NodeType | null {
    const nodeBranchHeight = node.branch_length ?? 0;
    const currentHeightV: number = currentHeight ?? 0;
    if ((currentHeightV + nodeBranchHeight) < cutHeight) {
      if (!bio.isLeaf(node)) {
        const res: bio.NodeType = Object.assign({}, node);
        res.children = node.children!
          .map((child) => {
            return this.treeCutAsTree(child, cutHeight, currentHeightV + nodeBranchHeight);
          })
          .filter((n) => n != null) as bio.NodeType[];
        return res;
      } else {
        return null;
      }
    } else {
      const res: bio.NodeCuttedType = Object.assign({}, node) as bio.NodeCuttedType;
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
    tree: bio.NodeType, grid: DG.Grid, leafColName: string,
    removeMissingDataRows: boolean = false
  ): [bio.NodeType, string[]] {
    console.debug('PhyloTreeViewer.setGridOrder() start');

    const dataDf: DG.DataFrame = grid.dataFrame;

    const missedDataNodeList: string[] = [];
    const missedTreeLeafList: string[] = [];

    // Build TREE node name dictionary
    const treeLeafList = this.getLeafList(tree); // ordered by tree
    const treeLeafDict: { [nodeName: string]: bio.NodeType } = {};
    for (const treeLeaf of treeLeafList) {
      if (treeLeaf.name in treeLeafDict)
        throw new Error('Non unique key tree leaf name');

      treeLeafDict[treeLeaf.name] = treeLeaf;
    }

    const dataNodeCol: DG.Column = dataDf.getCol(leafColName);
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

    const resTree: bio.NodeType = this.treeFilterByLeaves(tree, dataNodeDict)!;
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

    console.debug('PhyloTreeViewer.setGridOrder() ' + `resTreeLeafList.length = ${resTreeLeafList.length}`);

    console.debug('PhyloTreeViewer.setGridOrder() end');
    return [
      resTree,
      Array<string>().concat(missedDataNodeList, missedTreeLeafList),
    ];
  }

  /** Cuts tree at threshold (from root), returns array of subtrees, marks leafs for cluster in data */
  cutTreeToGrid(
    node: bio.NodeType, cutHeight: number, dataDf: DG.DataFrame, leafColName: string, clusterColName: string
  ): void {
    const clusterList: bio.NodeType[] = this.treeCutAsLeaves(node, cutHeight, 0);
    /* A leaf with cumulative height less than threshold
       will not be included in nor marked as cluster */

    const dataNodeDict: DataNodeDict = {};
    for (let dataRowI: number = 0; dataRowI < dataDf.rowCount; dataRowI++) {
      const dataNodeName = dataDf.get(leafColName, dataRowI);
      dataNodeDict[dataNodeName] = dataRowI;
    }

    function markCluster(node: bio.NodeType, cluster: number) {
      if (bio.isLeaf(node)) {
        const nodeName = node.name;
        const dataRowI: number = dataNodeDict[nodeName];
        dataDf.set(clusterColName, dataRowI, cluster);
      } else {
        for (const childNode of node.children!) {
          markCluster(childNode, cluster);
        }
      }
    }

    // const leafCol: DG.Column = dataDf.getCol(leafColName);

    for (let clusterI = 1; clusterI <= clusterList.length; clusterI++) {
      markCluster(clusterList[clusterI - 1], clusterI);
    }
  }

}

export class TreeToGridSyncer {
  private _th: TreeHelper = new TreeHelper();

  private readonly _nDiv: HTMLElement;
  private readonly _tree: bio.NodeType;
  private readonly _grid: DG.Grid;
  private readonly _dataDf: DG.DataFrame;

  private readonly _pc: PhylocanvasGL;
  private readonly _pcDiv: HTMLDivElement;

  private readonly _leafCol: DG.Column;

  /** order of tree leaves */
  private readonly _treeLeafList: bio.NodeType[];
  private readonly _treeLeafDict: TreeLeafDict;
  private readonly _dataNodeDict: DataNodeDict;

  private readonly _warnings: string[] = [];

  private readonly _missedDataNodeList: string[] = [];
  private readonly _missedTreeLeafList: string[] = [];

  get tree(): bio.NodeType { return this._tree; }

  get grid(): DG.Grid { return this._grid; }

  get dataDf(): DG.DataFrame { return this._dataDf; }

  get warnings(): string[] { return this._warnings; }

  constructor(
    nDiv: HTMLElement, tree: bio.NodeType, phylocanvas: PhylocanvasGL, grid: DG.Grid,
    leafColName?: string, fixDf: boolean = false
  ) {
    this._nDiv = nDiv;
    this._tree = tree;
    this._pc = phylocanvas;
    this._pcDiv = phylocanvas.view;
    this._grid = grid;
    this._dataDf = grid.dataFrame;

    this._leafCol = grid.dataFrame.getCol(leafColName!);

    if (fixDf) {
      [this._treeLeafList, this._treeLeafDict] = TreeToGridSyncer._getTreeLeafListAndDict(this._tree, this._th);
      this._dataDf = this.grid.dataFrame = TreeToGridSyncer._fixDataFrameForTree(
        grid.dataFrame, leafColName!, this._treeLeafList, this._treeLeafDict,
        (nn) => { this._missedDataNodeList.push(nn);});
    }

    // Recalculate after fix dataDf
    [this._treeLeafList, this._treeLeafDict] = TreeToGridSyncer._getTreeLeafListAndDict(this._tree, this._th);
    this._dataNodeDict = TreeToGridSyncer._getDataNodeDict(
      this._leafCol, this._treeLeafDict, (nn) => { this._missedDataNodeList.push(nn);});

    ui.onSizeChanged(nDiv).subscribe(this.nDivOnSizeChanged.bind(this));
    this._grid.dataFrame.onFilterChanged.subscribe(this.gridDataFrameOnFilterChanged.bind(this));
    this._grid.onBeforeDrawContent.subscribe(this.gridOnBeforeDrawContent.bind(this));
  }

  private static _getTreeLeafListAndDict(tree: bio.NodeType, th: bio.ITreeHelper): [bio.NodeType[], TreeLeafDict] {
    // Build TREE node name dictionary
    const treeLeafList = th.getLeafList(tree); // ordered by tree
    const treeLeafDict: TreeLeafDict = {};
    for (const treeLeaf of treeLeafList) {
      if (treeLeaf.name in treeLeafDict)
        throw new Error('Non unique key tree leaf name');

      treeLeafDict[treeLeaf.name] = treeLeaf;
    }
    return [treeLeafList, treeLeafDict];
  }

  /** Removes data rows missed in tree, and orders data rows in tree leaves order. */
  private static _fixDataFrameForTree(
    dataDf: DG.DataFrame, leafColName: string, treeLeafList: bio.NodeType[], treeLeafDict: TreeLeafDict,
    missedDataNodeCallback: NodeNameCallback
  ): DG.DataFrame {
    // Skip dataFrame fix for now
    // return dataNodeCol.dataFrame;

    const dataNodeDict: DataNodeDict = {};
    const dataNodeCol: DG.Column = dataDf.getCol(leafColName);
    for (let dataRowI = dataNodeCol.length - 1; dataRowI >= 0; dataRowI--) {
      const dataNodeName: string = dataNodeCol.get(dataRowI);
      if (!(dataNodeName in treeLeafDict)) {
        missedDataNodeCallback(dataNodeName);
        // dataDf.rows.removeAt(dataRowI, 1, false);
      } else {
        dataNodeDict[dataNodeName] = dataRowI;
      }
    }

    const colNameList: string[] = dataDf.columns.names();
    const resObjList: Object[] = treeLeafList.map((n) => {
      const oldDataRowI = dataNodeDict[n.name];
      return Object.assign({},
        ...colNameList.map((colName) => { return {[colName]: dataDf.get(colName, oldDataRowI)};}));
    });

    const resDataDf: DG.DataFrame = DG.DataFrame.fromObjects(resObjList)!;
    return resDataDf;
  }

  private static _getDataNodeDict(
    dataNodeCol: DG.Column, treeLeafDict: TreeLeafDict,
    missedDataNodeCallback: NodeNameCallback,
  ) {
    const dataDf: DG.DataFrame = dataNodeCol.dataFrame;

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
        missedDataNodeCallback(dataNodeName);

        // Hide grid data rows with node name not in tree leaf list
        if (dataNodeOldFilter) {
          dataDf.filter.set(dataRowI, false, false);
          dataDfFilterModified = true;
        }
      }
    }

    return dataNodeDict;
  }

  private sync(): void {
    const resTree: bio.NodeType = this._th.treeFilterByLeaves(this._tree, this._dataNodeDict)!;
    const resTreeLeafList = this._th.getLeafList(resTree);

    const order: number[] = new Array<number>(resTreeLeafList.length); // rowCount filtered for leaves
    for (let leafI = 0; leafI < resTreeLeafList.length; ++leafI) {
      const leafNodeName = resTreeLeafList[leafI].name;
      if (leafNodeName in this._dataNodeDict) {
        order[leafI] = this._dataNodeDict[leafNodeName];
      }
      // else {
      //   missedDataNodeList.push(leafNodeName);
      //   // TODO: consider to add empty row with leaf name
      // }
    }

    this._grid.setRowOrder(order);

    const source = resTree ? {type: 'biojs', data: resTree} :
      {type: 'biojs', data: {name: 'NONE', branch_length: 1, children: []}};

    this._pc.setProps({source: {type: 'biojs', data: resTree}});
  }

  // -- Layout --

  private calcLayout() {
    const leafCount = this._grid.dataFrame.filter.trueCount;

    const width: number = this._nDiv.clientWidth;
    let height: number = this._grid.props.rowHeight * leafCount;

    const firstVisibleRowIdx: number = Math.floor(this._grid.vertScroll.min);
    //const firstVisibleRowIdx: number = grid.vertScroll.min;

    console.debug('PhyloTreeViewer: injectTreeGridUI pcCalcSize() ' + `{ size: ${width} x ${height} }`);
    if (height == 0)
      height = 1;
    try {
      this._pc.setProps({size: {width, height}});
    } catch (ex) {}
    this._pcDiv.style.top = `${this._grid.colHeaderHeight - firstVisibleRowIdx * this._grid.props.rowHeight}px`;

    // if (cutSlider) {
    //   const cutDiv = cutSlider.root;
    //   cutDiv.style.top = `${0}px`;
    //   cutDiv.style.height = `${grid.colHeaderHeight}px`;
    //   cutDiv.style.width = `${width}px`;
    // }
  }

  // -- Handle events --

  private nDivOnSizeChanged(): void {
    this.calcLayout();
  }

  private gridOnBeforeDrawContent(): void {
    this.calcLayout();
  }

  private gridDataFrameOnFilterChanged(): void {
    // TODO: Filter newick tree
    console.debug('PhyloTreeViewer: injectTreeGridUI() grid.dataFrame.onFilterChanged()');

    // to prevent nested fire event in event handler
    window.setTimeout(() => {
      this.sync();
    }, 0 /* next event cycle */);
  }
}

