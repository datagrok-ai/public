//(\(*[^:]+:[^,]+\)*)+
// import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PhylocanvasGL, TreeTypes, Shapes} from '@phylocanvas/phylocanvas.gl';
import {TreeAnalyzer, PhylocanvasTreeNode} from './utils/tree-stats';

export class TreeBrowser {// extends DG.JsViewer {
  idColumnName: string = 'v id';
  title: string;
  phyloTreeViewer: PhylocanvasGL;
  networkViewer: DG.Viewer;
  dataFrame: DG.DataFrame;
  mlbView: DG.TableView;
  // network: DG.DataFrame;
  leaves: TreeLeaves = {};
  treeLeaves: TreeLeavesMap = {};
  idMappings: TreeLeavesMap = {};
  treeAnalyser: TreeAnalyzer;
  treeGrid: DG.Grid;

  protected _modifyTreeNodeIds(nwk: string): string {
    if (TreeAnalyzer.newickRegEx.test(nwk.trim()))
      return nwk.replaceAll(/([^|,:()]+)\|([^|,:()]+)\|([^|,:()]+)\|([^|,:()]+)/g, '$3');

    return nwk;
  }

  protected _takeTreeAt(index: number, treeColumnName = 'TREE'): string {
    return this.dataFrame.get(treeColumnName, index);
  }

  // constructor() {
  //   super();
  //   this.title = this.string('title', 'Phylogenetic tree');
  // }

  private _initGrid(treeCol = 'TREE', semType = 'newick'): DG.Grid {
    this.dataFrame.col(treeCol).semType = semType;

    const grid = DG.Viewer.grid(this.dataFrame);
    grid.setOptions({
      rowHeight: 160,
    });

    grid.onCellRender.subscribe(function(args) {
      if (args.cell.isTableCell && args.cell.tableColumn.semType == semType) {
        const nwk: string = args.cell.cell.value;
        const regexp = new RegExp(/^(\(*[^:]+:[^,]+\)*)+$/);

        if (regexp.test(nwk.trim())) {
          /*const ctx = args.g.canvas;//.getContext('webgl').canvas;
          const phTree = new PhylocanvasGL(ctx, {
            shape: Shapes.Dot,
            size: ctx.getBoundingClientRect(),
            source: nwk,
            type: TreeTypes.Rectangular,
          });
          console.warn([ctx.getBoundingClientRect(), phTree]);*/
          //console.warn(newickParser.parse_newick(nwk));
          //args.preventDefault();
        }
      }
    });

    return grid;
  }

  /**
   * Maps tree leaf {node id} to a {list of indices} of tree which it is found in.
   * @param {PhylocanvasTreeNode} node Node record.
   * @param {number} index Index of tree containing the node in trees column.
   * @return {PhylocanvasTreeNode} Modified node.
   */
  private _leavesProcessor(node: PhylocanvasTreeNode, index: number) {
    if (!node.isLeaf)
      return;

    const id = node.id;
    const cloneId = this.dataFrame.get('CLONE', index);

    if (!this.leaves[cloneId])
      this.leaves[cloneId] = {index: [index], items: []};

    this.leaves[cloneId].items.push(id);

    if (!this.leaves[id])
      this.leaves[id] = {index: [], items: []};

    this.leaves[id]['index'].push(index);

    if (!this.treeLeaves[id])
      this.treeLeaves[id] = [];

    this.treeLeaves[id].push(index);
    return node;
  }

  /**
   * Adds a grid with simple tree statistics.
   * @param {string[]} trees List of strings in Newick format.
   * @param {DG.Column[]} [baseColumns=[]] Columns to include into creating grid.
   * @return  {DG.Grid} Grid with statistics.
   */
  private _makeTreeGrid(trees: string[], baseColumns: DG.Column[] = []): DG.Grid {
    const _calcFilteredItems = function() {
      const filteredIndices = new Set(this.mlbView.dataFrame.filter.getSelectedIndexes());
      const isntFiltered = (x: number) => filteredIndices.has(x);
      return Object.keys(this.idMappings).filter((_, i) => isntFiltered(i));
    }.bind(this);

    const itemsToFind = _calcFilteredItems();

    this.treeAnalyser = new TreeAnalyzer(itemsToFind);
    this.treeAnalyser.addNodeInspector(this._leavesProcessor.bind(this));

    this.mlbView.dataFrame.onRowsFiltered.subscribe((args) => {
      this.treeAnalyser.items = _calcFilteredItems();
    });

    const stats = this.treeAnalyser.analyze(trees);
    const df = DG.DataFrame.fromColumns(baseColumns);

    for (const k of Object.keys(stats)) {
      const dType = k.toLowerCase().includes('length') ? 'double' : 'int';
      (df.columns as DG.ColumnList).add(DG.Column.fromList(dType, k, stats[k]));
    }

    // TODO: fix this ad hoc.
    // df.rows.removeAt(0);

    return DG.Viewer.grid(df);
  }

  /**
   * Modifies trees containing column of the given data frame to reduce leafs label.
   * @param {DG.DataFrame} df Target data frame.
   * @param {string} [treeColumnName='TREE'] Trees column name.
   * @return {DG.DataFrame} Modified data frame.
   */
  private _modifyTreeColumn(df: DG.DataFrame, treeColumnName: string = 'TREE'): DG.DataFrame {
    const treeCol = df.col(treeColumnName);
    const trees = treeCol.toList().map((v) => this._modifyTreeNodeIds(v));
    (df.columns as DG.ColumnList).replace(treeColumnName, DG.Column.fromStrings(treeColumnName, trees));
    return df;
  }

  /**
   * Fires when tree grid row is changed.
   * @param {*} args Callback arguments.
   */
  private _onTreeGridCurrentRowChanged(args: any) {
    const currentRowIdx = this.treeGrid.dataFrame.currentRowIdx;
    const cloneId = this.treeGrid.dataFrame.get('CLONE', currentRowIdx);

    if (this.leaves[cloneId])
      this.selectTree(this.leaves[cloneId]['index'][0]);
  }

  /**
   * Initializes tree browser.
   * @param {DG.DataFrame} df Table to convert.
   * CLONE ... TREE
   * 123 ... ((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930...
   * 123 ... (Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386)...
   * 456 ... (Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:...
   * @param {DG.TableView} mlbView
   */
  init(df: DG.DataFrame, mlbView: DG.TableView) {
    this.dataFrame = this._modifyTreeColumn(df);
    this.mlbView = mlbView;

    const treeCol = this.dataFrame.col('TREE');

    this.idMappings = this._collectMappings();
    this.treeGrid = this._makeTreeGrid(treeCol.toList(), ['CLONE'].map((v) => this.dataFrame.col(v)));

    const treeDiv = ui.div([]);
    const treeNode = mlbView.dockManager.dock(treeDiv, DG.DOCK_TYPE.DOWN);

    mlbView.dockManager.dock(this.treeGrid, DG.DOCK_TYPE.RIGHT, treeNode);
    this.treeGrid.dataFrame.currentRowIdx = 1;

    this.phyloTreeViewer = new PhylocanvasGL(treeDiv, {
      interactive: true,
      showLabels: true,
      showLeafLabels: true,
      shape: Shapes.Dot,
      size: treeDiv.getBoundingClientRect(),
      source: this._takeTreeAt(this.treeGrid.dataFrame.currentRowIdx),
      type: TreeTypes.Rectangular,
    });

    this.phyloTreeViewer.selectNode = this.selectNode.bind(this);
    this.treeGrid.dataFrame.onCurrentRowChanged.subscribe(this._onTreeGridCurrentRowChanged.bind(this));
    this._matchMappings();
  }

  /**
   * Makes mapping of {tree node id} taken from column value to {list of indices} which it is found in.
   * @return {TreeLeavesMap} Mapping.
   */
  private _collectMappings(): TreeLeavesMap {
    const col = this.mlbView.dataFrame.col(this.idColumnName);
    const mapper: TreeLeavesMap = {};

    for (let i = 0; i < col.length; ++i) {
      const id: string = col.get(i);

      if (!mapper[id])
        mapper[id] = [];

      mapper[id].push(i);
    }
    return mapper;
  }

  /**
   * Matches node ids from initial table found within trees.
   * Adds an auxiliary column containing clone ids of matched trees.
   * @param {string} [columnName='clones'] Auxiliary column name.
   */
  private _matchMappings(columnName = 'clones') {
    const df = this.mlbView.dataFrame;
    const col = (df.columns as DG.ColumnList).addNewString(columnName);

    for (const [id, indices] of Object.entries(this.idMappings)) {
      if (this.treeLeaves[id]) {
        for (const i of indices)
          col.set(i, [col.get(i), this.treeLeaves[id].join('|')].join('|'));
      }
    }
  }

  /**
   * Finds and selects tree node chosen.
   * @param {*} node Node to consider.
   */
  selectNode(node: any) {
    if (node) {
      if (node.label) {
        const nodeId: string = node?.label;
        const df = this.mlbView.dataFrame;
        const col = df.col(this.idColumnName);

        df.selection.init((i) => (col.get(i) as string).includes(nodeId));

        if (df.selection.trueCount > 0)
          df.currentRowIdx = df.selection.getSelectedIndexes()[0];
      }
    } else {
    }
  }

  /**
   * Selects a tree from the given index in the table.
   * @param {number} index Index to take tree from.
   */
  selectTree(index: number) {
    const maxIndex = this.dataFrame.rowCount;
    const currentIndex = index >= maxIndex ? maxIndex - 1 : (index < 0 ? 0 : index);

    this.dataFrame.currentRowIdx = currentIndex;
    this.phyloTreeViewer.setProps({source: this._takeTreeAt(currentIndex)});

    const treeItems = this.leaves[this.dataFrame.get('CLONE', currentIndex)].items;

    /**
     * Modifies node styles to mark intersected node ids.
     * @return {any}
     */
    const _modifyNodeStyles = function() {
      let styles = {};

      for (const item of treeItems) {
        const style = {};

        if (this.idMappings[item]) {
          const color = this.treeAnalyser.getItemsAsSet().has(item) ? '#0000ff' : '#ff0000';
          style[item] = {fillColour: color};
        } else
          style[item] = {shape: Shapes.Dot};
        styles = {...styles, ...style};
      }
      return styles;
    }.bind(this);

    this.phyloTreeViewer.setProps({styles: _modifyNodeStyles()});
  }

  // get root(): HTMLElement {
  //   const title = ui.h1(this.title, {style: {'align-self': 'center', 'alignContent': 'center'}});
  //   if (this.phyloTreeViewer && this.networkViewer) {
  //     [this.phyloTreeViewer.root, this.networkViewer.root].forEach((v) => v.style.width = 'auto');
  //     return ui.divV([title, ui.divH([this.phyloTreeViewer.root, this.networkViewer.root])]);
  //   }
  //   return title;
  // }
}

interface TreeMap {
  index: number[];
  items: string[]
}

type TreeLeaves = {[key: string]: TreeMap};
type TreeLeavesMap = {[key: string]: number[]};
