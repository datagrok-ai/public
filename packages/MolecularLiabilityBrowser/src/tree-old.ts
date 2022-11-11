//(\(*[^:]+:[^,]+\)*)+
// import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {TreeAnalyzer, mlbTreeNodeRe} from './utils/tree-stats';

export class TreeBrowserOld {// extends DG.JsViewer {
  static treeGridColumnsNameMapping: { [key: string]: { name: string, dType: string } } = {
    totalLeaves: {name: 'Leaves', dType: 'int'},
    leavesIntersected: {name: 'Intersected', dType: 'int'},
    totalSubtreeLength: {name: 'Height', dType: 'double'},
  };


  idColumnName: string = 'v id';
  treeSemanticType = 'newick';
  title: string;
  phyloTreeViewer: bio.PhylocanvasGL;
  dataFrame: DG.DataFrame;
  mlbView: DG.TableView;
  treeGrid: DG.Grid;

  cloneIndex: TreeLeaves = {};
  leavesIndex: TreeLeavesMap = {};
  vIdIndex: TreeLeavesMap = {};

  treeAnalyser: TreeAnalyzer;

  protected _modifyTreeNodeIds(nwk: string): string {
    // if (TreeAnalyzer.newickRegEx.test(nwk.trim()))
    //   return nwk.replaceAll(mlbTreeNodeRe, '$3');

    return nwk;
  }

  protected _takeTreeAt(index: number, columnName = 'TREE'): string {
    return this.dataFrame.get(columnName, index);
  }

  protected _getCloneIdAt(index: number, columnName = 'CLONE'): string {
    return this.dataFrame.get(columnName, index);
  }

  // constructor() {
  //   super();
  //   this.title = this.string('title', 'Phylogenetic tree');
  // }

  /**
   * Maps tree leaf {node id} or {clone id} to a {list of indices} of tree which it is found in.
   * @param {PhylocanvasTreeNode} node Node record.
   * @param {number} treeIndex Index of tree containing the node in trees column.
   * @return {PhylocanvasTreeNode} Modified node.
   */
  private _collectLeavesMapping(node: bio.PhylocanvasTreeNode, treeIndex: number) {
    if (node.isLeaf) {
      const id = node.id;
      const cloneId = this._getCloneIdAt(treeIndex);

      if (!this.cloneIndex[cloneId])
        this.cloneIndex[cloneId] = {index: [treeIndex], items: []};

      this.cloneIndex[cloneId].items.push(id);

      if (!this.leavesIndex[id])
        this.leavesIndex[id] = [];

      this.leavesIndex[id].push(treeIndex);
    }
    return node;
  }

  /**
   * Fires when tree grid cell is rendered.
   * @param {DG.GridCellRenderArgs} args Cell arguments.
   */
  private _onTreeGridCellRender(args: DG.GridCellRenderArgs) {
    if (args.cell.isTableCell && args.cell.tableColumn!.semType == this.treeSemanticType) {
      const nwk: string = args.cell.cell.value;

      if (TreeAnalyzer.newickRegEx.test(nwk.trim())) {
        //const ctx = args.g.canvas.getContext('webgl').canvas;
        /*const phTree = new PhylocanvasGL(ctx, {
          shape: Shapes.Dot,
          size: ctx.getBoundingClientRect(),
          source: nwk,
          type: TreeTypes.Rectangular,
        });*/
        //args.preventDefault();
      }
    }
  }

  /**
   * Returns filtered rows in MLB grid.
   * @return {string[]} Rows is left after filtering.
   */
  private _calcFilteredItems(): string[] {
    const filteredIndices = new Set(this.mlbView.dataFrame.filter.getSelectedIndexes());
    const isntFiltered = (x: number) => filteredIndices.has(x);
    return Object.keys(this.vIdIndex).filter((_, i) => isntFiltered(i));
  }

  /**
   * Fires when MLB grid row was filtered.
   * @param {*} args Unused.
   */
  private _onMLBGridRowsFiltered(args: any) {
    const treeDf = this.treeGrid.dataFrame;
    const sourceDf = this.mlbView.dataFrame;
    const selectedVIdIndices = new Set(sourceDf.filter.getSelectedIndexes());

    this.treeAnalyser.items = this._calcFilteredItems();

    treeDf.rows.filter((row: DG.Row) => {
      const cloneId = this._getCloneIdAt(row.idx);

      if (this.cloneIndex[cloneId]) {
        const cloneVIds = this.cloneIndex[cloneId].items;

        for (const vId of cloneVIds) {
          const indices = this.vIdIndex[vId];

          if (indices) {
            for (const i of indices) {
              if (selectedVIdIndices.has(i))
                return true;
            }
          }
        }
      }
      return false;
    });
    //treeDf.rows.requestFilter
  }

  onTreeGridCellTooltip(cell: DG.GridCell, x: number, y: number) {
    if (cell.isColHeader) {
      const msg = {
        [TreeBrowserOld.treeGridColumnsNameMapping.totalLeaves.name]: 'Total tree leaves',
        [TreeBrowserOld.treeGridColumnsNameMapping.leavesIntersected.name]: 'Tree leaves found also in the main taible',
        [TreeBrowserOld.treeGridColumnsNameMapping.totalSubtreeLength.name]: 'Tree height',
      };
      ui.tooltip.show(msg[cell.tableColumn!.name] ?? cell.tableColumn!.name, x, y);
      return true;
    }
  }

  /**
   * Adds a grid with simple tree statistics.
   * @param {string} treesColumnName Column containing trees in Newick format.
   * @param {string[]} [baseColumnNames=[]] Columns to include into creating grid.
   * @return  {DG.Grid} Grid with statistics.
   */
  private _makeTreeGrid(treesColumnName: string, baseColumnNames: string[] = []): DG.Grid {
    const itemsToFind = this._calcFilteredItems();

    this.treeAnalyser = new TreeAnalyzer(itemsToFind);
    this.treeAnalyser.addNodeInspector(this._collectLeavesMapping.bind(this));

    const treesColumn = this.dataFrame.getCol(treesColumnName);
    const stats = this.treeAnalyser.analyze(treesColumn.toList());
    const df: DG.DataFrame = DG.DataFrame.fromColumns(
      baseColumnNames.map((colName: string) => this.dataFrame.getCol(colName)));

    for (const k of Object.keys(stats)) {
      const col = TreeBrowserOld.treeGridColumnsNameMapping[k];
      (df.columns as DG.ColumnList).add(DG.Column.fromList(col.dType as DG.ColumnType, col.name, stats[k]));
    }

    // TODO: fix this ad hoc.
    // df.rows.removeAt(0);

    (df.columns as DG.ColumnList).add(treesColumn);
    df.getCol(treesColumnName).semType = this.treeSemanticType;

    const grid = DG.Viewer.grid(df, {
      //rowHeight: 160,
    });

    grid.col('TREE')!.visible = false; // this.treeColumnName

    grid.onCellRender.subscribe(this._onTreeGridCellRender.bind(this));
    grid.onCellTooltip(this.onTreeGridCellTooltip.bind(this));
    return grid;
  }

  /**
   * Modifies trees by reducing leafs label.
   * @param {DG.DataFrame} df Target data frame.
   * @param {string} [treeColumnName='TREE'] Trees column name.
   * @return {DG.DataFrame} Modified data frame.
   */
  private _modifyTreeColumn(df: DG.DataFrame, treeColumnName: string = 'TREE'): DG.DataFrame {
    const treeCol = df.getCol(treeColumnName);
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
    const cloneId = this._getCloneIdAt(currentRowIdx);

    if (this.cloneIndex[cloneId])
      this.selectTree(this.cloneIndex[cloneId]['index'][0]);
  }

  /**
   * Fires when MLB grid row is changed.
   * @param {*} args Callback arguments.
   */
  private _onMLBGridCurrentRowChanged(args: any) {
    const df = this.mlbView.dataFrame;
    const currentRowIdx = df.currentRowIdx;
    const clones: string = df.get('clones', currentRowIdx);

    if (!clones || clones.length == 0)
      return;

    const cloneId = clones.split('|')[0];

    if (this.cloneIndex[cloneId]) {
      const index = this.cloneIndex[cloneId]['index'][0];

      this.treeGrid.dataFrame.currentRowIdx = index;
      this.selectTree(index, df.get(this.idColumnName, currentRowIdx));
    }
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

    this.vIdIndex = this._collectMappings();
    this.treeGrid = this._makeTreeGrid('TREE', ['CLONE']);

    const treeDiv = ui.div([]);
    const treeNode = mlbView.dockManager.dock(treeDiv, DG.DOCK_TYPE.DOWN);

    mlbView.dockManager.dock(this.treeGrid, DG.DOCK_TYPE.RIGHT, treeNode);
    this.treeGrid.dataFrame.currentRowIdx = 1;

    this.phyloTreeViewer = new bio.PhylocanvasGL(treeDiv, {
      interactive: true,
      showLabels: true,
      showLeafLabels: true,
      shape: bio.Shapes.Dot,
      fontFamily: 'Roboto',
      fontSize: 13,
      size: treeDiv.getBoundingClientRect(),
      source: this._takeTreeAt(this.treeGrid.dataFrame.currentRowIdx),
      type: bio.TreeTypes.Rectangular,
    });

    this.phyloTreeViewer.selectNode = this.selectNode.bind(this);
    this.treeGrid.dataFrame.onCurrentRowChanged.subscribe(this._onTreeGridCurrentRowChanged.bind(this));
    this.mlbView.dataFrame.onCurrentRowChanged.subscribe(this._onMLBGridCurrentRowChanged.bind(this));
    this.mlbView.dataFrame.onRowsFiltered.subscribe(this._onMLBGridRowsFiltered.bind(this));
    this._matchMappings();
  }

  /**
   * Makes mapping of {tree node id} taken from column value to {list of indices} which it is found in.
   * @return {TreeLeavesMap} Mapping.
   */
  private _collectMappings(): TreeLeavesMap {
    const col = this.mlbView.dataFrame.getCol(this.idColumnName);
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

    for (const [vId, indices] of Object.entries(this.vIdIndex)) {
      if (this.leavesIndex[vId]) {
        for (const i of indices) {
          const values: string | null = col.get(i);
          const unique = new Set(values ? values.split('|').filter((v) => v.length > 0) : []);

          for (const v of this.leavesIndex[vId]) {
            const cloneId = this._getCloneIdAt(v);
            unique.add(cloneId);
          }

          col.set(i, Array.from(unique.values()).join('|'));
        }
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
        const col = df.getCol(this.idColumnName);

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
   * @param {string} [node] Node id to select.
   */
  selectTree(index: number, node?: string) {
    const cloneId: string = this._getCloneIdAt(index);
    const treeItems: string[] = this.cloneIndex[cloneId].items;

    this.treeGrid.dataFrame.currentRowIdx = index;
    this.phyloTreeViewer.setProps({source: this._takeTreeAt(index)});

    /**
     * Modifies node styles to mark intersected node ids.
     * @return {any}
     */
    const _modifyNodeStyles = () => {
      let styles: { [item: string]: { [property: string]: any } } = {};

      for (const item of treeItems) {
        const style: { [property: string]: any } = {};

        if (this.vIdIndex[item]) {
          const color = this.treeAnalyser.getItemsAsSet().has(item) ? '#0000ff' : '#ff0000';
          style[item] = {fillColour: color};
        } else {
          style[item] = {shape: bio.Shapes.Dot};
        }
        styles = {...styles, ...style};
      }
      return styles;
    };

    this.phyloTreeViewer.setProps({styles: _modifyNodeStyles(), selectedIds: node ? [node] : []});
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
  items: string[];
}

type TreeLeaves = { [key: string]: TreeMap };
type TreeLeavesMap = { [key: string]: number[] };
