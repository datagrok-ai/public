//(\(*[^:]+:[^,]+\)*)+
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {TreeAnalyzer, getVId} from './utils/tree-stats';
import {Subscription, Unsubscribable} from 'rxjs';
import {PickingInfo} from '@deck.gl/core/typed';
import {MjolnirPointerEvent} from 'mjolnir.js';

export class TreeBrowser extends DG.JsViewer {
  viewed: boolean = false;
  th: bio.ITreeHelper;
  newickHelper: bio.INewickHelper | null = null;

  static treeGridColumnsNameMapping: { [key: string]: { name: string, dType: string } } = {
    totalLeaves: {name: 'Leaves', dType: 'int'},
    leavesIntersected: {name: 'Intersected', dType: 'int'},
    totalSubtreeLength: {name: 'Height', dType: 'double'},
  };

  idColumnName: string = 'v id';
  treeSemanticType = 'newick';
  title: string;

  treeDiv?: HTMLDivElement;
  treeViewer?: bio.IPhylocanvasGlViewer;

  _treeGrid?: DG.Grid;
  _treeGridSubs: Unsubscribable[] = [];

  cloneIndex: TreeLeaves = {};
  leavesIndex: TreeLeavesMap = {};
  vIdIndex: TreeLeavesMap = {};

  treeAnalyser?: TreeAnalyzer;

  private _treeDf: DG.DataFrame = DG.DataFrame.fromCsv(['TREE', 'CLONE'].join(','));
  private _mlbDf: DG.DataFrame;

  get treeDf(): DG.DataFrame { return this._treeDf; }

  get mlbDf(): DG.DataFrame { return this._mlbDf; }

  private viewSubs: Subscription[] = [];

  /** Assigns data to viewer, calls destroyView() / buildView(). Use it instead of .dataFrame */
  public async setData(treeDf: DG.DataFrame, mlbDf: DG.DataFrame): Promise<void> {
    console.debug('MLB: TreeBrowser.setData()');

    if (!this.newickHelper)
      this.newickHelper = await grok.functions.call('PhyloTreeViewer:getNewickHelper') as bio.INewickHelper;

    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    this._mlbDf = mlbDf;
    this._treeDf = treeDf;

    if (!this.viewed) {
      await this.buildView();
      this.viewed = true;
    }
  }

  protected _takeTreeAt(index: number, columnName = 'TREE'): string {
    return this.treeDf.get(columnName, index);
  }

  protected _getCloneIdAt(index: number, columnName = 'CLONE'): string {
    return this.treeDf.get(columnName, index);
  }

  constructor() {
    super();

    this.subs.push(ui.onSizeChanged(this.root).subscribe(
      this.rootOnSizeChanged.bind(this)));
  }

  /**
   * Maps tree leaf {node id} or {clone id} to a {list of indices} of tree which it is found in.
   * @param {PhylocanvasTreeNode} node Node record.
   * @param {number} treeIndex Index of tree containing the node in trees column.
   * @return {PhylocanvasTreeNode} Modified node.
   */
  private _collectLeavesMapping(node: bio.PhylocanvasTreeNode, treeIndex: number) {
    if (node.isLeaf) {
      // console.debug(`MLB: TreeBrowser._collectLeavesMapping( node: ${node.name}, treeIndex: ${treeIndex} )`);
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

  // /**
  //  * Returns filtered rows in MLB grid.
  //  * @return {string[]} Rows is left after filtering.
  //  */
  // private _calcFilteredItems(): string[] {
  //   const filteredIndices = new Set(this.mlbDf.filter.getSelectedIndexes());
  //   const isntFiltered = (x: number) => filteredIndices.has(x);
  //   return Object.keys(this.vIdIndex).filter((_, i) => isntFiltered(i));
  // }

  /**
   * Makes mapping of {tree node id} taken from column value to {list of indices} which it is found in.
   * @return {TreeLeavesMap} Mapping.
   */
  private _collectMappings(): TreeLeavesMap {
    const mapper: TreeLeavesMap = {};
    if (this.mlbDf) {
      const col: DG.Column = this.mlbDf.getCol(this.idColumnName);
      for (let i = 0; i < col.length; ++i) {
        const id: string = col.get(i);

        if (!(id in mapper))
          mapper[id] = [];

        mapper[id].push(i);
      }
    }
    return mapper;
  }

  /**
   * Returns filtered rows in MLB grid.
   * @param {DG.DataFrame} mlbDf
   * @param {{}} vIdIndex
   * @return {string[]} Rows is left after filtering.
   */
  private static _calcFilteredItems(mlbDf: DG.DataFrame, vIdIndex: {}): string[] {
    let res: string [] = [];
    if (mlbDf) {
      const filteredIndices = new Set(mlbDf.filter.getSelectedIndexes());
      const isntFiltered = (x: number) => filteredIndices.has(x);
      res = Object.keys(vIdIndex).filter((_, i) => isntFiltered(i));
    }
    return res;
  }

  /**
   * Fires when MLB grid row was filtered.
   * @param {*} args Unused.
   */
  private _onMLBGridRowsFiltered(args: any) {
    const treeDf = this.treeDf;
    const sourceDf = this.mlbDf;
    const selectedVIdIndices = new Set(sourceDf.filter.getSelectedIndexes());

    this.treeAnalyser!.items = TreeBrowser._calcFilteredItems(sourceDf, this.vIdIndex);

    // Filtering of the data loaded for an antigen is not required.
    // treeDf.rows.filter((row: DG.Row) => {
    //   const cloneId = this._getCloneIdAt(row.idx);
    //
    //   if (this.cloneIndex[cloneId]) {
    //     const cloneVIds = this.cloneIndex[cloneId].items;
    //
    //     for (const vId of cloneVIds) {
    //       const indices = this.vIdIndex[vId];
    //
    //       if (indices) {
    //         for (const i of indices) {
    //           if (selectedVIdIndices.has(i))
    //             return true;
    //         }
    //       }
    //     }
    //   }
    //   return false;
    // });
    //treeDf.rows.requestFilter
  }

  onTreeGridCellTooltip(cell: DG.GridCell, x: number, y: number) {
    if (cell.isColHeader) {
      const msg = {
        [TreeBrowser.treeGridColumnsNameMapping.totalLeaves.name]: 'Total tree leaves',
        [TreeBrowser.treeGridColumnsNameMapping.leavesIntersected.name]: 'Tree leaves found also in the main taible',
        [TreeBrowser.treeGridColumnsNameMapping.totalSubtreeLength.name]: 'Tree height',
      };
      ui.tooltip.show(msg[cell.tableColumn!.name] ?? cell.tableColumn!.name, x, y);
      return true;
    }
  }

  private destroyTreeGrid() {
    this._treeGridSubs.forEach((s) => { s.unsubscribe(); });
    this._treeGridSubs = [];
    delete this._treeGrid;
  }

  /**
   * Adds a grid with simple tree statistics.
   * @param {string} treesColumnName Column containing trees in Newick format.
   * @param {string[]} [baseColumnNames=[]] Columns to include into creating grid.
   */
  private buildTreeGrid(treesColumnName: string, baseColumnNames: string[] = []): void {
    const treesColumn: DG.Column = this.treeDf.getCol(treesColumnName);
    const stats = this.treeAnalyser!.analyze(treesColumn.toList());
    const df = DG.DataFrame.fromColumns(baseColumnNames.map((colName: string) => this.treeDf.getCol(colName)));

    for (const k of Object.keys(stats)) {
      const col = TreeBrowser.treeGridColumnsNameMapping[k];
      (df.columns as DG.ColumnList).add(DG.Column.fromList(col.dType as DG.ColumnType, col.name, stats[k]));
    }

    // TODO: fix this ad hoc.
    // df.rows.removeAt(0);

    const treesCol: DG.Column = df.columns.add(treesColumn);
    treesCol.semType = this.treeSemanticType;

    const grid = DG.Viewer.grid(df, {
      rowHeight: 160,
    });

    grid.col('TREE')!.visible = false;

    this._treeGrid = grid;
    this._treeGridSubs = [];
    this._treeGridSubs.push(grid.onCellRender.subscribe(this._onTreeGridCellRender.bind(this)));
    this._treeGridSubs.push(grid.onCellTooltip(this.onTreeGridCellTooltip.bind(this)));
  }

  private rootOnSizeChanged(args: any): void {
    console.debug('MLB: TreeBrowser.rootOnSizeChanged()');
    this.calcSize();
  }

  calcSize(): void {
    console.debug('MLB: TreeBrowser.calcSize()');
    // this.treeDiv.innerText = `${this.root.clientWidth} x ${this.root.clientHeight}`;
    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;
    console.debug(`MLB: TreeBrowser.calcSize( ${cw.toString()} x ${ch.toString()} )`);

    if (this.treeDiv) {
      this.treeDiv.style.width = `${cw}px`;
      this.treeDiv.style.height = `${ch}px`;
      this.treeDiv.style.setProperty('overflow', 'hidden', 'important');
    }

    if (this.treeViewer) {
      //this.phyloTreeViewer.setProps({size: this.treeDiv.getBoundingClientRect()});
      this.treeViewer.root.style.width = `${cw}px`;
      this.treeViewer.root.style.height = `${ch}px`;
    }
  }

  public static cleanMlbNewick(nwkTxt: string): string {
    const treeClosePos: number = nwkTxt.substring(0, nwkTxt.length - 2).lastIndexOf(')');
    const treeTxtClean = nwkTxt.substring(1, treeClosePos + 1) + ';';
    return treeTxtClean;
  }

  private prepareNwkDf(nwkTxt: string, name: string): DG.DataFrame {
    // const nwkDf = await grok.functions.call('PhyloTreeViewer:_newickToDf',
    //   {newick: treeTxtClean, name: ''});
    const treeTxtClean = TreeBrowser.cleanMlbNewick(nwkTxt);
    const nwkDf: DG.DataFrame = this.newickHelper!.newickToDf(treeTxtClean, name);
    return nwkDf;
  }

  private async destroyView() {
    console.debug('MLB: TreeBrowser.destroyView() ');

    this.destroyTreeGrid();

    this.cloneIndex = {};
    this.leavesIndex = {};
    this.vIdIndex = {};

    delete this.treeAnalyser;

    this.viewSubs.forEach((s: Subscription) => { s.unsubscribe(); });
  }

  private async buildView() {
    console.debug('MLB: TreeBrowser.buildView() ');

    if (!this.th)
      this.th = await bio.getTreeHelper();

    this.vIdIndex = this._collectMappings();
    const itemsToFind = TreeBrowser._calcFilteredItems(this.mlbDf, this.vIdIndex);

    this.treeAnalyser = new TreeAnalyzer(itemsToFind);
    this.treeAnalyser.addNodeInspector(this._collectLeavesMapping.bind(this));

    // Rebuild treeGrid is required because of this.treeAnalyser.analyze() call required
    this.buildTreeGrid('TREE', ['CLONE']);

    //const color: string = `#bbff${Math.ceil(128 + Math.random() * 127).toString(16)}`;
    if (!this.treeDiv) {
      this.treeDiv = ui.div([], {
        style: {
          //backgroundColor: color,
          width: '100px',
          height: '100px',
          //margin: '10px',
          //backgroundColor: '#F0FFF0',
        }
      });
      this.root.appendChild(this.treeDiv);

      let nwkTxt: string = '((a:1),GERM);'; // side GERM will be removed in prepareNwkDf();
      const nwkDf: DG.DataFrame = this.prepareNwkDf(nwkTxt, '');
      const nodeShape: string = bio.Shapes.Circle;
      this.treeViewer = (await nwkDf.plot.fromType('PhylocanvasGL', {
        interactive: true,
        showLabels: true,
        showLeafLabels: true,
        nodeShape: nodeShape,
        nodeSize: 3,
        treeToCanvasRatio: 0.50,
        fontFamily: 'Roboto',
        fontSize: 13,
        treeType: bio.TreeTypesNames.Rectangular,
        // size: this.treeDiv.getBoundingClientRect(),
        // type: TreeTypes.Rectangular,
      })) as unknown as bio.IPhylocanvasGlViewer;

      // this.treeViewer.root.style.backgroundColor = '#FFF0F0';
      this.treeDiv.appendChild(this.treeViewer.root);
    }

    this.viewSubs.push(this.treeViewer!.onHover.subscribe(this.treeViewerOnHover.bind(this)));

    //this.phyloTreeViewer.selectNode = this.tvSelectNode.bind(this);
    //this.phyloTreeViewer.handleHover = this.tvHandleHover.bind(this);
    this.viewSubs.push(this.treeViewer!.nwkDf.onSelectionChanged.subscribe(() => {
      let k = 11;
      this.tvSelectNode.bind(this);
    }));

    // to fix blured image
    this.calcSize();

    this.viewSubs.push(
      this.treeDf.onCurrentRowChanged.subscribe(this._onTreeGridCurrentRowChanged.bind(this)));

    if (this.mlbDf) {
      this.viewSubs.push(
        this.mlbDf.onCurrentRowChanged.subscribe(this._onMLBGridCurrentRowChanged.bind(this)));
      this.viewSubs.push(
        this.mlbDf.onRowsFiltered.subscribe(this._onMLBGridRowsFiltered.bind(this)));
      this._matchMappings();
    }

    this.treeViewer!.root.style.visibility = this.treeDf.rowCount > 0 ? 'inherited' : 'hidden';
    if (this.treeDf.rowCount > 0)
      this.treeDf.currentRowIdx = 0;
  }

  // Do not handle onTableAttached() because data assigned with setData() method()
  // and it calls destroyView() / buildView() properly
  // override async onTableAttached() {
  //   console.debug(`MLB: TreeBrowser.onTableAttached( dataFrame = ${!this.dataFrame ? 'null' : 'value'} )`);
  //   super.onTableAttached();
  // }

  override async detach() {
    console.debug('MLB: TreeBrowser.detach() ' + `this.viewed = ${this.viewed}`);
    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    if (this.newickHelper)
      this.newickHelper = null;

    super.detach();
  }

  /**
   * Matches node ids from initial table found within trees.
   * Adds an auxiliary column containing clone ids of matched trees.
   * @param {string} [columnName='clones'] Auxiliary column name.
   */
  private _matchMappings(columnName = 'clones') {
    // Adding columns after tableView breaks display (no columns show)
    // const mlbDf = this.mlbDf;
    // let col: DG.Column<string> = mlbDf.columns.byName(columnName) as DG.Column<string> ||
    //   mlbDf.columns.addNewString(columnName);
    const clonesCol: DG.Column<string> = this.mlbDf.getCol('clones');

    for (const [vId, indices] of Object.entries(this.vIdIndex)) {
      if (this.leavesIndex[vId]) {
        for (const i of indices) {
          const values: string = clonesCol.get(i)!;
          const unique = new Set(values ? values.split('|').filter((v) => v.length > 0) : []);

          for (const v of this.leavesIndex[vId]) {
            const cloneId = this._getCloneIdAt(v);
            unique.add(cloneId);
          }

          clonesCol.set(i, Array.from(unique.values()).join('|'));
        }
      }
    }
  }

  treeViewerOnHover({info, event}: { info: PickingInfo, event: MjolnirPointerEvent }) {
    if (info.picked) {
      const args = event.srcEvent as MouseEvent;
      const tgtRect = this.treeDiv!.getBoundingClientRect();
      const obj = info.object;
      //@ts-ignore
      const centre = this.treeViewer.viewer.getCanvasCentrePoint();

      const leafList: bio.NodeType[] = this.th.getLeafList(obj as bio.NodeType);

      //@ts-ignore
      const toolTipText = leafList.map((l) => l.id).join('\n');
      ui.tooltip.show(ui.div(toolTipText),
        tgtRect.left + centre[0] + info.object.x, tgtRect.top + centre[1] + info.object.y + 12);


    } else {
      ui.tooltip.hide();
    }
  }

  /**
   * Finds and selects tree node chosen.
   * @param {*} node Node to consider.
   */
  tvSelectNode(node: any) {
    if (node) {
      if (node.label) {
        const nodeVId: string = getVId(node.id);
        const df = this.mlbDf;
        const col: DG.Column = df.getCol(this.idColumnName);

        df.selection.init((i) => (col.get(i) as string).includes(nodeVId));

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
    if (!this.treeViewer) return;

    const cloneId = this._getCloneIdAt(index);
    const treeItems = this.cloneIndex[cloneId].items;

    this.treeDf.currentRowIdx = index;
    const treeTxt: string = this._takeTreeAt(index);

    const nwkDf: DG.DataFrame = this.prepareNwkDf(treeTxt, cloneId);
    this.treeViewer.nwkDf = nwkDf;
    //this.phyloTreeViewer.setProps({source: treeTxtClean});

    /**
     * Modifies node styles to mark intersected node ids.
     * @return {any}
     */
    const _modifyNodeStyles = () => {
      let styles: { [nodeId: string]: { [propName: string]: any } } = {};

      for (const item of treeItems) {
        let style: bio.NodeStyleType = {};

        const nodeVId = getVId(item);
        if (nodeVId in this.vIdIndex) {
          const color = this.treeAnalyser!.getItemsAsSet().has(item) ? '#0000ff' : '#ff0000';
          style = {
            fillColour: color,
            shape: bio.Shapes.Circle,
            nodeSize: 5,
          };
        } else {
          style = {
            shape: bio.Shapes.Circle,
            nodeSize: 3,
          };
        }

        style['label'] = nodeVId;
        style['tooltip'] = 'hallo';

        styles[item] = style;
      }
      return styles;
    };

    const props: { [propName: string]: any } = {styles: _modifyNodeStyles()};
    if (!!node)
      props['selectedIds'] = [node];
    this.treeViewer.setProps(props);
    //this.treeViewer.nwkDf.setTag('styles', JSON.stringify(nodeStyles));
  }

  /**
   * Fires when tree grid row is changed.
   * @param {*} args Callback arguments.
   */
  private _onTreeGridCurrentRowChanged(args: any) {
    window.setTimeout(async () => {
      const currentRowIdx = this.treeDf.currentRowIdx;
      const cloneId = this._getCloneIdAt(currentRowIdx);

      if (cloneId in this.cloneIndex) {
        const treeIndex: number = this.cloneIndex[cloneId]['index'][0];
        await this.selectTree(treeIndex);
      }
    }, 0 /* next event cycle */);
  }

  /**
   * Fires when MLB grid row is changed.
   * @param {*} args Callback arguments.
   */
  private _onMLBGridCurrentRowChanged(args: any) {
    window.setTimeout(async () => {
      const df = this.mlbDf;
      const currentRowIdx = df.currentRowIdx;
      const currentCloneListStr: string = df.get('clones', currentRowIdx);
      if (!currentCloneListStr || currentCloneListStr.length == 0)
        return;

      const currentCloneList = currentCloneListStr.split('|');
      const currentCloneId = this.treeDf.get('CLONE', this.treeDf.currentRowIdx);
      const currentCloneIdStr = currentCloneId.toString();
      const cloneId = currentCloneList.includes(currentCloneIdStr) ? currentCloneId : currentCloneList[0];

      if (cloneId in this.cloneIndex) {
        const index = this.cloneIndex[cloneId]['index'][0];

        this.treeDf.currentRowIdx = index;
        await this.selectTree(index, df.get(this.idColumnName, currentRowIdx));
      }
    }, 0 /* next event cycle */);
  }
}

interface TreeMap {
  index: number[];
  items: string[];
}

type TreeLeaves = { [key: string]: TreeMap };
type TreeLeavesMap = { [key: string]: number[] };
