import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {newickToDf} from '../utils';
import {Unsubscribable} from 'rxjs';
import {IPhylocanvasGlViewer} from '@datagrok-libraries/bio';

export class GridWithTreeViewer extends DG.JsViewer {
  nodeNameColumnName: string | null;

  grid: DG.Grid | null = null;
  tree: IPhylocanvasGlViewer | DG.JsViewer | null = null;

  _newick: string | null = null;
  _nwkDf: DG.DataFrame | null = null;

  readonly divH: HTMLDivElement;
  readonly leftDiv: HTMLDivElement;
  readonly rightDiv: HTMLDivElement;

  readonly textDiv: HTMLDivElement;

  get newick(): string { return this._newick!; }

  get nwkDf(): DG.DataFrame { return this._nwkDf!; }

  constructor() {
    super();

    this.nodeNameColumnName = this.string('nodeNameColumnName', null);

    this.textDiv = ui.div('some text', {style: {backgroundColor: 'window'}});

    this.leftDiv = ui.div([], {style: {backgroundColor: '#F0FFF0'}});
    this.rightDiv = ui.div([], {style: {backgroundColor: '#FFF0F0', overflow: 'hidden'}});
    this.divH = ui.splitH([this.leftDiv, this.rightDiv], {}, true);

    this.root.append(this.divH);

    this.subs.push(ui.onSizeChanged(this.root).subscribe(() => { this.rootCalcSize(); }));
    this.subs.push(ui.onSizeChanged(this.leftDiv).subscribe(() => { this.leftCalcSize(); }));
    this.subs.push(ui.onSizeChanged(this.rightDiv).subscribe(() => { this.rightCalcSize(); }));
  }

  override onPropertyChanged(property: DG.Property | null) {
    super.onPropertyChanged(property);

    if (property) {
      switch (property.name) {
      case 'nodeNameColumnName':
        this.updateGrid();
        break;
      }
    }
  }

  override onTableAttached() {
    window.setTimeout(async () => {
      await this.destroyView();

      super.onTableAttached();
      this._newick = this.dataFrame.getTag('.newick');

      this._nwkDf = newickToDf(this.newick, 'nwkDf');
      const leafCol = this._nwkDf.getCol('leaf');
      this._nwkDf.filter.init((rowI: number) => { return leafCol.get(rowI); });

      await this.buildView();
    }, 0 /* next event cycle */);
  }

  //# region -- View --

  rootCalcSize() {
    const cw: number = this.root.clientWidth;
    const ch: number = this.root.clientHeight;
    this.divH.style.width = `${cw}px`;
    this.divH.style.height = `${ch}px`;
  }

  leftCalcSize() {
    if (this.tree) {
      const cw: number = this.leftDiv.clientWidth;
      const ch: number = this.leftDiv.clientHeight;

      this.tree.root.style.width = `${cw}px`;
      // this.tree.root.style.height = `${ch}px`;
    }
  }

  rightCalcSize() {
    if (this.grid) {
      const cw: number = this.rightDiv.clientWidth;
      const ch: number = this.rightDiv.clientHeight;
      this.grid.root.style.width = `${cw}px`;
      // this.grid.root.style.height = `100%`;
    }
  }

  viewSubs: Unsubscribable[] = [];

  async destroyView() {
    for (const sub of this.viewSubs)
      sub.unsubscribe();
    this.viewSubs = [];

    if (this.grid) {
      this.grid.close();
      this.grid.removeFromView();
      this.grid = null;
    }

    if (this.tree) {
      this.tree.close();
      this.tree.removeFromView();

      this.textDiv.remove();
    }
  }

  async buildView() {
    if (!this.grid) {
      this.grid = await this.dataFrame.plot.fromType(DG.VIEWER.GRID, {}) as DG.Grid;
      this.rightDiv.append(this.grid.root);
      this.grid.root.style.height = '100%';

      this.viewSubs.push(this.grid.onRowsResized.subscribe((value: any) => {
        console.debug('PhyloTreeViewer: GridWithTreeViewer.buildView/grid.onRowsResized()');
        this.updateTree();
      }));
    }

    if (!this.tree) {
      this.tree = (await this.nwkDf.plot.fromType('PhylocanvasGl', {
        interactive: false,
        alignLabels: true,
        showLabels: true,
        showLeafLabels: true,
        nodeSize: 1,
        padding: 10,
        treeToCanvasRatio: 1,
      })) as DG.JsViewer;
      this.leftDiv.append(this.tree.root, this.textDiv);

      this.tree.root.style.position = 'absolute';
      this.tree.root.style.top = `${this.grid.colHeaderHeight}px`;
      this.tree.root.style.border = 'solid 1px green';
    }

    // const styles: { [nodeName: string]: { [prop: string]: any } } = {};
    // const nodeCol: DG.Column = this.nwkDf.getCol('node');
    // for (const nodeName of nodeCol.values()) {
    //   styles[nodeName] = {label: ' ',};
    // }
    // this.tree.setStyles(styles);

    this.updateGrid();
    this.updateTree();
  }

  /** Update data grid (this.dataFrame) for order of leaves in newick */
  updateGrid() {
    if (!this.grid)
      return;

    // Build data node name dictionary
    const dataNodeIdx: { [nodeName: string]: number } = {};
    const dataNodeCol: DG.Column = this.dataFrame.getCol(this.nodeNameColumnName!);
    for (let dataRowI = 0; dataRowI < dataNodeCol.length; dataRowI++) {
      const dataNodeName: string = dataNodeCol.get(dataRowI);
      dataNodeIdx[dataNodeName] = dataRowI;
    }

    const order: number[] = new Array<number>(this.nwkDf.filter.trueCount); // rowCount filtered for leaves
    const nwkNodeCol: DG.Column = this.nwkDf.getCol('node');
    let leafI: number = 0;
    for (let nodeI = 0; nodeI < this.nwkDf.rowCount; ++nodeI) {
      if (this.nwkDf.filter.get(nodeI)) {
        const leafNodeName = nwkNodeCol.get(nodeI);
        if (!(leafNodeName in dataNodeIdx)) {
          throw new Error(`There is no data row for leaf node '${leafNodeName}'`);
          // TODO: consider to add empty row with leaf name
        }
        order[leafI] = dataNodeIdx[leafNodeName];
        ++leafI;
      }
    }

    this.grid.setRowOrder(order);
  }

  //# endregion -- View --

  // -- Handle events --

  updateTree() {
    if (!this.grid || !this.tree)
      return;

    /* this.dataFrame contains data (about leaves)
     * this.nwkDf.contains nodes from newick filtered for leaves only
     */

    const headerHeight = this.grid.colHeaderHeight; // this.grid.props.colHeaderHeight is null
    const rowHeight = this.grid.props.rowHeight;

    const leafCount = this.nwkDf.filter.trueCount;

    this.tree.root.style.height = `${rowHeight * (leafCount)}px`;
    this.tree.root.style.top = `${this.grid.colHeaderHeight}px`;

    let k = 11;
  }
}
