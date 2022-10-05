import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {newickToDf} from '../utils';
import {Shapes} from '@phylocanvas/phylocanvas.gl';
import {PhylocanvasGlViewer, TreeTypesNames} from './phylocanvas-gl-viewer';
import {DOCK_TYPE} from 'datagrok-api/dg';
import {Unsubscribable} from 'rxjs';

export class GridWithTreeViewer extends DG.JsViewer {

  grid: DG.Grid | null = null;
  tree: PhylocanvasGlViewer | null = null;

  _newick: string | null = null;
  _nwkDf: DG.DataFrame | null = null;

  get newick(): string { return this._newick!; }

  get nwkDf(): DG.DataFrame { return this._nwkDf!; }

  constructor() {
    super();

    const divLeft = ui.div([], {style: {backgroundColor: '#E0FFE0'}});
    const divRight = ui.div([], {style: {backgroundColor: '#FFE0E0'}});
    const divH = ui.splitH([divLeft, divRight], {}, true);

    this.root.append(divH);

    this.subs.push(ui.onSizeChanged(this.root).subscribe(this.rootOnSizeChanged.bind(this)));
  }

  rootOnSizeChanged() {
    this.calcSize();
  }

  override async onTableAttached(): Promise<void> {
    window.setTimeout(async () => {
      await this.destroyView();

      super.onTableAttached();
      this._newick = this.dataFrame.getTag('.newick');
      this._nwkDf = newickToDf(this.newick, 'nwkDf');

      await this.buildView();
    }, 0 /* next event cycle */);
  }

  //# region -- View --

  calcSize() {

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
    }
  }

  async buildView() {
    if (!this.grid) {
      this.grid = await this.dataFrame.plot.fromType(DG.VIEWER.GRID, {}) as DG.Grid;
    }

    if (!this.tree) {
      this.tree = await this.nwkDf.plot.fromType('PhylocanvasGl', {
        alignLabels: true,
        showLabels: true,
        showLeafLabels: true,
        nodeSize: 1,
      }) as PhylocanvasGlViewer;
    }

    const styles: { [nodeName: string]: { [prop: string]: any } } = {};
    for (const nodeName in this.dataFrame.getCol('node')) {
      styles[nodeName] = {label: ' ',};
    }
    this.tree.setStyles(styles);
  }

  //# endregion -- View --
}
