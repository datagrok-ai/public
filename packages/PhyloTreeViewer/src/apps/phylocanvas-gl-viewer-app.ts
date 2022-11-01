import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {newickToDf} from '../utils';
import {Unsubscribable} from 'rxjs';


export class PhylocanvasGlViewerApp {
  private viewed: boolean = false;
  private tv!: DG.TableView;

  // private ptv!: DG.Viewer; // PhyloTreeViewer

  treeHost!: HTMLDivElement | null;
  treeViewer!: bio.IPhylocanvasGlViewer | DG.JsViewer | null; // PhylocanvasGL
  treeDn!: DG.DockNode | null;

  _treeDf!: DG.DataFrame;

  get treeDf() { return this._treeDf; }

  constructor() {

  }

  async init(df?: DG.DataFrame): Promise<void> {
    await this.loadData(df);
  }

  async loadData(df?: DG.DataFrame): Promise<void> {
    if (df) {
      await this.setData(df);
    } else {
      const newickStr: string = await grok.dapi.files.readAsText('System:AppData/PhylotreeViewer/data/tree95.nwk');
      const treeDf: DG.DataFrame = newickToDf(newickStr, 'tree95');
      await this.setData(treeDf);
    }
  }

  async setData(treeDf: DG.DataFrame): Promise<void> {
    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    this._treeDf = treeDf;

    if (!this.viewed) {
      await this.buildView();
      this.viewed = true;
    }
  }

  //#region -- View --

  viewSubs: Unsubscribable[] = [];

  async destroyView(): Promise<void> {
    this.viewSubs.forEach((sub) => sub.unsubscribe());
    this.viewSubs = [];

    if (this.treeViewer) {
      this.treeViewer.close();
      this.treeViewer.removeFromView();
      this.treeViewer = null;
    }

    if (this.treeHost) {
      this.treeHost.remove();
      this.treeHost = null;
    }

    if (this.treeDn) {
      this.treeDn.detachFromParent();
      this.treeDn = null;
    }
  }

  async buildView(): Promise<void> {
    if (!this.tv) {
      // filter for leafs only, to align tree with grid
      const leafCol: DG.Column = this.treeDf.getCol('leaf');
      this.treeDf.filter.init((rowI: number) => { return leafCol.get(rowI); });

      this.tv = grok.shell.addTableView(this.treeDf, DG.DOCK_TYPE.FILL);
      this.tv.path = 'apps/PhyloTreeViewer/PhylocanvasGlViewer';

      this.viewSubs.push(this.tv.grid.onRowsResized.subscribe((args) => {
        if (!this.treeHost || !this.treeViewer) return;

        const cw: number = this.treeHost.clientWidth;
        const ch: number = this.treeHost.clientHeight;

        this.treeViewer.root.style.width = `${cw}px`;
        this.treeViewer.root.style.height = `${ch}px`;
      }));
    }

    // if (!this.treeHost) {
    //   this.treeHost = ui.div();
    // }
    //
    // if (!this.treeDn) {
    //   this.treeDn = this.tv.dockManager.dock(this.treeHost, DOCK_TYPE.LEFT);
    // }

    if (!this.treeViewer) {
      // TableView.addViewer() returns JsViewer (no access to viewer's attributes)
      // DataFrame.plot.fromType() return viewers type object (with attributes)
      this.treeViewer = (await this.treeDf.plot.fromType('PhylocanvasGl', {
        interactive: true,
        alignLabels: true,
        showLabels: true,
        showLeafLabels: true,
        padding: 0,
        treeToCanvasRatio: 1,
      })) as DG.JsViewer;
      this.treeDn = this.tv.dockManager.dock(this.treeViewer, DG.DOCK_TYPE.LEFT);
      let k = 11;
      //this.treeViewerDn = this.tv.dockManager.dock(this.treeViewer, DOCK_TYPE.LEFT);
    }
  }

  //#endregion -- View --
}