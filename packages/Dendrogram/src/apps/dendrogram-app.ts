import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import {_package} from '../../src/package';
import {TreeHelper} from '../utils/tree-helper';

export class DendrogramApp {
  private viewed: boolean = false;

  private _funcName: string = '';

  _treeDf!: DG.DataFrame;

  get treeDf() { return this._treeDf; }

  constructor() {}

  async init(df?: DG.DataFrame, funcName: string = 'dendrogramApp'): Promise<void> {
    this._funcName = funcName;
    await this.loadData(df);
  }

  async loadData(df?: DG.DataFrame): Promise<void> {
    const th = new TreeHelper();
    if (!df) {
      const newickStr: string = await _package.files.readAsText('data/tree95.nwk');
      df = th.newickToDf(newickStr, 'tree95');
    }
    await this.setData(df);
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

  // -- View --

  viewSubs: Unsubscribable[] = [];

  tv!: DG.TableView;
  treeHost!: HTMLDivElement | null;
  treeViewer!: DG.JsViewer | null; // Dendrogram
  treeDn!: DG.DockNode | null;


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

      this.tv = grok.shell.addTableView(this.treeDf, DG.DOCK_TYPE.FILL);
      this.tv.path = this.tv.basePath = `/func/${_package.name}.${this._funcName}`;
    }

    if (!this.treeViewer) {
      this.treeViewer = (await this.treeDf.plot.fromType('Dendrogram', {
        // interactive: true,
        // alignLabels: true,
        // showLabels: true,
        // showLeafLabels: true,
        // padding: 0,
        // treeToCanvasRatio: 1,
        nodeColumnName: 'node',
      })) as DG.JsViewer;
      this.treeDn = this.tv.dockManager.dock(this.treeViewer, DG.DOCK_TYPE.RIGHT);
    }
  }
}
