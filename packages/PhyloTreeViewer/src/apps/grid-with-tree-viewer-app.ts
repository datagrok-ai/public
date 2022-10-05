import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {newickToDf} from './utils';
import {Unsubscribable} from 'rxjs';
import {GridWithTreeViewer} from './grid-with-tree-viewer';

class AppView extends DG.ViewBase {

}

export class GridWithTreeViewerApp {
  private view: DG.TableView | null;

  viewer: GridWithTreeViewer | null;
  viewerDn: DG.DockNode | null;

  _dataDf: DG.DataFrame;

  get dataDf(): DG.DataFrame { return this._dataDf; }

  async init(): Promise<void> {
    await this.loadData();
  }

  async loadData(): Promise<void> {
    const csv = await grok.dapi.files.readAsText('System:AppData/PhyloTreeViewer/data/tree95df.csv');
    const newick = await grok.dapi.files.readAsText('System:AppData/PhyloTreeViewer/data/tree95.nwk');

    const dataDf = DG.DataFrame.fromCsv(csv);
    dataDf.setTag('.newick', newick);

    await this.setData(dataDf);
  }

  async setData(dataDf: DG.DataFrame): Promise<void> {
    await this.destroyView();

    this._dataDf = dataDf;

    await this.buildView();
  }

  //# region -- View --

  viewSubs: Unsubscribable[] = [];

  async destroyView(): Promise<void> {
    if (this.viewer) {
      this.viewer.close();
      this.viewer.removeFromView();
    }

    if (this.viewerDn) {
      this.viewerDn.detachFromParent();
      this.viewer = null;
    }

    if (this.view) {
      this.view.close();
      this.view = null;
    }
  }

  async buildView(): Promise<void> {
    if (!this.view) {
      this.view = grok.shell.addTableView(this.dataDf);
      this.view.path = 'apps/PhyloTreeViewer/GridWithTreeViewer';
    }

    if (!this.viewer) {
      this.viewer = (await this.dataDf.plot.fromType('GridWithTree', {})) as GridWithTreeViewer;
    }

    if (!this.viewerDn) {
      this.viewerDn = this.view?.dockManager.dock(this.viewer, DG.DOCK_TYPE.LEFT)!;
    }
  }

  //# endregion -- View --
}