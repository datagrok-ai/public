import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import {GridWithTreeViewer} from '../viewers/grid-with-tree-viewer';
import {injectTreeToGridUI} from '../viewers/inject-tree-to-grid';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import {_package} from '../package';

class AppView extends DG.ViewBase {

}

export class GridWithTreeApp {
  private view: DG.TableView | null;

  // viewer: GridWithTreeViewer | null;
  // viewerDn: DG.DockNode | null;
  tableView: DG.TableView | null = null;
  gridN: GridNeighbor | null = null;

  _dataDf: DG.DataFrame;
  _newickText: string;

  get dataDf(): DG.DataFrame { return this._dataDf; }

  get newickText(): string { return this._newickText; }

  async init(): Promise<void> {
    await this.loadData();
  }

  async loadData(): Promise<void> {
    const dataPath: string = 'System:AppData/PhyloTreeViewer/data';
    // const csv = await _package.files.readAsText(`data/220524_FcRn_MH_FcRn_pepclust.csv`);
    // const newick = await _package.files.readAsText(`data/220524_FcRn_MH_FcRn_model_clustering.height.nwk`);

    const csv = await _package.files.readAsText(`data/tree95df.csv`);
    const newick = await _package.files.readAsText(`data/tree95.nwk`);

    const dataDf = DG.DataFrame.fromCsv(csv);
    dataDf.setTag('.newick', newick);

    await this.setData(dataDf, newick);
  }

  async setData(dataDf: DG.DataFrame, newickText: string): Promise<void> {
    await this.destroyView();

    this._dataDf = dataDf;
    this._newickText = newickText;

    // const colClusterNum: DG.Column = this.dataDf.getCol('Cluster');
    // const colClusterStr: DG.Column = colClusterNum.convertTo('string', '#');
    // const colClusterIdx: number = this.dataDf.columns.toList().findIndex((col) => col == colClusterNum);
    // this.dataDf.columns.remove(colClusterNum.name);
    // this.dataDf.columns.insert(colClusterStr, colClusterIdx);

    await this.buildView();
  }

  // region -- View --

  viewSubs: Unsubscribable[] = [];

  async destroyView(): Promise<void> {
    if (this.gridN) {
      this.gridN.close();
      this.gridN = null;
    }

    if (this.tableView) {
      this.tableView.close();
      this.tableView = null;
    }
    // if (this.viewer) {
    //   this.viewer.close();
    //   this.viewer.removeFromView();
    // }
    //
    // if (this.viewerDn) {
    //   this.viewerDn.detachFromParent();
    //   this.viewer = null;
    // }
    //
    // if (this.view) {
    //   this.view.close();
    //   this.view = null;
    // }
  }

  async buildView(): Promise<void> {
    if (!this.tableView) {
      this.tableView = grok.shell.addTableView(this.dataDf);
      this.tableView.path = this.tableView.basePath = 'apps/PhyloTreeViewer/GridWithTree';

      this.dataDf.columns.addNewInt('Cluster').init((rowI: number) => { return 0;});

      this.gridN = injectTreeToGridUI(
        this.tableView.grid, this.newickText, 'id', 250,
        {min: 0, max: 1, clusterColName: 'Cluster'});

      // this.tableView.filters({
      //   filters: [
      //     {type: 'categorical', column: 'Cluster', label: 'Cluster filter',},
      //   ]
      // });

      this.tableView.filters({
        filters: [
          // {type: 'categorical', column: 'value', label: 'value cats'},
          {type: 'histogram', column: 'value', label: 'value hist'},
        ]
      });
    }
    // if (!this.view) {
    //   this.view = grok.shell.addTableView(this.dataDf);
    //   this.view.path = 'apps/PhyloTreeViewer/GridWithTreeViewer';
    // }
    //
    // if (!this.viewer) {
    //   this.viewer = (await this.dataDf.plot.fromType('GridWithTree', {
    //     nodeNameColumnName: 'id',
    //   })) as GridWithTreeViewer;
    // }
    //
    //
    // // if (!this.viewerDn) {
    // //   this.viewerDn = this.view?.dockManager.dock(this.viewer, DG.DOCK_TYPE.FILL)!;
    // // }
  }
}