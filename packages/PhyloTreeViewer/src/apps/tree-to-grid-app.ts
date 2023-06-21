import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import {GridWithTreeViewer} from '../viewers/grid-with-tree-viewer';
import {injectTreeToGridUI} from '../viewers/inject-tree-to-grid';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';
import {_package} from '../package';
import {TAGS as treeTAGS} from '@datagrok-libraries/bio/src/trees';

class AppView extends DG.ViewBase {

}

export class TreeToGridApp {
  private viewed: boolean = false;
  // viewer: GridWithTreeViewer | null;
  // viewerDn: DG.DockNode | null;
  tableView: DG.TableView | null = null;
  gridN: GridNeighbor | null = null;

  _dataDf: DG.DataFrame;
  _leafCol: DG.Column;
  _newickStr: string;

  get dataDf(): DG.DataFrame { return this._dataDf; }

  get leafCol(): DG.Column { return this._leafCol; }

  get newickStr(): string { return this._newickStr; }

  async init(): Promise<void> {
    await this.loadData();
  }

  async loadData(): Promise<void> {
    const dataPath: string = 'System:AppData/PhyloTreeViewer/data';
    // const csv = await _package.files.readAsText(`data/220524_FcRn_MH_FcRn_pepclust.csv`);
    // const newick = await _package.files.readAsText(`data/220524_FcRn_MH_FcRn_model_clustering.height.nwk`);

    const csv = await _package.files.readAsText(`data/tree95df.csv`);
    const newick = await _package.files.readAsText(`data/tree95.nwk`);
    const leafColName = 'id';

    // const csv = await _package.files.readAsText('data/tree-gen-100000.csv');
    // const newick = await _package.files.readAsText('data/tree-gen-100000.nwk');
    // const leafColName = 'Leaf';

    const dataDf = DG.DataFrame.fromCsv(csv);
    dataDf.setTag(treeTAGS.NEWICK, newick);
    dataDf.setTag(treeTAGS.NEWICK_LEAF_COL_NAME, leafColName);

    await this.setData(dataDf, newick);
  }

  async setData(dataDf: DG.DataFrame, newickStr: string): Promise<void> {
    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    this._dataDf = dataDf;
    const leafColName = dataDf.getTag(treeTAGS.NEWICK_LEAF_COL_NAME);
    if (!leafColName)
      throw new Error(`Specify leaf column name in DataFrame tag '${treeTAGS.NEWICK_LEAF_COL_NAME}'.`);
    this._leafCol = dataDf.getCol(leafColName);
    this._newickStr = newickStr;

    // const colClusterNum: DG.Column = this.dataDf.getCol('Cluster');
    // const colClusterStr: DG.Column = colClusterNum.convertTo('string', '#');
    // const colClusterIdx: number = this.dataDf.columns.toList().findIndex((col) => col == colClusterNum);
    // this.dataDf.columns.remove(colClusterNum.name);
    // this.dataDf.columns.insert(colClusterStr, colClusterIdx);

    if (!this.viewed) {
      await this.buildView();
      this.viewed = true;
    }
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
  }

  async buildView(): Promise<void> {
    if (!this.tableView) {
      this.tableView = grok.shell.addTableView(this.dataDf);
      //this.tableView.path = this.tableView.basePath = 'apps/PhyloTreeViewer/GridWithTree';

      // TODO: No cluster mark
      this.dataDf.columns.addNewInt('Cluster').init((rowI: number) => { return 0; });

      this.gridN = await injectTreeToGridUI(
        this.tableView.grid, this.newickStr, this.leafCol.name, 250,
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
