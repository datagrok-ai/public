import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';

import {_package} from '../package';
import {NodeType, parseNewick} from '@datagrok-libraries/bio';
import {TAGS as treeTAGS} from '@datagrok-libraries/bio/src/trees';

export class TreeCutAsTreeApp {
  private viewed: boolean = false;
  private tableView: DG.TableView | null;
  gridN: GridNeighbor | null;

  _dataDf: DG.DataFrame;
  _leafCol: DG.Column;
  _newickStr: string;
  _newickRoot: NodeType;

  get dataDf(): DG.DataFrame { return this._dataDf; }

  get leafCol(): DG.Column { return this._leafCol; }

  // get newickStr(): string { return this._newickStr; }

  get newickRoot(): NodeType { return this._newickRoot; }

  async init(): Promise<void> {
    await this.loadData();
  }

  async loadData(): Promise<void> {
    const csv = await _package.files.readAsText('data/tree-gen-10000.csv');
    const newick = await _package.files.readAsText('data/tree-gen-10000.nwk');
    const leafColName = 'Leaf';

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
    this._leafCol = dataDf.getCol(leafColName!);
    this._newickStr = newickStr;
    this._newickRoot = parseNewick(newickStr);

    if (!this.viewed) {
      await this.buildView();
      this.viewed = true;
    }
  }

  async destroyView() {
    if (this.gridN) {
      this.gridN.close();
      this.gridN = null;
    }

    if (this.tableView) {
      this.tableView.close();
      this.tableView = null;
    }
  }

  async buildView() {
    if (!this.tableView) {
      this.tableView = grok.shell.addTableView(this.dataDf);
      this.tableView.path = this.tableView.basePath = 'func/PhyloTreeViewer.treeForGridApp';

      this.dataDf.columns.addNewInt('Cluster').init((rowI) => { return null; });
      // this.gridN = injectTreeForGridUI(
      //   this.tableView.grid, this.newickRoot, this.leafCol.name, 250,
      //   {min: 0, max: 20, clusterColName: 'Cluster'});

      // const activityCol = this.dataDf.col('Activity');
      // if (activityCol) {
      //   this.tableView.filters({
      //     filters: [
      //       {type: 'histogram', column: 'Activity', label: 'Activity hist'}
      //     ]
      //   });
      // }
    }
  }
}
