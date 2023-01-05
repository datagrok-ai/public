import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';

import {_package} from '../package';
import {TAGS} from '../utils/tree-helper';
import {injectTreeForGridUI2} from '../viewers/inject-tree-for-grid2';
import {NodeType, parseNewick, TreeCutOptions} from '@datagrok-libraries/bio';

export class TreeForGridApp {
  private viewed: boolean = false;
  private tableView: DG.TableView | null = null;
  gridN: GridNeighbor | null = null;

  _dataDf!: DG.DataFrame;
  _leafCol!: DG.Column;
  _newickStr!: string;
  _newickRoot!: NodeType;

  get dataDf(): DG.DataFrame { return this._dataDf; }

  get leafCol(): DG.Column { return this._leafCol; }

  // get newickStr(): string { return this._newickStr; }

  get newickRoot(): NodeType { return this._newickRoot; }

  async init(): Promise<void> {
    await this.loadData();
  }

  async loadData(): Promise<void> {
    // const csv = await _package.files.readAsText('data/tree-gen-100000.csv');
    // const newick = await _package.files.readAsText('data/tree-gen-100000.nwk');
    // const leafColName = 'Leaf';

    const csv = await _package.files.readAsText('data/tree95df.csv');
    const newick = await _package.files.readAsText('data/tree95.nwk');
    const leafColName = 'id';

    const dataDf = DG.DataFrame.fromCsv(csv);
    dataDf.setTag(TAGS.DF_NEWICK, newick);
    dataDf.setTag(TAGS.DF_NEWICK_LEAF_COL_NAME, leafColName);
    // For debug purposes numbering rows with index
    dataDf.columns.addNewVirtual('[index]', (rowI: number) => { return rowI; }, DG.TYPE.INT);

    await this.setData(dataDf, newick);
  }

  async setData(dataDf: DG.DataFrame, newickStr: string): Promise<void> {
    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    this._dataDf = dataDf;
    const leafColName = dataDf.getTag(TAGS.DF_NEWICK_LEAF_COL_NAME);
    if (!leafColName)
      throw new Error(`Specify leaf column name in DataFrame tag '${TAGS.DF_NEWICK_LEAF_COL_NAME}'.`);
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
      const dataDf: DG.DataFrame = this.dataDf;
      dataDf.columns.addNewInt('Cluster').init((rowI) => { return null; });

      const clusterDf: DG.DataFrame = DG.DataFrame.create(0);
      clusterDf.columns.addNewInt('Cluster');
      clusterDf.columns.addNewString(this.leafCol.name);
      clusterDf.columns.addNewInt(`${this.leafCol.name}_Count`);

      this.tableView = grok.shell.addTableView(dataDf);
      this.tableView.path = this.tableView.basePath = `/func/${_package.name}.treeForGridApp`;

      const cutOpts: TreeCutOptions = {
        min: 0, max: 20, clusterColName: 'Cluster',
        dataDf: dataDf, clusterDf: clusterDf
      };
      this.gridN = injectTreeForGridUI2(
        this.tableView.grid, this.newickRoot, this.leafCol.name, 250
        /* cutOpts */);

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
