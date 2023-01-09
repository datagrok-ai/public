import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';

import {TAGS, TreeHelper} from '../utils/tree-helper';
import {injectTreeForGridUI2} from '../viewers/inject-tree-for-grid2';
import {ITreeHelper, NodeType, parseNewick} from '@datagrok-libraries/bio';
import {_package} from '../package';

export class TreeForGridFilterApp {
  private th!: ITreeHelper;

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
    this.th = new TreeHelper();
    await this.loadData();
  }

  async loadData(): Promise<void> {
    const tree = this.th.generateTree(5000);
    const newick = this.th.toNewick(tree);
    const leafColName = 'Leaf';

    const leafList = this.th.getLeafList(tree);
    const leafCol: DG.Column = DG.Column.fromList(DG.COLUMN_TYPE.STRING, leafColName,
      leafList.map((n) => n.name));
    const activityCol: DG.Column = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Activity',
      leafList.map((n) => Math.random()));
    const dataDf = DG.DataFrame.fromColumns([leafCol, activityCol]);

    // const csv = await _package.files.readAsText('data/tree-gen-100000.csv');
    // const newick = await _package.files.readAsText('data/tree-gen-100000.nwk');
    // const leafColName = 'Leaf';
    // const dataDf = DG.DataFrame.fromCsv(csv);

    dataDf.setTag(TAGS.DF_NEWICK, newick);
    dataDf.setTag(TAGS.DF_NEWICK_LEAF_COL_NAME, leafColName);

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

      this.tableView = grok.shell.addTableView(dataDf);
      this.tableView.path = this.tableView.basePath = `/func/${_package.name}.treeForGridFilterApp`;

      this.gridN = injectTreeForGridUI2(
        this.tableView.grid, this.newickRoot, this.leafCol.name, 300);

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
