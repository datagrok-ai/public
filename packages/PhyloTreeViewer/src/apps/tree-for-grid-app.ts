import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Unsubscribable} from 'rxjs';
import {GridNeighbor} from '@datagrok-libraries/gridext/src/ui/GridNeighbor';

import {_package} from '../package';
import {TAGS} from '../utils/tree-helper';
import {data} from 'datagrok-api/grok';
import {injectTreeToGridUI} from '../viewers/inject-tree-to-grid';
import {injectTreeForGridUI} from '../viewers/inject-tree-for-grid';

export class TreeForGridApp {
  private viewed: boolean = false;
  private tableView: DG.TableView | null;
  gridN: GridNeighbor | null;

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
    const csv = await _package.files.readAsText('data/tree-gen-1000000.csv');
    const newick = await _package.files.readAsText('data/tree-gen-1000000.nwk');
    const leafColName = 'Leaf';

    const dataDf = DG.DataFrame.fromCsv(csv);
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

      this.gridN = injectTreeForGridUI(
        this.tableView.grid, this.newickStr, this.leafCol.name, 250);

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