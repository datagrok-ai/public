import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {hierarchicalClusteringUI} from '../utils/hierarchical-clustering';

export class HierarchicalClusteringApp {
  private viewed: boolean = false;

  private _df!: DG.DataFrame;
  get df(): DG.DataFrame { return this._df; }

  constructor() {}

  async init(df?: DG.DataFrame): Promise<void> {
    await this.loadData(df);
  }

  async loadData(df?: DG.DataFrame): Promise<void> {
    if (!df) {
      const csv: string = await _package.files.readAsText('data/demog-short.csv');
      df = DG.DataFrame.fromCsv(csv);
    }

    await this.setData(df);
  }

  async setData(df: DG.DataFrame) {
    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    this._df = df;

    if (!this.viewed) {
      await this.buildView();
      this.viewed = true;
    }
  }

  // -- View --

  tv?: DG.TableView;

  private async destroyView(): Promise<void> {
    if (this.tv) {
      this.tv.close();
      delete this.tv;
    }
  }

  private async buildView(): Promise<void> {
    if (!this.tv)
      this.tv = grok.shell.addTableView(this.df, DG.DOCK_TYPE.FILL);

    this.tv.path = this.tv.basePath = `/func/${_package.name}.hierarchicalClusteringApp`;
    hierarchicalClusteringUI(this.df, ['HEIGHT'], 'euclidean', 'ward');
  }
}
