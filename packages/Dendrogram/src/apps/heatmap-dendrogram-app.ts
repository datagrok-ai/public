import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {hierarchicalClusteringUI} from '../utils/hierarchical-clustering';
import {DistanceMetric} from '@datagrok-libraries/bio/src/trees';

export class HeatmapDedndrogramApp {
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

    this.tv.path = this.tv.basePath = `/func/${_package.name}.heatmapDendrogramApp`;
    await hierarchicalClusteringUI(this.df, ['AGE'], DistanceMetric.Euclidean, 'ward');
    const adjustTreeHeight = () => {
      const rowCount = this.tv!.dataFrame.filter.trueCount;
      const rowsGridHeight: number = this.tv!.grid.root.clientHeight - this.tv!.grid.colHeaderHeight;
      this.tv!.grid.props.rowHeight = rowsGridHeight / (rowCount + 1);
    };
    this.tv!.grid.onBeforeDrawContent.subscribe(() => {
      adjustTreeHeight();
    });
    adjustTreeHeight();
    // tv.grid.props.showRowHeader = false;
    this.tv!.grid.props.isGrid = false;
    this.tv!.grid.props.isHeatmap = true;
    this.tv!.grid.props.showRowHeader = false;
    // this.tv!.grid.props.showAddNewRowIcon = false;
    setTimeout(() => {
      adjustTreeHeight();
      this.tv!.grid.invalidate();
    }, 10);
  }
}
