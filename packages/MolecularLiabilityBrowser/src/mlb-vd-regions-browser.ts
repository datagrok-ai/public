import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {IVdRegionsViewer, VdRegion} from '@datagrok-libraries/bio/src/vd-regions';
import {MlbDataFrame} from './types/dataframe';
import {Unsubscribable} from 'rxjs';

export class MlbVdRegionsBrowser {
  private viewed: boolean = false;

  private _mlbDf: MlbDataFrame;
  get mlbDf(): MlbDataFrame { return this._mlbDf; }

  private _regions: VdRegion[] = [];
  get regions(): VdRegion[] { return this._regions; }

  private viewer: IVdRegionsViewer;

  get root(): HTMLElement { return this.viewer!.root; }

  constructor() {}

  async init(): Promise<void> {
    this._mlbDf = MlbDataFrame.Empty;
    this.viewer = await this._mlbDf.plot.fromType('VdRegions', {
      skipEmptyPositions: true
    }) as unknown as IVdRegionsViewer;
  }

  async setData(mlbDf: MlbDataFrame, regions: VdRegion[]): Promise<void> {
    if (this.viewed) {
      await this.destroyView();
      this.viewed = false;
    }

    this._mlbDf = mlbDf;
    this._regions = regions;

    this.viewSubs = [];

    if (!this.viewed) {
      await this.buildView();
      this.viewed = true;
    }
  }

  // -- View --

  private viewSubs: Unsubscribable[] = [];

  async destroyView(): Promise<void> {
    await this.viewer.setData(MlbDataFrame.Empty, []);
  }

  async buildView(): Promise<void> {
    await this.viewer.setData(this.mlbDf, this.regions);
  }
}