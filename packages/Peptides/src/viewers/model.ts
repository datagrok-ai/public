import * as DG from 'datagrok-api/dg';

import {describe} from '../describe';
import { Subject } from 'rxjs';

class SARViewerModel {
  private viewerGrid: Subject<DG.Grid> = new Subject<DG.Grid>();
  private viewerVGrid: Subject<DG.Grid> = new Subject<DG.Grid>();
  private statsDf: Subject<DG.DataFrame> = new Subject<DG.DataFrame>();
  public viewerGrid$ = this.viewerGrid.asObservable();
  public viewerVGrid$ = this.viewerVGrid.asObservable();
  public statsDf$ = this.statsDf.asObservable();
  private dataFrame: DG.DataFrame | null;
  private activityColumn: string | null;
  private activityScaling: string | null;
  private sourceGrid: DG.Grid | null;
  private twoColorMode: boolean | null;
  private initialBitset: DG.BitSet | null;
  private isUpdating = false;

  constructor() {
    this.dataFrame = null;
    this.activityColumn = null;
    this.activityScaling = null;
    this.sourceGrid = null;
    this.twoColorMode = null;
    this.initialBitset = null;
  }

  async updateData(
      df: DG.DataFrame,
      activityCol: string,
      activityScaling: string,
      sourceGrid: DG.Grid,
      twoColorMode: boolean,
      initialBitset: DG.BitSet | null,
    ) {
    this.dataFrame = df;
    this.activityColumn = activityCol;
    this.activityScaling = activityScaling;
    this.sourceGrid = sourceGrid;
    this.twoColorMode = twoColorMode;
    this.initialBitset = initialBitset;
    await this.updateDefault();
  }

  async updateDefault() {
    if (this.dataFrame && this.activityColumn && this.activityScaling && this.sourceGrid && this.twoColorMode !== null && !this.isUpdating) {
      this.isUpdating = true;
      const [viewerGrid, viewerVGrid, statsDf] = await describe(
        this.dataFrame, this.activityColumn, this.activityScaling,
        this.sourceGrid, this.twoColorMode, this.initialBitset
      );
      this.viewerGrid.next(viewerGrid);
      this.viewerVGrid.next(viewerVGrid);
      this.statsDf.next(statsDf);
      this.isUpdating = false;
    }
  }
}

export let model = new SARViewerModel();
  