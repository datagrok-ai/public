import * as DG from 'datagrok-api/dg';

import {describe} from './describe';
import {Subject} from 'rxjs';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';

/**
 * Model class for SAR viewers that retrieves and stores data.
 *
 * @class SARViewerModel
 */
class SARViewerModel {
  private viewerGrid: Subject<DG.Grid> = new Subject<DG.Grid>();
  private viewerVGrid: Subject<DG.Grid> = new Subject<DG.Grid>();
  private statsDf: Subject<DG.DataFrame> = new Subject<DG.DataFrame>();
  private groupMapping: Subject<StringDictionary> = new Subject<StringDictionary>();
  public viewerGrid$;
  public viewerVGrid$;
  public statsDf$;
  public groupMapping$;
  private dataFrame: DG.DataFrame | null;
  private activityColumn: string | null;
  private activityScaling: string | null;
  private sourceGrid: DG.Grid | null;
  private twoColorMode: boolean | null;
  private initialBitset: DG.BitSet | null;
  private isUpdating = false;
  grouping: boolean;

  /**
   * Creates an instance of SARViewerModel.
   *
   * @memberof SARViewerModel
   */
  constructor() {
    this.dataFrame = null;
    this.activityColumn = null;
    this.activityScaling = null;
    this.sourceGrid = null;
    this.twoColorMode = null;
    this.initialBitset = null;
    this.grouping = false;
    this.viewerGrid$ = this.viewerGrid.asObservable();
    this.viewerVGrid$ = this.viewerVGrid.asObservable();
    this.statsDf$ = this.statsDf.asObservable();
    this.groupMapping$ = this.groupMapping.asObservable();
  }

  /**
   * Updates data with using specified parameters.
   *
   * @param {DG.DataFrame} df Working table.
   * @param {string} activityCol Activity column name.
   * @param {string} activityScaling Activity scaling method.
   * @param {DG.Grid} sourceGrid Working table grid.
   * @param {boolean} twoColorMode Bidirectional analysis enabled.
   * @param {(DG.BitSet | null)} initialBitset Initial bitset.
   * @param {boolean} grouping Grouping enabled.
   * @memberof SARViewerModel
   */
  async updateData(
    df: DG.DataFrame,
    activityCol: string,
    activityScaling: string,
    sourceGrid: DG.Grid,
    twoColorMode: boolean,
    initialBitset: DG.BitSet | null,
    grouping: boolean,
  ) {
    this.dataFrame = df;
    this.activityColumn = activityCol;
    this.activityScaling = activityScaling;
    this.sourceGrid = sourceGrid;
    this.twoColorMode = twoColorMode;
    this.initialBitset = initialBitset;
    this.grouping = grouping;
    await this.updateDefault();
  }

  /**
   * Update data using current parameters.
   *
   * @memberof SARViewerModel
   */
  async updateDefault() {
    if (
      this.dataFrame && this.activityColumn && this.activityScaling &&
        this.sourceGrid && this.twoColorMode !== null && !this.isUpdating
    ) {
      this.isUpdating = true;
      const [viewerGrid, viewerVGrid, statsDf, groupMapping] = await describe(
        this.dataFrame, this.activityColumn, this.activityScaling,
        this.sourceGrid, this.twoColorMode, this.initialBitset, this.grouping,
      );
      this.viewerGrid.next(viewerGrid);
      this.viewerVGrid.next(viewerVGrid);
      this.statsDf.next(statsDf);
      this.groupMapping.next(groupMapping);
      this.isUpdating = false;
    }
  }
}

export const model = new SARViewerModel();
