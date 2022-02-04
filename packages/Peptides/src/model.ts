import * as DG from 'datagrok-api/dg';

import {describe} from './describe';
import {Subject, Observable} from 'rxjs';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';

/**
 * Model class for SAR viewers that retrieves and stores data.
 *
 * @class SARViewerModel
 */
export class PeptidesModel {
  // private _viewerGrid: DG.Grid;
  // private viewerVGrid: DG.Grid;
  // private _statsDf: DG.DataFrame;
  // private groupMapping: StringDictionary;
  private dataFrame: DG.DataFrame | null;
  private activityColumn: string | null;
  private activityScaling: string | null;
  private sourceGrid: DG.Grid | null;
  private twoColorMode: boolean | null;
  private initialBitset: DG.BitSet | null;
  private isUpdating: boolean = false;
  private grouping: boolean = false;
  private statsDataFrameSubject = new Subject<DG.DataFrame>();
  private sarGridSubject = new Subject<DG.Grid>();
  private sarVGridSubject = new Subject<DG.Grid>();
  private groupMappingSubject = new Subject<StringDictionary>();
  private static _modelName = 'peptidesModel';

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

    // this._statsDf = DG.DataFrame.create();
    // this._viewerGrid = DG.Grid.create(this.statsDf);
    // this.viewerVGrid = DG.Grid.create(this.statsDf);
    // this.groupMapping = {};

    // this.statsDataFrameObservable = new Observable(subject => subject.next(this.statsDf));
    // this.sarGridObservable = new Observable(subject => subject.next(this.viewerGrid));
    // this.sarVGridObservable = new Observable(subject => subject.next(this.viewerVGrid));
    // this.groupMappingObservable = new Observable(subject => subject.next(this.groupMapping));
  }

  // get statsDf() {
  //   return this._statsDf;
  // }

  // get viewerGrid() {
  //   return this._viewerGrid;
  // }

  get onStatsDataFrameChanged(): Observable<DG.DataFrame> {
    return this.statsDataFrameSubject.asObservable();
  }

  get onSARGridChanged(): Observable<DG.Grid> {
    return this.sarGridSubject.asObservable();
  }

  get onSARVGridChanged(): Observable<DG.Grid> {
    return this.sarVGridSubject.asObservable();
  }

  get onGroupMappingChanged(): Observable<StringDictionary> {
    return this.groupMappingSubject.asObservable();
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
      df: DG.DataFrame, activityCol: string, activityScaling: string, sourceGrid: DG.Grid, twoColorMode: boolean,
      initialBitset: DG.BitSet | null, grouping: boolean) {
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
    if (this.dataFrame && this.activityColumn && this.activityScaling && this.sourceGrid &&
        this.twoColorMode !== null && !this.isUpdating) {
      this.isUpdating = true;
      const [viewerGrid, viewerVGrid, statsDf, groupMapping] = await describe(
        this.dataFrame, this.activityColumn, this.activityScaling,
        this.sourceGrid, this.twoColorMode, this.initialBitset, this.grouping,
      );
      this.statsDataFrameSubject.next(statsDf);
      this.groupMappingSubject.next(groupMapping);
      this.sarGridSubject.next(viewerGrid);
      this.sarVGridSubject.next(viewerVGrid);

      
      this.isUpdating = false;
    }
  }

  static get modelName() {
    return PeptidesModel._modelName;
  }

  static getOrInit(dataFrame: DG.DataFrame): PeptidesModel {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel();
    return dataFrame.temp[PeptidesModel.modelName];
  }
}
