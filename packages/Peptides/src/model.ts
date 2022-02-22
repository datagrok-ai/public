import * as DG from 'datagrok-api/dg';

import {describe} from './describe';
import {Subject, Observable} from 'rxjs';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {addViewerToHeader, StackedBarChart} from './viewers/stacked-barchart-viewer';

export class PeptidesModel {
  private _dataFrame: DG.DataFrame;
  private _activityColumn: string | null;
  private _activityScaling: string | null;
  private _sourceGrid: DG.Grid | null;
  private _twoColorMode: boolean | null;
  private _initialBitset: DG.BitSet | null;
  private _grouping: boolean = false;
  private _isUpdating: boolean = false;
  private _substFlag = false;
  private _statsDataFrameSubject = new Subject<DG.DataFrame>();
  private _sarGridSubject = new Subject<DG.Grid>();
  private _sarVGridSubject = new Subject<DG.Grid>();
  private _groupMappingSubject = new Subject<StringDictionary>();
  private _substFlagSubject = new Subject<boolean>();
  private static _modelName = 'peptidesModel';

  private constructor(dataFrame: DG.DataFrame) {
    this._dataFrame = dataFrame;
    this._activityColumn = null;
    this._activityScaling = null;
    this._sourceGrid = null;
    this._twoColorMode = null;
    this._initialBitset = null;
  }
    
  static getInstance(dataFrame: DG.DataFrame): PeptidesModel {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    return dataFrame.temp[PeptidesModel.modelName];
  }

  get dataFrame(): DG.DataFrame {
    return this._dataFrame;
  }

  get onStatsDataFrameChanged(): Observable<DG.DataFrame> {
    return this._statsDataFrameSubject.asObservable();
  }

  get onSARGridChanged(): Observable<DG.Grid> {
    return this._sarGridSubject.asObservable();
  }

  get onSARVGridChanged(): Observable<DG.Grid> {
    return this._sarVGridSubject.asObservable();
  }

  get onGroupMappingChanged(): Observable<StringDictionary> {
    return this._groupMappingSubject.asObservable();
  }

  get onSubstFlagChanged(): Observable<boolean> {
    return this._substFlagSubject.asObservable();
  }

  async updateData(activityCol: string | null, activityScaling: string | null, sourceGrid: DG.Grid | null,
    twoColorMode: boolean | null, initialBitset: DG.BitSet | null, grouping: boolean | null) {
    this._activityColumn = activityCol ?? this._activityColumn;
    this._activityScaling = activityScaling ?? this._activityScaling;
    this._sourceGrid = sourceGrid ?? this._sourceGrid;
    this._twoColorMode = twoColorMode ?? this._twoColorMode;
    this._initialBitset = initialBitset ?? this._initialBitset;
    this._grouping = grouping ?? this._grouping;
    await this.updateDefault();
  }

  async updateDefault() {
    if (
      this._activityColumn && this._activityScaling && this._sourceGrid && this._twoColorMode !== null &&
      !this._isUpdating
    ) {
      this._isUpdating = true;
      const [viewerGrid, viewerVGrid, statsDf, groupMapping] = await describe(
        this._dataFrame, this._activityColumn, this._activityScaling, this._sourceGrid, this._twoColorMode,
        this._initialBitset, this._grouping);
      this._statsDataFrameSubject.next(statsDf);
      this._groupMappingSubject.next(groupMapping);
      this._sarGridSubject.next(viewerGrid);
      this._sarVGridSubject.next(viewerVGrid);
      this._substFlag = !this._substFlag;
      this._substFlagSubject.next(this._substFlag);

      this._sourceGrid.invalidate();

      this._isUpdating = false;
    }

    await this.updateBarchart();
  }

  async updateBarchart() {
    const stackedBarchart = await this._dataFrame?.plot.fromType('StackedBarChartAA') as StackedBarChart;
    if (stackedBarchart && this._sourceGrid)
      addViewerToHeader(this._sourceGrid, stackedBarchart);
  }

  static get modelName() {
    return PeptidesModel._modelName;
  }
}
