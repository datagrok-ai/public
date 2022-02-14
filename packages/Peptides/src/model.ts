import * as DG from 'datagrok-api/dg';

import {describe} from './describe';
import {Subject, Observable} from 'rxjs';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {addViewerToHeader, StackedBarChart} from './viewers/stacked-barchart-viewer';

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
  private substFlag = false;
  private statsDataFrameSubject = new Subject<DG.DataFrame>();
  private sarGridSubject = new Subject<DG.Grid>();
  private sarVGridSubject = new Subject<DG.Grid>();
  private groupMappingSubject = new Subject<StringDictionary>();
  private substFlagSubject = new Subject<boolean>();
  private static _modelName = 'peptidesModel';

  private constructor(dataFrame: DG.DataFrame) {
    this.dataFrame = dataFrame;
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

  get onSubstFlagChanged(): Observable<boolean> {
    return this.substFlagSubject.asObservable();
  }

  async updateData(
    df: DG.DataFrame | null, activityCol: string | null, activityScaling: string | null, sourceGrid: DG.Grid | null,
    twoColorMode: boolean | null, initialBitset: DG.BitSet | null, grouping: boolean | null) {
    this.dataFrame = df ?? this.dataFrame;
    this.activityColumn = activityCol ?? this.activityColumn;
    this.activityScaling = activityScaling ?? this.activityScaling;
    this.sourceGrid = sourceGrid ?? this.sourceGrid;
    this.twoColorMode = twoColorMode ?? this.twoColorMode;
    this.initialBitset = initialBitset ?? this.initialBitset;
    this.grouping = grouping ?? this.grouping;
    await this.updateDefault();
  }

  async updateDefault() {
    if (this.dataFrame && this.activityColumn && this.activityScaling && this.sourceGrid &&
        this.twoColorMode !== null && !this.isUpdating) {
      this.isUpdating = true;
      const [viewerGrid, viewerVGrid, statsDf, groupMapping] = await describe(
        this.dataFrame, this.activityColumn, this.activityScaling, this.sourceGrid, this.twoColorMode,
        this.initialBitset, this.grouping);
      this.statsDataFrameSubject.next(statsDf);
      this.groupMappingSubject.next(groupMapping);
      this.sarGridSubject.next(viewerGrid);
      this.sarVGridSubject.next(viewerVGrid);
      this.substFlag = !this.substFlag;
      this.substFlagSubject.next(this.substFlag);

      this.sourceGrid.invalidate();

      this.isUpdating = false;
    }

    await this.updateBarchart();
  }

  async updateBarchart() {
    const stackedBarchart = await this.dataFrame?.plot.fromType('StackedBarChartAA') as StackedBarChart;
    if (stackedBarchart && this.sourceGrid)
      addViewerToHeader(this.sourceGrid, stackedBarchart);
  }

  static get modelName() {
    return PeptidesModel._modelName;
  }

  static getOrInit(dataFrame: DG.DataFrame): PeptidesModel {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    return dataFrame.temp[PeptidesModel.modelName];
  }
}
