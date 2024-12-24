/* eslint-disable max-len */
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';
import {ACTIVITY_TARGET, COLUMN_NAME, COLUMNS_NAMES} from '../utils/constants';
import wu from 'wu';
import $ from 'cash-dom';
import * as rxjs from 'rxjs';
import {filter} from 'rxjs/operators';
export const enum ClusterMaxActivityProps {
    CLUSTER_COLUMN = 'cluster',
    ACTIVITY_COLUMN = 'activity',
    CONNECTIVITY_COLUMN = 'connectivity',
    COLOR_COLUMN = 'color',
    CLUSTER_SIZE_THRESHOLD = 'clusterSizeThreshold',
    ACTIVITY_THRESHOLD = 'activityThreshold',
}

export interface IClusterMaxActivity {
    clusterColumnName: string;
    activityColumnName: string;
    connectivityColumnName?: string;
    colorColumnName?: string;
    activityTarget: ACTIVITY_TARGET;
    clusterSizeThreshold: number;
    activityThreshold: number;
}


export class ClusterMaxActivityViewer extends DG.JsViewer implements IClusterMaxActivity {
  _titleHost = ui.divText(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY, {id: 'pep-viewer-title', style: {marginRight: 'auto'}});
  _selsectIcon: HTMLElement = ui.div();
  clusterColumnName: string;
  activityColumnName: string;
  colorColumnName: string;
  connectivityColumnName?: string | undefined;
  activityTarget: ACTIVITY_TARGET = ACTIVITY_TARGET.HIGH;
  _scViewer?: DG.ScatterPlotViewer | null;
  viewerError: string = '';
  renderTimeout: NodeJS.Timeout | number | null = null;
  renderDebounceTime = 500;
  clusterSizeThreshold: number;
  activityThreshold: number;
  static clusterSizeColName = '~cluster.size' as const;
  static maxActivityInClusterColName = '~max.activity.for.cluster' as const;
  static maxConnectivityInClusterColName = '~max.connectivity.for.cluster' as const;
  static synSelectionColName = 'Syn Selection' as const;
  static maxActivityLabel = 'Max Activity' as const;
  static maxConnectivityLabel = 'Max Connectivity' as const;
  private scFilterQuery = `\$\{${ClusterMaxActivityViewer.maxActivityInClusterColName}\} == 1` as const;
  private selectionSubscription: rxjs.Subscription | null = null;
  private linesDrawSubscription: rxjs.Subscription | null = null;
  get scViewer(): DG.ScatterPlotViewer | null {
    if (!this._scViewer)
      this._scViewer = this.createSCViewer();
    return this._scViewer;
  }
  constructor() {
    super();
    this.clusterColumnName = this.column(ClusterMaxActivityProps.CLUSTER_COLUMN, {nullable: false});
    this.activityColumnName = this.column(ClusterMaxActivityProps.ACTIVITY_COLUMN, {nullable: false});
    this.colorColumnName = this.column(ClusterMaxActivityProps.COLOR_COLUMN,
      {nullable: true, defaultValue: null});
    this.activityTarget = this.string(
      'activityTarget', ACTIVITY_TARGET.HIGH, {choices: [ACTIVITY_TARGET.HIGH, ACTIVITY_TARGET.LOW]},
    ) as ACTIVITY_TARGET;
    this.clusterSizeThreshold = this.int(ClusterMaxActivityProps.CLUSTER_SIZE_THRESHOLD, 20);
    this.activityThreshold = this.int(ClusterMaxActivityProps.ACTIVITY_THRESHOLD, 1000);
    this.connectivityColumnName = this.column(ClusterMaxActivityProps.CONNECTIVITY_COLUMN, {nullable: true});
  }

  /**
   * Returns PeptidesModel instance that belongs to the attached dataframe.
   * @return - PeptidesModel instance.
   */
  get model(): PeptidesModel {
    return PeptidesModel.getInstance(this.dataFrame);
  }

  private createSCViewer(): DG.ScatterPlotViewer | null {
    // @ts-ignore TODO: fix after api update
    const scatterPlotProps: Partial<DG.IScatterPlotSettings> & Options = {
      showXAxis: true,
      showYAxis: true,
      showXSelector: false,
      showYSelector: false,
      showColorSelector: false,
      xAxisType: DG.AxisType.logarithmic,
      yAxisType: DG.AxisType.logarithmic,
      invertYAxis: this.activityTarget === ACTIVITY_TARGET.LOW,
      xColumnName: ClusterMaxActivityViewer.clusterSizeColName,
      markerType: 'circle',
      markerDefaultSize: 10,
      showSizeSelector: false,
    };
    if (this.clusterColumnName == null || this.activityColumnName == null ||
         !this.dataFrame.columns.contains(this.clusterColumnName) ||
         !this.dataFrame.columns.contains(this.activityColumnName)
    ) {
      this.viewerError = 'Please set valid cluster and activity columns';
      return null;
    }
    const activityCol = this.dataFrame.columns.byName(this.activityColumnName);
    const clusterCol = this.dataFrame.columns.byName(this.clusterColumnName);
    const connectivityCol = this.connectivityColumnName != null ?
      this.dataFrame.columns.byName(this.connectivityColumnName) : null;
    const numericColTypes: DG.ColumnType[] =
        [DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.BIG_INT, DG.COLUMN_TYPE.QNUM];
    if (!numericColTypes.includes(activityCol.type)) {
      this.viewerError = 'Activity column should be numeric';
      return null;
    }
    const clusterSizeCol = this.dataFrame.columns.getOrCreate(ClusterMaxActivityViewer.clusterSizeColName,
      DG.TYPE.INT);
    const clusterSizeMap: {[key: number | string]: number} = {};
    for (let i = 0; i < this.dataFrame.rowCount; i++) {
      const cluster: string | number = clusterCol.get(i);
      if (cluster == null || cluster == '-1')//fyi -1 == '-1' (with double ==)
        continue;
      clusterSizeMap[cluster] = (clusterSizeMap[cluster] ?? 0) + 1;
    }

    clusterSizeCol.init((i) => clusterCol.isNone(i) ? null : clusterSizeMap[clusterCol.get(i)] ?? null);

    if (activityCol.stats.min <= 0)
      scatterPlotProps.yAxisType = DG.AxisType.linear;

    // create a new column to store max activity for each cluster size
    const maxActivityIndexPerClusterMap: {[key: number]: number} = {};
    const maxConnectivityIndexPerClusterMap: {[key: number]: number} = {};
    for (let i = 0; i < this.dataFrame.rowCount; i++) {
      const cluster: number | null = clusterCol.get(i);
      if (cluster == null || clusterCol.isNone(i) || activityCol.isNone(i))
        continue;
      const activity: number = activityCol.get(i);
      const prevMaxActivityIndex = maxActivityIndexPerClusterMap[cluster];
      if (prevMaxActivityIndex == null || prevMaxActivityIndex == undefined)
        maxActivityIndexPerClusterMap[cluster] = i;
      else {
        if (activity > activityCol.get(prevMaxActivityIndex) && this.activityTarget === ACTIVITY_TARGET.HIGH)
          maxActivityIndexPerClusterMap[cluster] = i;
        else if (activity < activityCol.get(prevMaxActivityIndex) && this.activityTarget === ACTIVITY_TARGET.LOW)
          maxActivityIndexPerClusterMap[cluster] = i;
      }
      if (connectivityCol) {
        const connectivity: number = connectivityCol.get(i);
        const prevMaxConnectivityIndex = maxConnectivityIndexPerClusterMap[cluster];
        if (prevMaxConnectivityIndex == null || prevMaxConnectivityIndex == undefined)
          maxConnectivityIndexPerClusterMap[cluster] = i;
        else {
          if (connectivity > connectivityCol.get(prevMaxConnectivityIndex))
            maxConnectivityIndexPerClusterMap[cluster] = i;
        }
      }
    }

    const maxAtivityInClusterSizeCol = this.dataFrame.columns.getOrCreate(
      ClusterMaxActivityViewer.maxActivityInClusterColName, DG.COLUMN_TYPE.INT);
    maxAtivityInClusterSizeCol.init((i) => {
      if (clusterCol.isNone(i))
        return 0;
      return i === maxActivityIndexPerClusterMap[clusterCol.get(i)] ? 1 : 0;
    });

    const maxConnectivityInClusterSizeCol = this.dataFrame.columns.getOrCreate(
      ClusterMaxActivityViewer.maxConnectivityInClusterColName, DG.COLUMN_TYPE.INT);
    maxConnectivityInClusterSizeCol.init((i) => {
      if (clusterCol.isNone(i))
        return 0;
      return i === maxConnectivityIndexPerClusterMap[clusterCol.get(i)] ? 1 : 0;
    });

    const synSelectionCol = this.dataFrame.columns.getOrCreate(
      ClusterMaxActivityViewer.synSelectionColName, DG.TYPE.STRING);

    synSelectionCol.init((i) => {
      if (clusterCol.isNone(i))
        return null;
      let r: string | null = null;
      if (i === maxActivityIndexPerClusterMap[clusterCol.get(i)])
        r = ClusterMaxActivityViewer.maxActivityLabel;
      if (connectivityCol && i === maxConnectivityIndexPerClusterMap[clusterCol.get(i)])
        r = r ? `${r}, ${ClusterMaxActivityViewer.maxConnectivityLabel}` : ClusterMaxActivityViewer.maxConnectivityLabel;
      return r;
    });

    scatterPlotProps.xColumnName = ClusterMaxActivityViewer.clusterSizeColName;
    scatterPlotProps.yColumnName = this.activityColumnName;
    scatterPlotProps.filter = this.scFilterQuery;
    this.viewerError = '';
    const sc = DG.Viewer.scatterPlot(this.dataFrame, scatterPlotProps);

    if (this.selectionSubscription)
      this.selectionSubscription.unsubscribe();
    this.selectionSubscription = sc.onDataEvent.pipe(filter((e) => e.type == 'd4-select')).subscribe((e) => {
      const indexes = e.bitset?.getSelectedIndexes() ?? [];
      const currentSelection = this.dataFrame.selection;
      const filterBitset = DG.BitSet.create(this.dataFrame.rowCount, (i) => maxAtivityInClusterSizeCol.get(i) == 1);
      filterBitset.and(currentSelection);
      for (let i = 0; i < indexes.length; i++) {
        const index = indexes[i];
        const cluster = clusterCol.get(index);
        if (clusterCol.isNone(index))
          continue;
        filterBitset.set(index, true);
        if (maxConnectivityIndexPerClusterMap[cluster] != null)
          filterBitset.set(maxConnectivityIndexPerClusterMap[cluster], true);
      }
      filterBitset.fireChanged();
      const filteredIndexes = filterBitset.getSelectedIndexes();
      for (const i of filteredIndexes) {
        const cluster = clusterCol.get(i);
        if (cluster == null)
          continue;
        if (maxConnectivityIndexPerClusterMap[cluster] != null)
          filterBitset.set(maxConnectivityIndexPerClusterMap[cluster], true);
      }
      setTimeout(() => {
        this.dataFrame.selection.copyFrom(filterBitset, true);
        setTimeout(() => {
          if (this.model)
            this.model.createAccordion();
        }, 200);
      }, 200);
    });

    const selectTopQuadrants = (): void => {
      const selectionBitset = DG.BitSet.create(this.dataFrame.rowCount);
      Object.entries(maxActivityIndexPerClusterMap).forEach(([cluster, index]) => {
        if (cluster == null || index == null)
          return;
        const clusterInt = parseInt(cluster);
        const activity = activityCol.get(index);
        const clusterSize = clusterSizeMap[clusterInt] ?? clusterSizeMap[cluster];
        if (activity < this.activityThreshold && clusterSize < this.clusterSizeThreshold)
          return;
        selectionBitset.set(index, true, false);
        if (maxConnectivityIndexPerClusterMap[clusterInt] != null)
          selectionBitset.set(maxConnectivityIndexPerClusterMap[clusterInt], true, false);
      });
      selectionBitset.fireChanged();
      this.dataFrame.selection.copyFrom(selectionBitset, true);
    };

    this._selsectIcon = ui.iconSvg('select-all', () => {
      selectTopQuadrants();
    }, 'Select 3 Active quadrants');
    this._selsectIcon.style.cursor = 'pointer';
    this._selsectIcon.style.marginRight = '5px';
    selectTopQuadrants();
    if (this.linesDrawSubscription)
      this.linesDrawSubscription.unsubscribe();
    this.linesDrawSubscription = sc.onBeforeDrawScene.subscribe(() => {
      const canvas = sc.getInfo().canvas;
      const ctx: CanvasRenderingContext2D = canvas.getContext('2d');
      const viewPort = sc.viewport;
      const startPointHor = sc.worldToScreen(viewPort.x, this.activityThreshold);
      const startPointVer = sc.worldToScreen(this.clusterSizeThreshold, viewPort.y);
      const endPointHor = sc.worldToScreen(viewPort.x + viewPort.width, this.activityThreshold);
      const endPointVer = sc.worldToScreen(this.clusterSizeThreshold, viewPort.y + viewPort.height);
      ctx.beginPath();
      ctx.strokeStyle = 'rgb(0,0,0)';
      ctx.lineWidth = 1;
      ctx.moveTo(startPointHor.x, startPointHor.y);
      ctx.lineTo(endPointHor.x, endPointHor.y);
      ctx.moveTo(startPointVer.x, startPointVer.y);
      ctx.lineTo(endPointVer.x, endPointVer.y);
      ctx.stroke();
      ctx.closePath();
    });

    return sc;
  }

  onTableAttached(): void {
    super.onTableAttached();
    const activityCol: DG.Column | null = this.dataFrame?.col(COLUMNS_NAMES.ACTIVITY) ??
        wu(this.dataFrame?.columns.numerical).next()?.value;
    if (activityCol != null)
      this.getProperty(`${ClusterMaxActivityProps.ACTIVITY_COLUMN}${COLUMN_NAME}`)?.set(this, activityCol.name);

    const clusterCol: DG.Column | null = wu(this.dataFrame?.columns.categorical).next()?.value;
    if (clusterCol != null)
      this.getProperty(`${ClusterMaxActivityProps.CLUSTER_COLUMN}${COLUMN_NAME}`)?.set(this, clusterCol.name);

    const connectivityCol: DG.Column | null = wu(this.dataFrame?.columns.numerical).next()?.value;
    if (connectivityCol != null)
      this.getProperty(`${ClusterMaxActivityProps.CONNECTIVITY_COLUMN}${COLUMN_NAME}`)?.set(this, connectivityCol.name);

    this.render();

    this.dataFrame?.onDataChanged.subscribe(() => {
      this.render();
    });
  }

  render(): void {
    if (this.renderTimeout)
      clearTimeout(this.renderTimeout);
    this.renderTimeout = setTimeout(() => {
      if (!this.dataFrame)
        return;
      $(this.root).empty();
      const scViewer = this.scViewer;
      if (scViewer == null) {
        this.root.appendChild(ui.divText(this.viewerError ?? 'Error creating scatter plot'));
        return;
      }
      const clusterSizeLabel = ui.div('Cluster Size', {style: {
        alignSelf: 'center',
        color: 'var(--grey-6)',
        marginBottom: '5px'},
      });
      const maxActivityLabel = ui.div(
        this.activityTarget === ACTIVITY_TARGET.HIGH ? 'Maximum Activity' : 'Minimum Activity', {style: {
          color: 'var(--grey-6)',
          alignSelf: 'center',
          textOrientation: 'mixed',
          writingMode: 'tb',
          transform: 'rotate(180deg)',
          marginLeft: '5px',
        }});
      scViewer.props.colorColumnName = this.colorColumnName ?? null;

      this.root.appendChild(
        ui.divH([
          maxActivityLabel,
          ui.divV([
            ui.divH([this._titleHost, this._selsectIcon], {
              style: {
                alignSelf: 'center',
                lineHeight: 'normal',
                width: '100%',
              },
            }),
            scViewer.root,
            clusterSizeLabel,
          ], {style: {flexGrow: '1'},
          }),
        ]),
      );
      scViewer.root.style.width = '100%';
      //this.root.appendChild(scViewer.root);
      setTimeout(() => {
        scViewer.props.filter = this.scFilterQuery;
        scViewer.invalidateCanvas();
      }, 100);
    }, this.renderDebounceTime);
  }

  onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);
    if (property?.name !== `${ClusterMaxActivityProps.COLOR_COLUMN}${COLUMN_NAME}`)
      this._scViewer = null;
    this.render();
  }
}

