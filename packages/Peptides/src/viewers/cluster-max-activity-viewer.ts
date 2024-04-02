import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';
import {ACTIVITY_TARGET, COLUMN_NAME, COLUMNS_NAMES} from '../utils/constants';
import wu from 'wu';
import $ from 'cash-dom';
export const enum ClusterMaxActivityProps {
    CLUSTER_COLUMN = 'cluster',
    ACTIVITY_COLUMN = 'activity',
    COLOR_COLUMN = 'color',
}

export interface IClusterMaxActivity {
    clusterColumnName: string;
    activityColumnName: string;
    colorColumnName?: string;
    activityTarget: ACTIVITY_TARGET;
}


export class ClusterMaxActivityViewer extends DG.JsViewer implements IClusterMaxActivity {
  _titleHost = ui.divText(VIEWER_TYPE.CLUSTER_MAX_ACTIVITY, {id: 'pep-viewer-title'});
  clusterColumnName: string;
  activityColumnName: string;
  colorColumnName: string;
  activityTarget: ACTIVITY_TARGET = ACTIVITY_TARGET.HIGH;
  _scViewer?: DG.ScatterPlotViewer | null;
  viewerError: string = '';
  renderTimeout: NodeJS.Timeout | number | null = null;
  renderDebounceTime = 500;
  static clusterSizeColName = '~cluster.size' as const;
  static maxActivityInClusterSizeColName = '~max.activity.for.cluster.size' as const;
  private scFilterQuery = `\$\{${ClusterMaxActivityViewer.maxActivityInClusterSizeColName}\} == true` as const;
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
  }

  /**
   * Returns PeptidesModel instance that belongs to the attached dataframe.
   * @return - PeptidesModel instance.
   */
  get model(): PeptidesModel {
    return PeptidesModel.getInstance(this.dataFrame);
  }

  private createSCViewer(): DG.ScatterPlotViewer | null {
    const scatterPlotProps: Partial<DG.IScatterPlotLookSettings> & Options = {
      showXAxis: true,
      showYAxis: true,
      showXSelector: false,
      showYSelector: false,
      showColorSelector: false,
      xAxisType: 'logarithmic',
      yAxisType: 'logarithmic',
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
    const numericColTypes: DG.ColumnType[] =
        [DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.BIG_INT, DG.COLUMN_TYPE.QNUM];
    if (!numericColTypes.includes(activityCol.type)) {
      this.viewerError = 'Activity column should be numeric';
      return null;
    }
    const clusterSizeCol = this.dataFrame.columns.getOrCreate(ClusterMaxActivityViewer.clusterSizeColName,
      DG.TYPE.INT, this.dataFrame.rowCount);
    const clusterSizeMap: {[key: number | string]: number} = {};
    for (let i = 0; i < this.dataFrame.rowCount; i++) {
      const cluster: string | number = clusterCol.get(i);
      if (cluster == null)
        continue;
      clusterSizeMap[cluster] = (clusterSizeMap[cluster] ?? 0) + 1;
    }
    // for (let i = 0; i < this.dataFrame.rowCount; i++) {
    //   const cluster: string | number = clusterCol.get(i);
    //   if (clusterCol.isNone(i) || !clusterSizeMap[cluster])
    //     continue;
    //   clusterSizeCol.set(i, clusterSizeMap[cluster]);
    // }
    // clusterSizeCol.init((i) => {
    //   const cluster: string | number = clusterCol.get(i);
    //   if (clusterCol.isNone(i) || !clusterSizeMap[cluster])
    //     return null;
    //   return clusterSizeMap[cluster];
    // });
    clusterSizeCol.init((i) => clusterCol.isNone(i) ? null : clusterSizeMap[clusterCol.get(i)] ?? null);

    if (activityCol.stats.min <= 0)
      scatterPlotProps.yAxisType = 'linear';

    // create a new column to store max activity for each cluster size
    const maxActivityIndexPerClusterSizeMap: {[key: number]: number} = {};

    for (let i = 0; i < this.dataFrame.rowCount; i++) {
      const clusterSize: number | null = clusterSizeCol.get(i);
      if (clusterSize == null || clusterSizeCol.isNone(i) || activityCol.isNone(i))
        continue;
      const activity: number = activityCol.get(i);
      const prevMaxActivityIndex = maxActivityIndexPerClusterSizeMap[clusterSize];
      if (prevMaxActivityIndex == null || prevMaxActivityIndex == undefined)
        maxActivityIndexPerClusterSizeMap[clusterSize] = i;
      else {
        if (activity > activityCol.get(prevMaxActivityIndex) && this.activityTarget === ACTIVITY_TARGET.HIGH)
          maxActivityIndexPerClusterSizeMap[clusterSize] = i;
        else if (activity < activityCol.get(prevMaxActivityIndex) && this.activityTarget === ACTIVITY_TARGET.LOW)
          maxActivityIndexPerClusterSizeMap[clusterSize] = i;
      }
    }

    const maxAtivityInClusterSizeCol = this.dataFrame.columns.getOrCreate(
      ClusterMaxActivityViewer.maxActivityInClusterSizeColName, DG.COLUMN_TYPE.BOOL, this.dataFrame.rowCount);
    maxAtivityInClusterSizeCol.init((i) => {
      if (clusterSizeCol.isNone(i))
        return false;
      return i === maxActivityIndexPerClusterSizeMap[clusterSizeCol.get(i)];
    });
    scatterPlotProps.xColumnName = ClusterMaxActivityViewer.clusterSizeColName;
    scatterPlotProps.yColumnName = this.activityColumnName;
    scatterPlotProps.filter = this.scFilterQuery;
    this.viewerError = '';
    const sc = DG.Viewer.scatterPlot(this.dataFrame, scatterPlotProps);

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

    this.render();

    this.dataFrame.onDataChanged.subscribe(() => {
    //   this._scViewer = null;
    //   this.render();

      this.render();
    });
  }

  render(): void {
    if (this.renderTimeout)
      clearTimeout(this.renderTimeout);
    this.renderTimeout = setTimeout(() => {
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
            ui.divH([this._titleHost], {
              style: {
                alignSelf: 'center',
                lineHeight: 'normal',
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

