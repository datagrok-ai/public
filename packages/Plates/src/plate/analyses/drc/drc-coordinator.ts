/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Subscription} from 'rxjs';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {FitCellOutlierToggleArgs} from '@datagrok-libraries/statistics/src/fit/new-fit-API';
import {IFitPoint, FitMarkerType} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {Plate} from '../../plate';
import {PlateWidget} from '../../plate-widget/plate-widget';
import {getDoseResponseSeries} from './utils';

export class DrcAnalysisCoordinator {
  private subs: Subscription[] = [];
  private prevSelection: {
    seriesIndex: number;
    pointIndex: number;
    markerType: FitMarkerType;
    markerSize: number;
    markerColor: string;
    curvesGridCell?: DG.GridCell;
  } | null = null;

  constructor(
    private plate: Plate,
    private plateWidget: PlateWidget,
    private curvesGrid: DG.Grid,
    private curveCol: DG.Column,
    private seriesVals: Array<[string, any]>,
    private mappingOptions: {
      roleName: string,
      concentrationName: string,
      valueName: string,
    }
  ) {
    this.listenToPlateWidgetEvents();
    this.listenToGridEvents();
  }

  public regenerateCurves(): void {
    const newSeriesData = getDoseResponseSeries(this.plate, {
      value: this.mappingOptions.valueName,
      concentration: this.mappingOptions.concentrationName,
      groupBy: this.mappingOptions.roleName,
    });

    const newSeriesVals = Object.entries(newSeriesData);

    for (let i = 0; i < Math.min(newSeriesVals.length, this.curveCol.length); i++) {
      const [name, series] = newSeriesVals[i];
      const originalJson = JSON.parse(this.curveCol.get(i));

      originalJson.series[0].points = series.points;
      originalJson.chartOptions.title = `${name}`;

      this.curveCol.set(i, JSON.stringify(originalJson), true);
    }

    this.seriesVals.length = 0;
    Array.prototype.push.apply(this.seriesVals, newSeriesVals);

    this.curvesGrid.invalidate();
  }

  private listenToPlateWidgetEvents(): void {
    const gridSub = this.plateWidget.grid.onCurrentCellChanged.subscribe((gc) => {
      this.clearPreviousSelection();
      if (!gc.isTableCell || !gc.gridColumn || gc.gridColumn.idx === 0) return;

      const rowIdx = this.plate._idx(gc.gridRow, gc.gridColumn.idx - 1);
      if (rowIdx === undefined || rowIdx < 0) return;

      const catValue = this.plate.data.get(this.mappingOptions.roleName, rowIdx)?.toLowerCase();
      if (!catValue) return;

      const seriesIndex = this.seriesVals.findIndex(([serName, _]) => serName?.toLowerCase() === catValue);
      if (seriesIndex < 0) return;

      const concentration: number = this.plate.data.get(this.mappingOptions.concentrationName, rowIdx);
      const value: number = this.plate.data.get(this.mappingOptions.valueName, rowIdx);
      const pointInSeriesIndex: number = this.seriesVals[seriesIndex][1].points.findIndex((p: IFitPoint) => p.x === concentration && p.y === value);
      if (pointInSeriesIndex < 0) return;

      const currentPoints = this.seriesVals[seriesIndex][1].points;
      this.prevSelection = {
        seriesIndex: seriesIndex,
        pointIndex: pointInSeriesIndex,
        markerType: currentPoints[pointInSeriesIndex].marker ?? DG.MARKER_TYPE.CIRCLE,
        markerSize: currentPoints[pointInSeriesIndex].size ?? FitConstants.POINT_PX_SIZE,
        markerColor: currentPoints[pointInSeriesIndex].color ?? DG.Color.toHtml(DG.Color.getCategoricalColor(0))
      };

      const curveJSON = JSON.parse(this.curveCol.get(seriesIndex)!);
      const points: IFitPoint[] = curveJSON.series[0]?.points;
      if (points?.length > pointInSeriesIndex) {
        points[pointInSeriesIndex].marker = DG.MARKER_TYPE.SQUARE;
        points[pointInSeriesIndex].size = FitConstants.POINT_PX_SIZE * 2;
        points[pointInSeriesIndex].color = DG.Color.toHtml(DG.Color.green);
        this.curveCol.set(seriesIndex, JSON.stringify(curveJSON));
      }

      const curveGridRow = this.curvesGrid.tableRowToGrid(seriesIndex);
      this.prevSelection.curvesGridCell = this.curvesGrid.cell(this.curveCol.name, curveGridRow);
      this.prevSelection.curvesGridCell?.grid.scrollToCell(this.curveCol.name, curveGridRow);
    });
    this.subs.push(gridSub);
  }

  private listenToGridEvents(): void {
    const gridEventSub = grok.events.onCustomEvent('fit-cell-outlier-toggle').subscribe((args: FitCellOutlierToggleArgs) => {
      if (!args?.gridCell || !args.series || args.pointIdx == null || args.gridCell.cell.column !== this.curveCol)
        return;

      const point: IFitPoint = args.series.points[args.pointIdx];
      if (point.meta !== null) {
        const [row, col] = this.plate.rowIndexToExcel(point.meta);
        this.plate.markOutlier(row, col, !!point.outlier);
        this.plateWidget.grid.invalidate();
      }
    });
    this.subs.push(gridEventSub);
  }

  private clearPreviousSelection(): void {
    try {
      if (!this.prevSelection || (this.prevSelection.seriesIndex ?? -1) < 0) return;
      const s = this.curveCol.get(this.prevSelection.seriesIndex);
      if (!s) return;
      const parsed = JSON.parse(s);
      const points: IFitPoint[] = parsed.series[0]?.points;
      if (points && points.length > this.prevSelection.pointIndex) {
        points[this.prevSelection.pointIndex].marker = this.prevSelection.markerType;
        points[this.prevSelection.pointIndex].size = this.prevSelection.markerSize;
        points[this.prevSelection.pointIndex].color = this.prevSelection.markerColor;
        this.curveCol.set(this.prevSelection.seriesIndex, JSON.stringify(parsed));
      }
    } finally {
      this.prevSelection = null;
    }
  }

  public destroy(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }
}
