/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {merge} from 'rxjs';

import {IFitChartData, IFitPoint, LogOptions} from '@datagrok-libraries/statistics/src/fit/fit-curve';
import {getChartBounds} from '@datagrok-libraries/statistics/src/fit/fit-data';
import {getStatistic} from '@datagrok-libraries/statistics/src/fit/fit-engine';
import {FitConstants} from '@datagrok-libraries/statistics/src/fit/const';
import {Viewport} from '@datagrok-libraries/utils/src/transform';

import {isNativeFormat} from './curve-converter';
import {getOrCreateParsedChartData} from './fit-chart-data';
import {calculateSeriesFit, getChartDataAggrStats} from './fit-statistics';
import {areAxesLabelsShown, inflateScreenBounds, isTitleShown, layoutChart} from './fit-layout';
import {isLegendVisible, legendTooltip, legendTooltipElement} from './fit-legend';

/** What a user does to a curve in the grid: toggling an outlier, hovering a point, opening the editor. */

export function setOutlier(gridCell: DG.GridCell, p: IFitPoint, seriesIdx: number, pointIdx: number, _data?: IFitChartData) {
  p.outlier = !p.outlier;
  // outlier toggling works only for native JSON format (no converter)
  if (isNativeFormat(gridCell.cell.column)) {
    const gridCellValue = JSON.parse(gridCell.cell.value) as IFitChartData;
    gridCellValue.series![seriesIdx].points[pointIdx].outlier = p.outlier;
    gridCell.cell.column.set(gridCell.cell.rowIndex, JSON.stringify(gridCellValue), true);
    grok.events.fireCustomEvent('fit-cell-outlier-toggle', {
      gridCell: gridCell,
      series: gridCellValue.series![seriesIdx],
      seriesIdx: seriesIdx,
      pointIdx: pointIdx,
    });
  }
  const columns = gridCell.grid.dataFrame.columns.byTags({'.sourceColumn': gridCell.cell.column.name});
  if (columns) {
    _data ??= getOrCreateParsedChartData(gridCell.cell);
    for (const column of columns) {
      const chartLogOptions: LogOptions = {logX: _data.chartOptions?.logX, logY: _data.chartOptions?.logY};
      const statName = column.tags['.statistics'];
      // resolved by name - the legacy seven-key shape has no ic50/pIC50/maxY slot
      let value: number | undefined;
      if (column.tags['.seriesAggregation'] !== null)
        value = getChartDataAggrStats(_data, column.tags['.seriesAggregation'], gridCell.cell)[statName];
      else if (column.tags['.seriesNumber'] === seriesIdx.toString())
        value = getStatistic(calculateSeriesFit(_data.series![seriesIdx], seriesIdx, chartLogOptions, gridCell.cell), statName);
      else
        continue;
      column.set(gridCell.cell.rowIndex, value ?? null);
    }
  }
}

export function inspectCurve(gridCell: DG.GridCell, size?: Partial<DG.Size>, rerenderOnChange: boolean = false): void {
  if (!gridCell.cell.value)
    return;

  const gridCellWidget = DG.GridCellWidget.fromGridCell(gridCell);
  gridCellWidget.root.style.removeProperty('aspectRatio');
  gridCellWidget.root.style.height = '100%';
  gridCellWidget.canvas.style.removeProperty('height');
  gridCellWidget.canvas.style.removeProperty('width');
  gridCellWidget.canvas.style.left = '18px';
  gridCellWidget.canvas.style.top = '18px';

  const dlg = ui.dialog({title: 'Edit chart'})
    .add(gridCellWidget.root)
    .show({resizable: true, width: size?.width ?? 500, height: size?.height ?? 430});
  dlg.getButton('CANCEL').textContent = 'Close';
  if (rerenderOnChange)
    dlg.sub(merge(gridCell.grid.dataFrame.onDataChanged, gridCell.grid.dataFrame.onMetadataChanged).subscribe(() => gridCellWidget.render()));
}

export function hitTest(e: MouseEvent, point: IFitPoint, viewport: Viewport): boolean {
  const screenX = viewport.xToScreen(point.x);
  const screenY = viewport.yToScreen(point.y);
  const pxPerMarkerType = ((point.outlier ? FitConstants.OUTLIER_PX_SIZE : FitConstants.POINT_PX_SIZE) / 2) + FitConstants.OUTLIER_HITBOX_RADIUS;
  const pointRect = new DG.Rect(screenX - pxPerMarkerType, screenY - pxPerMarkerType,
    2 * pxPerMarkerType, 2 * pxPerMarkerType);
  return pointRect.containsPoint(new DG.Point(e.offsetX, e.offsetY));
}

function pointViewport(screenBounds: DG.Rect, data: IFitChartData): Viewport {
  const dataBox = layoutChart(screenBounds, areAxesLabelsShown(screenBounds, data), isTitleShown(screenBounds, data))[0];
  return new Viewport(getChartBounds(data), dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);
}

export function handleClick(gridCell: DG.GridCell, e: MouseEvent): void {
  if (!gridCell.cell.value || !gridCell.cell.column)
    return;

  const data = getOrCreateParsedChartData(gridCell.cell, false);

  for (const [message, condition] of Object.entries(FitConstants.CONDITION_MAP)) {
    if (condition(data.series)) {
      grok.shell.o = ui.divText(message, {style: {color: 'red'}});
      return;
    }
  }

  grok.shell.o = gridCell;

  const screenBounds = inflateScreenBounds(gridCell.bounds);
  const viewport = pointViewport(screenBounds, data);

  for (let i = 0; i < data.series?.length!; i++) {
    if (data.series![i].connectDots || !data.series![i].clickToToggle || data.series![i].showPoints !== 'points' ||
      screenBounds.width < FitConstants.MIN_AXES_CELL_PX_WIDTH || screenBounds.height < FitConstants.MIN_AXES_CELL_PX_HEIGHT)
      continue;
    for (let j = 0; j < data.series![i].points.length!; j++) {
      const p = data.series![i].points[j];
      if (p.outlier && !data.series![i].showOutliers)
        continue;
      if (hitTest(e, p, viewport)) {
        setOutlier(gridCell, p, i, j, data);
        return;
      }
    }
  }
}

/** `renderer` is the base type on purpose - taking the concrete one back would make this module and
 * the renderer import each other. */
export function handleMouseMove(gridCell: DG.GridCell, e: MouseEvent, renderer: DG.GridCellRenderer): void {
  if (!gridCell?.cell?.value)
    return;

  const screenBounds = inflateScreenBounds(gridCell.bounds);
  if (screenBounds.width < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH ||
    screenBounds.height < FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT) {
    const canvas = ui.canvas(300, 200);
    renderer.render(canvas.getContext('2d')!, 0, 0, 300, 200, gridCell, null as any);
    ui.tooltip.show(ui.divV([canvas]), e.x, e.y);
  }

  const data = getOrCreateParsedChartData(gridCell?.cell);

  if (screenBounds.width >= FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_WIDTH &&
    screenBounds.height >= FitConstants.MIN_POINTS_AND_STATS_VISIBILITY_PX_HEIGHT) {
    const viewport = pointViewport(screenBounds, data);

    // what the legend had no room to say: the whole of a shortened name, or everything it left out
    const legendRows = isLegendVisible(data, screenBounds) ?
      legendTooltip(data, viewport.screen, e.offsetX, e.offsetY) : null;
    if (legendRows !== null) {
      ui.tooltip.show(legendTooltipElement(legendRows), e.x + 16, e.y + 16);
      return;
    }

    for (let i = 0; i < data.series?.length!; i++) {
      if (data.series![i].showPoints !== 'points')
        continue;
      for (let j = 0; j < data.series![i].points.length!; j++) {
        const p = data.series![i].points[j];
        if (p.outlier && !data.series![i].showOutliers)
          continue;
        if (hitTest(e, p, viewport)) {
          ui.tooltip.show(ui.divV([ui.divText(`${data.chartOptions?.xAxisName ?? 'x'}: ${DG.format(p.x, !data.chartOptions?.logX ? '#0.000' : 'scientific')}`),
            ui.divText(`${data.chartOptions?.yAxisName ?? 'y'}: ${DG.format(p.y, !data.chartOptions?.logY ? '#0.000' : 'scientific')}`)]), e.x + 16, e.y + 16);
          if (!data.series![i].connectDots && data.series![i].clickToToggle && screenBounds.width >= FitConstants.MIN_AXES_CELL_PX_WIDTH &&
            screenBounds.height >= FitConstants.MIN_AXES_CELL_PX_HEIGHT)
            document.body.style.cursor = 'pointer';
          return;
        }
      }
    }
    ui.tooltip.hide();
  }
  document.body.style.cursor = 'default';
}
