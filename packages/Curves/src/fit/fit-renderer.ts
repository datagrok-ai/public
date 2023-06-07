import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {
  FitResult,
  fitResultProperties,
  getFittedCurve,
  SigmoidFunction,
} from '@datagrok-libraries/statistics/src/parameter-estimation/fit-curve';
import {StringUtils} from '@datagrok-libraries/utils/src/string-utils';

import {fitSeries, getChartData, getChartBounds, IFitChartData, IFitSeries,
  CONFIDENCE_INTERVAL_FILL_COLOR, CONFIDENCE_INTERVAL_STROKE_COLOR, CURVE_CONFIDENCE_INTERVAL_BOUNDS,
  TAG_FIT_CHART_FORMAT, TAG_FIT_CHART_FORMAT_3DX, getSeriesConfidenceInterval, getSeriesStatistics} from './fit-data';
import {convertXMLToIFitChartData} from './fit-parser';
import {Viewport} from './transform';
import {MultiCurveViewer} from './multi-curve-viewer';

/** Performs a chart layout, returning [viewport, xAxis, yAxis] */
function layoutChart(rect: DG.Rect): [DG.Rect, DG.Rect?, DG.Rect?] {
  if (rect.width < 100 || rect.height < 100)
    return [rect, undefined, undefined];
  return [
    rect.cutLeft(30).cutBottom(30),
    rect.getBottom(30).cutLeft(30),
    rect.getLeft(30).cutBottom(30)
  ];
}

/** Performs a curve confidence interval drawing */
function drawConfidenceInterval(g: CanvasRenderingContext2D, 
  confIntervals: {confidenceTop: (x: number) => number, confidenceBottom: (x: number) => number},
  screenBounds: DG.Rect, transform: Viewport, confidenceType: string): void {
  g.beginPath();
  for (let i = 0; i <= screenBounds.width; i++) {
    const x = screenBounds.x + i;
    const y = confidenceType === CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP ?
      transform.yToScreen(confIntervals.confidenceTop(transform.xToWorld(x))) :
      transform.yToScreen(confIntervals.confidenceBottom(transform.xToWorld(x)));
    if (i === 0)
      g.moveTo(x, y);
    else
      g.lineTo(x, y);
  }
  g.stroke();
}

/** Performs a curve confidence interval color filling */
function fillConfidenceInterval(g: CanvasRenderingContext2D, confIntervals: {confidenceTop: (x: number) => number, confidenceBottom: (x: number) => number},
  screenBounds: DG.Rect, transform: Viewport): void {
  g.beginPath();
  for (let i = 0; i <= screenBounds.width; i++) {
    const x = screenBounds.x + i;
    const y = transform.yToScreen(confIntervals.confidenceTop(transform.xToWorld(x)));
    if (i === 0)
      g.moveTo(x, y);
    else
      g.lineTo(x, y);
  }

  // reverse traverse to make a shape of confidence interval to fill it
  for (let i = screenBounds.width; i >= 0; i--) {
    const x = screenBounds.x + i;
    const y = transform.yToScreen(confIntervals.confidenceBottom(transform.xToWorld(x)));
    g.lineTo(x, y);
  }
  g.closePath();
  g.fill();
}

export class FitChartCellRenderer extends DG.GridCellRenderer {
  get name() { return 'fit'; }

  get cellType() { return 'fit'; }

  getDefaultSize(gridColumn: DG.GridColumn): {width?: number | null, height?: number | null} {
    return {width: 160, height: 100};
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    grok.shell.o = gridCell;

    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    const screenBounds = gridCell.bounds.inflate(-6, -6);
    const dataBox = layoutChart(screenBounds)[0];
    const dataBounds = getChartBounds(data);
    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);

    for (let i = 0; i < data.series?.length!; i++) {
      if (!data.series![i].clickToToggle || data.series![i].showBoxPlot)
        continue;
      for (let j = 0; j < data.series![i].points.length!; j++) {
        const p = data.series![i].points[j];
        const screenX = viewport.xToScreen(p.x);
        const screenY = viewport.yToScreen(p.y);
        const pxPerMarkerType = (p.outlier ? 6 : 4) / 2;
        if (e.offsetX >= screenX - pxPerMarkerType && e.offsetX <= screenX + pxPerMarkerType &&
          e.offsetY >= screenY - pxPerMarkerType && e.offsetY <= screenY + pxPerMarkerType) {
            p.outlier = !p.outlier;
            // temporarily works only for JSON structure
            if (gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) !== TAG_FIT_CHART_FORMAT_3DX) {
              const gridCellValue = JSON.parse(gridCell.cell.value) as IFitChartData;
              gridCellValue.series![i].points[j].outlier = p.outlier;
              gridCell.cell.value = JSON.stringify(gridCellValue);
            }
            return;
          }
      }
    }   
  }

  onDoubleClick(gridCell: DG.GridCell, e: MouseEvent): void {
    ui.dialog({title: 'Edit chart'})
      .add(MultiCurveViewer.fromChartData(getChartData(gridCell)).root)
      .show();
  }

  renderCurves(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, data: IFitChartData): void {
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();

    const screenBounds = new DG.Rect(x, y, w, h).inflate(-6, -6);
    const [dataBox, xAxisBox, yAxisBox] = layoutChart(screenBounds);

    const dataBounds = getChartBounds(data);
    const viewport = new Viewport(dataBounds, dataBox, data.chartOptions?.logX ?? false, data.chartOptions?.logY ?? false);
    const minSize = Math.min(dataBox.width, dataBox.height);
    const ratio = minSize > 100 ? 1 : 0.2 + (minSize / 100) * 0.8;

    viewport.drawCoordinateGrid(g, xAxisBox, yAxisBox);

    for (const series of data.series!) {
      series.points.sort((a, b) => a.x - b.x);
      let curve: (x: number) => number;
      if (series.parameters)
        curve = getFittedCurve(new SigmoidFunction().y, series.parameters);
      else {
        const fitResult = fitSeries(series);
        curve = fitResult.fittedCurve;
        series.parameters = fitResult.parameters;
      }

      if (series.showPoints ?? true) {
        g.strokeStyle = series.pointColor ?? '0xFF40699c';
        for (let i = 0, candleStart = null; i < series.points.length!; i++) {
          const p = series.points[i];
          const nextSame = i + 1 < series.points.length && series.points[i + 1].x === p.x;
          if (!candleStart && nextSame)
            candleStart = i;
          else if ((series.showBoxPlot ?? false) && candleStart !== null && !nextSame) {
            let minY = series.points[candleStart].y;
            let maxY = minY;
            for (let j = candleStart; j < i; j++) {
              minY = Math.min(minY, series.points[j].y);
              maxY = Math.max(minY, series.points[j].y);
            }

            g.beginPath();
            g.moveTo(viewport.xToScreen(p.x), viewport.yToScreen(minY));
            g.lineTo(viewport.xToScreen(p.x), viewport.yToScreen(maxY));
            g.stroke();

            candleStart = null;
          }
          else if (!candleStart || !series.showBoxPlot) {
            DG.Paint.marker(g,
              p.outlier ? DG.MARKER_TYPE.OUTLIER : DG.MARKER_TYPE.CIRCLE,
              viewport.xToScreen(p.x), viewport.yToScreen(p.y),
              series.pointColor ? DG.Color.fromHtml(series.pointColor) : DG.Color.scatterPlotMarker,
              (p.outlier ? 6 : 4) * ratio);
          }
        }
      }

      if (series.showFitLine ?? true) {
        g.strokeStyle = series.fitLineColor ?? 'black';
        g.lineWidth = 2 * ratio;

        g.beginPath();
        for (let i = 0; i <= screenBounds.width; i++) {
          const x = screenBounds.x + i;
          const y = viewport.yToScreen(curve(viewport.xToWorld(x)));
          if (i === 0)
            g.moveTo(x, y);
          else
            g.lineTo(x, y);
        }
        g.stroke();
      }

      if ((series.showFitLine ?? true) && (series.showCurveConfidenceInterval ?? true)) {
        g.strokeStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_STROKE_COLOR;
        g.fillStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_FILL_COLOR;

        const confidenceIntervals = getSeriesConfidenceInterval(series);
        drawConfidenceInterval(g, confidenceIntervals, screenBounds, viewport, CURVE_CONFIDENCE_INTERVAL_BOUNDS.TOP);
        drawConfidenceInterval(g, confidenceIntervals, screenBounds, viewport, CURVE_CONFIDENCE_INTERVAL_BOUNDS.BOTTOM);
        fillConfidenceInterval(g, confidenceIntervals, screenBounds, viewport);
      }


      if (data.chartOptions?.showStatistics) {
        const statistics = getSeriesStatistics(series);
        for (let i = 0; i < data.chartOptions.showStatistics.length; i++) {
          const statName = data.chartOptions.showStatistics[i];
          const prop = fitResultProperties.find(p => p.name === statName);
          if (prop) {
            const s = StringUtils.formatNumber(prop.get(statistics));
            g.fillStyle = series.fitLineColor ?? 'black';
            g.textAlign = 'left';
            g.fillText(prop.name + ': ' + s, dataBox.x + 5, dataBox.y + 20 + 20 * i);
          }
        }
      }
    }
    g.restore();
  }

  render(g: CanvasRenderingContext2D,
         x: number, y: number, w: number, h: number,
         gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    if (w < 20 || h < 10) return;

    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_FORMAT) === TAG_FIT_CHART_FORMAT_3DX
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    this.renderCurves(g, x, y, w, h, data);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (gridCell.bounds.width < 50) {
      const canvas = ui.canvas(300, 200);
      this.render(canvas.getContext('2d')!, 0, 0, 300, 200, gridCell, null as any);
      const content = ui.divV([canvas]);
      ui.tooltip.show(content, e.x, e.y);
    }
  }
}
