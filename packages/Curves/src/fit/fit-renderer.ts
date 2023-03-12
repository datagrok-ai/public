import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {GridColumn, Paint} from 'datagrok-api/dg';

import {fitResultProperties} from "@datagrok-libraries/statistics/src/parameter-estimation/fit-curve";
import {StringUtils} from "@datagrok-libraries/utils/src/string-utils";

import {fitSeries, getChartData, getChartBounds, getFittedCurve, getConfidenceIntrevals,
  CONFIDENCE_INTERVAL_FILL_COLOR, CONFIDENCE_INTERVAL_STROKE_COLOR, IFitChartData,
  TAG_FIT_CHART_NAME, TAG_FIT_CHART_FORMAT} from './fit-data';
import {convertXMLToIFitChartData} from './fit-parser';


interface ITransform {
  xToScreen(world: number): number;
  yToScreen(world: number): number;
}

class Transform {
  static linear(world: DG.Rect, screen: DG.Rect): ITransform {
    return {
      xToScreen(worldX: number): number {
        return screen.left + screen.width * (worldX - world.left) / world.width;
      },
      yToScreen(worldY: number): number {
        return screen.bottom - screen.height * (worldY - world.top) / world.height;
      }
    };
  }
}

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

export class FitChartCellRenderer extends DG.GridCellRenderer {
  get name() { return 'fit'; }

  get cellType() { return 'fit'; }

  getDefaultSize(gridColumn: GridColumn): {width?: number | null, height?: number | null} {
    return {width: 160, height: 100};
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    grok.shell.o = gridCell;
  }

  renderCurves(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, data: IFitChartData) {
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();

    const screenBounds = new DG.Rect(x, y, w, h).inflate(-6, -6);
    const [dataBox, xAxisBox, yAxisBox] = layoutChart(screenBounds);

    const dataBounds = getChartBounds(data);
    const transform = Transform.linear(dataBounds, dataBox);
    const minSize = Math.min(dataBox.width, dataBox.height);
    const ratio = minSize > 100 ? 1 : 0.2 + (minSize / 100) * 0.8;

    DG.Paint.coordinateGrid(g, dataBounds, xAxisBox, yAxisBox, dataBox);

    for (const series of data.series!) {
      if (series.showPoints ?? true) {
        g.strokeStyle = series.pointColor ?? '0xFF40699c';
        for (let i = 0, candleStart = null; i < series.points.length!; i++) {
          const p = series.points[i];
          const nextSame = i + 1 < series.points.length && series.points[i + 1].x == p.x;
          if (!candleStart && nextSame)
            candleStart = i;
          else if (candleStart != null && !nextSame) {
            let minY = series.points[candleStart].y;
            let maxY = minY;
            for (let j = candleStart; j < i; j++) {
              minY = Math.min(minY, series.points[j].y);
              maxY = Math.max(minY, series.points[j].y);
            }

            g.beginPath();
            g.moveTo(transform.xToScreen(p.x), transform.yToScreen(minY));
            g.lineTo(transform.xToScreen(p.x), transform.yToScreen(maxY));
            g.stroke();

            candleStart = null;
          }
          else if (!candleStart) {
            DG.Paint.marker(g,
              p.outlier ? DG.MARKER_TYPE.OUTLIER : DG.MARKER_TYPE.CIRCLE,
              transform.xToScreen(p.x), transform.yToScreen(p.y),
              series.pointColor ? DG.Color.fromHtml(series.pointColor) : DG.Color.scatterPlotMarker,
              (p.outlier ? 6 : 4) * ratio);
          }
        }
      }

      if (series.showFitLine ?? true ) {
        g.strokeStyle = series.fitLineColor ?? 'black';
        g.lineWidth = 2 * ratio;
        const curve = getFittedCurve(series);

        g.beginPath();
        for (let i = 0; i < series.points.length!; i++) {
          const x = transform.xToScreen(series.points[i].x);
          const y = transform.yToScreen(curve(series.points[i].x));
          if (i == 0)
            g.moveTo(x, y);
          else
            g.lineTo(x, y);
        }
        g.stroke();
      }

      if (series.showCurveConfidenceInterval ?? true) {
        g.strokeStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_STROKE_COLOR;
        g.fillStyle = series.confidenceIntervalColor ?? CONFIDENCE_INTERVAL_FILL_COLOR;
        const confidence = getConfidenceIntrevals(series);

        // draw confidence top
        g.beginPath();
        for (let i = 0; i < series.points.length!; i++) {
          const x = transform.xToScreen(series.points[i].x);
          const y = transform.yToScreen(confidence.top(series.points[i].x));
          if (i == 0)
            g.moveTo(x, y);
          else
            g.lineTo(x, y);
        }
        g.stroke();

        // draw confidence bottom
        g.beginPath();
        for (let i = 0; i < series.points.length!; i++) {
          const x = transform.xToScreen(series.points[i].x);
          const y = transform.yToScreen(confidence.bottom(series.points[i].x));
          if (i == 0)
            g.moveTo(x, y);
          else
            g.lineTo(x, y);
        }
        g.stroke();

        // fill the interval
        g.beginPath();
        for (let i = 0; i < series.points.length!; i++) {
          const x = transform.xToScreen(series.points[i].x);
          const y = transform.yToScreen(confidence.top(series.points[i].x));
          if (i == 0)
            g.moveTo(x, y);
          else
            g.lineTo(x, y);
        }

        // reverse traverse to make a shape of confidence interval to fill it
        for (let i = series.points.length! - 1; i >= 0; i--) {
          const x = transform.xToScreen(series.points[i].x);
          const y = transform.yToScreen(confidence.bottom(series.points[i].x));
          g.lineTo(x, y);
        }
        g.closePath();
        g.fill();
      }


      if (data.chartOptions?.showStatistics) {
        let fitResult = fitSeries(series, true);
        for (let i = 0; i < data.chartOptions.showStatistics.length; i++) {
          const statName = data.chartOptions.showStatistics[i];
          const prop = fitResultProperties.find(p => p.name == statName);
          if (prop) {
            const s = StringUtils.formatNumber(prop.get(fitResult));
            g.fillStyle = series.fitLineColor ?? 'black';
            g.textAlign = "left";
            g.fillText(prop.name + ': ' + s, dataBox.x + 5, dataBox.y + 20 + 20 * i);
          }
        }
      }
    }

    g.restore();
  }

  render(g: CanvasRenderingContext2D,
         x: number, y: number, w: number, h: number,
         gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    if (w < 20 || h < 10) return;

    const data = gridCell.cell.column.getTag(TAG_FIT_CHART_NAME) === TAG_FIT_CHART_FORMAT
      ? convertXMLToIFitChartData(gridCell.cell.value)
      : getChartData(gridCell);

    this.renderCurves(g, x, y, w, h, data);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent) {
    if (gridCell.bounds.width < 50) {
      const canvas = ui.canvas(300, 200);
      this.render(canvas.getContext("2d")!, 0, 0, 300, 200, gridCell, null as any);
      const content = ui.divV([canvas]);
      ui.tooltip.show(content, e.x, e.y);
    }
  }
}
