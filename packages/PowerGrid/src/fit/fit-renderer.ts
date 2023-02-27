import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {GridColumn} from "datagrok-api/dg";
import {fitSeries, getChartData, getChartBounds, getFittedCurve} from "./fit-data";

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
    }
  }
}

export class FitChartCellRenderer extends DG.GridCellRenderer {
  get name() { return 'fit'; }

  get cellType() { return 'fit'; }

  getDefaultSize(gridColumn: GridColumn): {width?: number | null, height?: number | null} {
    return { width: 160, height: 100};
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    grok.shell.o = gridCell;
  }

  render(g: CanvasRenderingContext2D,
         x: number, y: number, w: number, h: number,
         gridCell: DG.GridCell, cellStyle: DG.GridCellStyle)
  {
    g.save();
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();

    if (w < 20 || h < 10) return;
    const screenBounds = new DG.Rect(x, y, w, h).inflate(-6, -6);

    const data = getChartData(gridCell);
    let dataBounds = getChartBounds(data);
    let transform = Transform.linear(dataBounds, screenBounds);

    for (let series of data.series!) {
      if (series.showPoints ?? true) {
        for (let i = 0; i < series.points.length!; i++) {
          let p = series.points[i];
          DG.Paint.marker(g,
            p.outlier ? DG.MARKER_TYPE.OUTLIER : DG.MARKER_TYPE.CIRCLE,
            transform.xToScreen(p.x), transform.yToScreen(p.y),
            series.pointColor ? DG.Color.fromHtml(series.pointColor) : DG.Color.scatterPlotMarker,
            p.outlier ? 6 : 4);
        }
      }

      if (series.showFitLine ?? true ) {
        g.strokeStyle = series.fitLineColor ?? 'black';
        let curve = getFittedCurve(series);

        g.beginPath();
        for (let i = 0; i < series.points.length!; i++) {
          let x = transform.xToScreen(series.points[i].x);
          let y = transform.yToScreen(curve(series.points[i].x));
          if (i == 0)
            g.moveTo(x, y);
          else
            g.lineTo(x, y);
        }
        g.stroke();
      }
    }

    g.restore();
  }
}
