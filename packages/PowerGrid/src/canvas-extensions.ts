import * as DG from 'datagrok-api/dg';
import {Point} from "datagrok-api/src/grid";
import {ColorType} from "datagrok-api/src/const";

declare global {

  interface CanvasRenderingContext2D {

    setFillStyle(fill: string | CanvasGradient | CanvasPattern): CanvasRenderingContext2D;

    setStrokeStyle(stroke: string | CanvasGradient | CanvasPattern): CanvasRenderingContext2D;

    line(x1: number, y1: number, x2: number, y2: number, color: ColorType): CanvasRenderingContext2D;

    /**
     * Use stroke() or fill() after.
     * @param pa: Array of points
     */
    polygon(pa: Point[]): CanvasRenderingContext2D;
  }
}

CanvasRenderingContext2D.prototype.setFillStyle = function (fill) {
  this.fillStyle = fill;
  return this;
}

CanvasRenderingContext2D.prototype.setStrokeStyle = function (stroke) {
  this.strokeStyle = stroke;
  return this;
}

CanvasRenderingContext2D.prototype.line = function (x1, y1, x2, y2, color) {

  this.beginPath();
  this.strokeStyle = DG.Color.toRgb(color);
  this.moveTo(x1, y1);
  this.lineTo(x2, y2);
  this.stroke();

  return this;
}

CanvasRenderingContext2D.prototype.polygon = function (pa) {

  this.beginPath();

  const last_p = pa[pa.length - 1]
  this.moveTo(last_p.x, last_p.y);

  for (let p of pa) {
    this.lineTo(p.x, p.y);
  }
  this.closePath();

  return this;
}
