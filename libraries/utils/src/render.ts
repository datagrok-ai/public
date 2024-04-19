import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

declare global {
  interface CanvasRenderingContext2D {
    drawCross(r: DG.Rect, strokeStyle?: string, lineWidth?: number): void;
  }
}

CanvasRenderingContext2D.prototype.drawCross = function drawCross(
  r: DG.Rect, strokeStyle: string = '#FF0000', lineWidth: number = 2
): void {
  const q = r.fitSquare();
  this.clearRect(r.left, r.top, r.width, r.height);

  this.strokeStyle = strokeStyle;
  this.lineWidth = lineWidth;
  this.beginPath();
  this.moveTo(q.left, q.top);
  this.lineTo(q.right, q.bottom);
  this.moveTo(q.right, q.top);
  this.lineTo(q.left, q.bottom);
  this.stroke();
};
