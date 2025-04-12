import {Point} from "datagrok-api/src/grid";

export function drawLines(g: CanvasRenderingContext2D, points: Iterable<Point>) {
  g.beginPath();

  let first = true;
  for (const point of points) {
    if (first) {
      g.moveTo(point.x, point.y);
      first = false;
    } else {
      g.lineTo(point.x, point.y);
    }
  }

  g.stroke();
  return this;
}