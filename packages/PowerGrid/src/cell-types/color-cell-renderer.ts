import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';


const COLOR_CELL_PADDING_PX = 4;
const COLOR_CELL_TYPE = 'Color';

/** Validates any CSS color string (hex, rgb/rgba, hsl/hsla, named) without rendering. */
export function isValidColor(color: string): boolean {
  const s = new Option().style;
  s.color = color;
  return s.color !== '';
}

@grok.decorators.cellRenderer({
  name: 'Color',
  cellType: 'Color',
})
export class ColorCellRenderer extends DG.GridCellRenderer {
  get name() { return COLOR_CELL_TYPE; }
  get cellType() { return COLOR_CELL_TYPE; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    if (gridCell.cell.isNone())
      return;

    const color = gridCell.cell.valueString.trim();
    if (!color || !isValidColor(color))
      return;

    const p = Math.min(COLOR_CELL_PADDING_PX, Math.floor(Math.min(w, h) / 5));
    const rx = x + p, ry = y + p, rw = w - 2 * p, rh = h - 2 * p;
    g.fillStyle = color;
    g.fillRect(rx, ry, rw, rh);
    g.strokeStyle = 'rgba(0, 0, 0, 0.15)';
    g.lineWidth = 1;
    g.strokeRect(rx + 0.5, ry + 0.5, rw - 1, rh - 1);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (gridCell.cell.isNone())
      return;
    const color = gridCell.cell.valueString.trim();
    if (color)
      ui.tooltip.show(color, e.x + 16, e.y + 16);
  }

  onMouseLeave(gridCell: DG.GridCell, e: MouseEvent): void {
    ui.tooltip.hide();
  }
}
