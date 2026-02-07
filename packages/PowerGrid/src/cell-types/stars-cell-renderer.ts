import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {GridCellRenderer} from "datagrok-api/dg";


const STAR_PADDING = 2;
const MIN_STAR_SIZE = 8;
const FILLED_COLOR = '#FFB400';
const EMPTY_COLOR = '#D0D0D0';

function getMaxStars(gridCell: DG.GridCell): number {
  const tag = gridCell.tableColumn?.tags?.['.maxStars'];
  if (tag) {
    const n = parseInt(tag);
    if (n > 0)
      return n;
  }
  return 5;
}

function drawStar(g: CanvasRenderingContext2D, cx: number, cy: number, outerR: number) {
  const innerR = outerR * 0.4;
  const points = 5;
  g.beginPath();
  for (let i = 0; i < points * 2; i++) {
    const r = i % 2 === 0 ? outerR : innerR;
    const angle = (i * Math.PI / points) - Math.PI / 2;
    const px = cx + r * Math.cos(angle);
    const py = cy + r * Math.sin(angle);
    if (i === 0)
      g.moveTo(px, py);
    else
      g.lineTo(px, py);
  }
  g.closePath();
}

@grok.decorators.cellRenderer({
  name: 'Stars',
  cellType: 'Stars',
})
export class StarsCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Stars'; }

  get cellType() { return 'Stars'; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    if (gridCell.cell.isNone())
      return;

    const maxStars = getMaxStars(gridCell);
    const value = Math.round(gridCell.cell.value);
    const starSize = Math.min((w - STAR_PADDING) / maxStars - STAR_PADDING, h - STAR_PADDING * 2);

    if (starSize < MIN_STAR_SIZE || (!gridCell.cell.isNone() && (value < 0 || value > maxStars))) {
      GridCellRenderer.byName('number')?.render(g, x, y, w, h, gridCell, cellStyle);
      return;
    }

    const outerR = starSize / 2;
    const totalWidth = maxStars * (starSize + STAR_PADDING) - STAR_PADDING;
    const startX = x + (w - totalWidth) / 2;
    const cy = y + h / 2;

    for (let i = 0; i < maxStars; i++) {
      const cx = startX + i * (starSize + STAR_PADDING) + outerR;
      drawStar(g, cx, cy, outerR);
      g.fillStyle = i < value ? FILLED_COLOR : EMPTY_COLOR;
      g.fill();
    }
  }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.gridColumn.editable)
      return;

    const maxStars = getMaxStars(gridCell);
    const w = gridCell.bounds.width;
    const h = gridCell.bounds.height;
    const starSize = Math.min((w - STAR_PADDING) / maxStars - STAR_PADDING, h - STAR_PADDING * 2);

    if (starSize < MIN_STAR_SIZE)
      return;

    const totalWidth = maxStars * (starSize + STAR_PADDING) - STAR_PADDING;
    const startX = (w - totalWidth) / 2;
    const localX = e.offsetX - gridCell.bounds.left;
    const idx = Math.floor((localX - startX) / (starSize + STAR_PADDING));

    if (idx < 0 || idx >= maxStars)
      return;

    const newValue = idx + 1;
    gridCell.cell.value = newValue === gridCell.cell.value ? 0 : newValue;
    gridCell.grid.invalidate();
  }
}
