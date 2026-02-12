import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {PADDING_Y} from './combined-measurements-cell-renderer';

const Y_AXIS_CELL_TYPE = 'y-axis';
const TICK_LENGTH = 4;
const FONT_SIZE = 9;

@grok.decorators.cellRenderer({
  name: 'yAxisRenderer',
  cellType: 'y-axis',
})
export class YAxisCellRenderer extends DG.GridCellRenderer {
  get name(): string {
    return 'yAxisRenderer';
  }

  get cellType(): string {
    return Y_AXIS_CELL_TYPE;
  }

  get defaultWidth(): number {
    return 50;
  }

  get defaultHeight(): number {
    return 60;
  }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle,
  ): void {
    const value = gridCell.cell.value;
    if (!value)
      return;

    let data: {min: number, max: number};
    try {
      data = typeof value === 'string' ? JSON.parse(value) : value;
    }
    catch {
      return;
    }

    if (data.min === data.max && data.min === 0)
      return;

    const paddingY = PADDING_Y - 3;
    const lineX = x + w - 8;
    const topY = y + paddingY;
    const bottomY = y + h - paddingY;

    g.strokeStyle = '#888';
    g.lineWidth = 1;

    g.beginPath();
    g.moveTo(lineX, topY);
    g.lineTo(lineX, bottomY);
    g.stroke();

    g.beginPath();
    g.moveTo(lineX - TICK_LENGTH, topY);
    g.lineTo(lineX, topY);
    g.stroke();

    g.beginPath();
    g.moveTo(lineX - TICK_LENGTH, bottomY);
    g.lineTo(lineX, bottomY);
    g.stroke();

    g.fillStyle = '#888';
    g.font = `${FONT_SIZE}px sans-serif`;
    g.textAlign = 'right';

    const maxLabel = formatNumber(data.max);
    const minLabel = formatNumber(data.min);

    g.textBaseline = 'top';
    g.fillText(maxLabel, lineX - TICK_LENGTH - 2, topY + 1);
    g.textBaseline = 'middle';
    g.fillText(minLabel, lineX - TICK_LENGTH - 2, bottomY - 1);
  }
}

function formatNumber(v: number): string {
  if (Number.isInteger(v))
    return v.toString();
  return v.toPrecision(3);
}
