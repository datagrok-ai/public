import * as DG from 'datagrok-api/dg';
import {TYPE} from 'datagrok-api/dg';

interface BarCellSettings {
  color: string;
  radius: number;
}

const BarCellSettingsProperties = {
  color: DG.Property.js('color', TYPE.STRING),
  radius: DG.Property.js('radius', TYPE.STRING),
};

export class BarCellRenderer extends DG.GridCellRenderer {
  get name() {return 'bar';}

  get cellType() {return 'bar';}

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle,
  ) {
    if (gridCell.cell.isNone())
      return;

    const ratio = gridCell.cell.value / gridCell.cell.column.max;
    g.fillStyle = '#006400';
    g.roundRect(x + 2, y + 4, (w - 4) * ratio, h - 8, 4).fill();
  }
}
