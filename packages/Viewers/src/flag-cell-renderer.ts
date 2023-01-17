import * as DG from 'datagrok-api/dg';

export class FlagCellRenderer extends DG.GridCellRenderer {
  get name() {
    return 'Flag cell renderer';
  }

  get cellType() {
    return 'flag';
  }

  get defaultWidth() {
    return 50;
  }

  get defaultHeight() {
    return 50;
  }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    g.fillStyle = 'black';
    g.fillText('flag', x, y);
  }
}
