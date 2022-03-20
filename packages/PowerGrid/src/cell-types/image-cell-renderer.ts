import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class ImageCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Image'; }
  get cellType() { return 'Image'; }

  get defaultWidth(): number | null { return 200; }
  get defaultHeight(): number | null { return 100; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    g.fillText(`image ${gridCell.cell.value}`, x + 3, y + 3);
  }
}