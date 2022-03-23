import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class HyperlinkCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Hyperlink'; }
  get cellType() { return 'Hyperlink'; }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    g.fillText(`link ${gridCell.cell.value}`, x + 3, y + 3);
  }
}