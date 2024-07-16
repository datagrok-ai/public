import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


@grok.decorators.cellRenderer({
  name: 'hyperlinkCellRenderer',
  cellType: 'Hyperlink',
})
export class HyperlinkCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Hyperlink'; }

  get cellType() { return 'Hyperlink'; }

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    g.fillText(`link ${gridCell.cell.value}`, x + 3, y + 3);
  }
}
