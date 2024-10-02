import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


@grok.decorators.cellRenderer({
  name: 'binaryImageCellRenderer',
  cellType: 'BinaryImage',
})
export class BinaryImageCellRenderer extends DG.GridCellRenderer {
  get name() { return 'BinaryImage'; }

  get cellType() { return 'BinaryImage'; }

  get defaultWidth(): number | null { return 200; }

  get defaultHeight(): number | null { return 100; }

  static images: DG.LruCache = new DG.LruCache();

  render(
    g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ) {
    if (w < 5 || h < 5) return;
    const bytes = gridCell.cell.value;
    if (bytes?.constructor === Uint8Array)
      DG.Paint.pngImage(g, new DG.Rect(x, y, w, h), bytes);
  }

  onDoubleClick(gridCell: DG.GridCell, e: MouseEvent): void {
    DG.Utils.openFileBytes({
      accept: 'image/jpg, image/jpeg, image/png',
      open: (bytes) => {
        gridCell.cell.value = bytes;
      },
    });
  }
}
