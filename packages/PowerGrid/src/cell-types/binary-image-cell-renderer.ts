import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {GridCell, Rect, Utils} from "datagrok-api/dg";
import {Paint} from "datagrok-api/src/utils";

// Put on hold - waiting for the binary column to be implemented.

export class BinaryImageCellRenderer extends DG.GridCellRenderer {
  get name() { return 'Binary Image'; }
  get cellType() { return 'binary image'; }

  get defaultWidth(): number | null { return 200; }
  get defaultHeight(): number | null { return 100; }

  static images: DG.LruCache = new DG.LruCache();

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    if (w < 5 || h < 5) return;
    var bytes = gridCell.cell.value;
    if (bytes?.constructor === Uint8Array)
      Paint.pngImage(g, new Rect(x, y, w, h), bytes);
  }

  onDoubleClick(gridCell: GridCell, e: MouseEvent): void {
    Utils.openFileBytes({
      accept: 'image/png',
      open: (bytes) => {
        gridCell.cell.value = bytes;
      }
    });
  }
}