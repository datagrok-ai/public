
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export class GasteigerPngRenderer extends DG.GridCellRenderer {
  get name(): string {return 'customGasteigerPNG';}

  get cellType(): string {return 'customGasteigerPNG';}

  get defaultHeight(): number {return 200;}

  get defaultWidth(): number {return 200;}

  /**
     * Cell renderer function.
     *
     * @param {CanvasRenderingContext2D} g Canvas rendering context.
     * @param {number} x x coordinate on the canvas.
     * @param {number} y y coordinate on the canvas.
     * @param {number} w width of the cell.
     * @param {number} h height of the cell.
     * @param {DG.GridCell} gridCell Grid cell.
     * @param {DG.GridCellStyle} _cellStyle Cell style.
     */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    _cellStyle: DG.GridCellStyle,
  ): void {
    if (gridCell.gridRow < 0 || !gridCell.cell.value) return;

    const pngString: string = gridCell.cell.value;
    const img = new Image(w-2, h-2);

    img.src = 'data:image/png;base64,' + pngString;
    img.onload = function() {
      g.drawImage(img, x+2, y+2, w-2, h-2);
    };
  }
}
