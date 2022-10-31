import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import * as bio from '@datagrok-libraries/bio';

export class PdbRenderer extends DG.GridCellRenderer {
  get name(): string { return 'protein'; }

  get cellType(): string { return 'protein'; }

  get defaultHeight(): number { return 30; }

  get defaultWidth(): number { return 230; }

  onClick(gridCell: DG.GridCell, e: MouseEvent): void {
    gridCell.cell.column.temp['current-word'] = gridCell.cell.value;
    gridCell.grid.invalidate();
  }

  /**
   * Cell renderer function.
   *
   * @param {CanvasRenderingContext2D} g Canvas rendering context.
   * @param {number} x x coordinate on the canvas.
   * @param {number} y y coordinate on the canvas.
   * @param {number} w width of the cell.
   * @param {number} h height of the cell.
   * @param {DG.GridCell} gridCell Grid cell.
   * @param {DG.GridCellStyle} cellStyle Cell style.
   * @memberof AlignedSequenceCellRenderer
   */
  render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell,
    cellStyle: DG.GridCellStyle
  ): void {

    const indexStart = gridCell.cell.value.indexOf('TITLE');
    const indexEnd = gridCell.cell.value.indexOf('\n', indexStart + 1);
    const name = gridCell.cell.value.substring(indexStart + 10, indexEnd);
    
    g.fillStyle = '#0095B6';
    g.fillText(name, x, y);

    g.restore();
    return;
  }
}
