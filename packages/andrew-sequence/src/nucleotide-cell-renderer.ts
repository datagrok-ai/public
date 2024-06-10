/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as nu from './nucleotide-utils';

export default class NucleotideBoxCellRenderer extends DG.GridCellRenderer {
  get name() {return 'Nucleotide cell renderer';}
  get cellType() {return nu.NUCLEOTIDE_SEMTYPE;}
  get defaultWidth() {return 200;}
  get defaultHeight() {return 30;}

  render(ctx: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle): void {
    const sequence: string[] = gridCell.cell.value;
    if (!ctx) return;

    // ctx.font = '11px courier';
    ctx.font = cellStyle.font;

    const paddingX = 3;
    const paddingY = 3;
    const charOffsetY = 13;
    const x2 = x + w - paddingX - 2;

    ctx.save();
    ctx.beginPath();
    ctx.rect(x + paddingX, y + paddingY, w - paddingX - 1, h - paddingY * 2);
    ctx.clip();

    let curX = x + paddingX;
    let curY = y + paddingY + charOffsetY;
    for (const nucleotideRaw of sequence) {
      const nucleotide = nucleotideRaw.toUpperCase() as nu.Nucleotide;
      ctx.fillStyle = nu.NUCLEOTIDE_COLORS[nucleotide] ?? 'black';
      ctx.fillText(nucleotideRaw, curX, curY);
      curX += ctx.measureText(nucleotideRaw).width + 1;
      if (curX >= x2) {
        curX = x + paddingX;
        curY += charOffsetY;
      }
    }

    ctx.restore();
  }
}
