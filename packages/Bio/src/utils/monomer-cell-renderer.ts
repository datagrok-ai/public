import * as DG from 'datagrok-api/dg';
import {GridCell} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ALPHABET, monomerToShort} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {GAP_SYMBOL, TAGS as bioTAGS,} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {MONOMER_RENDERER_TAGS} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import {getGridCellColTemp} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';

import {CellRendererWithMonomerLibBackBase} from './monomer-cell-renderer-base';
import * as C from './constants';
import {undefinedColor} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {HelmTypes} from '@datagrok-libraries/js-draw-lite/src/types/org';

const Tags = new class {
  tooltipHandlerTemp = 'tooltip-handler.Monomer';
}();

const svgMolOptions = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};

export class MonomerCellRendererBack extends CellRendererWithMonomerLibBackBase {
  constructor(gridCol: DG.GridColumn | null, tableCol: DG.Column) {
    super(gridCol, tableCol);
  }

  render(g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ): void {
    g.save();
    try {
      if (gridCell.gridRow < 0) return;
      const applyToBackground = gridCell.cell?.column && gridCell.cell.column.getTag(MONOMER_RENDERER_TAGS.applyToBackground) === 'true';

      g.font = `12px monospace`;
      g.textBaseline = 'middle';
      g.textAlign = 'center';

      const symbol: string = gridCell.cell.value;
      if (!symbol || symbol == GAP_SYMBOL) return;

      let color = undefinedColor;
      if (this.monomerLib) {
        const alphabet = this.tableCol.getTag(bioTAGS.alphabet);
        const biotype = alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA;
        color = this.monomerLib.getMonomerTextColor(biotype, symbol);
      }

      //cell width of monomer should dictate how many characters can be displayed
      // for width 40, 6 characters can be displayed (0.15 is 6 / 40)
      const maxChars = Math.max(2, Math.floor(w * 0.15));
      g.fillStyle = color;
      if (applyToBackground) {
        g.fillRect(x, y, w, h);
        g.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(DG.Color.fromHtml(color)));
      }
      g.fillText(monomerToShort(symbol, maxChars), x + (w / 2), y + (h / 2), w);
    } finally {
      g.restore();
    }
  }

  override onMouseMove(gridCell: GridCell, e: MouseEvent) {
    if (
      gridCell.grid.dart != this.gridCol?.grid.dart || gridCell.gridColumn.dart != this.gridCol?.dart ||
      !gridCell.tableColumn || !gridCell.isTableCell
    ) return false;

    const alphabet = gridCell.tableColumn.getTag(bioTAGS.alphabet) as ALPHABET;
    const monomerName = gridCell.cell.value;
    const canvasClientRect = gridCell.grid.canvas.getBoundingClientRect();
    const x1 = gridCell.bounds.right + canvasClientRect.left - 4;
    const y1 = gridCell.bounds.bottom + canvasClientRect.top - 4;

    if (monomerName == GAP_SYMBOL) {
      ui.tooltip.show(ui.divText('gap'), x1, y1);
      return true;
    }

    if (!this.monomerLib) {
      ui.tooltip.show(ui.divText('Monomer library is not available.'), x1, y1);
      return true;
    }

    const biotype = alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA;
    const tooltipEl = this.monomerLib.getTooltip(biotype, monomerName);
    ui.tooltip.show(tooltipEl, x1, y1);

    return true; // To prevent default tooltip behaviour
  }

  override async awaitRendered(timeout: number = 10000, reason: string = `${timeout} timeout`): Promise<void> {
    return Promise.resolve();
  }

  static getOrCreate(gridCell: DG.GridCell): MonomerCellRendererBack {
    const [gridCol, tableCol, temp] =
      getGridCellColTemp<string, MonomerCellRendererBack>(gridCell);

    let res: MonomerCellRendererBack = temp.rendererBack;
    if (!res) res = temp.rendererBack = new MonomerCellRendererBack(gridCol, tableCol);
    return res;
  }
}

export class MonomerCellRenderer extends DG.GridCellRenderer {
  get name(): string { return C.SEM_TYPES.MONOMER; }

  get cellType(): string { return C.SEM_TYPES.MONOMER; }

  get defaultHeight(): number { return 15; }

  get defaultWidth(): number { return 40; }

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
    _cellStyle: DG.GridCellStyle
  ): void {
    const back = MonomerCellRendererBack.getOrCreate(gridCell);
    back.render(g, x, y, w, h, gridCell, _cellStyle);
  }

  onMouseMove(gridCell: GridCell, e: MouseEvent) {
    const back = MonomerCellRendererBack.getOrCreate(gridCell);
    back.onMouseMove(gridCell, e);
  }
}
