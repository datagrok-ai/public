/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {GridCell} from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {ALPHABET, monomerToShort} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {BioTags, GAP_SYMBOL, MONOMER_MOTIF_SPLITTER, TAGS as bioTAGS,} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {MONOMER_RENDERER_TAGS} from '@datagrok-libraries/bio/src/utils/cell-renderer';
import {getGridCellColTemp} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';

import {CellRendererWithMonomerLibBackBase} from './monomer-cell-renderer-base';
import * as C from './constants';
import {undefinedColor} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {HelmType, HelmTypes, PolymerType, PolymerTypes} from '@datagrok-libraries/js-draw-lite/src/types/org';
import {polymerTypeToHelmType} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';

const Tags = new class {
  tooltipHandlerTemp = 'tooltip-handler.Monomer';
}();

const DASH_GAP_SYMBOL = '-';

export class MonomerCellRendererBack extends CellRendererWithMonomerLibBackBase {
  constructor(gridCol: DG.GridColumn | null, tableCol: DG.Column) {
    super(gridCol, tableCol);
  }

  render(g: CanvasRenderingContext2D,
    x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ): void {
    g.save();
    // clip the cell
    g.beginPath();
    g.rect(x, y, w, h);
    g.clip();
    try {
      if (!gridCell.isTableCell) return;
      const applyToBackground = gridCell.cell?.column && gridCell.cell.column?.dart && gridCell.cell.column.getTag(MONOMER_RENDERER_TAGS.applyToBackground) === 'true';

      g.font = `12px monospace`;
      g.textBaseline = 'middle';
      g.textAlign = 'left';

      let value: string = gridCell.cell.value;
      if (!value || value === GAP_SYMBOL)
        value = DASH_GAP_SYMBOL;
      const symbols = value.split(MONOMER_MOTIF_SPLITTER).map((s) => !s || s === GAP_SYMBOL ? DASH_GAP_SYMBOL : s.trim());
      //cell width of monomer should dictate how many characters can be displayed
      // for width 40, 6 characters can be displayed (0.15 is 6 / 40)
      const shortSymbols = symbols.map((s) => monomerToShort(s, Math.max(2, Math.floor(w * 0.15 / symbols.length))));
      const symbolWidths = shortSymbols.map((s) => g.measureText(s).width);
      const totalWidth = symbolWidths.reduce((a, b) => a + b, 0);
      const xOffset = (w - totalWidth) / 2;
      let xPos = x + xOffset;
      const alphabet = this.tableCol.getTag(bioTAGS.alphabet);
      const biotype = alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA;
      for (let i = 0; i < shortSymbols.length; i++) {
        const symbol: string = symbols[i];
        const actBioType: HelmType = this.getHelmType(gridCell, biotype);

        let textcolor = undefinedColor;
        let backgroundcolor = 'rgb(255, 255, 255)';
        if (this.monomerLib) {
          if (applyToBackground) {
            const colors = this.monomerLib.getMonomerColors(actBioType, symbol);
            textcolor = colors?.textcolor ?? textcolor;
            backgroundcolor = colors?.backgroundcolor ?? backgroundcolor;
          } else
            textcolor = this.monomerLib.getMonomerTextColor(actBioType, symbol);
        }
        if (applyToBackground && symbols.length == 1) {
          g.fillStyle = backgroundcolor;
          g.fillRect(x, y, w, h);
        }
        g.fillStyle = textcolor;
        g.fillText(shortSymbols[i], xPos, y + (h / 2), w);
        xPos += symbolWidths[i];
      }
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
    const monomerName: string = gridCell.cell.value;
    const canvasClientRect = gridCell.grid.canvas.getBoundingClientRect();
    const x1 = gridCell.bounds.right + canvasClientRect.left - 4;
    const y1 = gridCell.bounds.bottom + canvasClientRect.top - 4;

    if (!monomerName || monomerName == GAP_SYMBOL || monomerName == DASH_GAP_SYMBOL) {
      ui.tooltip.show(ui.divText('gap'), x1, y1);
      return true;
    }

    if (!this.monomerLib) {
      ui.tooltip.show(ui.divText('Monomer library is not available.'), x1, y1);
      return true;
    }

    const biotype = alphabet === ALPHABET.RNA || alphabet === ALPHABET.DNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA;
    const tooltipEls = monomerName.split(MONOMER_MOTIF_SPLITTER)
      .map((s) => {
        if (!s || s === GAP_SYMBOL || s === DASH_GAP_SYMBOL)
          return ui.divText('gap');
        const actBioType: HelmType = this.getHelmType(gridCell, biotype);
        return this.monomerLib!.getTooltip(actBioType, s);
      });
    const tooltipEl = ui.divH(tooltipEls, {style: {alignItems: 'top'}});
    // tooltip max width is 600px, so we need to shrink the canvases a bit if needed. by default, it is 250px
    const canvases = Array.from(tooltipEl.querySelectorAll('canvas'));
    if (canvases.length > 2) {
      const side = Math.floor(550 / canvases.length);
      canvases.forEach((c) => {
        c.style.setProperty('max-width', `${side}px`, 'important');
        c.style.setProperty('max-height', `${side}px`, 'important');
      });
    }
    ui.tooltip.show(tooltipEl, x1, y1);

    return true; // To prevent default tooltip behaviour
  }

  private getHelmType(gridCell: GridCell, defaultType: HelmType): HelmType {
    let biotype = defaultType;
    if ((gridCell.tableRowIndex ?? -1) > -1 && gridCell.tableColumn?.getTag(BioTags.polymerTypeColumnName)) {
      const ptColName = gridCell.tableColumn.getTag(BioTags.polymerTypeColumnName);
      const ptCol = gridCell.tableColumn.dataFrame?.col(ptColName);
      if (ptCol) {
        const ptrString = ptCol.get(gridCell.tableRowIndex!);
        if (ptrString && [PolymerTypes.BLOB, PolymerTypes.CHEM, PolymerTypes.G, PolymerTypes.PEPTIDE, PolymerTypes.RNA].includes(ptrString))
          biotype = polymerTypeToHelmType(ptrString as PolymerType);
      }
    }
    return biotype;
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
    e.preventDefault();
  }
}
