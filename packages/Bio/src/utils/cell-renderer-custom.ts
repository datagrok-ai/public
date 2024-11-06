/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {getGridCellColTemp, CellRendererBackBase} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';

import {_package} from '../package';

/// NB! UNUSED LEGACY CODE
export class MacromoleculeCustomCellRenderer extends DG.GridCellRenderer {
  private readonly seqHelper: ISeqHelper;

  get name(): string { return 'sequence'; }

  get cellType(): string { return 'sequence'; }

  get defaultHeight(): number | null { return 30; }

  get defaultWidth(): number | null { return 230; }

  constructor() {
    super();
    this.seqHelper = _package.seqHelper;
  }

  getRendererBack(gridCell: DG.GridCell): CellRendererBackBase<string> | null {
    const [gridCol, tableCol, _temp] = getGridCellColTemp<string, any>(gridCell);
    let back: CellRendererBackBase<string> | null = null;
    if (this.seqHelper) {
      const sh = this.seqHelper.getSeqHandler(tableCol);
      if (sh.notation !== NOTATION.CUSTOM)
        throw new Error(`Unexpected notation: '${sh.notation}'.`);
      back = sh.getRendererBack(gridCol, tableCol);
    }
    return back;
  }

  override onMouseEnter(gridCell: DG.GridCell, e: MouseEvent) {
    const back = this.getRendererBack(gridCell);
    back?.onMouseEnter(gridCell, e);
  }

  override onMouseLeave(_gridCell: DG.GridCell, _e: MouseEvent) {
    // TODO: We get gridCell from another column here, so we can not get back object from the column rendered.
    ui.tooltip.hide();
    // const back = this.getRendererBack(gridCell);
    // back.onMouseLeave(gridCell, e);
  }

  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent) {
    const back = this.getRendererBack(gridCell);
    back?.onMouseMove(gridCell, e);
  }

  override onClick(gridCell: DG.GridCell, e: MouseEvent) {
    const back = this.getRendererBack(gridCell);
    back?.onClick(gridCell, e);
  }

  override onDoubleClick(gridCell: DG.GridCell, e: MouseEvent) {
    const back = this.getRendererBack(gridCell);
    back?.onDoubleClick(gridCell, e);
  }

  override onKeyDown(gridCell: DG.GridCell, e: KeyboardEvent) {
    const back = this.getRendererBack(gridCell);
    back?.onKeyDown(gridCell, e);
  }

  override onKeyPress(gridCell: DG.GridCell, e: KeyboardEvent) {
    const back = this.getRendererBack(gridCell);
    back?.onKeyPress(gridCell, e);
  }

  override render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, gridCell: DG.GridCell, cellStyle: DG.GridCellStyle) {
    const back = this.getRendererBack(gridCell);
    back?.render(g, x, y, w, h, gridCell, cellStyle);
  }
}
