import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {
  CellRendererBackAsyncBase, RenderServiceBase
} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';
import {HelmProps} from '@datagrok-libraries/bio/src/viewers/helm-service';

import {IEditorMol, ISeqMonomer} from '../helm-monomer-placer';
import {getHoveredMonomerFromEditorMol} from './get-hovered';
import {_getHelmService} from '../package-utils';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {getGridCellRendererBack} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types/index';

import {_package, getMonomerLib} from '../package';

class WrapLogger implements ILogger {
  constructor(
    private readonly base: ILogger
  ) {}

  error(message: any, params?: object | undefined, stackTrace?: string | undefined): void {
    this.base.error(message, params, stackTrace);
  }

  warning(message: string, params?: object | undefined): void {
    this.base.warning(message, params);
  }

  info(message: string, params?: object | undefined): void {
    // this.base.info(message, params);
  }

  debug(message: string, params?: object | undefined): void {
    // this.base.debug(message, params);
  }
}

class HelmGridCellRendererBack extends CellRendererBackAsyncBase<HelmProps> {
  private readonly monomerLib: IMonomerLib;

  constructor(
    gridCol: DG.GridColumn | null,
    tableCol: DG.Column<string>,
  ) {
    super(gridCol, tableCol, new WrapLogger(_package.logger) /* _package.logger */, true);
    this.monomerLib = getMonomerLib();
    this.subs.push(this.monomerLib.onChanged.subscribe(this.monomerLibOnChanged.bind(this)));
  }

  protected override reset(): void {
    super.reset();
    this._editorMolList = new Array<IEditorMol | null>(this.tableCol.length).fill(null);
    if (this.gridCol && this.gridCol.dart) this.gridCol.grid?.invalidate();
  }

  protected override getRenderService(): RenderServiceBase<HelmProps> {
    return _getHelmService();
  }

  protected override getRenderTaskProps(gridCell: DG.GridCell, dpr: number): HelmProps {
    return new HelmProps(
      gridCell,
      gridCell.grid.props.backColor,
      gridCell.gridColumn.width * dpr,
      gridCell.grid.props.rowHeight * dpr);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    /* Can not do anything without tableColumn containing temp */
    let tableCol: DG.Column | null = null;
    try { tableCol = gridCell.tableColumn; } catch { }
    if (!tableCol) return;

    /** {@link gridCell}.bounds */ const gcb = gridCell.bounds;
    const argsX = e.offsetX - gcb.x;
    const argsY = e.offsetY - gcb.y;

    const editorMol: IEditorMol | null = this.getEditorMol(gridCell.tableRowIndex!);
    if (!editorMol) return; // The gridCell is not rendered yet
    let seqMonomer: ISeqMonomer | null;
    seqMonomer = getHoveredMonomerFromEditorMol(argsX, argsY, gridCell, editorMol);

    if (seqMonomer) {
      ui.tooltip.show(ui.divV([
        ui.divText(`Monomer '${seqMonomer.symbol}' not found.`),
        ui.divText('Open the Context Panel, then expand Manage Libraries'),
      ]), e.x + 16, e.y + 16);
    } else
      ui.tooltip.hide();
  }

  // -- Handle events --
  private monomerLibOnChanged(_value: any): void {
    this.reset();
    this.invalidateGrid();
  }


  static getOrCreate(gridCell: DG.GridCell): HelmGridCellRendererBack {
    const [gridCol, tableCol, temp] =
      getGridCellRendererBack<string, HelmGridCellRendererBack>(gridCell);

    let res: HelmGridCellRendererBack = temp['rendererBack'];
    if (!res) res = temp['rendererBack'] = new HelmGridCellRendererBack(gridCol, tableCol);
    return res;
  }

  // -- Mol --

  private _editorMolList: (IEditorMol | null)[];

  private setEditorMol(tableRowIndex: number, editorMol: IEditorMol): void {
    this._editorMolList[tableRowIndex] = editorMol;
  }

  private getEditorMol(tableRowIndex: number): IEditorMol | null {
    return this._editorMolList[tableRowIndex];
  }
}

export class HelmGridCellRenderer extends DG.GridCellRenderer {
  get name() { return 'helm'; }

  get cellType() { return 'helm'; }

  get defaultWidth(): number | null { return 400; }

  get defaultHeight(): number | null { return 100; }

  override onMouseMove(gridCell: DG.GridCell, e: MouseEvent) {
    const logPrefix = `Helm: HelmGridCellRenderer.onMouseMove()`;
    try {
      HelmGridCellRendererBack.getOrCreate(gridCell).onMouseMove(gridCell, e);
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
    } finally {
      e.preventDefault();
      e.stopPropagation();
    }
  }

  override render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  ): void {
    HelmGridCellRendererBack.getOrCreate(gridCell).render(g, x, y, w, h, gridCell, cellStyle);
  }
}
