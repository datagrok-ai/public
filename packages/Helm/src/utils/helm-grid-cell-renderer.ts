import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {
  CellRendererBackAsyncBase, RenderServiceBase
} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';
import {IEditorMol} from '@datagrok-libraries/bio/src/types/helm-web-editor';
import {HelmAux, HelmProps} from '@datagrok-libraries/bio/src/viewers/helm-service';

import {ISeqMonomer} from '../helm-monomer-placer';
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

export class HelmGridCellRendererBack extends CellRendererBackAsyncBase<HelmProps, HelmAux> {
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

  protected override getRenderService(): RenderServiceBase<HelmProps, HelmAux> {
    return _getHelmService();
  }

  protected override getRenderTaskProps(gridCell: DG.GridCell, dpr: number): HelmProps {
    return new HelmProps(
      gridCell,
      gridCell.grid.props.backColor,
      gridCell.gridColumn.width * dpr - 2,
      gridCell.grid.props.rowHeight * dpr - 2);
  }

  protected override storeAux(gridCell: DG.GridCell, aux: HelmAux): void {
    if (gridCell.tableRowIndex !== null)
      this.setEditorMol(gridCell.tableRowIndex, aux.mol);
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    const logPrefix = `${this.toLog()}.onMouseMove()`;
    const dpr = window.devicePixelRatio;

    /** {@link gridCell}.bounds */ const gcb = gridCell.bounds;
    const argsX = (e.offsetX - gcb.x) * dpr;
    const argsY = (e.offsetY - gcb.y) * dpr;

    const editorMol: IEditorMol | null = this.getEditorMol(gridCell.tableRowIndex!);
    if (!editorMol) {
      this.logger.warning(`${logPrefix}, editorMol of the cell not found.`);
      return; // The gridCell is not rendered yet
    }
    const seqMonomer: ISeqMonomer | null = getHoveredMonomerFromEditorMol(argsX, argsY, gridCell, editorMol);

    if (seqMonomer) {
      const tooltipEl = this.monomerLib.getTooltip(seqMonomer.polymerType, seqMonomer.symbol);
      ui.tooltip.show(tooltipEl, e.x + 16, e.y + 16);
    } else {
      // Tooltip for missing monomers
      ui.tooltip.hide();
    }
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
