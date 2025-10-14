import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {HelmType, IHelmBio, ISeqMonomer, Mol} from '@datagrok-libraries/bio/src/helm/types';
import {
  CellRendererBackAsyncBase, RenderServiceBase
} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';
import {HelmAux, HelmProps, HelmServiceBase} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {ILogger} from '@datagrok-libraries/bio/src/utils/logger';
import {getGridCellColTemp} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/index';
import {execMonomerHoverLinks} from '@datagrok-libraries/bio/src/monomer-works/monomer-hover';
import {MmcrTemps} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {getHoveredMonomerFromEditorMol, getSeqMonomerFromHelmAtom} from './get-hovered';

import {_package, getHelmService} from '../package';

export class HelmGridCellRendererBack extends CellRendererBackAsyncBase<HelmProps, HelmAux> {
  private _auxList: Map<string, HelmAux | null>;

  private sysMonomerLib: IMonomerLibBase | null = null;
  private helmHelper: IHelmHelper | null = null;
  private helmRenderService: HelmServiceBase | null = null;

  // eslint-disable-next-line max-len
  private readonly uuid: string = wu.repeat(1).map(() => Math.floor((Math.random() * 36)).toString(36)).take(4).toArray().join('');

  constructor(
    gridCol: DG.GridColumn | null,
    tableCol: DG.Column<string>,
  ) {
    super(gridCol, tableCol, _package.logger, true);
  }

  public async init(): Promise<void> {
    await Promise.all([
      (async () => {
        const libHelper = await getMonomerLibHelper();
        this.sysMonomerLib = libHelper.getMonomerLib();
      })(),
      (async () => {
        this.helmHelper = await getHelmHelper();
      })(),
      (async () => {
        this.helmRenderService = await getHelmService();
      })(),
    ]);

    this.subs.push(this.sysMonomerLib!.onChanged.subscribe(() => {
      this.dirty = true;
      this.invalidateGrid();
    }));

    this.dirty = true;
    this.invalidateGrid();
  }

  protected getMonomerLib(): IMonomerLibBase {
    return this.tableCol.temp[MmcrTemps.overriddenLibrary] ?? this.sysMonomerLib;
  }

  protected override reset(): void {
    const logPrefix = `${this.toLog()}.reset()`;
    this.logger.debug(`${logPrefix}, start`);
    super.reset();
    this._auxList = new Map<string, HelmAux | null>();
    this.invalidateGrid();
    this.logger.debug(`${logPrefix}, end`);
  }

  protected override getRenderService(): RenderServiceBase<HelmProps, HelmAux> | null {
    return this.helmRenderService;
  }

  protected override getRenderTaskProps(
    gridCell: DG.GridCell, backColor: number, width: number, height: number,
  ): HelmProps {
    const monomerLib = this.getMonomerLib();
    return new HelmProps(gridCell.cell.value, monomerLib, backColor, width, height);
  }

  protected override storeAux(gridCell: DG.GridCell, aux: HelmAux): void {
    const logPrefix = `${this.toLog()}.storeAux()`;
    this.logger.debug(`${logPrefix}, start`);
    if (!!gridCell.cell?.value)
      this._auxList.set(gridCell.cell.value, aux);
  }

  /** Renders cell from image data (cache), returns true to update the cell by service.
   * @param gridCellBounds {DG.Rect} Grid cell bounds
   */
  protected override renderCellImageData(
    gridCtx: CanvasRenderingContext2D, gridCellBounds: DG.Rect, gridCell: DG.GridCell, cellImageData: ImageData,
  ): boolean {
    // super.renderCellImageData(gridCtx, gridCellBounds, gridCell, cellImageData);
    // return true;

    const gcb = gridCellBounds;
    const dpr = window.devicePixelRatio;
    if (!gridCell.cell?.value)
      return false;
    const aux = this._auxList.get(gridCell.cell.value);
    if (!aux)
      return true;

    const [cellWidth, cellHeight] = [gcb.width * dpr - 2, gcb.height * dpr - 2];
    const cellDScale = Math.min(0.95 * cellWidth / aux.dBox.width, 0.95 * cellHeight / aux.dBox.height);
    const [cellDWidth, cellDHeight] = [cellDScale * aux.dBox.width, cellDScale * aux.dBox.height];
    const cellDBox = new DG.Rect((cellWidth - cellDWidth) / 2, (cellHeight - cellDHeight) / 2, cellDWidth, cellDHeight);


    // Draw cell image data to scale it with drawImage() while transform()
    const cellCanvas = ui.canvas(cellImageData.width, cellImageData.height);
    const cellCtx = cellCanvas.getContext('2d')!;
    cellCtx.putImageData(cellImageData, 0, 0);

    // // Get bbox canvas
    // const bBoxImageData = cellCtx.getImageData(aux.bBox.x, aux.bBox.y, aux.bBox.width, aux.bBox.height);
    // const bBoxCanvas = ui.canvas(aux.bBox.width, aux.bBox.height);
    // const bBoxCtx = bBoxCanvas.getContext('2d')!;
    // bBoxCtx.putImageData(bBoxImageData, 0, 0);
    //
    const fitCanvasWidth = gridCellBounds.width * dpr - 2;
    const fitCanvasHeight = gridCellBounds.height * dpr - 2;
    //
    // const bBoxRatio = Math.max(aux.bBox.width / aux.cBox.width, aux.bBox.height / aux.cBox.height);
    // const bBoxFullWidth = aux.bBox.width / bBoxRatio;
    // const bBoxFullHeight = aux.bBox.height / bBoxRatio;
    //
    // const fitScale = Math.min(fitCanvasWidth / bBoxFullWidth, fitCanvasHeight / bBoxFullHeight);
    //
    // // Relative shift of bbox center to cbox center
    // const bBoxShiftHR = (aux.bBox.left + aux.bBox.width / 2 - aux.cBox.width / 2) / aux.cBox.width;
    // const bBoxShiftVR = (aux.bBox.top + aux.bBox.height / 2 - aux.cBox.height / 2) / aux.cBox.height;
    //
    // const bBoxFitLeft = fitCanvasWidth / 2 - bBoxCanvas.width * fitScale * (1 - 2 * bBoxShiftHR) / 2;
    // const bBoxFitTop = fitCanvasHeight / 2 - bBoxCanvas.height * fitScale * (1 - 2 * bBoxShiftVR) / 2;
    //
    const fitCanvas = ui.canvas(fitCanvasWidth, fitCanvasHeight);
    const fitCtx = fitCanvas.getContext('2d')!;
    // fitCtx.fillStyle = '#FFFFA0';
    // fitCtx.fillRect(0, 0, fitCanvasWidth, fitCanvasHeight);
    // fitCtx.transform(fitScale, 0, 0, fitScale, bBoxFitLeft, bBoxFitTop);
    // fitCtx.drawImage(bBoxCanvas, 0, 0); // draw with scale transform
    fitCtx.drawImage(cellCanvas,
      aux.dBox.x, aux.dBox.y, aux.dBox.width, aux.dBox.height,
      cellDBox.x, cellDBox.y, cellDBox.width, cellDBox.height);

    const fitCanvasData = fitCtx.getImageData(0, 0, fitCanvasWidth, fitCanvasHeight);
    this.renderOnGrid(gridCtx, gridCellBounds, gridCell, fitCanvasData);

    return aux.cBox.width != cellWidth || aux.cBox.height != cellHeight; // request rendering
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.cell?.value || !this._auxList || !!e.buttons) return;
    const aux = this._auxList.get(gridCell.cell.value);
    if (!aux) return;

    const logPrefix = `${this.toLog()}.onMouseMove()`;
    const dpr = window.devicePixelRatio;

    /** {@link gridCell}.bounds */ const gcb = gridCell.bounds;
    const cX = (e.offsetX - gcb.x) * dpr;
    const cY = (e.offsetY - gcb.y) * dpr;

    const bX = aux.bBox.x + (cX - aux.dBox.x) * aux.bBox.width / aux.dBox.width;
    const bY = aux.bBox.y + (cY - aux.dBox.y) * aux.bBox.height / aux.dBox.height;

    const editorMol: Mol<HelmType, IHelmBio> | null = aux.mol; // in bBox world
    if (!editorMol) {
      this.logger.warning(`${logPrefix}, editorMol of the cell not found.`);
      return; // The gridCell is not rendered yet
    }
    const hoveredAtom = getHoveredMonomerFromEditorMol(bX, bY, editorMol, gridCell.bounds.height);

    let seqMonomer: ISeqMonomer | null = null;
    if (hoveredAtom && this.helmHelper && gridCell.tableRowIndex != null) {
      const seqHandler = this.helmHelper.seqHelper.getSeqHandler(this.tableCol);
      const seqValue = seqHandler.getValue(gridCell.tableRowIndex);
      seqMonomer = getSeqMonomerFromHelmAtom(seqValue, hoveredAtom);
      const monomerLib = this.tableCol.temp[MmcrTemps.overriddenLibrary] ?? _package.monomerLib;
      const tooltipEl = monomerLib ? monomerLib.getTooltip(seqMonomer.biotype, seqMonomer.symbol) :
        ui.divText('Monomer library is not available');
      ui.tooltip.show(tooltipEl, e.x + 16, e.y + 16);
    } else {
      // Tooltip for missing monomers
      ui.tooltip.hide();
    }
    execMonomerHoverLinks(gridCell, seqMonomer);

    e.preventDefault();
    e.stopPropagation();

  }

  onMouseLeave(gridCell: DG.GridCell, e: MouseEvent): void {
    execMonomerHoverLinks(gridCell, null);
  }

  // public override render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
  //   gridCell: DG.GridCell, cellStyle: DG.GridCellStyle
  // ) {
  //   super.render(g, x, y, w, h, gridCell, cellStyle);
  // }

  // -- Handle events --

  static getOrCreate(gridCell: DG.GridCell): HelmGridCellRendererBack {
    const [gridCol, tableCol, temp] =
      getGridCellColTemp<string, HelmGridCellRendererBack>(gridCell);

    let res: HelmGridCellRendererBack = temp.rendererBack;
    if (!res) {
      res = temp.rendererBack = new HelmGridCellRendererBack(gridCol, tableCol);
      res.init().then(() => {});
    }
    return res;
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
    }
  }

  override onMouseLeave(gridCell: DG.GridCell, e: MouseEvent) {
    const logPrefix = `Helm: HelmGridCellRenderer.onMouseLeave()`;
    try {
      HelmGridCellRendererBack.getOrCreate(gridCell).onMouseLeave(gridCell, e);
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
