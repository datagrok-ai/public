import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {HelmType, IHelmBio, ISeqMonomer, HelmMol, Mol} from '@datagrok-libraries/bio/src/helm/types';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {CellRendererBackBase, getGridCellColTemp} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types/monomer-library';
import {execMonomerHoverLinks} from '@datagrok-libraries/bio/src/monomer-works/monomer-hover';
import {MmcrTemps} from '@datagrok-libraries/bio/src/utils/cell-renderer-consts';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {
  DEFAULT_MACROMOLECULE_HIGHLIGHT_FILL, DEFAULT_MACROMOLECULE_HIGHLIGHT_STROKE,
  MACROMOLECULE_HIGHLIGHT_EVENT_ID, MACROMOLECULE_HIGHLIGHT_TEMP,
  MacromoleculeHighlightEntry, MacromoleculeHighlightEventArgs, macromoleculeHighlightColorToCss,
} from '@datagrok-libraries/bio/src/utils/macromolecule-highlight';

// Phase 4 (hwe migration): the cell renderer paints synchronously straight into
// the grid canvas via the standalone `@datagrok-libraries/hwe` library — no
// async SVG→image / ImageData round-trip. `renderToCanvas` fits + centers +
// clips the molecule into the cell bounds and returns the laid-out mol plus the
// model/draw boxes used for hover + highlight overlays.
import {HelmService as HweHelmService, bridgeMonomerLib, buildHelmMolView} from '@datagrok-libraries/hwe';
import type {CanvasBounds, IMonomerLibBaseLike} from '@datagrok-libraries/hwe';

import {getHoveredMonomerFromEditorMol, getSeqMonomerFromHelmAtom} from './get-hovered';

import {_package} from '../package';

/** Per-cell render bookkeeping for hover hit-testing + the highlight overlay. */
type HelmCellAux = {
  /** Laid-out molecule as a legacy-shaped view (atoms[].p in model space). */
  molView: HelmMol;
  /** Tight model-space bbox of the drawn content. */
  bBox: DG.Rect;
  /** Absolute (device-pixel) draw box the mol was painted into. */
  dBox: DG.Rect;
  /** Cell height in device pixels (single-atom hover threshold). */
  cellHDev: number;
};

export class HelmGridCellRendererBack extends CellRendererBackBase<string> {
  /**
   * Compat shim for the `measureCellRenderer` benchmark: hwe renders
   * synchronously with an internal LRU parse/layout cache, so this is always
   * `true` (the benchmark then waits on `grid.onAfterDrawContent`).
   */
  public readonly cacheEnabled: boolean = true;

  private sysMonomerLib: IMonomerLibBase | null = null;
  private helmHelper: IHelmHelper | null = null;

  /** hwe service per monomer-library source (a cell may carry an overridden lib). */
  private readonly serviceByLib: Map<string, HweHelmService> = new Map();
  /** Per-cell-value render bookkeeping (hover + highlight). */
  private auxList: Map<string, HelmCellAux> = new Map();

  constructor(gridCol: DG.GridColumn | null, tableCol: DG.Column<string>) {
    super(gridCol, tableCol, _package.logger);
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
    ]);

    this.subs.push(this.sysMonomerLib!.onChanged.subscribe(() => {
      // A library change invalidates every cached service + render.
      for (const svc of this.serviceByLib.values()) svc.destroy();
      this.serviceByLib.clear();
      this.auxList.clear();
      this.dirty = true;
      this.invalidateGrid();
    }));

    this.subs.push(grok.events.onCustomEvent(MACROMOLECULE_HIGHLIGHT_EVENT_ID)
      .subscribe((args: MacromoleculeHighlightEventArgs) => this.handleHighlightEvent(args)));

    this.dirty = true;
    this.invalidateGrid();
  }

  protected override destroy(): void {
    for (const svc of this.serviceByLib.values()) svc.destroy();
    this.serviceByLib.clear();
    this.auxList.clear();
    super.destroy();
  }

  private handleHighlightEvent(args: MacromoleculeHighlightEventArgs): void {
    if (!args || !this.tableCol || !this.tableCol.dataFrame) return;
    if (args.columnName !== this.tableCol.name) return;
    const df = this.tableCol.dataFrame;
    if (args.tableId && df.id !== args.tableId) return;
    if (args.tableName && df.name !== args.tableName) return;

    const map = (this.tableCol.temp[MACROMOLECULE_HIGHLIGHT_TEMP] ??=
      new Map<number, MacromoleculeHighlightEntry>()) as Map<number, MacromoleculeHighlightEntry>;
    if (args.monomers == null || args.monomers.length === 0)
      map.delete(args.rowIdx);
    else {
      map.set(args.rowIdx, {
        monomers: args.monomers.slice(),
        fillColor: args.fillColor,
        strokeColor: args.strokeColor,
      });
    }
    this.invalidateGrid();
  }

  private getHighlightForRow(rowIdx: number | null): MacromoleculeHighlightEntry | null {
    if (rowIdx == null) return null;
    const map = this.tableCol.temp[MACROMOLECULE_HIGHLIGHT_TEMP] as
      Map<number, MacromoleculeHighlightEntry> | undefined;
    return map?.get(rowIdx) ?? null;
  }

  protected getMonomerLib(): IMonomerLibBase {
    return this.tableCol.temp[MmcrTemps.overriddenLibrary] ?? this.sysMonomerLib;
  }

  /** Get (or lazily build) the hwe service bound to `lib`, keyed by its source. */
  private getService(lib: IMonomerLibBase): HweHelmService {
    const key = lib.source;
    let svc = this.serviceByLib.get(key);
    if (svc === undefined) {
      svc = new HweHelmService({monomerLib: bridgeMonomerLib(lib as unknown as IMonomerLibBaseLike)});
      this.serviceByLib.set(key, svc);
    }
    return svc;
  }

  protected override reset(): void {
    super.reset();
    this.auxList = new Map<string, HelmCellAux>();
    this.invalidateGrid();
  }

  /** Synchronous paint straight into the grid canvas via hwe `renderToCanvas`. */
  public render(
    g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
    gridCell: DG.GridCell, _cellStyle: DG.GridCellStyle,
  ): void {
    const helm = gridCell.cell?.value;
    if (!helm || this.sysMonomerLib === null) return;
    if (w < 5 || h < 5) return;

    const dpr = window.devicePixelRatio;
    const svc = this.getService(this.getMonomerLib());

    g.save();
    try {
      // The grid canvas works in device pixels after resetTransform(); use the
      // cell's own bounds when painting onto the grid's main canvas (robust to
      // scroll), otherwise the passed x/y (off-screen tooltip canvas, etc.).
      g.resetTransform();
      const onGrid = !!(this.gridCol?.dart && this.gridCol.grid?.canvas === g.canvas);
      const bx = (onGrid ? gridCell.bounds.x : x) * dpr;
      const by = (onGrid ? gridCell.bounds.y : y) * dpr;
      const bounds: CanvasBounds = {x: bx, y: by, width: w * dpr, height: h * dpr};

      // No backColor → the grid keeps its own cell background (selection, etc.).
      // maxScale: Infinity → the drawing scales UP to fill the cell (hwe's
      // default caps at native size, leaving small molecules tiny in big cells).
      // showIds:false → no continuous-id numbers (keeps the grid uncluttered).
      // showAnnotations:true → still draw the 5'/3' direction + ss/as strand-type
      // markers; those are gated separately from the id numbers.
      const res = svc.renderToCanvas(g, helm, bounds, {
        maxScale: Infinity, showIds: false, showAnnotations: true,
      });

      const aux: HelmCellAux = {
        molView: buildHelmMolView(res.mol) as unknown as HelmMol,
        bBox: new DG.Rect(res.bbox.x, res.bbox.y, res.bbox.width, res.bbox.height),
        dBox: new DG.Rect(res.drawBox.x, res.drawBox.y, res.drawBox.width, res.drawBox.height),
        cellHDev: h * dpr,
      };
      this.auxList.set(helm, aux);

      // drawMolToContext restores its own transform, so we are back in
      // device-pixel space here for the overlay.
      this.drawHighlightOverlay(g, aux, gridCell.tableRowIndex);
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      this.logger.error(errMsg, undefined, errStack);
    } finally {
      g.restore();
    }
    this._onRendered.next();
  }

  // Overlays a translucent ring around each highlighted monomer on top of the
  // already-drawn helm. Highlight data comes from
  // `tableCol.temp[MACROMOLECULE_HIGHLIGHT_TEMP]` for the given row. Atom
  // positions (model space) map to device pixels via dBox / bBox.
  private drawHighlightOverlay(
    ctx: CanvasRenderingContext2D, aux: HelmCellAux, rowIdx: number | null,
  ): void {
    const entry = this.getHighlightForRow(rowIdx);
    if (!entry || !entry.monomers || entry.monomers.length === 0) return;
    const atoms = aux.molView?.atoms;
    if (!atoms || atoms.length === 0) return;
    if (aux.bBox.width <= 0 || aux.bBox.height <= 0) return;

    const sx = aux.dBox.width / aux.bBox.width;
    const sy = aux.dBox.height / aux.bBox.height;
    const sAvg = (sx + sy) / 2;

    const svgMonomerR = this.estimateMonomerRadiusModel(aux);
    const markerR = Math.max(6, svgMonomerR * sAvg * 0.7);
    const lineWidth = Math.max(1.5, markerR * 0.35);

    const fillCss = macromoleculeHighlightColorToCss(entry.fillColor ?? DEFAULT_MACROMOLECULE_HIGHLIGHT_FILL, 0.2);
    const strokeCss = macromoleculeHighlightColorToCss(entry.strokeColor ?? DEFAULT_MACROMOLECULE_HIGHLIGHT_STROKE);

    ctx.save();
    try {
      ctx.lineWidth = lineWidth;
      ctx.strokeStyle = strokeCss;
      ctx.fillStyle = fillCss;
      for (const idx of entry.monomers) {
        const a = atoms[idx];
        if (!a || !a.p) continue;
        const fx = aux.dBox.x + (a.p.x - aux.bBox.x) * sx;
        const fy = aux.dBox.y + (a.p.y - aux.bBox.y) * sy;
        ctx.beginPath();
        ctx.arc(fx, fy, markerR, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
      }
    } finally {
      ctx.restore();
    }
  }

  // Best-effort monomer glyph radius in model (bBox) coordinates from the median
  // bond length (glyph width ≈ bond length), else an area-density estimate.
  private estimateMonomerRadiusModel(aux: HelmCellAux): number {
    const bonds = (aux.molView as any)?.bonds as
      Array<{a1?: {p?: {x: number, y: number}}, a2?: {p?: {x: number, y: number}}}> | undefined;
    const lens: number[] = [];
    if (bonds && bonds.length > 0) {
      for (const b of bonds) {
        const p1 = b?.a1?.p; const p2 = b?.a2?.p;
        if (!p1 || !p2) continue;
        const d = Math.hypot(p1.x - p2.x, p1.y - p2.y);
        if (d > 0) lens.push(d);
      }
    }
    if (lens.length > 0) {
      lens.sort((p, q) => p - q);
      return lens[Math.floor(lens.length / 2)] * 0.55;
    }
    const n = Math.max(1, aux.molView?.atoms?.length ?? 1);
    return Math.sqrt((aux.bBox.width * aux.bBox.height) / n) * 0.45;
  }

  onMouseMove(gridCell: DG.GridCell, e: MouseEvent): void {
    if (!gridCell.cell?.value || !!e.buttons) return;
    const aux = this.auxList.get(gridCell.cell.value);
    if (!aux) return;

    const logPrefix = `${this.toLog()}.onMouseMove()`;
    const dpr = window.devicePixelRatio;

    /** {@link gridCell}.bounds */ const gcb = gridCell.bounds;
    const cX = (e.offsetX - gcb.x) * dpr;
    const cY = (e.offsetY - gcb.y) * dpr;

    // Map the cursor (device px, relative to the cell) into model space using
    // the absolute draw box. dBox is absolute, so subtract the cell origin.
    const dboxRelX = aux.dBox.x - gcb.x * dpr;
    const dboxRelY = aux.dBox.y - gcb.y * dpr;
    const mX = aux.bBox.x + (cX - dboxRelX) * aux.bBox.width / aux.dBox.width;
    const mY = aux.bBox.y + (cY - dboxRelY) * aux.bBox.height / aux.dBox.height;
    // Single-atom hover threshold (0.35·cellHeight in the legacy) mapped to model space.
    const cellHModel = aux.cellHDev * aux.bBox.height / aux.dBox.height;

    const editorMol: Mol<HelmType, IHelmBio> | null = aux.molView as unknown as Mol<HelmType, IHelmBio>;
    if (!editorMol) {
      this.logger.warning(`${logPrefix}, editorMol of the cell not found.`);
      return;
    }
    const hoveredAtom = getHoveredMonomerFromEditorMol(mX, mY, aux.molView, cellHModel);

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
      ui.tooltip.hide();
    }
    execMonomerHoverLinks(gridCell, seqMonomer);

    e.preventDefault();
    e.stopPropagation();
  }

  onMouseLeave(gridCell: DG.GridCell, e: MouseEvent): void {
    execMonomerHoverLinks(gridCell, null);
  }

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
    try {
      HelmGridCellRendererBack.getOrCreate(gridCell).onMouseMove(gridCell, e);
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
    }
  }

  override onMouseLeave(gridCell: DG.GridCell, e: MouseEvent) {
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
