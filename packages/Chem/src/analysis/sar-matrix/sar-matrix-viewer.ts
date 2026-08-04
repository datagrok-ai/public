import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {merge} from 'rxjs';

import '../../../css/sar-matrix.css';
import {drawMoleculeToCanvas, drawRdKitMoleculeToOffscreenCanvas, getRdKitModule, getRdKitService}
  from '../../utils/chem-common-rdkit';
import {getMolSafe} from '../../utils/mol-creation_rdkit';
import {renderMolecule} from '../../rendering/render-molecule';
import {SCALING_METHODS} from '../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {rankMatrices, SarRankScheme} from './sar-matrix-ranking';
import {computeAllTransfers, Transfer, TransferSide, transferStats} from './sar-matrix-transfer';
import {runSarMatrix, SarMatrixParams} from './sar-matrix-run';
import {SarMatrix, SarMatrixCell} from './sar-matrix-types';

/** Transparent (alpha 0) so a drawn core has no white box — it blends with the card/pane. */
const CORE_BG_ARGB = 0x00000000;
const HEADER_ARGB = 0xFFF7F7F9;
/** Linear scheme for potency, matching how the platform colors activity/confidence scores:
 *  red at the low end, green at the high end of the scaled activity. */
const ACTIVITY_SCHEME = [DG.Color.red, DG.Color.green];
const CELL_W = 104;
const CELL_H = 76;
const HEADER_W = 78;
const HEADER_H = 46;
const CORE_W = 132;
const CORE_H = 60;
/** Small inline thumbnail of the transfer's benefiting substituent, in the statistics block. */
const BENEFIT_MOL_W = 62;
const BENEFIT_MOL_H = 34;
/** Cells grow to fill the pane, but never past this — beyond it the structures just float in space. */
const CELL_W_MAX = 210;
/** Navigator width + paddings + border-spacing, subtracted when fitting cells to the pane. */
const NAV_W = 260;
const TABLE_CHROME = 60;
/** Cell-tint alpha (0-255): solid for observed compounds, fainter for virtual predictions, and
 *  fainter still when a prediction rests on few observations. */
const REAL_ALPHA = 102;
const VIRTUAL_ALPHA = 46;
const VIRTUAL_ALPHA_MIN = 16;
/** Support (min of backing row/column observations) at/above which a virtual cell is at full alpha. */
const FULL_SUPPORT = 3;
/** Method label carried on every exported/predicted analog. */
const FREE_WILSON_METHOD = 'local Free-Wilson (row + column effects)';
/** A real cell's observed value "matches" its additive-model fit when the residual is within this
 *  fraction of the matrix's activity range; beyond it the cell is flagged as non-additive. */
const MATCH_FRACTION = 0.25;
/** Name of the running make-list table single-analog "Generate" appends to. */
const MAKELIST_NAME = 'SAR virtual analogs';

type AnalogPanelBuilder = () => HTMLElement;
/** Column-annotation modes offered in the matrix pane — they label the columns, never reorder them. */
const COLSORT_FREQUENCY = 'None';
const COLSORT_POTENCY = 'Potency';
const COLSORT_MW = 'Molecular weight';
const COLUMN_SORTS = [COLSORT_FREQUENCY, COLSORT_POTENCY, COLSORT_MW];

/** Properties that only reorder/recolor the already-assembled matrices. Changing one must NOT re-run
 *  fragmentation and decomposition — those cost seconds of RDKit worker time and produce identical
 *  matrices. Every other property (columns, scaling, cutoffs, prediction) does change the assembly. */
const RERANK_ONLY_PROPS = ['rankScheme', 'activityDirection'];

/** Which end of the activity scale is "more potent". Auto derives it from scaling (only −lg is
 *  higher-is-better); the explicit options cover pre-computed pIC50/pKi/%-inhibition left on `none`. */
const DIR_AUTO = 'Auto (from scaling)';
const DIR_HIGHER = 'Higher is better';
const DIR_LOWER = 'Lower is better';
const ACTIVITY_DIRECTIONS = [DIR_AUTO, DIR_HIGHER, DIR_LOWER];

/** Average molecular weight of a substituent, capping its `[*:n]` attachment points with H so the
 *  weight is that of the capped fragment. Cached (alignment is deterministic); Infinity on failure so
 *  unparseable substituents sort last. */
const mwCache = new Map<string, number>();
const MW_CACHE_MAX = 4000;

function substituentMW(smiles: string): number {
  const cached = mwCache.get(smiles);
  if (cached !== undefined)
    return cached;
  const rdkit = getRdKitModule();
  let mol = null;
  let mw = Number.POSITIVE_INFINITY;
  try {
    mol = rdkit.get_mol(smiles.replace(/\[\*:\d+\]/g, '[H]'));
    if (mol && mol.is_valid()) {
      const amw = JSON.parse(mol.get_descriptors()).amw;
      if (typeof amw === 'number' && Number.isFinite(amw))
        mw = amw;
    }
  } catch {
    // leave Infinity — the substituent sorts last
  } finally {
    mol?.delete();
  }
  if (mwCache.size >= MW_CACHE_MAX)
    mwCache.clear();
  mwCache.set(smiles, mw);
  return mw;
}


const ALIGN_OPTS = JSON.stringify({useCoordGen: true, allowRGroups: true, acceptFailure: true, alignOnly: true});

/** Aligned-molblock cache keyed by template + molecule. Alignment is deterministic, so a resize
 *  re-render reuses the layout instead of re-parsing every core and cell in RDKit. Cleared wholesale
 *  when it grows past the cap so it can't leak across many datasets in one session. */
const alignCache = new Map<string, string>();
const ALIGN_CACHE_MAX = 4000;
/** Per-core alignment templates, memoized by core SMILES. The layout is deterministic, so re-renders
 *  (selection, resize) reuse it instead of re-parsing the core in RDKit — and it keeps the alignCache
 *  keys stable. `null` (RDKit couldn't build a template) is cached too, so failures aren't retried. */
const templateCache = new Map<string, string | null>();
const TEMPLATE_CACHE_MAX = 4000;

function argbToRgba(argb: number): [number, number, number, number] {
  return [DG.Color.r(argb) / 255, DG.Color.g(argb) / 255, DG.Color.b(argb) / 255, DG.Color.a(argb) / 255];
}

/**
 * A molblock template, laid out canonically, that every cell/core is aligned to so the shared core
 * points the same way across the whole matrix. The `[*:n]` dummies are swapped for hydrogens so the
 * template is a plain substructure (following the MMP viewer's alignment).
 */
function buildAlignmentTemplate(coreSmiles: string): string | null {
  const cached = templateCache.get(coreSmiles);
  if (cached !== undefined)
    return cached;
  const rdkit = getRdKitModule();
  let mol = null;
  let template: string | null = null;
  try {
    mol = rdkit.get_mol(coreSmiles.replace(/\[\*:\d+\]/g, '[H]'));
    if (mol && mol.is_valid()) {
      mol.set_new_coords();
      template = mol.get_molblock();
    }
  } catch {
    template = null;
  } finally {
    mol?.delete();
  }
  if (templateCache.size >= TEMPLATE_CACHE_MAX)
    templateCache.clear();
  templateCache.set(coreSmiles, template);
  return template;
}

/** Aligned molblock for `molStr` oriented to the template, or the input unchanged if it can't align. */
function alignToTemplate(molStr: string, templateMolblock: string): string {
  const key = `${templateMolblock}${molStr}`;
  const cached = alignCache.get(key);
  if (cached !== undefined)
    return cached;
  const rdkit = getRdKitModule();
  let template = null;
  let mol = null;
  let result = molStr;
  try {
    template = rdkit.get_mol(templateMolblock);
    mol = rdkit.get_mol(molStr);
    if (template && mol && mol.is_valid()) {
      mol.generate_aligned_coords(template, ALIGN_OPTS);
      result = mol.get_molblock() || molStr;
    }
  } catch {
    result = molStr;
  } finally {
    template?.delete();
    mol?.delete();
  }
  if (alignCache.size >= ALIGN_CACHE_MAX)
    alignCache.clear();
  alignCache.set(key, result);
  return result;
}

/** Render a molecule onto a fresh canvas with the given ARGB background baked directly into the
 *  RDKit draw call — avoids relying on the grid cell renderer, whose putImageData overwrites any
 *  cell background. Straightens the depiction for a tidy standalone layout (used for substituents). */
function renderMoleculeOnColor(smiles: string, w: number, h: number, argb: number): HTMLElement {
  const canvas = ui.canvas(w, h);
  if (!smiles)
    return canvas;
  try {
    drawMoleculeToCanvas(0, 0, w, h, canvas, smiles, null,
      {normalizeDepiction: true, straightenDepiction: true}, null,
      {clearBackground: true, backgroundColour: argbToRgba(argb)});
  } catch (e) {
    // leave the canvas blank on a malformed structure
  }
  return canvas;
}

/**
 * Render a molecule aligned to `templateMolblock`, preserving the aligned coordinates. Goes through
 * `drawRdKitMoleculeToOffscreenCanvas` rather than `drawMoleculeToCanvas` because the latter always
 * re-straightens molecules that already carry coordinates, which would undo the alignment.
 */
function renderAlignedOnColor(molStr: string, w: number, h: number, argb: number,
  templateMolblock: string | null): HTMLElement {
  if (!molStr || !templateMolblock)
    return renderMoleculeOnColor(molStr, w, h, argb);
  const canvas = ui.canvas(w, h);
  const r = window.devicePixelRatio;
  const nW = Math.floor(w * r);
  const nH = Math.floor(h * r);
  canvas.width = nW;
  canvas.height = nH;
  canvas.style.width = `${w}px`;
  canvas.style.height = `${h}px`;
  let molCtx = null;
  try {
    molCtx = getMolSafe(alignToTemplate(molStr, templateMolblock), {}, getRdKitModule());
    if (!molCtx.mol)
      return canvas;
    const offscreen = new OffscreenCanvas(nW, nH);
    drawRdKitMoleculeToOffscreenCanvas(molCtx, nW, nH, offscreen, null,
      {clearBackground: true, backgroundColour: argbToRgba(argb)});
    // Blit the offscreen across directly. getImageData/putImageData would walk every pixel through
    // an unpremultiplied byte array on the CPU, which is both slower and lossy on anti-aliased edges.
    canvas.getContext('2d')!.drawImage(offscreen, 0, 0);
  } catch (e) {
    // leave the canvas blank on a malformed structure
  } finally {
    molCtx?.mol?.delete();
  }
  return canvas;
}

/** Rendered-depiction cache for the grid path, keyed by molecule + alignment template + rounded CSS
 *  size + baked background. `onCellRender` fires per visible cell on every repaint (scroll, selection,
 *  invalidate), so a scrolled/re-highlighted matrix reuses the produced canvas instead of re-running
 *  RDKit. Bounded by BYTES, not entries: the canvases are device-pixel sized, so at devicePixelRatio 2
 *  one body cell already costs a quarter of a megabyte and an entry cap would let the cache reach
 *  hundreds of megabytes. Least-recently-used entries are evicted first (Map insertion order, refreshed
 *  on every hit). */
const cellCanvasCache = new Map<string, HTMLCanvasElement>();
const CELL_CANVAS_CACHE_MAX_BYTES = 64 * 1024 * 1024;
let cellCanvasCacheBytes = 0;

/** Backing-store size of a cached canvas: 4 bytes (RGBA) per device pixel. */
function canvasBytes(canvas: HTMLCanvasElement): number {
  return canvas.width * canvas.height * 4;
}

/** Drop every rendered depiction and the aligned-molblock layouts behind them. Called when the
 *  matrices are rebuilt: none of the cached keys can be hit again, and the canvases hold tens of
 *  megabytes of device-pixel bitmaps until they are released. */
export function clearDepictionCaches(): void {
  cellCanvasCache.clear();
  cellCanvasCacheBytes = 0;
  alignCache.clear();
}

/** Depiction canvas for a grid cell/header. A null template gives a straightened standalone depiction
 *  (substituent headers); a template aligns to the shared core (cells, cores). The returned canvas is
 *  intrinsically device-pixel sized (the render helpers multiply by devicePixelRatio), so drawing it
 *  into the grid context at CSS size stays crisp under that context's own DPR scale. */
function cachedCellCanvas(smiles: string, template: string | null, w: number, h: number,
  argb: number): HTMLCanvasElement {
  const key = `${smiles}|${template ?? ''}|${Math.round(w)}x${Math.round(h)}|${argb}`;
  const cached = cellCanvasCache.get(key);
  if (cached !== undefined) {
    // Re-insert so this key moves to the end of the Map's insertion order: eviction takes the front,
    // which is then the least recently used entry.
    cellCanvasCache.delete(key);
    cellCanvasCache.set(key, cached);
    return cached;
  }
  const canvas = renderAlignedOnColor(smiles, w, h, argb, template) as HTMLCanvasElement;
  cellCanvasCache.set(key, canvas);
  cellCanvasCacheBytes += canvasBytes(canvas);
  // Never evict the entry just inserted — the caller is about to draw it.
  while (cellCanvasCacheBytes > CELL_CANVAS_CACHE_MAX_BYTES && cellCanvasCache.size > 1) {
    const oldest = cellCanvasCache.keys().next();
    if (oldest.done)
      break;
    const evicted = cellCanvasCache.get(oldest.value);
    cellCanvasCache.delete(oldest.value);
    if (evicted !== undefined)
      cellCanvasCacheBytes -= canvasBytes(evicted);
  }
  return canvas;
}

/** Resolve a Datagrok CSS palette variable (e.g. `--grey-6`) to a color string for canvas painting,
 *  memoized by name so per-cell paints don't re-run `getComputedStyle`. Falls back to `fallback` when
 *  the variable is unset. */
const cssColorCache = new Map<string, string>();
function cssColor(root: HTMLElement, name: string, fallback: string): string {
  const cached = cssColorCache.get(name);
  if (cached !== undefined)
    return cached;
  let value = '';
  try {
    value = getComputedStyle(root).getPropertyValue(name).trim();
  } catch {
    value = '';
  }
  const resolved = value || fallback;
  cssColorCache.set(name, resolved);
  return resolved;
}

const GRID_FONT = 'Roboto, "Segoe UI", sans-serif';

/** Canonical SMILES for display. The cell's structure is the raw molecule-column value, which for a
 *  molblock/SDF-format column is a V2000 block — this converts it so the "SMILES" row never shows raw
 *  molblock text. Returns '' when the structure can't be parsed. */
function toDisplaySmiles(molStr: string): string {
  let molCtx = null;
  try {
    molCtx = getMolSafe(molStr, {}, getRdKitModule());
    return molCtx.mol ? molCtx.mol.get_smiles() : '';
  } catch {
    return '';
  } finally {
    molCtx?.mol?.delete();
  }
}

/** Identity of a transfer's source core within its section (same- and cross-series kept apart). The
 *  nav grouping and the pane's sibling lookup both key on this, so the card "+N" and the pane dropdown
 *  can never disagree about which targets a source reaches. */
function transferSourceKey(t: Transfer): string {
  return `${t.crossMatrix}:${t.a.matrixIndex}:${t.a.rowIndex}`;
}

/** Leave-one-out R² as text: two decimals, or "<0" when the additive model is worse than the mean. */
function formatR2(r2: number): string {
  return r2 >= 0 ? r2.toFixed(2) : '<0';
}

/** The rendered matrix grid plus the per-render header state. `onCellRender` fires for every visible
 *  header cell on every repaint, so the visible-column list, the group boundaries and the column
 *  captions are computed once when the grid is built instead of per cell per repaint. */
interface MatrixGridState {
  grid: DG.Grid;
  colKeyToCi: Map<string, number>;
  matrix: SarMatrix;
  /** Visible column indices, in display order. */
  colIdxs: number[];
  /** Visible column indices that begin an R-position group — only those draw the position label. */
  firstOfGroup: Set<number>;
  /** `columnSortCaption` per visible column index; recomputed whenever the grid is rebuilt, which is
   *  what the Label control does, so a caption can never outlive the metric it was computed for. */
  captions: Map<number, string>;
}

export class SarMatrixViewer extends DG.JsViewer {
  moleculesColumnName: string;
  activityColumnName: string;
  scaling: string;
  activityDirection: string;
  fragmentCutoff: number;
  threshold: number;
  predictVirtual: boolean;
  rankScheme: string;

  private matrices: SarMatrix[] = [];
  /** SAR-transfer entries, listed as their own navigator section (the mockup's "SAR TRANSFER").
   *  Global: every correlated core pair across all matrices, including cross-series pairs. */
  private transfers: Transfer[] = [];
  /** Current navigator selection: either a matrix or a transfer entry. */
  private selKind: 'matrix' | 'transfer' = 'matrix';
  private selIndex = 0;
  /** "Vary" filter: show only this R-position's column group, or all when empty. */
  private varyPosition = '';
  /** Which metric annotates each substituent column (None / mean potency / molecular weight). Columns
   *  keep their as-assembled order — this only controls the caption shown under each, never the order. */
  private sortColumnsBy = COLSORT_FREQUENCY;
  /** The virtual-analog cell under the last right-click, so the context menu can offer a per-cell
   *  make-list add. Null when the right-click wasn't on an assembled virtual cell. */
  private contextCell: {matrix: SarMatrix, ri: number, ci: number} | null = null;
  /** Per-SMILES builders for this viewer's gated "SAR analysis" info panel, so a clicked analog shows
   *  its SAR context (prediction, support, decomposition, "Add to make-list") alongside the native
   *  Molecule panels. Registered on click; cleared on recompute and detach. Per-instance so one
   *  viewer closing can't wipe another's panels. */
  private readonly analogPanels = new Map<string, AnalogPanelBuilder>();
  /** Real cells of the currently-rendered matrix, keyed by parent-df molecule index, so host-grid
   *  selection / current-row changes can highlight the matching cells without a full re-render. */
  private cellByMolIdx = new Map<number, HTMLElement>();
  /** The cell currently carrying the "current" ring (real or virtual) — only one at a time. */
  private currentCellEl: HTMLElement | null = null;
  private readonly host = ui.divH([], 'chem-sar-matrix');
  private computing = false;
  /** Set when a recompute is requested while one is already running, so it is re-queued after. */
  private dirty = false;
  private computeTimer = 0;
  /** Set in `detach`, so an in-flight compute can't render into a closed viewer. */
  private detached = false;
  /** Cell width for the current render, fitted to the pane by `fitCellWidth`. */
  private cellW = CELL_W;
  /** Molblock template for the current matrix; cells and cores are aligned to it so the shared
   *  core points the same way everywhere. Set per row by the grid painters, and per side by the
   *  transfer pane. */
  private alignTemplate: string | null = null;
  /** The virtualized matrix grid currently on screen (null while a transfer pane is shown). Held so
   *  selection / current-row changes can repaint it with a cheap `invalidate()` — the highlight is
   *  drawn per-cell in `paintBodyCell`, so no DOM rebuild is needed. */
  private matrixGrid: MatrixGridState | null = null;
  /** Subscriptions owned by the current matrix grid (render/click/tooltip/overlay), unsubscribed and
   *  rebuilt whenever a new grid is created so repeated renders don't leak into detached grids. */
  private matrixGridSubs: {unsubscribe(): void}[] = [];
  /** Last pointer event over the matrix grid overlay, so a cell click can honor ctrl/shift for the
   *  host-grid selection extend (`onCellClick` carries no DOM event). */
  private lastGridMouseEvent: MouseEvent | null = null;

  constructor() {
    super();
    this.moleculesColumnName = this.string('moleculesColumnName', '', {userEditable: false});
    this.activityColumnName = this.string('activityColumnName', '', {userEditable: false});
    this.scaling = this.string('scaling', SCALING_METHODS.MINUS_LG, {choices: Object.values(SCALING_METHODS)});
    this.activityDirection = this.string('activityDirection', DIR_AUTO, {choices: ACTIVITY_DIRECTIONS});
    this.fragmentCutoff = this.float('fragmentCutoff', 0.4);
    this.threshold = this.float('threshold', 0.5);
    this.predictVirtual = this.bool('predictVirtual', true);
    // "SAR transfer" is a navigator section, not a ranking mode, so it is not offered as a rank scheme.
    this.rankScheme = this.string('rankScheme', SarRankScheme.Potency,
      {choices: [SarRankScheme.Potency, SarRankScheme.Discontinuity, SarRankScheme.Preferred]});
    this.host.style.height = '100%';
    this.root.appendChild(this.host);
  }

  onTableAttached(): void {
    this.detached = false; // the flag tracks the CURRENT attachment, not whether one ever ended
    // Cells are fitted to the pane width, so a resize needs a re-render (not a recompute).
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 200).subscribe(() => {
      if (!this.computing && this.matrices.length)
        this.render();
    }));
    this.subs.push(this.onContextMenu.subscribe((menu) => this.buildContextMenu(menu)));
    // Two-way link with the host grid: a selection or current-row change there highlights the
    // matching matrix cells (subscribed in onTableAttached, per the chem-search-base-viewer idiom).
    // Both streams feed one debounced subscription: as separate streams they cannot coalesce, so a
    // single click that moves the current row AND changes the selection repaints the matrix twice.
    this.subs.push(DG.debounce(
      merge(this.dataFrame.selection.onChanged, this.dataFrame.onCurrentRowChanged), 50)
      .subscribe(() => this.syncSelection()));
    // Capture-phase reset runs before a cell's own bubbling handler, so contextCell reflects only a
    // right-click that actually landed on a virtual cell (stale otherwise).
    this.host.addEventListener('contextmenu', () => this.contextCell = null, true);
    // Surface a clicked analog's SAR context at the top of its context panel. Scoped to this viewer
    // (this.subs, cleared on detach) so it can't inject panes platform-wide after the viewer closes.
    this.subs.push(grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
      const context = acc.context;
      const smiles = context instanceof DG.SemanticValue ? String(context.value) :
        (typeof context === 'string' ? context : null);
      if (!smiles || !this.analogPanels.has(smiles) || acc.getPane('SAR analysis'))
        return;
      const build = this.analogPanels.get(smiles)!;
      acc.addPane('SAR analysis', () => build(), true, acc.panes.length ? acc.panes[0] : null);
    }));
    // Start the RDKit workers now so their spawn + WASM instantiation overlaps the compute debounce
    // instead of being paid serially once the first fragmentation call arrives. A failure here is not
    // actionable — compute() re-enters the same init and reports it properly.
    getRdKitService().catch(() => {});
    this.scheduleCompute();
  }

  /** Base `detach` unsubscribes `this.subs`; also stop a pending compute and drop the analog-panel
   *  builders so a closed viewer leaves no timer firing or stale entries injecting panes. */
  detach(): void {
    this.detached = true;
    window.clearTimeout(this.computeTimer);
    this.analogPanels.clear();
    this.releaseMatrixGrid();
    super.detach();
  }

  /** Let go of the matrix grid currently on screen. Dropping the DOM is not enough: the
   *  render/click/tooltip subscriptions keep the grid — and through it its scaffold DataFrame and the
   *  whole SarMatrix — reachable, and a Dart-backed grid is only released when it is closed. */
  private releaseMatrixGrid(): void {
    this.matrixGridSubs.forEach((s) => s.unsubscribe());
    this.matrixGridSubs = [];
    try {
      this.matrixGrid?.grid?.close?.();
    } catch (e) {
      // A standalone (view-less) grid may not support close; dropping the reference is enough.
    }
    this.matrixGrid = null;
  }

  /** Move the single "current" ring to `td` (real or virtual), clearing it from the previous cell. */
  private setCurrentCell(td: HTMLElement | null): void {
    if (this.currentCellEl === td)
      return;
    this.currentCellEl?.classList.remove('chem-sar-cell-current');
    this.currentCellEl = td;
    td?.classList.add('chem-sar-cell-current');
  }

  /** Reflect the host grid's selection and current row onto the rendered matrix cells. Selection is a
   *  set (any number of real cells); the current ring is a single cell. A virtual cell has no grid
   *  row, so a virtual-cell click parks the grid's current row at -1 and owns its ring directly —
   *  here we only re-point the ring when the grid actually has a current row. */
  private syncSelection(): void {
    if (!this.dataFrame)
      return;
    // The virtualized matrix grid draws the selection/current ring per cell in `paintBodyCell`, so a
    // host-grid change only needs a repaint of the visible cells.
    this.matrixGrid?.grid.invalidate();
    const selection = this.dataFrame.selection;
    this.cellByMolIdx.forEach((td, molIdx) => td.classList.toggle('chem-sar-cell-selected', selection.get(molIdx)));
    const current = this.dataFrame.currentRowIdx;
    if (current >= 0)
      this.setCurrentCell(this.cellByMolIdx.get(current) ?? null);
  }

  // ---- Export virtual analogs to a make-list --------------------------------------------------

  /** Right-click menu: export this matrix's (or every matrix's) predicted analogs as a new table. */
  private buildContextMenu(menu: DG.Menu): void {
    if (!this.matrices.length)
      return;
    const group = menu.group('SAR Matrix');
    if (this.contextCell) {
      const {matrix, ri, ci} = this.contextCell;
      group.item('Add this analog to make-list', () => this.addAnalogToMakeList(matrix, ri, ci));
    }
    const current = this.selKind === 'transfer' ?
      this.matrices[this.transfers[this.selIndex]?.a.matrixIndex ?? 0] :
      this.matrices[Math.min(this.selIndex, this.matrices.length - 1)];
    if (current?.virtualCount)
      group.item(`Export virtual analogs (${current.label})`, () => this.exportVirtualAnalogs([current]));
    if (this.matrices.reduce((n, m) => n + m.virtualCount, 0) > (current?.virtualCount ?? 0)) {
      group.item(`Export virtual analogs (all ${this.matrices.length} matrices)`,
        () => this.exportVirtualAnalogs(this.matrices));
    }
  }

  /** Every assembled virtual cell in the given matrices, as (matrix, row, column) references. */
  private analogCells(matrices: SarMatrix[]): {matrix: SarMatrix, ri: number, ci: number}[] {
    const out: {matrix: SarMatrix, ri: number, ci: number}[] = [];
    for (const matrix of matrices) {
      for (let ri = 0; ri < matrix.rows.length; ri++) {
        for (let ci = 0; ci < matrix.columns.length; ci++) {
          const cell = matrix.cells[ri][ci];
          if (cell.kind === 'virtual' && cell.smiles !== null && cell.value !== null)
            out.push({matrix, ri, ci});
        }
      }
    }
    return out;
  }

  /**
   * Build a make-list molecule table from a set of virtual cells. Each row carries its predicted
   * activity, the support behind that prediction, and full provenance (series, core, R-position,
   * substituent, method) so it stands on its own once detached from the viewer.
   */
  private buildAnalogTable(cells: {matrix: SarMatrix, ri: number, ci: number}[], name: string): DG.DataFrame {
    const molCol = (name: string, values: string[]): DG.Column => {
      const col = DG.Column.fromStrings(name, values);
      col.semType = DG.SEMTYPE.MOLECULE;
      return col;
    };
    const cell = (c: {matrix: SarMatrix, ri: number, ci: number}): SarMatrixCell => c.matrix.cells[c.ri][c.ci];
    const df = DG.DataFrame.fromColumns([
      molCol('Structure', cells.map((c) => cell(c).smiles!)),
      DG.Column.fromList('double', `Predicted activity (${this.scaling})`, cells.map((c) => cell(c).value!)),
      DG.Column.fromList('int', 'Support (n)', cells.map((c) => cell(c).support ?? 0)),
      DG.Column.fromStrings('Series', cells.map((c) => c.matrix.label)),
      molCol('Core', cells.map((c) => c.matrix.rows[c.ri].coreSmiles)),
      DG.Column.fromStrings('Position', cells.map((c) => c.matrix.columns[c.ci].position)),
      molCol('Substituent', cells.map((c) => c.matrix.columns[c.ci].substSmiles)),
      DG.Column.fromStrings('Method', cells.map(() => FREE_WILSON_METHOD)),
    ]);
    df.name = name;
    return df;
  }

  /** A table name not currently open, so a fresh export never collides with an existing table. */
  private unusedTableName(base: string): string {
    if (!grok.shell.tableByName(base))
      return base;
    for (let i = 2; ; i++) {
      const name = `${base} (${i})`;
      if (!grok.shell.tableByName(name))
        return name;
    }
  }

  /** Bulk export: every assembled virtual analog in the given matrices to a new (uniquely-named) table. */
  private exportVirtualAnalogs(matrices: SarMatrix[]): void {
    const cells = this.analogCells(matrices);
    if (!cells.length) {
      grok.shell.info('No assembled virtual analogs to export.');
      return;
    }
    grok.shell.addTableView(this.buildAnalogTable(cells, this.unusedTableName(MAKELIST_NAME)));
    grok.shell.info(`Exported ${cells.length} virtual analog${cells.length === 1 ? '' : 's'}.`);
  }

  /** Per-cell design action: append one virtual analog to the running make-list, creating it (and
   *  its view) on first use, or reusing it if it is still open. */
  private addAnalogToMakeList(matrix: SarMatrix, ri: number, ci: number): void {
    const cell = matrix.cells[ri][ci];
    if (cell.kind !== 'virtual' || cell.smiles === null || cell.value === null) {
      grok.shell.info('This analog has no assembled structure to add.');
      return;
    }
    // Build a one-row frame so the column schema lives only in buildAnalogTable; append its row to
    // the running make-list, or open it as the make-list on first use.
    const analog = this.buildAnalogTable([{matrix, ri, ci}], MAKELIST_NAME);
    const existing = grok.shell.tableByName(MAKELIST_NAME);
    if (existing) {
      existing.rows.addNew(analog.columns.toList().map((c) => c.get(0)));
      grok.shell.info(`Added to make-list (${existing.rowCount} total).`);
    } else {
      grok.shell.addTableView(analog);
      grok.shell.info('Started a make-list with this analog.');
    }
  }

  onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);
    // Ranking and potency direction don't change the fragmentation, so they must not trigger a full
    // rebuild — re-fragmenting to reorder cards or flip the color direction costs seconds of RDKit
    // worker time for a result the already-assembled matrices can produce directly.
    if (property !== null && RERANK_ONLY_PROPS.includes(property.name)) {
      this.reRank();
      return;
    }
    this.scheduleCompute();
  }

  /** Re-rank the assembled matrices and redraw, without re-fragmenting. Used by the properties that
   *  only affect ordering/direction, and by the navigator's "Rank by" control. */
  private reRank(): void {
    if (!this.matrices.length) {
      this.scheduleCompute(); // nothing assembled yet — the first build still has to run
      return;
    }
    this.matrices = rankMatrices(this.matrices, this.rankScheme as SarRankScheme, this.higherIsBetter);
    this.computeTransfers();
    this.selKind = 'matrix';
    this.selIndex = 0;
    this.render();
  }

  /** Coalesce rapid property changes (e.g. from setOptions) into a single compute. */
  private scheduleCompute(): void {
    window.clearTimeout(this.computeTimer);
    this.computeTimer = window.setTimeout(() => this.compute(), 150);
  }

  private async compute(): Promise<void> {
    if (!this.dataFrame || !this.moleculesColumnName || !this.activityColumnName)
      return;
    const molecules = this.dataFrame.col(this.moleculesColumnName);
    const activity = this.dataFrame.col(this.activityColumnName);
    if (!molecules || !activity)
      return;
    // A change arrived mid-compute: mark dirty and let the running compute re-queue when it finishes,
    // so the final view always reflects the latest property values.
    if (this.computing) {
      this.dirty = true;
      return;
    }

    this.computing = true;
    this.analogPanels.clear(); // drop stale analog-panel closures from the previous matrices
    clearDepictionCaches(); // the previous matrices' depictions can never be hit again
    // Emptying the host detaches the grid's DOM but does not release the Dart-backed grid or its
    // subscriptions; without this the failure path below would leave one live and repainting forever.
    this.releaseMatrixGrid();
    ui.empty(this.host);
    this.host.appendChild(ui.loader());
    const progress = DG.TaskBarProgressIndicator.create('Building SAR matrices...');
    try {
      const params: SarMatrixParams = {
        scaling: this.scaling as SCALING_METHODS,
        fragmentCutoff: this.fragmentCutoff,
        predictVirtual: this.predictVirtual,
        threshold: this.threshold,
        rankScheme: this.rankScheme as SarRankScheme,
      };
      const matrices = await runSarMatrix(molecules, activity as DG.Column<number>, params);
      // The viewer may have been closed while the workers were running; rendering into it now would
      // build a grid nothing will ever release.
      if (this.detached)
        return;
      this.matrices = matrices;
      this.computeTransfers();
      this.selKind = 'matrix';
      this.selIndex = 0;
      this.render();
    } catch (e) {
      if (this.detached)
        return; // a viewer the user already closed must not report its own teardown as a failure
      const message = e instanceof Error ? e.message : String(e);
      ui.empty(this.host);
      this.host.appendChild(ui.divText(`SAR Matrix failed: ${message}`));
      grok.shell.error(`SAR Matrix: ${message}`);
    } finally {
      progress.close();
      this.computing = false;
      if (this.dirty) {
        this.dirty = false;
        if (!this.detached)
          this.scheduleCompute();
      }
    }
  }

  /** Whether a higher scaled activity means a more potent compound — drives coloring, ranking, the
   *  best-trend, and the transfer benefiting cell. The explicit Activity direction wins; on Auto only
   *  `-lg` (which turns a lower-is-better IC50 into a higher-is-better value) is higher-is-better,
   *  while raw (`none`) and `lg` keep the assay's native "lower is more potent" direction. */
  private get higherIsBetter(): boolean {
    if (this.activityDirection === DIR_HIGHER)
      return true;
    if (this.activityDirection === DIR_LOWER)
      return false;
    return this.scaling === SCALING_METHODS.MINUS_LG;
  }

  /** Rebuild the SAR-transfer list: every correlated core pair across all matrices — within a matrix
   *  AND across matrices (differently-scaffolded series) — strongest first. Recomputed whenever the
   *  matrix set or order changes. */
  private computeTransfers(): void {
    this.transfers = computeAllTransfers(this.matrices);
  }

  /** Potency tint at an explicit alpha (0-255), green = more potent. `DG.Color.scaleColor`'s own alpha
   *  arg can't be used (ignored at the scale ends, transparent mid-range), so take its RGB and set
   *  alpha here. For a lower-is-better activity the value is mirrored so the smallest maps to green. */
  private tint(matrix: SarMatrix, value: number, alpha: number): number {
    // A single observed value (or all-equal values) gives a zero-width range; scaleColor would
    // divide by zero. A lone value is trivially the most potent, so paint it the green end.
    if (matrix.maxActivity === matrix.minActivity)
      return DG.Color.argb(alpha, DG.Color.r(DG.Color.green), DG.Color.g(DG.Color.green), DG.Color.b(DG.Color.green));
    const potency = this.higherIsBetter ? value : matrix.minActivity + matrix.maxActivity - value;
    const base = DG.Color.scaleColor(potency, matrix.minActivity, matrix.maxActivity, undefined, ACTIVITY_SCHEME);
    return DG.Color.argb(alpha, DG.Color.r(base), DG.Color.g(base), DG.Color.b(base));
  }

  // ---- Navigator (left pane) ----------------------------------------------------------------

  /** Compact leave-one-out fit quality (e.g. " · R² 0.87"), the headline trust signal for the virtual
   *  predictions; empty when there were too few observations to cross-validate. */
  private confidenceShort(matrix: SarMatrix): string {
    const conf = matrix.confidence;
    return conf ? ` · R² ${formatR2(conf.r2)}` : '';
  }

  /** Human-readable name of the scaling applied to the activity, for labelling the cells' axis. */
  private get scalingLabel(): string {
    if (this.scaling === SCALING_METHODS.MINUS_LG)
      return '−log₁₀';
    if (this.scaling === SCALING_METHODS.LG)
      return 'log₁₀';
    return 'raw';
  }

  /** Mean activity (on the selected scale) over the matrix's observed cells, or null when none. */
  private meanRealActivity(matrix: SarMatrix): number | null {
    let sum = 0;
    let n = 0;
    for (const row of matrix.cells) {
      for (const cell of row) {
        if (cell.kind === 'real' && cell.value !== null) {
          sum += cell.value;
          n++;
        }
      }
    }
    return n ? sum / n : null;
  }

  /** The navigator card's rank score, made legible: best AND mean potency shown on the selected
   *  scale (matching the cells), plus a hover explanation naming the activity and scaling. */
  private cardScore(matrix: SarMatrix): {lines: {value: string, label: string}[], tip: string} {
    const scheme = this.rankScheme;
    const unit = this.activityColumnName || 'activity';
    if (scheme === SarRankScheme.Potency) {
      const best = this.higherIsBetter ? matrix.maxActivity : matrix.minActivity;
      const lines = [{value: this.formatActivity(best), label: 'best'}];
      const mean = this.meanRealActivity(matrix);
      if (mean !== null)
        lines.push({value: this.formatActivity(mean), label: 'mean'});
      return {lines, tip: `Best and mean potency of the series (${unit}), on the ${this.scalingLabel} ` +
        'scale — the same representation the cells show.'};
    }
    if (scheme === SarRankScheme.Preferred) {
      const raw = matrix.scores[scheme] ?? 0;
      return {lines: [{value: this.formatActivity(this.higherIsBetter ? raw : -raw), label: 'best R'}],
        tip: `Best mean potency of any single substituent across the cores (${unit}, ${this.scalingLabel} scale).`};
    }
    // The only remaining scheme is Discontinuity (the rank choices are Potency/Preferred/Discontinuity).
    return {lines: [{value: (matrix.scores[SarRankScheme.Discontinuity] ?? 0).toFixed(2), label: 'spread'}],
      tip: `Largest activity spread within a single core, on the ${this.scalingLabel} scale — ` +
        'how discontinuous the SAR is.'};
  }

  /** One selectable matrix card: the aligned core drawn (so the matrix is identified by its
   *  scaffold, not just "Series A"), a descriptor line, and the rank score. */
  private buildCard(matrix: SarMatrix, index: number): HTMLElement {
    const desc = `${matrix.rows.length} cores · ${matrix.positions.join('/')} · ${matrix.realCount} cpd` +
      (matrix.virtualCount ? ` · ${matrix.virtualCount} virtual` : '') + this.confidenceShort(matrix);
    const core = renderMoleculeOnColor(matrix.rows[0]?.coreSmiles ?? '', 78, 44, CORE_BG_ARGB);
    core.classList.add('chem-sar-card-core');
    const body = ui.divV([
      ui.divText(matrix.label, 'chem-sar-card-title'),
      ui.divText(desc, 'chem-sar-card-desc'),
    ], 'chem-sar-card-body');
    const sc = this.cardScore(matrix);
    const scoreBox = ui.divV(sc.lines.map((ln) => ui.divH([
      ui.divText(ln.value, 'chem-sar-card-score'),
      ui.divText(ln.label, 'chem-sar-card-score-cap'),
    ], 'chem-sar-card-score-line')), 'chem-sar-card-score-box');
    ui.tooltip.bind(scoreBox, () => sc.tip);
    const card = ui.divH([core, body, scoreBox], 'chem-sar-card');
    if (this.selKind === 'matrix' && index === this.selIndex)
      card.classList.add('selected');
    card.onclick = () => {
      this.selKind = 'matrix';
      this.selIndex = index;
      this.varyPosition = ''; // the Vary filter is per-matrix; don't carry it to a different one
      this.render();
    };
    return card;
  }

  /** Core label: "Series B · Core 1" with the series prefix, or just "Core 1" without it. */
  private coreLabel(side: TransferSide, withSeries: boolean): string {
    const matrix = this.matrices[side.matrixIndex];
    const core = matrix.rows[side.rowIndex].label;
    return withSeries ? `${matrix.label} · ${core}` : core;
  }

  /** Name a transfer side, prefixing the series whenever there is more than one matrix. */
  private sideLabel(side: TransferSide): string {
    return this.coreLabel(side, this.matrices.length > 1);
  }

  /** A transfer's target ("→" side) label with its R-position. The series prefix is dropped when the
   *  target shares the source's matrix, so it reads "Core 2 · R1" not "Series A · Core 2 · R1". */
  private transferTargetLabel(transfer: Transfer): string {
    return `${this.coreLabel(transfer.b, transfer.b.matrixIndex !== transfer.a.matrixIndex)} · ${transfer.a.position}`;
  }

  /** Group transfers by source core (best-correlation order preserved), so a source reaching several
   *  targets collapses into one card whose pane dropdown switches between them. */
  private groupTransfersBySource(transfers: Transfer[]): Transfer[][] {
    const groups = new Map<string, Transfer[]>();
    for (const t of transfers) {
      const key = transferSourceKey(t);
      let group = groups.get(key);
      if (!group) {
        group = [];
        groups.set(key, group);
      }
      group.push(t);
    }
    return [...groups.values()];
  }

  /** A transfer navigator section (same-series or cross-series): a header, then one card per source
   *  core — with a target dropdown when that source transfers to more than one core. */
  private buildTransferSection(title: string, transfers: Transfer[]): HTMLElement[] {
    if (!transfers.length)
      return [];
    const cards = this.groupTransfersBySource(transfers).map((group) => this.buildTransferCard(group));
    return [ui.divText(title, 'chem-sar-nav-section'), ...cards];
  }

  /** Select a specific transfer (by its identity in this.transfers) and re-open the trend view. */
  private selectTransfer(transfer: Transfer): void {
    this.selKind = 'transfer';
    this.selIndex = this.transfers.indexOf(transfer);
    this.varyPosition = '';
    this.render();
  }

  /** Every transfer sharing this one's source core — the alternatives the pane's target dropdown
   *  switches between. Keys on the same sourceKey as the nav, so card "+N" and dropdown always agree. */
  private transferSiblings(transfer: Transfer): Transfer[] {
    const key = transferSourceKey(transfer);
    return this.transfers.filter((t) => transferSourceKey(t) === key);
  }

  /** One selectable card for a group of transfers sharing a source core: the source and its (shown)
   *  target named in text, a "+N" badge when the source reaches more cores (switched in the pane),
   *  and the shown transfer's correlation. */
  private buildTransferCard(group: Transfer[]): HTMLElement {
    const selected = this.selKind === 'transfer' ? this.transfers[this.selIndex] : null;
    // Show the selected transfer when it belongs to this group, otherwise the strongest (first).
    const shown = (selected && group.includes(selected)) ? selected : group[0];

    const desc: HTMLElement[] = [
      ui.divText('→', 'chem-sar-card-arrow'),
      ui.divText(this.transferTargetLabel(shown), 'chem-sar-card-target-text'),
    ];
    if (group.length > 1) {
      const more = ui.divText(`+${group.length - 1}`, 'chem-sar-nav-badge');
      ui.tooltip.bind(more, () =>
        `${group.length - 1} more target core${group.length - 1 === 1 ? '' : 's'} — switch in the transfer panel.`);
      desc.push(more);
    }
    const body = ui.divV([
      ui.divText(this.sideLabel(shown.a), 'chem-sar-card-title'),
      ui.divH(desc, 'chem-sar-card-desc'),
    ], 'chem-sar-card-body');
    const scoreBox = ui.divV([ui.divH([
      ui.divText(`r ${shown.correlation.toFixed(2)}`, 'chem-sar-card-score'),
      ui.divText('transfer', 'chem-sar-card-score-cap'),
    ], 'chem-sar-card-score-line')], 'chem-sar-card-score-box');
    ui.tooltip.bind(scoreBox, () => shown.crossMatrix ?
      'Correlation of two different chemotypes’ potency trends across shared substituents — the SAR ' +
      'learned on one series should carry to the other. 1.00 = perfectly parallel.' :
      'Correlation of the two cores’ potency trends (SAR transfer): 1.00 = perfectly parallel.');
    const card = ui.divH([body, scoreBox], 'chem-sar-card');
    if (selected && group.includes(selected))
      card.classList.add('selected');
    card.onclick = () => this.selectTransfer(shown);
    return card;
  }

  private buildNavigator(): HTMLElement {
    const rankValue = this.rankScheme;
    const rankInput = ui.input.choice('Rank by', {
      value: rankValue,
      items: [SarRankScheme.Potency, SarRankScheme.Discontinuity, SarRankScheme.Preferred],
      // A field assignment raises no onPropertyChanged — the property has no accessor on the viewer —
      // so re-rank explicitly here. The property-panel path reaches the same method via that event.
      onValueChanged: (value) => {
        this.rankScheme = value!;
        this.reRank();
      },
    });

    const title = ui.divH([
      ui.divText('Select a SAR matrix'),
      ui.divText('MMP grouping', 'chem-sar-nav-badge'),
    ], 'chem-sar-nav-title');
    const sub = ui.divText(
      `${this.matrices.length} matrices · related cores grouped by similarity`, 'chem-sar-nav-sub');
    // The one place the scale is spelled out: every value (cards and cells) and the cell coloring use it.
    const units = ui.divText(
      `Activity: ${this.activityColumnName} · values on the ${this.scalingLabel} scale`,
      'chem-sar-nav-units');
    const header = ui.divV([title, sub, units, rankInput.root], 'chem-sar-nav-header');

    const list = ui.div([], 'chem-sar-nav-list');
    list.appendChild(ui.divText(rankValue.toUpperCase(), 'chem-sar-nav-section'));
    this.matrices.forEach((matrix, i) => list.appendChild(this.buildCard(matrix, i)));
    const sameSeries = this.transfers.filter((t) => !t.crossMatrix);
    const crossSeries = this.transfers.filter((t) => t.crossMatrix);
    for (const el of this.buildTransferSection('SAME-SERIES SAR TRANSFER', sameSeries))
      list.appendChild(el);
    for (const el of this.buildTransferSection('CROSS-SERIES SAR TRANSFER', crossSeries))
      list.appendChild(el);

    return ui.divV([header, list], 'chem-sar-nav');
  }

  // ---- Matrix table (right pane) ------------------------------------------------------------

  /**
   * Grow cells to fill the pane so a narrow matrix doesn't leave half the width empty, but never
   * shrink below CELL_W — with many columns the table scrolls horizontally instead.
   */
  private fitCellWidth(nCols: number): number {
    if (nCols <= 0)
      return CELL_W;
    const avail = (this.root.clientWidth || 900) - NAV_W - CORE_W - TABLE_CHROME - nCols * 6;
    return Math.max(CELL_W, Math.min(CELL_W_MAX, Math.floor(avail / nCols)));
  }

  /** Column indices to show: all, or only the "Vary" position's group when one is chosen. Columns
   *  keep their as-assembled (frequency) order — the "Label" control annotates them, never reorders. */
  private visibleColIdxs(matrix: SarMatrix): number[] {
    const all = matrix.columns.map((_, ci) => ci);
    return this.varyPosition && matrix.positions.includes(this.varyPosition) ?
      all.filter((ci) => matrix.columns[ci].position === this.varyPosition) : all;
  }

  /** Mean observed potency of a column (real cells only), or null when it has no observation. */
  private columnMeanPotency(matrix: SarMatrix, colIdx: number): number | null {
    let sum = 0;
    let n = 0;
    for (let ri = 0; ri < matrix.rows.length; ri++) {
      const cell = matrix.cells[ri][colIdx];
      if (cell.kind === 'real' && cell.value !== null) {
        sum += cell.value;
        n++;
      }
    }
    return n ? sum / n : null;
  }

  /** Caption annotating a column with the chosen metric — its mean potency (scaled) or the
   *  substituent MW — shown under the substituent header. Empty when labelling is off. */
  private columnSortCaption(matrix: SarMatrix, colIdx: number): string {
    if (this.sortColumnsBy === COLSORT_MW) {
      const mw = substituentMW(matrix.columns[colIdx].substSmiles);
      return Number.isFinite(mw) ? `MW ${mw.toFixed(0)}` : '';
    }
    if (this.sortColumnsBy === COLSORT_POTENCY) {
      const mean = this.columnMeanPotency(matrix, colIdx);
      return mean === null ? '' : `μ ${this.formatActivity(mean)}`;
    }
    return '';
  }

  /** Contiguous runs of the given columns sharing an R-position, in order. */
  private groupColumns(matrix: SarMatrix, colIdxs: number[]): {position: string, colIdxs: number[]}[] {
    const groups: {position: string, colIdxs: number[]}[] = [];
    for (const ci of colIdxs) {
      const position = matrix.columns[ci].position;
      const last = groups[groups.length - 1];
      if (last?.position === position)
        last.colIdxs.push(ci);
      else
        groups.push({position, colIdxs: [ci]});
    }
    return groups;
  }

  private buildMatrixCell(matrix: SarMatrix, cell: SarMatrixCell, ri: number, ci: number): HTMLElement {
    const td = ui.element('td') as HTMLTableCellElement;
    td.className = `chem-sar-cell chem-sar-cell-${cell.kind}`;
    td.style.width = `${this.cellW}px`;
    if (cell.kind === 'empty' || cell.value === null)
      return td;
    const value = cell.value;
    const isVirtual = cell.kind === 'virtual';
    const support = cell.support ?? 0;
    // Thin predictions (few observations behind the Free-Wilson estimate) are drawn fainter, so a
    // heavily-extrapolated cell reads as less certain than a well-supported one. Real cells are solid.
    const alpha = isVirtual ?
      Math.round(VIRTUAL_ALPHA_MIN + (VIRTUAL_ALPHA - VIRTUAL_ALPHA_MIN) * Math.min(1, support / FULL_SUPPORT)) :
      REAL_ALPHA;
    // Tint the whole cell (not just the drawn canvas) so no white gap shows around the structure;
    // the molecule is then drawn on a TRANSPARENT canvas over it, avoiding a double-tint seam.
    const bg = this.tint(matrix, value, alpha);
    td.style.backgroundColor = DG.Color.toHtml(bg);
    if (cell.smiles !== null)
      td.appendChild(renderAlignedOnColor(cell.smiles, this.cellW, CELL_H, CORE_BG_ARGB, this.alignTemplate));
    // '~' marks a predicted (virtual) value.
    const chip = ui.divText(`${isVirtual ? '~' : ''}${this.formatActivity(value)}`, 'chem-sar-chip');
    if (isVirtual && support <= 1)
      chip.classList.add('chem-sar-chip-weak');
    td.appendChild(chip);
    // Compare an observed cell with the additive model's fitted value: a large residual means the
    // Free-Wilson assumption fails at this core × substituent (a non-additive, cliff-like cell).
    const range = matrix.maxActivity - matrix.minActivity;
    const fitMatches = (fit: number): boolean => range === 0 || Math.abs(value - fit) <= MATCH_FRACTION * range;
    if (!isVirtual && cell.fit !== undefined && !fitMatches(cell.fit))
      td.classList.add('chem-sar-cell-deviant');
    ui.tooltip.bind(td, () => {
      if (isVirtual) {
        return `Predicted ${this.formatActivity(value)} · local Free-Wilson · support n=${support}` +
          (support <= 1 ? ' (low)' : '');
      }
      if (cell.fit === undefined)
        return `Observed ${this.formatActivity(value)}`;
      return `Observed ${this.formatActivity(value)} · Free-Wilson fit ${this.formatActivity(cell.fit)} · ` +
        (fitMatches(cell.fit) ? '✓ matches' : 'non-additive');
    });
    td.onclick = (event) => {
      grok.shell.windows.showContextPanel = true;
      // The "current" ring is a single cell — move it to whatever was clicked (real or virtual), so
      // clicking a virtual analog clears the previously-clicked cell's ring.
      this.setCurrentCell(td);
      // Link back to the host grid: a real compound becomes the current row (grid scrolls to it) and
      // ctrl/shift-click extends the grid selection (modifiedSelectOnly keeps a plain click from
      // wiping an existing selection). A virtual analog has no grid row, so park the current row at
      // -1 — otherwise the grid (and syncSelection) would keep the previous real cell ringed.
      if (cell.molIdx !== null) {
        this.dataFrame.currentRowIdx = cell.molIdx;
        this.dataFrame.selection.handleClick((i) => i === cell.molIdx, event, true);
      } else
        this.dataFrame.currentRowIdx = -1;
      // Any cell with an assembled structure — real OR virtual — opens the platform's full Molecule
      // context (properties, drug-likeness, toxicity, structural alerts, identifiers, 3D). We also
      // register its SAR context (potency, decomposition, and — for a prediction — the make-list
      // action) so the gated `SAR analysis` info panel shows it alongside the native panels. Only an
      // unassembled virtual (no structure) falls back to the standalone SAR panel.
      if (cell.smiles) {
        this.analogPanels.set(cell.smiles, () => this.buildCellPanel(matrix, ri, ci));
        grok.shell.o = DG.SemanticValue.fromValueType(cell.smiles, DG.SEMTYPE.MOLECULE);
      } else
        grok.shell.o = this.buildCellPanel(matrix, ri, ci);
    };
    // Track the right-clicked virtual analog so the context menu can offer a per-cell make-list add.
    if (isVirtual && cell.smiles)
      td.addEventListener('contextmenu', () => this.contextCell = {matrix, ri, ci});
    if (cell.molIdx !== null)
      this.cellByMolIdx.set(cell.molIdx, td); // for host-grid selection/current-row highlighting
    return td;
  }

  /**
   * Virtualized matrix render: a scaffold DataFrame (one row per core, one string column per visible
   * substituent plus a pinned 'Core' column) backs a DG.Grid whose every cell — body, core, and
   * R-group header — is hand-painted in `onCellRender`. Only viewport cells draw, so a large matrix no
   * longer renders every core×substituent up front, and selection/current-row changes repaint via
   * `invalidate` instead of rebuilding the DOM.
   */
  private buildMatrixGrid(matrix: SarMatrix): HTMLElement {
    const colIdxs = this.visibleColIdxs(matrix);
    this.cellW = this.fitCellWidth(colIdxs.length);
    // Header state is computed here, once per grid: `onCellRender` runs for every visible header cell
    // on every repaint, and rebuilding the column list (or rescanning rows for a column's mean
    // potency) there is quadratic in the column count for a result that cannot change until rebuild.
    const firstOfGroup = new Set(this.groupColumns(matrix, colIdxs).map((group) => group.colIdxs[0]));
    const captions = new Map<number, string>();
    for (const ci of colIdxs)
      captions.set(ci, this.columnSortCaption(matrix, ci));

    const df = DG.DataFrame.create(matrix.rows.length);
    df.columns.addNewString('Core');
    const colKeyToCi = new Map<string, number>();
    // Stable string keys (never the grid column idx, which pinning and the hidden row header shift).
    for (const ci of colIdxs) {
      const key = `c${ci}`;
      df.columns.addNewString(key);
      colKeyToCi.set(key, ci);
    }

    const grid = DG.Grid.create(df);
    // A fixed header height fits the position band + R-group depiction + sort caption; the built-in
    // row-number column is hidden because the cores live in the pinned 'Core' column.
    grid.setOptions({colHeaderHeight: HEADER_H + 36, rowHeight: CELL_H, showRowHeader: false});
    grid.col('Core')!.width = CORE_W;
    colKeyToCi.forEach((_ci, key) => grid.col(key)!.width = this.cellW);
    grid.col('Core')!.pin();
    const state: MatrixGridState = {grid, colKeyToCi, matrix, colIdxs, firstOfGroup, captions};

    // Owned by this grid instance; unsubscribed when the next grid replaces it (or on detach) so
    // repeated renders can't leak render/click/tooltip handlers into detached grids.
    this.matrixGridSubs.forEach((s) => s.unsubscribe());
    this.matrixGridSubs = [];

    this.matrixGridSubs.push(grid.onCellRender.subscribe((args) => {
      const c = args.cell;
      const isColHeader = c.isColHeader;
      if (!isColHeader && !c.isTableCell)
        return;
      const b = args.bounds;
      if (!b || b.width < 1 || b.height < 1)
        return;
      const g = args.g;
      // The grid's render context already carries a devicePixelRatio scale, so paint in CSS
      // coordinates; save/restore isolates the fillStyle/font/lineDash changes from other cells.
      g.save();
      try {
        const name = c.gridColumn.name;
        if (isColHeader)
          this.paintHeader(g, b, name, state);
        else {
          const ri = grid.gridRowToTable(c.gridRow);
          if (ri >= 0 && ri < matrix.rows.length) {
            if (name === 'Core')
              this.paintCore(g, b, matrix, ri);
            else {
              const ci = colKeyToCi.get(name);
              if (ci !== undefined)
                this.paintBodyCell(g, b, matrix, ri, ci);
            }
          }
        }
      } finally {
        g.restore();
      }
      args.preventDefault();
    }));

    this.matrixGridSubs.push(grid.onCellClick.subscribe((c) => this.onGridCellClick(grid, matrix, colKeyToCi, c)));
    this.matrixGridSubs.push(grid.onCellTooltip((c, x, y) =>
      this.onGridCellTooltip(grid, matrix, colKeyToCi, c, x, y)));

    // Track the pointer for ctrl/shift selection extend, and decode right-clicks to the assembled
    // virtual cell so the viewer's context menu can offer a per-cell make-list add.
    const overlay = grid.overlay;
    const onMouseDown = (e: MouseEvent): void => {this.lastGridMouseEvent = e;};
    const onContextMenu = (e: MouseEvent): void => {
      this.contextCell = null;
      const rect = overlay.getBoundingClientRect();
      const hit = grid.hitTest(e.clientX - rect.left, e.clientY - rect.top);
      if (!hit || !hit.isTableCell)
        return;
      const name = hit.gridColumn.name;
      if (name === 'Core')
        return;
      const ci = colKeyToCi.get(name);
      if (ci === undefined)
        return;
      const ri = grid.gridRowToTable(hit.gridRow);
      if (ri < 0 || ri >= matrix.rows.length)
        return;
      const cell = matrix.cells[ri][ci];
      if (cell.kind === 'virtual' && cell.smiles)
        this.contextCell = {matrix, ri, ci};
    };
    overlay.addEventListener('mousedown', onMouseDown);
    overlay.addEventListener('contextmenu', onContextMenu);
    this.matrixGridSubs.push({unsubscribe: () => overlay.removeEventListener('mousedown', onMouseDown)});
    this.matrixGridSubs.push({unsubscribe: () => overlay.removeEventListener('contextmenu', onContextMenu)});

    this.matrixGrid = state;

    // The grid virtualizes off its own height, so give it a plain flex host that fills the pane (not
    // ui.box, which would pin a fixed pixel width and stop the matrix growing with the pane).
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    return ui.div([grid.root], 'chem-sar-grid-host');
  }

  /** Click on a grid body cell: open the platform Molecule context, set the host grid's current row /
   *  selection, and register the cell's SAR panel — mirrors `buildMatrixCell`'s click handler. */
  private onGridCellClick(grid: DG.Grid, matrix: SarMatrix, colKeyToCi: Map<string, number>,
    c: DG.GridCell): void {
    if (!c.isTableCell)
      return;
    const name = c.gridColumn.name;
    if (name === 'Core')
      return;
    const ci = colKeyToCi.get(name);
    if (ci === undefined)
      return;
    const ri = grid.gridRowToTable(c.gridRow);
    if (ri < 0 || ri >= matrix.rows.length)
      return;
    const cell = matrix.cells[ri][ci];
    if (cell.kind === 'empty' || cell.value === null)
      return;
    grok.shell.windows.showContextPanel = true;
    // A real compound becomes the host grid's current row and (ctrl/shift) extends its selection; a
    // virtual analog has no row, so park the current row at -1 to drop any previous real-cell ring.
    if (cell.molIdx !== null) {
      this.dataFrame.currentRowIdx = cell.molIdx;
      const event = this.lastGridMouseEvent ?? new MouseEvent('click');
      this.dataFrame.selection.handleClick((i) => i === cell.molIdx, event, true);
    } else
      this.dataFrame.currentRowIdx = -1;
    // Any assembled cell opens the platform Molecule context; its SAR context is registered for the
    // gated "SAR analysis" pane. An unassembled virtual falls back to the standalone SAR panel.
    if (cell.smiles) {
      this.analogPanels.set(cell.smiles, () => this.buildCellPanel(matrix, ri, ci));
      grok.shell.o = DG.SemanticValue.fromValueType(cell.smiles, DG.SEMTYPE.MOLECULE);
    } else
      grok.shell.o = this.buildCellPanel(matrix, ri, ci);
    // Repaint so the freshly-set current/selection ring shows immediately.
    grid.invalidate();
  }

  /** Hover text for a grid cell: substituent/position for an R-group header, potency detail for a body
   *  cell. Returns true to suppress the default tooltip. */
  private onGridCellTooltip(grid: DG.Grid, matrix: SarMatrix, colKeyToCi: Map<string, number>,
    c: DG.GridCell, x: number, y: number): boolean {
    const name = c.gridColumn.name;
    if (c.isColHeader) {
      if (name === 'Core')
        return true;
      const ci = colKeyToCi.get(name);
      if (ci === undefined)
        return false;
      const col = matrix.columns[ci];
      const otherRefs = matrix.positions.filter((p) => p !== col.position);
      const parts = [`${col.position}: ${col.substSmiles}`];
      if (otherRefs.length)
        parts.push(`${otherRefs.join(', ')} at ref`);
      ui.tooltip.show(ui.divText(parts.join(' · ')), x, y);
      return true;
    }
    if (!c.isTableCell || name === 'Core')
      return false;
    const ci = colKeyToCi.get(name);
    if (ci === undefined)
      return false;
    const ri = grid.gridRowToTable(c.gridRow);
    if (ri < 0 || ri >= matrix.rows.length)
      return false;
    const cell = matrix.cells[ri][ci];
    if (cell.kind === 'empty' || cell.value === null)
      return false;
    ui.tooltip.show(ui.divText(this.cellTooltipText(matrix, cell)), x, y);
    return true;
  }

  /** Cell hover text — mirrors the HTML matrix: predicted (+ support) for a virtual analog, observed
   *  (+ additive-fit check) for a real one. */
  private cellTooltipText(matrix: SarMatrix, cell: SarMatrixCell): string {
    const value = cell.value!;
    const support = cell.support ?? 0;
    if (cell.kind === 'virtual') {
      return `Predicted ${this.formatActivity(value)} · local Free-Wilson · support n=${support}` +
        (support <= 1 ? ' (low)' : '');
    }
    if (cell.fit === undefined)
      return `Observed ${this.formatActivity(value)}`;
    const range = matrix.maxActivity - matrix.minActivity;
    const matches = range === 0 || Math.abs(value - cell.fit) <= MATCH_FRACTION * range;
    return `Observed ${this.formatActivity(value)} · Free-Wilson fit ${this.formatActivity(cell.fit)} · ` +
      (matches ? '✓ matches' : 'non-additive');
  }

  /** Paint an R-group column header (position band + straightened substituent depiction + sort
   *  caption) or the 'Core' column's "Aligned core" label. The header band is pinned on top natively. */
  private paintHeader(g: CanvasRenderingContext2D, b: DG.Rect, name: string, state: MatrixGridState): void {
    const grey6 = cssColor(this.root, '--grey-6', '#4a4a4a');
    const grey5 = cssColor(this.root, '--grey-5', '#7d7d7d');
    if (name === 'Core') {
      g.fillStyle = grey5;
      g.font = `italic 11px ${GRID_FONT}`;
      g.textAlign = 'center';
      g.textBaseline = 'middle';
      g.fillText('Aligned core', b.x + b.width / 2, b.y + b.height / 2, b.width - 6);
      return;
    }
    const ci = state.colKeyToCi.get(name);
    if (ci === undefined)
      return;
    const matrix = state.matrix;
    g.fillStyle = DG.Color.toHtml(HEADER_ARGB);
    g.fillRect(b.x, b.y, b.width, b.height);

    // Draw the position label only on the first column of a group — per-cell clipping can't paint a
    // true cross-column spanner, so this reads as the group header (the "at ref" detail is on hover).
    const posBandH = 16;
    if (state.firstOfGroup.has(ci)) {
      g.fillStyle = grey6;
      g.font = `600 11px ${GRID_FONT}`;
      g.textAlign = 'left';
      g.textBaseline = 'top';
      g.fillText(matrix.columns[ci].position, b.x + 5, b.y + 2, b.width - 10);
    }

    const caption = state.captions.get(ci) ?? '';
    const captionBandH = caption ? 14 : 0;
    const depH = Math.max(1, b.height - posBandH - captionBandH);
    const depW = Math.min(b.width, HEADER_W * 2);
    const depCanvas = cachedCellCanvas(matrix.columns[ci].substSmiles, null, depW, depH, HEADER_ARGB);
    g.drawImage(depCanvas, b.x + (b.width - depW) / 2, b.y + posBandH, depW, depH);

    if (caption) {
      g.fillStyle = grey5;
      g.font = `10px ${GRID_FONT}`;
      g.textAlign = 'center';
      g.textBaseline = 'top';
      g.fillText(caption, b.x + b.width / 2, b.y + posBandH + depH, b.width);
    }
  }

  /** Paint the pinned-left core cell: the core aligned to its own template + the row label beneath. */
  private paintCore(g: CanvasRenderingContext2D, b: DG.Rect, matrix: SarMatrix, ri: number): void {
    const row = matrix.rows[ri];
    // Each core aligns to its OWN template so every cell in the row shows that core the same way.
    this.alignTemplate = buildAlignmentTemplate(row.coreSmiles);
    const labelH = 16;
    const molH = Math.max(1, b.height - labelH);
    const canvas = cachedCellCanvas(row.coreSmiles, this.alignTemplate, b.width, molH, CORE_BG_ARGB);
    g.drawImage(canvas, b.x, b.y, b.width, molH);
    g.fillStyle = cssColor(this.root, '--grey-6', '#4a4a4a');
    g.font = `600 11px ${GRID_FONT}`;
    g.textAlign = 'center';
    g.textBaseline = 'top';
    g.fillText(row.label, b.x + b.width / 2, b.y + molH + 1, b.width);
  }

  /** Paint one core×substituent cell: potency tint (over white), aligned assembled molecule, value
   *  chip, virtual/deviant markers, and the host-grid selection/current ring — the grid analogue of
   *  `buildMatrixCell`. */
  private paintBodyCell(g: CanvasRenderingContext2D, b: DG.Rect, matrix: SarMatrix, ri: number,
    ci: number): void {
    const cell = matrix.cells[ri][ci];
    if (cell.kind === 'empty' || cell.value === null) {
      g.fillStyle = cssColor(this.root, '--white', '#ffffff');
      g.fillRect(b.x, b.y, b.width, b.height);
      return;
    }
    const value = cell.value;
    const isVirtual = cell.kind === 'virtual';
    const support = cell.support ?? 0;
    // Thin predictions read fainter than well-supported ones; real cells are solid.
    const alpha = isVirtual ?
      Math.round(VIRTUAL_ALPHA_MIN + (VIRTUAL_ALPHA - VIRTUAL_ALPHA_MIN) * Math.min(1, support / FULL_SUPPORT)) :
      REAL_ALPHA;
    // White base + the semi-transparent tint reproduces the HTML cell's "tint over the pane" without
    // depending on how the grid clears cells; the molecule then draws on a transparent canvas over it.
    g.fillStyle = cssColor(this.root, '--white', '#ffffff');
    g.fillRect(b.x, b.y, b.width, b.height);
    g.fillStyle = DG.Color.toHtml(this.tint(matrix, value, alpha));
    g.fillRect(b.x, b.y, b.width, b.height);
    if (cell.smiles !== null) {
      this.alignTemplate = buildAlignmentTemplate(matrix.rows[ri].coreSmiles);
      const canvas = cachedCellCanvas(cell.smiles, this.alignTemplate, b.width, b.height, CORE_BG_ARGB);
      g.drawImage(canvas, b.x, b.y, b.width, b.height);
    }

    // '~value' chip, top-left, fainter for a single-observation prediction.
    this.paintChip(g, b, `${isVirtual ? '~' : ''}${this.formatActivity(value)}`, isVirtual,
      isVirtual && support <= 1);

    // A real cell whose observed value departs from the additive fit gets an amber corner dot.
    const range = matrix.maxActivity - matrix.minActivity;
    const fitMatches = (fit: number): boolean => range === 0 || Math.abs(value - fit) <= MATCH_FRACTION * range;
    if (!isVirtual && cell.fit !== undefined && !fitMatches(cell.fit)) {
      g.fillStyle = cssColor(this.root, '--orange-2', '#e8853a');
      g.beginPath();
      g.arc(b.x + b.width - 8, b.y + 8, 3.5, 0, Math.PI * 2);
      g.fill();
    }

    // Cell outline: dashed for a virtual analog, a light solid frame otherwise.
    g.lineWidth = 1;
    if (isVirtual) {
      g.setLineDash([3, 2]);
      g.strokeStyle = cssColor(this.root, '--grey-4', '#b5b5b5');
    } else {
      g.setLineDash([]);
      g.strokeStyle = cssColor(this.root, '--grey-2', '#e0e0e0');
    }
    g.strokeRect(b.x + 0.5, b.y + 0.5, b.width - 1, b.height - 1);
    g.setLineDash([]);

    // Host-grid link: a selected row's cell gets a blue ring; the current row a stronger one.
    if (cell.molIdx !== null && this.dataFrame.selection.get(cell.molIdx)) {
      g.lineWidth = 2;
      g.strokeStyle = cssColor(this.root, '--blue-1', '#2083d5');
      g.strokeRect(b.x + 1, b.y + 1, b.width - 2, b.height - 2);
    }
    if (cell.molIdx !== null && cell.molIdx === this.dataFrame.currentRowIdx) {
      g.lineWidth = 3;
      g.strokeStyle = cssColor(this.root, '--blue-3', '#0d5ba6');
      g.strokeRect(b.x + 1.5, b.y + 1.5, b.width - 3, b.height - 3);
    }
  }

  /** Draw the potency chip (a translucent white pill + text) at the cell's top-left. */
  private paintChip(g: CanvasRenderingContext2D, b: DG.Rect, text: string, isVirtual: boolean,
    weak: boolean): void {
    g.font = `${isVirtual ? 'italic ' : ''}600 10px ${GRID_FONT}`;
    g.textAlign = 'left';
    g.textBaseline = 'top';
    const padX = 3;
    const chipW = g.measureText(text).width + padX * 2;
    g.globalAlpha = weak ? 0.65 : 1;
    g.fillStyle = 'rgba(255, 255, 255, 0.82)';
    g.fillRect(b.x + 3, b.y + 3, chipW, 13);
    g.fillStyle = cssColor(this.root, isVirtual ? '--grey-5' : '--grey-6', isVirtual ? '#7d7d7d' : '#4a4a4a');
    g.fillText(text, b.x + 3 + padX, b.y + 5, b.width - 6);
    g.globalAlpha = 1;
  }

  /** Format an activity value: no decimals at ≥100, one decimal below. */
  private formatActivity(value: number): string {
    return value >= 100 ? value.toFixed(0) : value.toFixed(1);
  }

  /** A context-panel section header with a small provenance badge (mirrors the mockup). */
  private cpSection(title: string, badge: string): HTMLElement {
    return ui.divH([
      ui.divText(title, 'chem-sar-cp-section-title'),
      ui.divText(badge, 'chem-sar-cp-badge'),
    ], 'chem-sar-cp-section');
  }

  /** A label / value row in the context panel. */
  private cpRow(label: string, value: HTMLElement | string): HTMLElement {
    const v = typeof value === 'string' ? ui.divText(value, 'chem-sar-cp-rv') : value;
    return ui.divH([ui.divText(label, 'chem-sar-cp-rl'), v], 'chem-sar-cp-row');
  }

  /** A small framed fragment (core or substituent) with a caption. */
  private cpFragment(smiles: string, label: string): HTMLElement {
    return ui.divV([
      ui.div([renderMolecule(smiles, {width: 78, height: 52, popupMenu: false})], 'chem-sar-cp-frag-box'),
      ui.divText(label, 'chem-sar-cp-frag-label'),
    ], 'chem-sar-cp-frag');
  }

  /**
   * Rich SAR context for the selected cell, mirroring the mockup: header (+ "not synthesized" badge
   * for a prediction), framed structure, SMILES + core, the observed or predicted potency (with the
   * additive-fit check for real cells), the decomposition, and — for a virtual analog — a design
   * action to add it to the make-list.
   */
  private buildCellPanel(matrix: SarMatrix, rowIdx: number, colIdx: number): HTMLElement {
    const cell = matrix.cells[rowIdx][colIdx];
    const row = matrix.rows[rowIdx];
    const column = matrix.columns[colIdx];
    const panel = ui.divV([], 'chem-sar-matrix-cp');

    if (cell.kind === 'empty' || cell.value === null) {
      panel.appendChild(ui.divText('No compound and no prediction for this combination.'));
      return panel;
    }
    const isVirtual = cell.kind === 'virtual';

    const header = ui.divH([ui.h2(isVirtual ? 'Virtual analog' : 'Compound')], 'chem-sar-cp-head');
    if (isVirtual)
      header.appendChild(ui.divText('not synthesized', 'chem-sar-cp-notsynth'));
    panel.appendChild(header);

    if (cell.smiles) {
      panel.appendChild(ui.div([renderMolecule(cell.smiles, {width: 240, height: 130, popupMenu: false})],
        'chem-sar-cp-structbox'));
    }

    const displaySmiles = cell.smiles ? toDisplaySmiles(cell.smiles) : '';
    if (displaySmiles)
      panel.appendChild(this.cpRow('SMILES', ui.divText(displaySmiles, 'chem-sar-cp-smiles')));
    const substs = matrix.positions
      .map((p) => (p === column.position ? column.substSmiles : matrix.refValues[p])).filter((v) => v);
    panel.appendChild(this.cpRow(isVirtual ? 'Core × R' : 'Core',
      isVirtual ? `${row.label} × ${substs.join(' / ')}` : `${row.label} · ${matrix.rows.length} cores`));

    panel.appendChild(this.cpSection(isVirtual ? 'Predicted potency' : 'Potency',
      isVirtual ? 'Free-Wilson' : 'observed'));
    panel.appendChild(this.cpRow(isVirtual ? 'Predicted' : 'Observed',
      ui.divText(this.formatActivity(cell.value), 'chem-sar-cp-value')));
    if (isVirtual) {
      panel.appendChild(this.cpRow('Method', FREE_WILSON_METHOD));
      panel.appendChild(this.cpRow('Neighbours', 'row + column analogs'));
      panel.appendChild(this.cpRow('Support', `n = ${cell.support ?? 0}${(cell.support ?? 0) <= 1 ? ' (low)' : ''}`));
    } else if (cell.fit !== undefined) {
      const range = matrix.maxActivity - matrix.minActivity;
      const matches = range === 0 || Math.abs(cell.value - cell.fit) <= MATCH_FRACTION * range;
      panel.appendChild(this.cpRow('Free-Wilson fit', ui.divH([
        ui.divText(this.formatActivity(cell.fit)),
        ui.divText(matches ? '✓ matches' : 'non-additive',
          matches ? 'chem-sar-cp-matches' : 'chem-sar-cp-mismatch'),
      ], 'chem-sar-cp-fit')));
    }

    panel.appendChild(this.cpSection('Decomposition', 'R-group'));
    const parts = [this.cpFragment(row.coreSmiles, 'core')];
    matrix.positions.forEach((position) => {
      const v = position === column.position ? column.substSmiles : matrix.refValues[position];
      if (v)
        parts.push(this.cpFragment(v, `${position} = ${v}`));
    });
    panel.appendChild(ui.divH(parts, 'chem-sar-cp-decomp'));

    if (isVirtual && cell.smiles) {
      panel.appendChild(this.cpSection('Design action', 'new'));
      panel.appendChild(ui.divText(
        'This core × substituent is not in the dataset. Add it to the make-list for synthesis triage.',
        'chem-sar-cp-hint'));
      const btn = ui.button('Add to make-list', () => this.addAnalogToMakeList(matrix, rowIdx, colIdx));
      btn.classList.add('chem-sar-cp-generate');
      panel.appendChild(btn);
    }
    return panel;
  }

  /** One core's row in the transfer view (from its own matrix), followed by its trend row. */
  private appendTransferSide(tbody: HTMLElement, side: TransferSide, cols: number[],
    template?: string | null): void {
    const matrix = this.matrices[side.matrixIndex];
    const row = matrix.rows[side.rowIndex];
    // Align to this side's own core; the caller may pass a template it already built for this core.
    this.alignTemplate = template ?? buildAlignmentTemplate(row.coreSmiles) ?? this.alignTemplate;
    const tr = ui.element('tr') as HTMLTableRowElement;
    const rhd = ui.element('td') as HTMLTableCellElement;
    rhd.className = 'chem-sar-rhd';
    rhd.appendChild(renderAlignedOnColor(row.coreSmiles, CORE_W, CORE_H, CORE_BG_ARGB, this.alignTemplate));
    rhd.appendChild(ui.divText(this.sideLabel(side), 'chem-sar-row-label'));
    tr.appendChild(rhd);
    for (const ci of cols)
      tr.appendChild(this.buildMatrixCell(matrix, matrix.cells[side.rowIndex][ci], side.rowIndex, ci));
    tbody.appendChild(tr);

    // Trend row: the potency progression along this core, so a parallel trend in the other core
    // is readable at a glance.
    const trendTr = ui.element('tr') as HTMLTableRowElement;
    trendTr.appendChild(ui.element('td'));
    const trendTd = ui.element('td') as HTMLTableCellElement;
    trendTd.colSpan = cols.length;
    const trend = ui.divH([], 'chem-sar-trend');
    const values = cols.map((ci) => matrix.cells[side.rowIndex][ci]);
    const observed = values.filter((c) => c.value !== null).map((c) => c.value as number);
    const bestValue = observed.length ?
      (this.higherIsBetter ? Math.max(...observed) : Math.min(...observed)) : null;
    values.forEach((cell, i) => {
      if (i > 0)
        trend.appendChild(ui.divText('→', 'chem-sar-trend-arrow'));
      if (cell.value === null) {
        trend.appendChild(ui.divText('·', 'chem-sar-trend-empty'));
        return;
      }
      const text = `${cell.kind === 'virtual' ? '~' : ''}${this.formatActivity(cell.value)}`;
      const el = ui.divText(text, 'chem-sar-trend-value');
      if (cell.kind === 'virtual')
        el.classList.add('chem-sar-trend-virtual');
      if (bestValue !== null && cell.value === bestValue && cell.kind === 'real')
        el.classList.add('chem-sar-trend-best');
      trend.appendChild(el);
    });
    trendTd.appendChild(trend);
    trendTr.appendChild(trendTd);
    tbody.appendChild(trendTr);
  }

  /**
   * The SAR transfer view: two cores whose potency trends run in parallel across the substituents
   * they share — within one matrix, or across two matrices (a cross-series transfer between different
   * chemotypes) — with a trend row under each so the parallel is visible directly.
   */
  private buildTransferPane(transfer: Transfer): HTMLElement {
    const matrixA = this.matrices[transfer.a.matrixIndex];
    const from = this.sideLabel(transfer.a);
    const to = this.sideLabel(transfer.b);
    const bar = ui.divH([
      ui.divText(transfer.crossMatrix ? 'Cross-series SAR transfer' : `${matrixA.label} · SAR transfer`,
        'chem-sar-main-title'),
      ui.divText(`${from} → ${to} · across ${transfer.a.position} · r = ${transfer.correlation.toFixed(2)}`),
    ], 'chem-sar-main-bar');

    // When this source core transfers to more than one target, offer a dropdown to switch between them.
    const siblings = this.transferSiblings(transfer);
    let controlBar: HTMLElement | null = null;
    if (siblings.length > 1) {
      const targetInput = ui.input.choice('Transfer to', {
        value: this.transferTargetLabel(transfer),
        items: siblings.map((t) => this.transferTargetLabel(t)),
        onValueChanged: (value) => {
          const pick = siblings.find((t) => this.transferTargetLabel(t) === value) ?? transfer;
          this.selectTransfer(pick);
        },
      });
      ui.tooltip.bind(targetInput.root, () =>
        `Other cores ${from} transfers to — the SAR carries to whichever you pick.`);
      controlBar = ui.divH([targetInput.root], 'chem-sar-control-bar');
    }

    const table = ui.element('table', 'chem-sar-table') as HTMLTableElement;
    const tbody = ui.element('tbody') as HTMLTableSectionElement;
    table.appendChild(tbody);
    this.cellW = this.fitCellWidth(transfer.substituents.length);
    // Build core a's alignment template once — the corner and side a's row both align to it.
    const aTemplate = buildAlignmentTemplate(matrixA.rows[transfer.a.rowIndex].coreSmiles);
    this.alignTemplate = aTemplate;

    const headRow = ui.element('tr') as HTMLTableRowElement;
    const corner = ui.element('td') as HTMLTableCellElement;
    corner.className = 'chem-sar-corner';
    corner.appendChild(ui.divText('Aligned core', 'chem-sar-corner-label'));
    corner.appendChild(renderAlignedOnColor(matrixA.rows[transfer.a.rowIndex].coreSmiles,
      CORE_W, CORE_H, CORE_BG_ARGB, aTemplate));
    headRow.appendChild(corner);
    for (const substSmiles of transfer.substituents) {
      const th = ui.element('th') as HTMLTableCellElement;
      th.className = 'chem-sar-chd';
      th.style.width = `${this.cellW}px`;
      th.appendChild(renderMoleculeOnColor(substSmiles, Math.min(this.cellW, HEADER_W * 2), HEADER_H, HEADER_ARGB));
      th.appendChild(ui.divText(transfer.a.position, 'chem-sar-chd-label'));
      headRow.appendChild(th);
    }
    tbody.appendChild(headRow);

    this.appendTransferSide(tbody, transfer.a, transfer.aCols, aTemplate);
    this.appendTransferSide(tbody, transfer.b, transfer.bCols);

    const noteText = transfer.crossMatrix ?
      `${from} and ${to} are different chemotypes, yet their potencies track together across these ` +
        `shared substituents — so a substituent change learned on ${from} should transfer to ${to}.` :
      `Both cores follow the same trend across ${transfer.a.position}, so a change learned on ${from} ` +
        `is expected to carry over to ${to}.`;
    const note = ui.divText(noteText, 'chem-sar-transfer-note');

    // A plain div, not ui.box: ui.box pins an explicit pixel width, which stops the matrix from
    // growing with the pane.
    const scroll = ui.div([table, this.buildTransferStats(transfer), note], 'chem-sar-main-scroll');
    return ui.divV(controlBar ? [bar, controlBar, scroll] : [bar, scroll], 'chem-sar-main');
  }

  /** The "Transfer statistics" block: rank correlation, fold-change match, and the benefiting cell —
   *  the follower core's untested analog the transferred rule predicts. */
  private buildTransferStats(transfer: Transfer): HTMLElement {
    const stats = transferStats(transfer, this.higherIsBetter);
    const row = (label: string, value: HTMLElement): HTMLElement =>
      ui.divH([ui.divText(label, 'chem-sar-xfer-stat-l'), value], 'chem-sar-xfer-stat');
    const text = (s: string): HTMLElement => ui.divText(s, 'chem-sar-xfer-stat-v');

    const rows = [
      row('Correlation (Pearson r)', text(`r = ${stats.correlation.toFixed(2)}`)),
      row('Fold-change match', text(stats.foldMatch === null ? '—' : `${Math.round(stats.foldMatch * 100)}%`)),
    ];
    if (stats.benefiting) {
      const side = stats.benefiting.side === 'a' ? transfer.a : transfer.b;
      rows.push(row('Benefiting cell', ui.divH([
        text(`${this.sideLabel(side)} ×`),
        renderMoleculeOnColor(transfer.substituents[stats.benefiting.substIndex],
          BENEFIT_MOL_W, BENEFIT_MOL_H, HEADER_ARGB),
      ], 'chem-sar-xfer-ben')));
    }
    return ui.divV([ui.divText('Transfer statistics', 'chem-sar-xfer-title'), ...rows], 'chem-sar-xfer-stats');
  }

  /** Compact summary chips for the current matrix — short labels, full detail on hover: compound
   *  count, cores×substituents, potency range, virtual count, Free-Wilson fit R², and transfer r. */
  private buildChips(matrix: SarMatrix): HTMLElement {
    const chip = (text: string, tip: string, cls = ''): HTMLElement => {
      const el = ui.divText(text, `chem-sar-chip-badge ${cls}`.trim());
      ui.tooltip.bind(el, () => tip);
      return el;
    };
    const items = [
      chip(`${matrix.realCount} cpd`, `${matrix.realCount} observed compounds`),
      chip(`${matrix.rows.length}×${matrix.columns.length}`,
        `${matrix.rows.length} cores × ${matrix.columns.length} substituents`),
    ];
    if (matrix.realCount) {
      items.push(chip(`${this.formatActivity(matrix.minActivity)}–${this.formatActivity(matrix.maxActivity)}`,
        `Activity range across the matrix, on the ${this.scalingLabel} scale`));
    }
    if (matrix.virtualCount)
      items.push(chip(`${matrix.virtualCount} virtual`, `${matrix.virtualCount} predicted (virtual) analog(s)`));
    if (matrix.virtualCount) {
      const conf = matrix.confidence;
      items.push(chip(conf ? `R² ${formatR2(conf.r2)}` : 'R² —',
        conf ? `Free-Wilson leave-one-out fit — R² ${conf.r2.toFixed(2)}, RMSE ${conf.rmse.toFixed(2)}; ` +
          `${conf.n} of ${conf.total} observed cells cross-validated` :
          'Too few observations to cross-validate the Free-Wilson fit'));
    }
    const idx = this.matrices.indexOf(matrix);
    // this.transfers is sorted by correlation (computeAllTransfers), so the first match is the strongest.
    const involving = this.transfers.filter((t) => t.a.matrixIndex === idx || t.b.matrixIndex === idx);
    if (involving.length) {
      const best = involving[0];
      const crossing = involving.some((t) => t.crossMatrix);
      items.push(chip(`r ${best.correlation.toFixed(2)}`,
        crossing ? 'Strongest SAR-transfer correlation for this series, incl. cross-series pairs' :
          'Best row-to-row potency-trend correlation (SAR transfer)', 'chem-sar-chip-transfer'));
    }
    return ui.divH(items, 'chem-sar-chips');
  }

  private buildMatrixPane(matrix: SarMatrix): HTMLElement {
    // Row 1 — title + compact chips (the Free-Wilson fit now rides along as an "R²" chip).
    const infoBar = ui.divH([ui.divText(matrix.label, 'chem-sar-main-title'), this.buildChips(matrix)],
      'chem-sar-main-bar');

    // Row 2 — controls, kept off the info row so neither crowds the other.
    const controls: HTMLElement[] = [];
    if (matrix.positions.length >= 2) {
      const varyValue = this.varyPosition && matrix.positions.includes(this.varyPosition) ?
        this.varyPosition : 'All';
      const varyInput = ui.input.choice('Vary', {
        value: varyValue,
        items: ['All', ...matrix.positions],
        onValueChanged: (value) => {
          this.varyPosition = value === 'All' ? '' : value!;
          this.render();
        },
      });
      controls.push(varyInput.root);
    }
    const labelInput = ui.input.choice('Label', {
      value: this.sortColumnsBy,
      items: COLUMN_SORTS,
      onValueChanged: (value) => {
        this.sortColumnsBy = value!;
        this.render();
      },
    });
    ui.tooltip.bind(labelInput.root, () => 'Annotate each substituent column with a metric — its mean ' +
      'potency (μ) or molecular weight (MW). Columns keep their order; only the caption is added.');
    controls.push(labelInput.root);
    const controlBar = ui.divH(controls, 'chem-sar-control-bar');

    // The grid scrolls and virtualizes internally, so it goes in a flex host that fills the pane (no
    // outer .chem-sar-main-scroll, which would add a second scroll container and its own padding).
    return ui.divV([infoBar, controlBar, this.buildMatrixGrid(matrix)], 'chem-sar-main');
  }

  private render(): void {
    // Keep the navigator scrolled where it was — selecting a card lower down must not jump it to the top.
    const prevNav = this.host.querySelector('.chem-sar-nav-list');
    const navScroll = prevNav instanceof HTMLElement ? prevNav.scrollTop : 0;
    // Before the DOM goes: a transfer pane leaves no grid behind, and buildMatrixGrid installs a new
    // one, so the previous grid (and everything its handlers close over) has to be let go here.
    this.releaseMatrixGrid();
    ui.empty(this.host);
    // Palette values are memoized by variable name only, so a theme switch would otherwise keep
    // painting the previous theme's colors for the life of the session.
    cssColorCache.clear();
    this.cellByMolIdx.clear(); // rebuilt as this render's cells are created
    this.currentCellEl = null; // its DOM element is gone; syncSelection re-points it from currentRowIdx
    if (this.matrices.length === 0) {
      this.host.appendChild(ui.divText(
        'No SAR matrices found. Try lowering the clustering threshold or raising the fragment cutoff.'));
      return;
    }
    this.host.appendChild(this.buildNavigator());
    const newNav = this.host.querySelector('.chem-sar-nav-list');
    if (newNav instanceof HTMLElement)
      newNav.scrollTop = navScroll;
    // A transfer entry opens the dedicated trend view; a matrix entry opens the full matrix.
    if (this.selKind === 'transfer' && this.transfers[this.selIndex])
      this.host.appendChild(this.buildTransferPane(this.transfers[this.selIndex]));
    else {
      const matrix = this.matrices[Math.min(this.selIndex, this.matrices.length - 1)];
      this.host.appendChild(this.buildMatrixPane(matrix));
    }
    this.syncSelection(); // apply the host grid's current selection/row to the freshly-built cells
  }
}
