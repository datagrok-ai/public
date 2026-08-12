/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../../css/sar-matrix.css';
import {drawMoleculeToCanvas, drawRdKitMoleculeToOffscreenCanvas, getRdKitModule, getRdKitService}
  from '../../utils/chem-common-rdkit';
import {getMolSafe} from '../../utils/mol-creation_rdkit';
import {renderMolecule} from '../../rendering/render-molecule';
import * as TextUtils from '@datagrok-libraries/gridext/src/utils/TextUtils';
import {SCALING_METHODS} from '../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {nestByContainment, rankMatrices, SarRankScheme} from './sar-matrix-ranking';
import {computeAllTransfers, DEFAULT_TRANSFER_SIMILARITY, Transfer, TransferSide, transferStats}
  from './sar-matrix-transfer';
import {MAX_SERIES_LEVELS, runSarMatrix, SarGrouping, SarMatrixParams} from './sar-matrix-run';
import {SarMatrix, SarMatrixCell} from './sar-matrix-types';

/** Transparent (alpha 0) so a drawn core has no white box — it blends with the card/pane. */
const CORE_BG_ARGB = 0x00000000;
const HEADER_ARGB = 0xFFF7F7F9;
const WHITE_ARGB = 0xFFFFFFFF;
/** Linear scheme for potency, matching how the platform colors activity/confidence scores:
 *  red at the low end, green at the high end of the scaled activity. */
const ACTIVITY_SCHEME = [DG.Color.red, DG.Color.green];
const CELL_W = 104;
const CELL_H = 76;
/** Cell chip plate: height, inner text padding, and inset from the cell's edge. */
const CHIP_H = 13;
const CHIP_PAD = 3;
const CHIP_MARGIN = 3;
const HEADER_W = 78;
const HEADER_H = 46;
/** Grid column-header height: the R-position band, the substituent depiction and the metric caption. */
const COL_HEADER_H = HEADER_H + 36;
const CORE_W = 132;
/** Room for the grid's own horizontal scrollbar, added when a grid is sized to its rows instead of
 *  being stretched to fill the pane. */
const GRID_SCROLLBAR_H = 18;
/** Small inline thumbnail of the transfer's benefiting substituent, in the statistics block. */
const BENEFIT_MOL_W = 62;
const BENEFIT_MOL_H = 34;
/** Core thumbnail on a navigator card. Shared by the matrix and transfer lists so a card reads the
 *  same in both and the description column is left the same width in each. */
const CARD_CORE_W = 78;
const CARD_CORE_H = 44;
/** Cells grow to fill the pane, but never past this — beyond it the structures just float in space. */
const CELL_W_MAX = 210;
/** Navigator width + paddings + border-spacing, subtracted when fitting cells to the pane. Must track
 *  the `.chem-sar-nav` width in the stylesheet, or the cells are fitted against the wrong pane. */
const NAV_W = 320;
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
/** Name of the running make-list table single-analog "Generate" appends to. */
const MAKELIST_NAME = 'SAR virtual analogs';

type AnalogPanelBuilder = () => HTMLElement;

const COLSORT_POTENCY = 'Potency';
const COLSORT_MW = 'Molecular weight';
/** No "None": sitting beside the potency threshold, a control reading "None" was read as a filter that
 *  was switched off rather than as a column caption, so it always annotates with something. */
const COLUMN_SORTS = [COLSORT_POTENCY, COLSORT_MW];

/** Properties that only reorder/recolor the already-assembled matrices. Changing one must NOT re-run
 *  fragmentation and decomposition — those cost seconds of RDKit worker time and produce identical
 *  matrices. Every other property (columns, scaling, cutoffs, prediction) does change the assembly. */
const RERANK_ONLY_PROPS = ['rankScheme', 'activityDirection'];

/** Properties that change nothing but what is drawn — no re-fragmentation, and no re-ranking either,
 *  which would reorder the navigator for what is only a change of column caption. */
const RENDER_ONLY_PROPS = ['columnCaption', 'idColumnName'];

/** Properties only the transfer scan reads. The matrices themselves are untouched, so re-fragmenting
 *  for one would throw away seconds of RDKit work to recompute a list the current matrices still support. */
const TRANSFER_ONLY_PROPS = ['transferSimilarity'];

/** The viewer's two tabs (MMP-style). The matrix tab is the default; the transfer tab computes its
 *  quadratic core-vs-core comparison lazily, the first time it is opened for the current matrices. */
const TAB_MATRIX = 'SAR Matrix';
const TAB_TRANSFER = 'SAR Transfer';

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

/** Draw a molecule onto an existing canvas with the given ARGB background baked directly into the
 *  RDKit draw call: bond edges are anti-aliased against whatever background the draw is handed, so a
 *  colour applied afterwards would leave a pale fringe around every bond. Straightens the depiction
 *  for a tidy standalone layout (used for substituents). */
function paintMoleculeOnColor(canvas: HTMLCanvasElement, smiles: string, w: number, h: number,
  argb: number): void {
  try {
    drawMoleculeToCanvas(0, 0, w, h, canvas, smiles, null,
      {normalizeDepiction: true, straightenDepiction: true}, null,
      {clearBackground: true, backgroundColour: argbToRgba(argb)});
  } catch (e) {
    // leave the canvas blank on a malformed structure
  }
}

/** {@link paintMoleculeOnColor} onto a fresh canvas, eagerly. */
function renderMoleculeOnColor(smiles: string, w: number, h: number, argb: number): HTMLElement {
  const canvas = ui.canvas(w, h);
  if (smiles)
    paintMoleculeOnColor(canvas as HTMLCanvasElement, smiles, w, h, argb);
  return canvas;
}

/** Drop the aligned-molblock layouts. Called when the matrices are rebuilt, since none of the cached
 *  keys can be hit again. */
export function clearDepictionCaches(): void {
  alignCache.clear();
}

/**
 * Label every attachment point in a molblock as a plain "R". A matrix varies one position, so the map
 * number distinguishes nothing and only invites reading it against the column header — and which
 * position the decomposition happened to call R1 is arbitrary on a symmetric core. RDKit draws a bare
 * dummy as `*` and appends any map number, so the map is dropped upstream and an atom alias supplies
 * the letter. Returns the input unchanged when there is no attachment point to label.
 */
function labelAttachmentPoints(molblock: string): string {
  const lines = molblock.split('\n');
  const atomCount = Number.parseInt((lines[3] ?? '').trim().split(/\s+/)[0], 10);
  const endIdx = lines.findIndex((l) => l.trim() === 'M  END');
  if (!Number.isFinite(atomCount) || endIdx < 0)
    return molblock;
  const aliases: string[] = [];
  for (let i = 0; i < atomCount; i++) {
    // Columns 31-34 of a V2000 atom line hold the element symbol.
    if ((lines[4 + i] ?? '').slice(31, 34).trim() === '*')
      aliases.push(`A  ${String(i + 1).padStart(3, ' ')}`, 'R');
  }
  return aliases.length === 0 ? molblock :
    [...lines.slice(0, endIdx), ...aliases, ...lines.slice(endIdx)].join('\n');
}

/** Depiction as an offscreen bitmap sized in device pixels. A null template gives a standalone
 *  depiction laid out on its own (substituent headers); a template aligns to the shared core. */
function depictionCanvas(molStr: string, template: string | null, dw: number, dh: number,
  argb: number): OffscreenCanvas | null {
  if (!molStr || dw < 1 || dh < 1)
    return null;
  let molCtx = null;
  let labelled = null;
  try {
    // Stripped before alignment, since a molblock encodes the map in a fixed atom-line column where
    // this substitution could not reach it.
    const plain = molStr.replace(/\[\*:\d+\]/g, '[*]');
    molCtx = getMolSafe(template ? alignToTemplate(plain, template) : plain, {}, getRdKitModule());
    if (!molCtx.mol)
      return null;
    if (plain !== molStr) {
      const aliased = labelAttachmentPoints(molCtx.mol.get_molblock());
      labelled = getMolSafe(aliased, {}, getRdKitModule());
      if (labelled.mol) {
        molCtx.mol.delete();
        molCtx = labelled;
        labelled = null;
      }
    }
    if (!molCtx.mol)
      return null;
    if (!template)
      molCtx.mol.set_new_coords(); // no core to align to, so give the fragment a tidy layout of its own
    const offscreen = new OffscreenCanvas(dw, dh);
    drawRdKitMoleculeToOffscreenCanvas(molCtx, dw, dh, offscreen, null,
      {clearBackground: true, backgroundColour: argbToRgba(argb)});
    return offscreen;
  } catch (e) {
    return null; // malformed structure — leave the cell to its background
  } finally {
    labelled?.mol?.delete(); // only set when the aliased re-parse failed and was discarded
    molCtx?.mol?.delete();
  }
}

/**
 * Draw a depiction onto the grid canvas in device pixels: the bitmap is rendered at the cell's device
 * size and blitted 1:1 at a whole-pixel offset, so it is never resampled and bond lines stay crisp —
 * drawing it through the grid's scaled transform instead lands the rect on fractional device pixels
 * and bilinear-filters every bond.
 *
 * The device rect comes from the context's own transform rather than from `devicePixelRatio`, which
 * is only the same number when the grid has installed a pure scale; any translate in the transform
 * would displace every cell by a constant.
 *
 * The transform is reset for the blit but the clip region deliberately is not, which is why this
 * cannot use `putImageData` — that ignores the clip as well, letting a cell scrolled under the pinned
 * core column, or one running past the grid's right edge, paint outside the grid's own bounds.
 */
function drawDepiction(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
  molStr: string, template: string | null, argb: number): void {
  const m = g.getTransform();
  const canvas = depictionCanvas(molStr, template, Math.round(m.a * w), Math.round(m.d * h), argb);
  if (canvas === null)
    return;
  g.save();
  g.setTransform(1, 0, 0, 1, 0, 0);
  g.drawImage(canvas, Math.round(m.a * x + m.e), Math.round(m.d * y + m.f));
  g.restore();
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

/** Context-panel structure sizing: the drawn width is clamped to this range and the height follows.
 *  Below the minimum a structure is unreadable; above the maximum it just wastes panel height. */
const CP_STRUCT_MIN_W = 200;
const CP_STRUCT_MAX_W = 520;
const CP_STRUCT_ASPECT = 0.55;
/** Horizontal chrome between a measured element's client width and the width a structure can occupy:
 *  the box's own padding, and for the panel additionally the panel body's padding and the box border. */
const BOX_CHROME = 12;
const PANEL_CHROME = 30;

const GRID_FONT = 'Roboto, "Segoe UI", sans-serif';

/** Column captions of the per-cell frame the filter runs over. They are what the filter widgets are
 *  labelled with, so they name the two structure axes as the matrix itself does. */
const STRUCT_CORE = 'Core';
const STRUCT_R = 'R';
const STRUCT_POTENCY = 'Potency';
const STRUCT_REFS = 'Reference points';
const STRUCT_MW = 'MW';
/** Column captions of the per-series frame the navigator filter runs over. */
const NAV_SERIES = 'Series';
const NAV_CORE = 'Core';
const NAV_BEST = 'Best';
const NAV_MEAN = 'Mean';
const NAV_SPREAD = 'Spread';
const NAV_BEST_R = 'Best R';
const NAV_COMPOUNDS = 'Compounds';
const NAV_CORES = 'Cores';
const NAV_LEVEL = 'Level';

/** How one of a cell's corner chips is drawn. The two corners carry different kinds of fact — the
 *  potency the matrix is colored by, and the identity of the compound it came from — so they are kept
 *  on opposite diagonals and never collide however wide either one gets. */
interface ChipStyle {
  corner: 'top-left' | 'bottom-right';
  color: string;
  /** Set for a predicted value, to keep it visually separable from a measured one. */
  italic?: boolean;
  /** Dimmed, for a prediction too thin to read as a firm number. */
  faint?: boolean;
}

/** Identity of a transfer's source core within its section (same- and cross-series kept apart). The
 *  nav grouping and the pane's sibling lookup both key on this, so the card "+N" and the pane dropdown
 *  can never disagree about which targets a source reaches. */
function transferSourceKey(t: Transfer): string {
  return `${t.a.matrixIndex}:${t.a.rowIndex}`;
}

/** One displayed grid row: which matrix row it is, and which of that matrix's columns map to the
 *  displayed columns. Each row carries its OWN matrix, because a transfer's two sides can come from
 *  different matrices — the core, the cells, the potency range and the context panel all resolve
 *  through this descriptor rather than through one pane-wide matrix. */
interface PaneRow {
  matrix: SarMatrix;
  rowIndex: number;
  /** One entry per displayed column, indexing into `matrix.columns` / `matrix.cells[rowIndex]`. */
  colIdxs: number[];
  /** Drawn under the core. */
  label: string;
  /** Highlight the row's most potent observed cell. Set for a transfer's two sides, where reading the
   *  trend along each row against the other is the whole point of the view. */
  markBest?: boolean;
}

/** One displayed grid column. */
interface PaneColumn {
  substSmiles: string;
  position: string;
  /** The metric caption under the depiction; empty when the Label control is set to None. */
  caption: string;
}

/** The rendered pane grid plus the per-render row/column state. `onCellRender` fires for every
 *  visible cell on every repaint, so the displayed row/column descriptors, the group boundaries and
 *  the column captions are computed once when the grid is built instead of per cell per repaint. */
interface MatrixGridState {
  grid: DG.Grid;
  /** Scaffold frame backing the grid — one row per displayed row, one string column per displayed
   *  column plus the pinned 'Core'. Held so the grid keeps a live reference for its lifetime. */
  df: DG.DataFrame;
  rows: PaneRow[];
  columns: PaneColumn[];
  /** Grid column NAME (never the grid column index, which pinning shifts) -> index into `columns`. */
  colKeyToIdx: Map<string, number>;
  /** Indices into `columns` that begin an R-position group — only those draw the position label. */
  firstOfGroup: Set<number>;
  /** Structure drawn in the pinned column's header, naming the series the rows belong to. Null when
   *  the rows span more than one matrix, where no single structure describes all of them. */
  headerCore: string | null;
  /** One alignment template for the whole pane, from the shared core. The header, every row key and
   *  every cell are drawn against it, so an attachment point sits in the same place everywhere — laid
   *  out per depiction instead, the same position appears at a different spot in each and the R labels
   *  read as though they were swapped. Null when rows span several matrices and share no core. */
  paneTemplate: string | null;
  /** Per displayed row, the molblock its cells align to. Filled on first paint and kept for the life
   *  of the grid: every repaint would otherwise re-run the alignment for every visible cell. */
  rowTemplates: (string | null)[];
}

/** One live pane grid together with the subscriptions it owns. The matrix and the SAR-transfer panel
 *  hold one each: both can be on screen at the same time, so a single shared slot would have each
 *  rebuild silently unsubscribe the other's painters and drop its Dart-backed grid unreleased. */
interface PaneGridSlot {
  state: MatrixGridState | null;
  subs: {unsubscribe(): void}[];
}

export class SarMatrixViewer extends DG.JsViewer {
  moleculesColumnName: string;
  activityColumnName: string;
  /** Optional column captioning each observed cell — empty leaves cells uncaptioned. */
  idColumnName: string;
  scaling: string;
  activityDirection: string;
  fragmentCutoff: number;
  grouping: string;
  fragmentationLevels: number;
  threshold: number;
  transferSimilarity: number;
  predictVirtual: boolean;
  rankScheme: string;

  private matrices: SarMatrix[] = [];
  /** SAR-transfer entries: every correlated core pair across all matrices, including cross-series
   *  pairs. Computed lazily — empty until the SAR Transfer tab is first opened (see
   *  {@link activateTransferTab}), and cleared whenever the matrices change under it. */
  private transfers: Transfer[] = [];
  /** Whether `transfers` reflects the current `matrices`. False until the transfer tab runs its
   *  lazy computation, and reset by {@link invalidateTransfers} on every rebuild or re-rank. */
  private transfersComputed = false;
  /** Bumped every time the transfer list is dropped. A scan that started before the bump is holding
   *  indices into matrices that have since been rebuilt or reordered, so it must not publish. */
  private transferGeneration = 0;
  /** A lazy transfer computation is already scheduled — a second tab activation must not stack one. */
  private transfersComputing = false;
  /** Index into `matrices` of the series the matrix pane is showing. */
  private selIndex = 0;
  /** Index into `transfers` of the entry the SAR-transfer tab is showing. Deliberately separate from
   *  `selIndex`: the two tabs are navigated independently, so a shared index would make picking a
   *  transfer move the matrix out from under the user. */
  private transferIndex = 0;
  /** "Vary" filter: show only this R-position's column group, or all when empty. */
  private varyPosition = '';
  /** Which metric annotates each substituent column (mean potency or molecular weight). Columns keep
   *  their as-assembled order — this only controls the caption shown under each, never the order. */
  columnCaption: string;
  /** The virtual-analog cell under the last right-click, so the context menu can offer a per-cell
   *  make-list add. Null when the right-click wasn't on an assembled virtual cell. */
  private contextCell: {matrix: SarMatrix, ri: number, ci: number} | null = null;
  /** Per-SMILES builders for this viewer's gated "SAR analysis" info panel, so a clicked analog shows
   *  its SAR context (prediction, support, decomposition, "Add to make-list") alongside the native
   *  Molecule panels. Registered on click; cleared on recompute and detach. Per-instance so one
   *  viewer closing can't wipe another's panels. */
  private readonly analogPanels = new Map<string, AnalogPanelBuilder>();
  private readonly host = ui.divH([], 'chem-sar-matrix');
  /** The SAR Transfer tab's content — sibling pane to `host` in `tabs`, filled lazily on first open. */
  private readonly transferHost = ui.divH([], 'chem-sar-xfer-panel');
  /** MMP-style tab control holding the matrix pane (default) and the SAR-transfer pane. */
  private readonly tabs: DG.TabControl;
  private computing = false;
  /** Set when a recompute is requested while one is already running, so it is re-queued after. */
  private dirty = false;
  private computeTimer = 0;
  /** Set in `detach`, so an in-flight compute can't render into a closed viewer. */
  private detached = false;
  /** Cell width for the current render, fitted to the pane by `fitCellWidth`. */
  private cellW = CELL_W;
  /** The matrix pane's virtualized grid, and the SAR-transfer panel's. Held so selection / current-row
   *  changes can repaint with a cheap `invalidate()` — the highlight is drawn per-cell in
   *  `paintBodyCell`, so no DOM rebuild is needed. */
  private readonly matrixSlot: PaneGridSlot = {state: null, subs: []};
  private readonly transferSlot: PaneGridSlot = {state: null, subs: []};
  /** Last pointer event over the matrix grid overlay, so a cell click can honor ctrl/shift for the
   *  host-grid selection extend (`onCellClick` carries no DOM event). */
  private lastGridMouseEvent: MouseEvent | null = null;
  /** Size observer for the context panel structure currently on screen; replaced when a new cell is
   *  opened so observers do not accumulate one per click. */
  private cpStructureSub: {unsubscribe(): void} | null = null;
  /** Generation of the context-panel structure awaiting layout, so a panel opened while an earlier one
   *  is still waiting cannot have its observer overwritten when that earlier wait resolves. */
  private cpStructureToken = 0;
  /** Matrices whose folded-in children are hidden, by matrix id. Kept on the viewer rather than
   *  rebuilt per render so re-ranking or resizing does not reopen everything the user closed. */
  private readonly collapsed = new Set<string>();
  /** The navigator's card per matrix (index-aligned with `matrices`) and the parent chain they were
   *  built against. Built ONCE per navigator build; selection, collapse and the threshold filter
   *  mutate these cards in place (a class swap / display toggles) instead of rebuilding the list —
   *  with hundreds of series, recreating every card and re-rasterizing every structure on each
   *  click is what made the navigator feel slow. */
  private navCards: HTMLElement[] = [];
  private navParents: number[] = [];
  /** Rasterizes a card's core structure only when the card first becomes visible (initial viewport,
   *  scroll, expand, or filter reveal). Recreated with each navigator build; `navPendingCores`
   *  holds each canvas's deferred draw until it fires. */
  private navCoreObserver: IntersectionObserver | null = null;
  private readonly navPendingCores: Map<Element, () => void> = new Map();
  /** Potency a cell must reach to stay drawn in the matrix pane; null leaves every cell alone. Rows and
   *  columns left with nothing drawn are dropped rather than shown as empty bands. */
  private cellPass: Set<string> | null = null;
  /** One row per series, so the navigator can be filtered by the platform’s own filters: a
   *  substructure sketcher over the cores, histograms over the numbers the cards print. */
  private navFrame: DG.DataFrame | null = null;
  private navFilters: DG.FilterGroup | null = null;
  private navSub: {unsubscribe(): void} | null = null;
  /** Matrix ids the navigator filter admits, or null while it admits everything. */
  private navPass: Set<string> | null = null;
  /** Frame row index -> matrix id, so a filtered frame maps back onto the cards. */
  private navKeys: string[] = [];
  /** Refreshes the navigator match count. Held because the filter group updates the list outside the
   *  render that created it. */
  private navMatchCount: (() => void) | null = null;
  /** Frame row index -> the cell it describes, so a filtered frame maps back onto the matrices. */
  private cellKeys: string[] = [];
  /** One row per distinct (core, substituent) the matrices contain, both typed as molecules, so the
   *  platform's own substructure filter applies to them. Built on first use rather than with the
   *  matrices: nothing needs it until the filter is opened, and building it costs a pass over every
   *  cell. Null until then. */
  private structFrame: DG.DataFrame | null = null;
  private structFilters: DG.FilterGroup | null = null;
  private structSub: {unsubscribe(): void} | null = null;
  /** Parts of the matrix pane the potency threshold updates in place, so moving it repaints rather
   *  than rebuilds. Null while no matrix pane is on screen. */
  private paneGridHost: HTMLElement | null = null;
  private paneEmptyNote: HTMLElement | null = null;
  private paneDimsChip: HTMLElement | null = null;
  /** Cleared when a new set of matrices arrives, so the tree closes to its roots once per analysis
   *  rather than snapping shut again every time the list is redrawn. */
  private collapseSeeded = false;
  get helpUrl() {
    return 'https://raw.githubusercontent.com/datagrok-ai/public/refs/heads/master/help/datagrok/solutions/domains/chem/chem.md#sar-matrix';
  }
  constructor() {
    super();
    // Data properties rather than hidden strings, so both columns are pickable from the property panel
    // like any other viewer's. `column()` appends "ColumnName" to the stem, so the field names and
    // everything that sets them through `setOptions` are unchanged.
    //
    // No options: the helper already registers these as string-typed data properties, and passing a
    // `type` or `semType` overwrites that descriptor and leaves the platform unable to resolve the
    // property at all ("Property type not found" at construction, which yields a viewer with no
    // properties whatsoever).
    this.moleculesColumnName = this.column('molecules');
    this.activityColumnName = this.column('activity');
    // Optional: when set, each observed cell captions its structure with that column's value, so a
    // compound in the matrix can be named against the source table. Empty leaves the cells as they are.
    this.idColumnName = this.column('id');
    this.scaling = this.string('scaling', SCALING_METHODS.MINUS_LG, {choices: Object.values(SCALING_METHODS)});
    this.activityDirection = this.string('activityDirection', DIR_AUTO, {choices: ACTIVITY_DIRECTIONS});
    this.fragmentCutoff = this.float('fragmentCutoff', 0.4);
    this.grouping = this.string('grouping', SarGrouping.Site, {choices: Object.values(SarGrouping)});
    // Named for what the navigator shows rather than for the fragmentation passes underneath: the
    // cards are series badged L1, L2, L3, and this is how many of those tiers get built.
    this.fragmentationLevels = this.int('fragmentationLevels', 3,
      {min: 1, max: MAX_SERIES_LEVELS, friendlyName: 'Series levels'});
    // Only read under Similarity grouping; Site grouping matches exactly and has nothing to tune.
    this.threshold = this.float('threshold', 0.5);
    // How alike the two scaffolds must be, once the R-group is held identical. Read only by the
    // transfer scan, so it retunes the tab without re-fragmenting anything.
    this.transferSimilarity = this.float('transferSimilarity', DEFAULT_TRANSFER_SIMILARITY,
      {min: 0.1, max: 1, friendlyName: 'Transfer similarity',
        description: 'How alike two compounds carrying the same R-group must be for a transfer between their scaffolds to count'});
    this.predictVirtual = this.bool('predictVirtual', true);
    // "SAR transfer" is a tab, not a ranking mode, so it is not offered as a rank scheme.
    this.rankScheme = this.string('rankScheme', SarRankScheme.Potency,
      {choices: [SarRankScheme.Potency, SarRankScheme.Discontinuity, SarRankScheme.Preferred]});
    // Annotates the substituent columns; changing it repaints and never re-fragments.
    this.columnCaption = this.string('columnCaption', COLSORT_POTENCY, {choices: COLUMN_SORTS});
    this.host.style.height = '100%';
    this.transferHost.style.height = '100%';

    // MMP-style layout: the matrices and the transfer view are tabs of one control rather than the
    // transfer being docked beneath the viewer. Detecting transfers compares every core against every
    // other one — quadratic in the total row count — so that tab computes on first open, not up front.
    this.tabs = ui.tabControl(null, false);
    const matrixPane = this.tabs.addPane(TAB_MATRIX, () => this.host);
    ui.tooltip.bind(matrixPane.header, 'Core × substituent potency matrices, one per series');
    const transferPane = this.tabs.addPane(TAB_TRANSFER, () => this.transferHost);
    ui.tooltip.bind(transferPane.header, 'Pairs of cores whose potency trends run in parallel across the ' +
      'R-groups they have both explored — detected when the tab is first opened');
    // Lives with the tab control itself (not this.subs): tab switching must keep working across a
    // detach/re-attach cycle, and the subscription dies with the viewer's DOM anyway.
    this.tabs.onTabChanged.subscribe(() => {
      if (this.tabs.currentPane?.name === TAB_TRANSFER)
        this.activateTransferTab();
      else {
        // A pane rebuilt while its tab was hidden sat in a display:none host; repaint it now that
        // its canvas has real dimensions.
        this.matrixSlot.state?.grid.invalidate();
      }
    });
    this.tabs.root.style.width = '100%';
    this.tabs.root.style.height = '100%';
    this.root.appendChild(this.tabs.root);
  }

  onTableAttached(): void {
    this.detached = false; // the flag tracks the CURRENT attachment, not whether one ever ended
    // Cells are fitted to the pane width, so a resize re-fits the columns of the grid already on
    // screen. Rebuilding the pane instead would replace that grid, and a new grid starts at the first
    // row — so every step of a splitter drag would throw away where the user had scrolled to.
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 200).subscribe(() => {
      if (!this.computing && this.matrices.length)
        this.refitColumns();
    }));
    this.subs.push(this.onContextMenu.subscribe((menu) => this.buildContextMenu(menu)));
    // Two-way link with the host grid: a selection change there re-rings the matching matrix cells.
    // The current row is deliberately not watched — nothing is painted from it, and a repaint
    // re-rasterizes every visible structure, so following it would cost a full redraw per click for
    // no visible change.
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50)
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
    this.cpStructureSub?.unsubscribe(); // the panel outlives the viewer; its observer must not
    this.cpStructureSub = null;
    this.releaseMatrixGrid();
    // The transfer tab's DOM goes with the root, but its Dart-backed grid has to be released by hand.
    this.releaseSlot(this.transferSlot);
    this.navCoreObserver?.disconnect();
    this.navCoreObserver = null;
    this.navPendingCores.clear();
    // The pair frame is held by a headless view of its own, so nothing else drops this subscription.
    this.structSub?.unsubscribe();
    this.structSub = null;
    this.navSub?.unsubscribe();
    this.navSub = null;
    this.navMatchCount = null;
    super.detach();
  }

  /** The matrix pane's grid — what the cell filter, the dimensions chip and `visibleDims` act on.
   *  Named separately from the transfer panel's so a filter can never reach across to it. */
  private get matrixGrid(): MatrixGridState | null {
    return this.matrixSlot.state;
  }

  /** Let go of one pane grid. Dropping the DOM is not enough: the render/click/tooltip subscriptions
   *  keep the grid — and through it its scaffold DataFrame and the whole SarMatrix — reachable, and a
   *  Dart-backed grid is only released when it is closed. */
  private releaseSlot(slot: PaneGridSlot): void {
    slot.subs.forEach((s) => s.unsubscribe());
    slot.subs = [];
    try {
      slot.state?.grid?.close?.();
    } catch (e) {
      // A standalone (view-less) grid may not support close; dropping the reference is enough.
    }
    slot.state = null;
  }

  /** Release the matrix pane's grid and the pane elements that point into it. */
  private releaseMatrixGrid(): void {
    this.releaseSlot(this.matrixSlot);
    // The pane those pointed into is about to be replaced; holding them would let the threshold write
    // into detached elements.
    this.paneGridHost = null;
    this.paneEmptyNote = null;
    this.paneDimsChip = null;
  }

  /** Reflect the host grid's selection and current row onto the rendered cells. The pane grid draws
   *  both rings per cell in `paintBodyCell` and reads them straight off the dataframe, so a host-grid
   *  change only needs a repaint of the visible cells. */
  private syncSelection(): void {
    if (!this.dataFrame)
      return;
    // Both panes paint the host grid's selection ring, so both have to be repainted for a click in one
    // to show up in the other.
    this.matrixSlot.state?.grid.invalidate();
    this.transferSlot.state?.grid.invalidate();
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
    const current = this.matrices[Math.min(this.selIndex, this.matrices.length - 1)];
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
      molCol('Core', cells.map((c) => c.matrix.rows[c.ri].keySmiles)),
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
    if (property !== null && RENDER_ONLY_PROPS.includes(property.name)) {
      // Captions and the id column show up in the pane's cells only, so the navigator stays as-is.
      this.renderMatrixPane();
      return;
    }
    if (property !== null && RERANK_ONLY_PROPS.includes(property.name)) {
      this.reRank();
      return;
    }
    if (property !== null && TRANSFER_ONLY_PROPS.includes(property.name)) {
      // Drop the list and let the tab rebuild it; if the tab is not showing, nothing is recomputed
      // at all until it next is.
      this.invalidateTransfers();
      if (this.transferTabActive)
        this.activateTransferTab();
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
    // Transfers index into `matrices` by position, so a reorder invalidates them; the transfer tab
    // recomputes lazily the next time it is shown (render() below refreshes it if it is up now).
    this.invalidateTransfers();
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
    // The transfers on screen index into the matrices about to be replaced, so the transfer tab is
    // cleared now rather than being left pointing at rows that no longer exist if this compute fails.
    this.invalidateTransfers();
    ui.empty(this.host);
    this.host.appendChild(ui.loader());
    const progress = DG.TaskBarProgressIndicator.create('Building SAR matrices...');
    try {
      const params: SarMatrixParams = {
        scaling: this.scaling as SCALING_METHODS,
        fragmentCutoff: this.fragmentCutoff,
        predictVirtual: this.predictVirtual,
        grouping: this.grouping as SarGrouping,
        fragmentationLevels: this.fragmentationLevels,
        higherIsBetter: this.higherIsBetter,
        threshold: this.threshold,
        rankScheme: this.rankScheme as SarRankScheme,
      };
      const matrices = await runSarMatrix(molecules, activity as DG.Column<number>, params);
      // The viewer may have been closed while the workers were running; rendering into it now would
      // build a grid nothing will ever release.
      if (this.detached)
        return;
      this.matrices = matrices;
      this.collapseSeeded = false;
      this.selIndex = 0;
      this.render();
    } catch (e) {
      if (this.detached)
        return; // a viewer the user already closed must not report its own teardown as a failure
      const message = e instanceof Error ? e.message : String(e);
      ui.empty(this.host);
      this.host.appendChild(ui.divText(`SAR Matrix failed: ${message}`));
      // The transfer tab may still be showing "Building SAR matrices..." from a mid-compute visit;
      // there is no render() on this path to replace it.
      ui.empty(this.transferHost);
      this.transferHost.appendChild(ui.divText(`SAR Matrix failed: ${message}`, 'chem-sar-empty-note'));
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
   *  AND across matrices (differently-scaffolded series) — strongest first. Runs only from the
   *  transfer tab's lazy activation, so the quadratic core-vs-core pass is never paid by a session
   *  that only ever reads the matrices.
   *
   * Left empty until then, which is what keeps it out of the rest of the viewer: the per-matrix
   * correlation chip and the trend view all read this list, so none of them has to test whether the
   * tab was ever opened. */
  /** Drop the transfer list and the tab's content: the matrices it indexed into are gone (rebuilt or
   *  reordered). The next visit to the transfer tab recomputes against the current matrices. */
  private invalidateTransfers(): void {
    this.transferGeneration++;
    this.releaseSlot(this.transferSlot);
    ui.empty(this.transferHost);
    this.transfers = [];
    this.transferIndex = 0;
    this.transfersComputed = false;
  }

  /** The potency tint composited over white and returned opaque, for use as a depiction's background.
   *  RDKit anti-aliases bond edges against whatever background the draw is handed, so a translucent
   *  tint would let every bond blend toward white and leave a pale fringe once the cell is painted. */
  private flatTint(matrix: SarMatrix, value: number, alpha: number): number {
    const c = this.tint(matrix, value, alpha);
    const a = alpha / 255;
    const overWhite = (channel: number): number => Math.round(channel * a + 255 * (1 - a));
    return DG.Color.argb(255, overWhite(DG.Color.r(c)), overWhite(DG.Color.g(c)), overWhite(DG.Color.b(c)));
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

  /** Rows and columns of a matrix that survive the current Vary and cell filters — what the pane
   *  actually draws, as opposed to what the matrix holds. */
  private visibleDims(matrix: SarMatrix): {rows: number, cols: number} {
    const all = this.visibleColIdxs(matrix);
    if (!this.cellFilterActive)
      return {rows: matrix.rows.length, cols: all.length};
    const cols = all.filter((ci) => matrix.cells.some((_row, ri) => this.cellVisible(matrix, ri, ci)));
    const rows = matrix.rows.filter((_row, ri) =>
      cols.some((ci) => this.cellVisible(matrix, ri, ci))).length;
    return {rows, cols: cols.length};
  }

  /** True while the filter is narrowing anything, so the untouched case skips the work entirely. */
  private get cellFilterActive(): boolean {
    return this.cellPass !== null;
  }

  /** Whether the filter admits this cell. Everything passes while no filter is set, which is what
   *  keeps the lookup off the painting path until the user asks for it. */
  private cellVisible(matrix: SarMatrix, ri: number, ci: number): boolean {
    return this.cellPass === null || this.cellPass.has(`${matrix.id}:${ri}:${ci}`);
  }


  /**
   * One row per matrix cell, carrying everything the filter can ask about: the row's core and the
   * column's substituent as molecules, and the cell's potency, reference points and substituent
   * weight as numbers. Typed this way the platform supplies every filter itself — sketchers for the
   * two structures, histograms for the three numbers — which is the same set MMP gets over its
   * fragment table, rather than controls hand-built to look like them.
   *
   * Per cell rather than per distinct pair, because potency and reference count belong to a cell: a
   * pair repeated across matrices has a different value in each, and there would be no single row to
   * put on the histogram.
   */
  private buildStructFrame(): DG.DataFrame {
    const cores: string[] = [];
    const subs: string[] = [];
    const potency: number[] = [];
    const refs: number[] = [];
    const weights: number[] = [];
    this.cellKeys = [];
    for (const matrix of this.matrices) {
      for (let ri = 0; ri < matrix.rows.length; ri++) {
        for (let ci = 0; ci < matrix.columns.length; ci++) {
          const cell = matrix.cells[ri][ci];
          if (cell.kind === 'empty' || cell.value === null)
            continue; // nothing to filter on, and it already draws blank
          cores.push(matrix.rows[ri].coreSmiles);
          subs.push(matrix.columns[ci].substSmiles);
          potency.push(cell.value);
          // A measured compound has no reference points to rest on, so it is given the count of the
          // measured neighbours it shares a row and column with — the same quantity, defined for it —
          // rather than a blank the histogram would have to guess at.
          refs.push(cell.references ?? this.observedNeighbours(matrix, ri, ci));
          const mw = substituentMW(matrix.columns[ci].substSmiles);
          weights.push(Number.isFinite(mw) ? mw : 0);
          this.cellKeys.push(`${matrix.id}:${ri}:${ci}`);
        }
      }
    }
    const coreCol = DG.Column.fromStrings(STRUCT_CORE, cores);
    const subCol = DG.Column.fromStrings(STRUCT_R, subs);
    coreCol.semType = DG.SEMTYPE.MOLECULE;
    subCol.semType = DG.SEMTYPE.MOLECULE;
    return DG.DataFrame.fromColumns([
      coreCol,
      subCol,
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, STRUCT_POTENCY, potency),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, STRUCT_REFS, refs),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, STRUCT_MW, weights),
    ]);
  }

  /** Measured compounds sharing this cell's row or column — what a prediction here would rest on. */
  private observedNeighbours(matrix: SarMatrix, ri: number, ci: number): number {
    let n = 0;
    for (let c = 0; c < matrix.columns.length; c++) {
      if (c !== ci && matrix.cells[ri][c].kind === 'real')
        n++;
    }
    for (let r = 0; r < matrix.rows.length; r++) {
      if (r !== ri && matrix.cells[r][ci].kind === 'real')
        n++;
    }
    return n;
  }

  /** Read the frame's filter back into the passing-cell set, then repaint. Left null while everything
   *  passes, so an untouched filter costs nothing on the per-cell path. */
  private syncStructFilter(): void {
    const frame = this.structFrame;
    if (frame === null)
      return;
    if (frame.filter.trueCount === frame.rowCount)
      this.cellPass = null;
    else {
      const pass = new Set<string>();
      for (let i = 0; i < frame.rowCount; i++) {
        if (frame.filter.get(i))
          pass.add(this.cellKeys[i]);
      }
      this.cellPass = pass;
    }
    this.applyCellFilter();
  }


  /** The pair frame's filter group, built on first use. A filter group belongs to a view, so a
   *  headless one is created for the frame and never shown — only its filters are, in a popup. */
  private structureFilterRoot(): HTMLElement {
    if (this.structFilters === null) {
      this.structFrame = this.buildStructFrame();
      const view = DG.TableView.create(this.structFrame, false);
      // The default set, not a hand-picked one: the platform already gives a molecule column its
      // sketcher filter and a numeric column its histogram, along with the per-filter and
      // whole-group on/off switches. Naming the filters explicitly only risks diverging from that.
      this.structFilters = view.getFiltersGroup();
      this.structSub = DG.debounce(this.structFrame.onFilterChanged, 300).subscribe(() => this.syncStructFilter());
    }
    return this.structFilters.root;
  }


  /**
   * One row per series, carrying its drawn core as a molecule and the numbers its card prints. Same
   * mechanism as the matrix pane's filter: the platform supplies a sketcher for the core and a
   * histogram for each number, so finding a series by substructure and by potency is one filter group
   * rather than a bespoke control per quantity.
   */
  private buildNavFrame(): DG.DataFrame {
    const series: string[] = [];
    const cores: string[] = [];
    const best: number[] = [];
    const mean: number[] = [];
    const spread: number[] = [];
    const bestR: number[] = [];
    const compounds: number[] = [];
    const coreCount: number[] = [];
    const level: number[] = [];
    this.navKeys = [];
    for (const matrix of this.matrices) {
      series.push(matrix.label);
      cores.push(this.matrixCore(matrix));
      // Direction-adjusted, so "more potent" is always the larger number on the histogram whichever
      // way the assay runs — otherwise the same drag means opposite things on -lg and raw scales.
      const top = matrix.realCount ? (this.higherIsBetter ? matrix.maxActivity : -matrix.minActivity) : NaN;
      best.push(top);
      const m = this.meanRealActivity(matrix);
      mean.push(m === null ? NaN : (this.higherIsBetter ? m : -m));
      spread.push(matrix.scores[SarRankScheme.Discontinuity] ?? NaN);
      const pref = matrix.scores[SarRankScheme.Preferred];
      bestR.push(pref === undefined ? NaN : (this.higherIsBetter ? pref : -pref));
      compounds.push(matrix.realCount);
      coreCount.push(matrix.rows.length);
      level.push(matrix.level - 1);
      this.navKeys.push(matrix.id);
    }
    const coreCol = DG.Column.fromStrings(NAV_CORE, cores);
    coreCol.semType = DG.SEMTYPE.MOLECULE;
    return DG.DataFrame.fromColumns([
      // First, and a plain string: the platform gives a string column a category list with a search
      // box, which is how a named series is actually looked for — typed, not picked off a histogram.
      DG.Column.fromStrings(NAV_SERIES, series),
      coreCol,
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, NAV_BEST, best),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, NAV_MEAN, mean),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, NAV_SPREAD, spread),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, NAV_BEST_R, bestR),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, NAV_COMPOUNDS, compounds),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, NAV_CORES, coreCount),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, NAV_LEVEL, level),
    ]);
  }

  /** Read the series frame's filter back into the passing-matrix set and reapply it to the cards. */
  private syncNavFilter(): void {
    const frame = this.navFrame;
    if (frame === null)
      return;
    if (frame.filter.trueCount === frame.rowCount)
      this.navPass = null;
    else {
      const pass = new Set<string>();
      for (let i = 0; i < frame.rowCount; i++) {
        if (frame.filter.get(i))
          pass.add(this.navKeys[i]);
      }
      this.navPass = pass;
    }
    this.updateNavVisibility();
    this.navMatchCount?.();
  }

  /** The series filter group, built on first use, mounted in a popup off the navigator's filter icon. */
  private navFilterRoot(): HTMLElement {
    if (this.navFilters === null) {
      this.navFrame = this.buildNavFrame();
      const view = DG.TableView.create(this.navFrame, false);
      this.navFilters = view.getFiltersGroup();
      // Asked for by name, because the default set leaves it out: every series name is distinct, so a
      // category list would hold one entry per row and the platform declines to build one. It is still
      // how a named series is looked for — the list has a search box, which is what turns "3" into
      // every series carrying a 3 — so the filter is added explicitly.
      this.navFilters.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: NAV_SERIES}, false);
      this.navSub = DG.debounce(this.navFrame.onFilterChanged, 300).subscribe(() => this.syncNavFilter());
    }
    return this.navFilters.root;
  }

  /** Whether a matrix clears the series filter. Unfiltered until the filter group narrows the frame. */
  private passesFilter(matrix: SarMatrix): boolean {
    return this.navPass === null || this.navPass.has(matrix.id);
  }

  /** The score block on the right of a navigator card: one big value per line with its caption under
   *  it. Shared by the series and transfer cards so a one-line score sits exactly where a two-line one
   *  does and the two lists read as one column of scores. */
  private cardScoreBox(lines: {value: string, label: string}[], tip: () => string): HTMLElement {
    const box = ui.divV(lines.map((ln) => ui.divH([
      ui.divText(ln.value, 'chem-sar-card-score'),
      ui.divText(ln.label, 'chem-sar-card-score-cap'),
    ], 'chem-sar-card-score-line')), 'chem-sar-card-score-box');
    ui.tooltip.bind(box, tip);
    return box;
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

  /**
   * The structure a matrix is keyed on, in a form the depiction can draw — what every one of its rows
   * has in common and, for a coarser matrix, what its children agree on a cut deeper. Drawing the
   * first row's core instead would show the same picture for a parent and its children, since a row
   * of the parent is a core of a child.
   *
   * Each fragmentation level marks the sites it inherits with its own isotope, which would render as
   * numbered stars; renumbering them as attachment points draws them as R labels like every other core.
   */
  private matrixCore(matrix: SarMatrix): string {
    if (!matrix.siteKey)
      return matrix.rows[0]?.coreSmiles ?? '';
    let n = 0;
    return matrix.siteKey.replace(/\[\d+\*\]|\[\*:\d+\]/g, () => `[*:${++n}]`);
  }

  /**
   * Parent index per matrix, `-1` for a root, against the current (ranked) order. A matrix built at a
   * coarser fragmentation level owns the ones folded into it; where no coarser level was asked for,
   * compound containment stands in, so the navigator still groups related views rather than listing
   * every matrix flat.
   */
  private matrixParents(): number[] {
    const byId = new Map<string, number>();
    this.matrices.forEach((matrix, i) => byId.set(matrix.id, i));
    const linked = this.matrices.map((matrix) =>
      matrix.parentId !== undefined ? byId.get(matrix.parentId) ?? -1 : -1);
    return linked.some((p) => p >= 0) ? linked : nestByContainment(this.matrices);
  }

  /** One selectable matrix card: the aligned core drawn (so the matrix is identified by its
   *  scaffold, not just "Series A"), a descriptor line, and the rank score. `depth` indents a matrix
   *  folded into the one above it; `folded` is how many are folded into this one, which gives the card
   *  its expander. */
  private buildCard(matrix: SarMatrix, index: number, depth = 0, folded = 0,
    onToggle?: () => void): HTMLElement {
    // Held to what the list alone has to answer — how big a series is, and how much sits under it.
    // The description column is barely 110px wide once the structure and the score are placed, so
    // every extra fact costs a wrapped line on every card; the virtual count and the Free-Wilson R²
    // are on the matrix pane's own chips the moment a card is opened.
    const desc = `${matrix.rows.length} cores · ${matrix.positions.join('/')} · ${matrix.realCount} cpd` +
      (folded ? ` · ${folded} inside` : '');
    const core = this.lazyCoreDepiction(this.matrixCore(matrix), CARD_CORE_W, CARD_CORE_H);
    core.classList.add('chem-sar-card-core');
    // Numbered from the matrices, not from the fragmentation passes: the first pass produces series,
    // which are the ROWS of a matrix rather than anything the navigator lists, so calling the first
    // listable thing "L2" invited the question of where L1 had gone.
    // Branches also bottom out at different depths — a key already reduced to a bare ring has nothing
    // left to strip — so indentation alone does not say which level a card sits at. The badge does.
    const shown = matrix.level - 1;
    const level = ui.divText(`L${shown}`, 'chem-sar-card-level');
    ui.tooltip.bind(level, () => shown === 1 ?
      'Level 1 — one matrix per site. Each of its rows is a series: one core with its substituents.' :
      `Level ${shown} — a coarser matrix, holding the compounds of the level-${shown - 1} matrices ` +
      'whose cores agree one further cut deeper.');
    const body = ui.divV([
      ui.divH([ui.divText(matrix.label, 'chem-sar-card-name'), level], 'chem-sar-card-title'),
      ui.divText(desc, 'chem-sar-card-desc'),
    ], 'chem-sar-card-body');
    const sc = this.cardScore(matrix);
    const scoreBox = this.cardScoreBox(sc.lines, () => sc.tip);
    // The expander is its own click target: the rest of the card selects the matrix, and collapsing
    // must not also switch the pane to it. Leaves get a spacer so every core lines up regardless.
    let twisty: HTMLElement;
    if (folded > 0 && onToggle) {
      twisty = ui.iconFA(!this.collapsed.has(matrix.id) ? 'chevron-down' : 'chevron-right', (e: MouseEvent) => {
        e.stopPropagation();
        onToggle();
        // The card is not rebuilt on toggle (the list mutates in place), so flip the icon here.
        const open = !this.collapsed.has(matrix.id);
        twisty.classList.toggle('fa-chevron-down', open);
        twisty.classList.toggle('fa-chevron-right', !open);
      });
      twisty.classList.add('chem-sar-card-twisty');
      ui.tooltip.bind(twisty, () => this.collapsed.has(matrix.id) ?
        `Show ${folded} matrices folded into this one` : 'Hide the matrices folded into this one');
    } else
      twisty = ui.div([], 'chem-sar-card-twisty');

    const card = ui.divH([twisty, core, body, scoreBox], 'chem-sar-card');
    if (depth > 0) {
      // Capped: past a couple of levels the indent would eat the card rather than show the nesting.
      card.style.marginLeft = `${Math.min(depth, 3) * 10}px`;
      card.classList.add('chem-sar-card-nested');
    }
    if (index === this.selIndex)
      card.classList.add('selected');
    card.onclick = () => this.selectMatrix(index);
    return card;
  }

  /** Card click: move the highlight and rebuild only the matrix pane. The navigator list — hundreds
   *  of cards and their depictions — is left untouched; rebuilding it per selection is exactly the
   *  cost this viewer used to pay on every click. */
  private selectMatrix(index: number): void {
    if (index === this.selIndex)
      return;
    this.navCards[this.selIndex]?.classList.remove('selected');
    this.selIndex = index;
    this.navCards[index]?.classList.add('selected');
    this.varyPosition = ''; // the Vary filter is per-matrix; don't carry it to a different one
    this.renderMatrixPane();
  }

  /** A card-core canvas whose molecule is drawn on first visibility (via `navCoreObserver`): with
   *  hundreds of series, rasterizing every structure up front is most of what made the navigator
   *  expensive to build. Eager fallback when no observer is active. */
  private lazyCoreDepiction(smiles: string, w: number, h: number): HTMLElement {
    if (this.navCoreObserver === null)
      return renderMoleculeOnColor(smiles, w, h, CORE_BG_ARGB);
    const canvas = ui.canvas(w, h);
    if (smiles) {
      this.navPendingCores.set(canvas, () => paintMoleculeOnColor(canvas, smiles, w, h, CORE_BG_ARGB));
      this.navCoreObserver.observe(canvas);
    }
    return canvas;
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
  /** Select a specific transfer (by its identity in this.transfers) and redraw the trend view.
   *  Only the transfer tab is rebuilt — the matrix tab is untouched and must not move. */
  private selectTransfer(transfer: Transfer): void {
    this.transferIndex = Math.max(0, this.transfers.indexOf(transfer));
    this.renderTransferPanel();
  }

  /** Whether the SAR Transfer tab is the one on screen. */
  private get transferTabActive(): boolean {
    return this.tabs.currentPane?.name === TAB_TRANSFER;
  }

  /**
   * Bring the SAR Transfer tab up to date, computing the transfer list on first open. The quadratic
   * core-vs-core comparison is exactly what the lazy tab defers: it runs here — behind a spinner and
   * a timeout so the tab switch paints first — and never anywhere else. Called on every switch to the
   * tab and from `render()` while the tab is up; both are no-ops when the list is already current.
   */
  private activateTransferTab(): void {
    if (this.transfersComputing)
      return;
    if (this.computing) {
      // The matrices themselves are still being built; compute() ends in render(), which re-enters
      // here once there is something to compare.
      ui.empty(this.transferHost);
      this.transferHost.appendChild(ui.divText('Building SAR matrices...', 'chem-sar-empty-note'));
      return;
    }
    if (this.transfersComputed) {
      // The list is current; rebuild the pane only if it is not on screen (first show after a render
      // that happened while this tab was hidden), otherwise just repaint the live grid.
      if (this.transferHost.childElementCount === 0)
        this.renderTransferPanel();
      else
        this.transferSlot.state?.grid.invalidate();
      return;
    }
    this.transfersComputing = true;
    ui.empty(this.transferHost);
    ui.setUpdateIndicator(this.transferHost, true, 'Detecting SAR transfers...');
    // Timeout, MMP-style: let the tab switch and the indicator paint before the synchronous
    // quadratic pass freezes the frame.
    setTimeout(async () => {
      // A matrix rebuild that started meanwhile is about to replace `matrices`; computing now would
      // mark stale transfers as current. Its final render() re-enters the activation instead.
      if (this.detached || this.computing) {
        this.transfersComputing = false;
        ui.setUpdateIndicator(this.transferHost, false);
        return;
      }
      // Held across the await, not just up to it: fingerprinting goes out to the RDKit workers, and a
      // tab bounce or a threshold change in that window would otherwise start a second full scan.
      const generation = this.transferGeneration;
      let detected: Transfer[] | null = null;
      try {
        detected = await computeAllTransfers(this.matrices, this.transferSimilarity, this.higherIsBetter);
      } catch (e) {
        console.error('SAR Matrix | transfer detection failed', e);
      }
      this.transfersComputing = false;
      ui.setUpdateIndicator(this.transferHost, false);
      // The await handed control back to the event loop. Anything that dropped the list while it was
      // out there bumped the generation, and these indices no longer address the matrices on screen —
      // publishing them would label a card with a core that is not the one it was computed from.
      if (this.detached || generation !== this.transferGeneration)
        return;
      if (detected === null) {
        ui.empty(this.transferHost);
        this.transferHost.appendChild(ui.divText('SAR transfer detection failed — see the browser console.',
          'chem-sar-empty-note'));
        return;
      }
      this.transfers = detected;
      this.transfersComputed = true;
      this.transferIndex = 0;
      this.renderTransferPanel();
    }, 100);
  }

  /**
   * Rebuild the SAR-transfer tab: the list of detected transfers alongside the trend view of the
   * selected one, or a note when the current matrices hold no transfer worth showing.
   */
  private renderTransferPanel(): void {
    // The tab's own grid, not the matrix's: releasing here is what keeps switching transfers from
    // stacking up Dart-backed grids, and it cannot touch the matrix tab's slot.
    this.releaseSlot(this.transferSlot);
    ui.empty(this.transferHost);
    if (this.transfers.length === 0) {
      this.transferHost.appendChild(ui.divText(this.matrices.length === 0 ?
        'No SAR matrices to compare — build matrices first on the SAR Matrix tab.' :
        'No SAR transfers detected: no two cores share enough R-groups, on scaffolds alike enough, ' +
        'with correlated potency trends. Lower "Transfer similarity" to accept less alike scaffolds.',
      'chem-sar-empty-note'));
      return;
    }
    if (this.transferIndex >= this.transfers.length)
      this.transferIndex = 0;
    const list = ui.div([], 'chem-sar-xfer-list');
    for (const group of this.groupTransfersBySource(this.transfers))
      list.appendChild(this.buildTransferCard(group));
    this.transferHost.appendChild(list);
    this.transferHost.appendChild(this.buildTransferPane(this.transfers[this.transferIndex]));
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
    const selected = this.transfers[this.transferIndex] ?? null;
    // Show the selected transfer when it belongs to this group, otherwise the strongest (first).
    const shown = (selected && group.includes(selected)) ? selected : group[0];

    // Only the source core is named. Which target it carries to is a choice, not an identity, and the
    // pane's own dropdown is where it is made — repeating it here cost a wrapped line per card.
    const body = ui.divV([
      ui.divText(this.sideLabel(shown.a), 'chem-sar-card-title'),
    ], 'chem-sar-card-body');
    // The source core the card is named after, drawn the way a matrix card draws its series core.
    // Painted up front rather than through the navigator's deferred observer: that observer belongs to
    // the matrix list and is disconnected whenever it refills, which would leave these blank.
    const sourceRow = this.matrices[shown.a.matrixIndex].rows[shown.a.rowIndex];
    const core = renderMoleculeOnColor(sourceRow.keySmiles, CARD_CORE_W, CARD_CORE_H, CORE_BG_ARGB);
    core.classList.add('chem-sar-card-core');
    const scoreBox = this.cardScoreBox([{value: `r ${shown.correlation.toFixed(2)}`, label: 'transfer'}],
      () => `${this.sideLabel(shown.a)} → ${this.transferTargetLabel(shown)}. Correlation of two ` +
        'different chemotypes’ potency trends across the R-groups they have both explored — the SAR ' +
        'learned on one should carry to the other. 1.00 = perfectly parallel.');
    const card = ui.divH([core, body, scoreBox], 'chem-sar-card');
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
    // The three schemes answer different questions and the names alone do not say which, so the whole
    // set is spelled out rather than only the one in force — choosing between them is the point.
    ui.tooltip.bind(rankInput.root, () => {
      const unit = this.activityColumnName || 'activity';
      const mark = (scheme: string): string => scheme === this.rankScheme ? '▸ ' : '   ';
      return ui.divV([
        ui.divText(`${mark(SarRankScheme.Potency)}Potent compounds — the single most potent ` +
          `compound in the matrix (${unit}). Where the best chemistry already is.`),
        ui.divText(`${mark(SarRankScheme.Discontinuity)}SAR discontinuity — the largest activity ` +
          'spread within one core. Where a small substituent change flips potency, so the SAR is ' +
          'steep and worth reading.'),
        ui.divText(`${mark(SarRankScheme.Preferred)}Preferred substituent — the best mean potency of ` +
          'any one substituent across the cores. A substituent that pays off on every scaffold, ' +
          'rather than one lucky compound.'),
      ], 'chem-sar-rank-tip');
    });

    const parents = this.matrixParents();
    const roots = parents.filter((p) => p < 0).length;
    const deepest = this.matrices.reduce((m, matrix) => Math.max(m, matrix.level), 2) - 1;
    const sub = ui.divText(`${this.matrices.length} matrices in ${roots} famil${roots === 1 ? 'y' : 'ies'} · ` +
      `${deepest} level${deepest === 1 ? '' : 's'}`, 'chem-sar-nav-sub');
    const list = ui.div([], 'chem-sar-nav-list');
    const matchCount = ui.divText('', 'chem-sar-nav-matches');
    // The cards are built once (fillNavList below); a filter change only shows/hides them, so
    // dragging the threshold never recreates a card or re-rasterizes a structure.
    const refill = (): void => {
      this.updateNavVisibility();
      const hits = this.matrices.filter((matrix) => this.passesFilter(matrix)).length;
      const filtered = this.navPass !== null;
      matchCount.innerText = filtered ? `${hits} of ${this.matrices.length} match` : '';
    };

    // One form rather than a loose input: it aligns the label column for us, which hand-placing does
    // not, and it keeps the spacing consistent with every other panel in the platform.
    const controls = ui.form([rankInput]);

    // Every filter lives in the platform's own filter group: a substructure sketcher over the cores
    // plus a histogram per printed number, which no pair of header inputs can express.
    const navIcon = ui.icons.filter(() => {
      ui.showPopup(ui.div(this.navFilterRoot(), 'chem-sar-struct-filters'), navIcon, {vertical: true});
    }, 'Filter series by core structure, potency, SAR spread, size');
    navIcon.classList.add('chem-sar-struct-icon');
    this.navMatchCount = refill;

    const header = ui.divV([sub, ui.divH([controls, navIcon], 'chem-sar-nav-controls'), matchCount],
      'chem-sar-nav-header');
    this.fillNavList(list, parents);
    refill(); // initial visibility + match count (a threshold can survive a navigator rebuild)
    return ui.divV([header, list], 'chem-sar-nav');
  }

  /** The navigator's child lists, recomputed from the parent chain the cards were built against.
   *  A coarser matrix owns the finer ones folded into it, so each family is one branch of the tree;
   *  ranked order decides the order of the families and of the siblings within each. */
  private navChildren(parents: number[]): Map<number, number[]> {
    const children = new Map<number, number[]>();
    parents.forEach((p, i) => {
      if (p >= 0)
        children.set(p, [...(children.get(p) ?? []), i]);
    });
    return children;
  }

  /**
   * Build the navigator's cards, ONCE per navigator build: one card per matrix, each parent
   * immediately followed by the subtree it owns, all of them in the DOM. What is *shown* is decided
   * separately by {@link updateNavVisibility} — selection, collapse toggles and the threshold filter
   * all mutate the existing cards (a class swap / display toggles) instead of coming back here.
   * The structures themselves rasterize lazily on first visibility ({@link lazyCoreDepiction}), so
   * building every card up front costs DOM only.
   */
  private fillNavList(list: HTMLElement, parents: number[]): void {
    ui.empty(list);
    this.navCards = [];
    this.navParents = parents;
    this.navCoreObserver?.disconnect();
    this.navPendingCores.clear();
    this.navCoreObserver = new IntersectionObserver((entries) => {
      for (const entry of entries) {
        if (!entry.isIntersecting)
          continue;
        this.navCoreObserver?.unobserve(entry.target);
        const draw = this.navPendingCores.get(entry.target);
        if (draw) {
          this.navPendingCores.delete(entry.target);
          draw();
        }
      }
    }, {root: list, rootMargin: '300px'}); // pre-draw a little past the viewport so scrolling reads smooth

    const children = this.navChildren(parents);
    const descendants = (i: number): number =>
      (children.get(i) ?? []).reduce((n, child) => n + 1 + descendants(child), 0);

    // A fresh analysis opens at its roots: with a few hundred matrices the fully expanded list is
    // unreadable, and the coarsest matrix is the one worth reading first anyway. Seeded once, so a
    // redraw after a toggle keeps whatever the user has opened since.
    if (!this.collapseSeeded) {
      this.collapsed.clear();
      for (const parent of children.keys())
        this.collapsed.add(this.matrices[parent].id);
      // The selection rises to the root that holds it rather than that root being opened to reveal it,
      // or the family containing whatever ranked first would always arrive expanded. The pane is built
      // after the navigator, so it picks this up and opens on the broadest matrix of that family.
      let root = Math.min(this.selIndex, this.matrices.length - 1);
      while (root >= 0 && parents[root] >= 0)
        root = parents[root];
      if (root >= 0)
        this.selIndex = root;
      this.collapseSeeded = true;
    }

    const emit = (i: number, depth: number): void => {
      const id = this.matrices[i].id;
      const card = this.buildCard(this.matrices[i], i, depth, descendants(i), () => {
        if (!this.collapsed.delete(id))
          this.collapsed.add(id);
        this.updateNavVisibility();
      });
      this.navCards[i] = card;
      list.appendChild(card);
      for (const child of children.get(i) ?? [])
        emit(child, depth + 1);
    };
    parents.forEach((p, root) => {
      if (p < 0)
        emit(root, 0);
    });
  }

  /**
   * Show/hide the already-built cards per the collapse state and the threshold filter — pure display
   * toggles, so neither dragging the threshold nor toggling a branch ever recreates a card.
   *
   * A parent that fails the filter is still shown when something under it passes, or the match would
   * be unreachable. While a filter is on, collapse is ignored down those paths for the same reason —
   * filtering is a search, and a hit hidden inside a closed node is not a hit the user can see.
   */
  private updateNavVisibility(): void {
    const parents = this.navParents;
    const filtering = this.navPass !== null;
    const children = this.navChildren(parents);
    const hits = new Map<number, boolean>();
    const anyHit = (i: number): boolean => {
      const cached = hits.get(i);
      if (cached !== undefined)
        return cached;
      hits.set(i, true); // guards against a malformed parent chain looping back on itself
      let ok = this.passesFilter(this.matrices[i]);
      for (const child of children.get(i) ?? [])
        ok = anyHit(child) || ok;
      hits.set(i, ok);
      return ok;
    };
    const hiddenByCollapse = (i: number): boolean => {
      for (let p = parents[i]; p >= 0; p = parents[p]) {
        if (this.collapsed.has(this.matrices[p].id))
          return true;
      }
      return false;
    };
    this.navCards.forEach((card, i) => {
      if (!card)
        return;
      const visible = filtering ? anyHit(i) : !hiddenByCollapse(i);
      card.style.display = visible ? '' : 'none';
    });
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

  /**
   * Apply the potency threshold to the grid already on screen instead of rebuilding it.
   *
   * Rows go through the grid's own row filter, columns through `GridColumn.visible`, and the blanking
   * of surviving-but-failing cells is just a repaint — the painter reads `cellMin` directly. A full
   * `render()` would rebuild the navigator, re-run the containment nesting, recreate the DataFrame and
   * grid, and recompute every row's alignment template, none of which the threshold can change.
   *
   * The painter maps a grid row back through `gridRowToTable`, so a filtered row set keeps resolving
   * to the right `PaneRow` without any index bookkeeping here.
   */
  private applyCellFilter(): void {
    const state = this.matrixGrid;
    if (state === null)
      return;
    const rowPasses = (i: number): boolean => {
      const row = state.rows[i];
      return row === undefined || !this.cellFilterActive ||
        row.colIdxs.some((ci) => this.cellVisible(row.matrix, row.rowIndex, ci));
    };
    state.df.filter.init(rowPasses);
    state.colKeyToIdx.forEach((idx, key) => {
      const column = state.grid.col(key);
      if (column === null)
        return;
      column.visible = !this.cellFilterActive || state.rows.some((row, i) =>
        rowPasses(i) && this.cellVisible(row.matrix, row.rowIndex, row.colIdxs[idx]));
    });
    this.refitColumns();
    state.grid.invalidate();

    const matrix = state.rows.length ? state.rows[0].matrix : null;
    if (matrix !== null && this.paneDimsChip !== null) {
      const {rows, cols} = this.visibleDims(matrix);
      this.paneDimsChip.innerText = `${rows}×${cols}`;
      const full = `${matrix.rows.length} cores × ${matrix.columns.length} substituents`;
      ui.tooltip.bind(this.paneDimsChip, () =>
        rows === matrix.rows.length && cols === matrix.columns.length ? full :
          `${rows} cores × ${cols} substituents shown, filtered from ${full}`);
    }
    const emptied = state.df.filter.trueCount === 0 || this.visibleGridCols(state) === 0;
    if (this.paneGridHost !== null)
      this.paneGridHost.style.display = emptied ? 'none' : '';
    if (this.paneEmptyNote !== null)
      this.paneEmptyNote.style.display = emptied ? '' : 'none';
  }

  /** How many of the pane's columns are currently shown — hidden ones must not claim width. */
  private visibleGridCols(state: MatrixGridState): number {
    let n = 0;
    state.colKeyToIdx.forEach((_idx, key) => {
      if (state.grid.col(key)?.visible !== false)
        n++;
    });
    return n;
  }

  /**
   * Re-fit the live grid's columns to the current pane width. Cell geometry reaches the painters as
   * the grid's own cell bounds rather than through `cellW`, so a width change needs no rebuild — which
   * is what lets the grid keep its scroll position across a resize. Widths are clamped to whole
   * pixels, so a drag that does not cross a pixel boundary is skipped rather than repainting the
   * viewport on every debounce tick.
   */
  private refitColumns(): void {
    const state = this.matrixGrid;
    if (state === null)
      return;
    const cellW = this.fitCellWidth(this.visibleGridCols(state));
    if (cellW === this.cellW)
      return;
    this.cellW = cellW;
    state.colKeyToIdx.forEach((_i, key) => {
      const column = state.grid.col(key);
      if (column !== null)
        column.width = cellW;
    });
  }

  /** Column indices to show: all, or only the "Vary" position's group when one is chosen. Columns
   *  keep their as-assembled (frequency) order — the "Label" control annotates them, never reorders. */
  private visibleColIdxs(matrix: SarMatrix): number[] {
    const all = matrix.columns.map((_, ci) => ci);
    // The potency threshold is NOT applied here: it hides columns through `GridColumn.visible` on the
    // grid that is already up, so the pane can be built once and filtered without a rebuild.
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
    if (this.columnCaption === COLSORT_MW) {
      const mw = substituentMW(matrix.columns[colIdx].substSmiles);
      return Number.isFinite(mw) ? `MW ${mw.toFixed(0)}` : '';
    }
    if (this.columnCaption === COLSORT_POTENCY) {
      const mean = this.columnMeanPotency(matrix, colIdx);
      return mean === null ? '' : `μ ${this.formatActivity(mean)}`;
    }
    return '';
  }

  /**
   * Virtualized render of a row/column spec: a scaffold DataFrame (one row per displayed row, one
   * string column per displayed column plus a pinned 'Core' column) backs a DG.Grid whose every cell
   * — body, core, and R-group header — is hand-painted in `onCellRender`. Only viewport cells draw,
   * so a large matrix no longer renders every core×substituent up front, and selection/current-row
   * changes repaint via `invalidate` instead of rebuilding the DOM. Both the matrix pane and the
   * transfer pane build through here; each displayed row resolves against its own matrix, which is
   * what lets a cross-series transfer show two rows drawn from two different matrices.
   */
  private buildPaneGrid(rows: PaneRow[], columns: PaneColumn[], slot: PaneGridSlot): HTMLElement {
    this.cellW = this.fitCellWidth(columns.length);
    // The group boundaries are computed here, once per grid: `onCellRender` runs for every visible
    // header cell on every repaint, and rescanning the columns there is quadratic in the column count
    // for a result that cannot change until the grid is rebuilt.
    const firstOfGroup = new Set<number>();
    columns.forEach((column, i) => {
      if (i === 0 || columns[i - 1].position !== column.position)
        firstOfGroup.add(i);
    });

    const df = DG.DataFrame.create(rows.length);
    df.columns.addNewString('Core');
    const colKeyToIdx = new Map<string, number>();
    // Stable string keys (never the grid column idx, which pinning and the hidden row header shift).
    columns.forEach((_column, i) => {
      const key = `c${i}`;
      df.columns.addNewString(key);
      colKeyToIdx.set(key, i);
    });

    const grid = DG.Grid.create(df);
    // A fixed header height fits the position band + R-group depiction + sort caption; the built-in
    // row-number column is hidden because the cores live in the pinned 'Core' column.
    grid.setOptions({colHeaderHeight: COL_HEADER_H, rowHeight: CELL_H, showRowHeader: false,
      currentRowColor: DG.Color.white});
    grid.col('Core')!.width = CORE_W;
    colKeyToIdx.forEach((_i, key) => grid.col(key)!.width = this.cellW);
    grid.col('Core')!.pin();
    // Resolved once per grid: `onCellRender` runs for the pinned header on every repaint, and the
    // answer cannot change until the grid is rebuilt. A cross-series transfer draws its two rows from
    // different matrices, so there is no one series structure to name and the header stays textual.
    // A single template only serves rows that genuinely share a core. A cluster built on a generic MCS
    // anchor gives each row a different concrete core, and aligning those to one another's template
    // silently fails — `generate_aligned_coords` leaves the molecule with its own layout — so every
    // such row would be drawn in its own orientation. Those panes align per row instead.
    const firstCore = rows.length > 0 ? rows[0].matrix.rows[rows[0].rowIndex]?.coreSmiles ?? null : null;
    const sharedCore = firstCore !== null &&
      rows.every((row) => row.matrix.rows[row.rowIndex]?.coreSmiles === firstCore) ? firstCore : null;
    // Attachment points are capped off the header. It names the scaffold every row shares, while an
    // open position only means something on a row key, where exactly one is left open; carrying both
    // here labels a site the rows have already filled and invites reading them as disagreeing with the
    // column header. Capping leaves the atoms implicit, so nothing is drawn in their place.
    const headerCore = sharedCore === null ? null : sharedCore.replace(/\[\*:\d+\]/g, '[H]');
    const paneTemplate = headerCore !== null ? buildAlignmentTemplate(headerCore) : null;
    const state: MatrixGridState = {grid, df, rows, columns, colKeyToIdx, firstOfGroup, headerCore, paneTemplate,
      rowTemplates: new Array(rows.length).fill(null)};

    // Owned by this grid instance; unsubscribed when the next grid replaces it (or on detach) so
    // repeated renders can't leak render/click/tooltip handlers into detached grids.
    slot.subs.forEach((s) => s.unsubscribe());
    slot.subs = [];

    slot.subs.push(grid.onCellRender.subscribe((args) => {
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
          if (ri >= 0 && ri < rows.length) {
            if (name === 'Core')
              this.paintCore(g, b, rows[ri], ri, state);
            else {
              const idx = colKeyToIdx.get(name);
              if (idx !== undefined)
                this.paintBodyCell(g, b, rows[ri], ri, idx, state);
            }
          }
        }
      } finally {
        g.restore();
      }
      args.preventDefault();
    }));

    slot.subs.push(grid.onCellClick.subscribe((c) => this.onGridCellClick(state, c)));
    slot.subs.push(grid.onCellTooltip((c, x, y) => this.onGridCellTooltip(state, c, x, y)));

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
      const resolved = this.resolveGridCell(state, hit);
      if (!resolved)
        return;
      const {paneRow, ci} = resolved;
      const cell = paneRow.matrix.cells[paneRow.rowIndex][ci];
      if (cell.kind === 'virtual' && cell.smiles)
        this.contextCell = {matrix: paneRow.matrix, ri: paneRow.rowIndex, ci};
    };
    overlay.addEventListener('mousedown', onMouseDown);
    overlay.addEventListener('contextmenu', onContextMenu);
    slot.subs.push({unsubscribe: () => overlay.removeEventListener('mousedown', onMouseDown)});
    slot.subs.push({unsubscribe: () => overlay.removeEventListener('contextmenu', onContextMenu)});

    slot.state = state;

    // The grid virtualizes off its own height, so give it a plain flex host that fills the pane (not
    // ui.box, which would pin a fixed pixel width and stop the matrix growing with the pane).
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    return ui.div([grid.root], 'chem-sar-grid-host');
  }

  /** Resolve a grid body cell to the displayed row descriptor and the index of the matrix column
   *  behind it. Null for the pinned 'Core' column and for anything outside the built row/column set. */
  private resolveGridCell(state: MatrixGridState, c: DG.GridCell): {paneRow: PaneRow, ci: number} | null {
    const name = c.gridColumn.name;
    if (name === 'Core')
      return null;
    const idx = state.colKeyToIdx.get(name);
    if (idx === undefined)
      return null;
    const ri = state.grid.gridRowToTable(c.gridRow);
    if (ri < 0 || ri >= state.rows.length)
      return null;
    const paneRow = state.rows[ri];
    return {paneRow, ci: paneRow.colIdxs[idx]};
  }

  /** Click on a grid body cell: open the platform Molecule context, set the host grid's current row /
   *  selection, and register the cell's SAR panel. The row descriptor carries the matrix, so a
   *  cross-series transfer opens the panel of whichever matrix the clicked side belongs to. */
  private onGridCellClick(state: MatrixGridState, c: DG.GridCell): void {
    if (!c.isTableCell)
      return;
    const resolved = this.resolveGridCell(state, c);
    if (!resolved)
      return;
    const {paneRow, ci} = resolved;
    const matrix = paneRow.matrix;
    const ri = paneRow.rowIndex;
    const cell = matrix.cells[ri][ci];
    // A cell the threshold has blanked is not selectable: it draws as empty, so opening the compound
    // behind it would answer a click the user never appeared to make, and the context panel would
    // describe a molecule that is not on screen.
    if (cell.kind === 'empty' || cell.value === null || !this.cellVisible(matrix, ri, ci))
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
    state.grid.invalidate();
  }

  /** Hover text for a grid cell: substituent/position for an R-group header, potency detail for a body
   *  cell. Returns true to suppress the default tooltip. */
  private onGridCellTooltip(state: MatrixGridState, c: DG.GridCell, x: number, y: number): boolean {
    if (c.isColHeader) {
      const name = c.gridColumn.name;
      if (name === 'Core')
        return true;
      const idx = state.colKeyToIdx.get(name);
      if (idx === undefined)
        return false;
      const column = state.columns[idx];
      // Which positions sit at their reference value is a property of the matrix the columns were
      // taken from — the first displayed row's, since the columns follow that row's position.
      const matrix = state.rows.length ? state.rows[0].matrix : null;
      const otherRefs = matrix ? matrix.positions.filter((p) => p !== column.position) : [];
      const parts = [`${column.position}: ${column.substSmiles}`];
      if (otherRefs.length)
        parts.push(`${otherRefs.join(', ')} at ref`);
      ui.tooltip.show(ui.divText(parts.join(' · ')), x, y);
      return true;
    }
    if (!c.isTableCell)
      return false;
    const resolved = this.resolveGridCell(state, c);
    if (!resolved)
      return false;
    const {paneRow, ci} = resolved;
    const cell = paneRow.matrix.cells[paneRow.rowIndex][ci];
    // Blanked by the threshold reads as empty, so it must hover as empty too — otherwise the value the
    // filter just hid comes straight back under the cursor.
    if (cell.kind === 'empty' || cell.value === null || !this.cellVisible(paneRow.matrix, paneRow.rowIndex, ci))
      return false;
    ui.tooltip.show(ui.divText(this.cellTooltipText(cell)), x, y);
    return true;
  }

  /**
   * The id chip's text for an observed cell, or null when there is nothing to label it with.
   *
   * Virtual analogs are skipped on purpose: a prediction has no row in the source table, so there is
   * no id to show and captioning it with anything would invent an identity for a compound nobody made.
   */
  private cellIdText(cell: SarMatrixCell): string | null {
    if (!this.idColumnName || cell.kind !== 'real' || cell.molIdx === null)
      return null;
    const column = this.dataFrame.col(this.idColumnName);
    if (column === null || column.isNone(cell.molIdx))
      return null;
    return String(column.get(cell.molIdx));
  }

  /** Cell hover text: predicted (+ support) for a virtual analog, observed for a real one. */
  private cellTooltipText(cell: SarMatrixCell): string {
    const value = cell.value!;
    const support = cell.support ?? 0;
    if (cell.kind === 'virtual')
      return `Predicted ${this.formatActivity(value)} · support n=${support}${support <= 1 ? ' (low)' : ''}`;
    return `Observed ${this.formatActivity(value)}`;
  }

  /** Paint an R-group column header (position band + straightened substituent depiction + sort
   *  caption) or the 'Core' column's "Aligned core" label. The header band is pinned on top natively. */
  private paintHeader(g: CanvasRenderingContext2D, b: DG.Rect, name: string, state: MatrixGridState): void {
    const grey6 = cssColor(this.root, '--grey-6', '#4a4a4a');
    const grey5 = cssColor(this.root, '--grey-5', '#7d7d7d');
    if (name === 'Core') {
      // Every header cell paints its own background because the grid's default paint is suppressed;
      // without this fill the pinned header keeps whatever pixels the previous frame happened to leave
      // there, so it changes shade as the matrix scrolls.
      g.fillStyle = DG.Color.toHtml(HEADER_ARGB);
      g.fillRect(b.x, b.y, b.width, b.height);
      // The series structure rides in the header rather than only in the navigator: the header is
      // pinned both ways, so it names what the rows are variations of even once row 1 has scrolled off.
      const captionH = 14;
      if (state.headerCore !== null) {
        // Aligned to its own template, the way each row's core is, so the header points the same way
        // as the column beneath it.
        drawDepiction(g, b.x, b.y, b.width, Math.max(1, b.height - captionH), state.headerCore,
          state.paneTemplate, HEADER_ARGB);
      }
      g.fillStyle = grey5;
      g.font = `italic 11px ${GRID_FONT}`;
      g.textAlign = 'center';
      g.textBaseline = state.headerCore !== null ? 'top' : 'middle';
      const captionY = state.headerCore !== null ? b.y + b.height - captionH : b.y + b.height / 2;
      g.fillText('Aligned core', b.x + b.width / 2, captionY, b.width - 6);
      return;
    }
    const idx = state.colKeyToIdx.get(name);
    if (idx === undefined)
      return;
    const column = state.columns[idx];
    g.fillStyle = DG.Color.toHtml(HEADER_ARGB);
    g.fillRect(b.x, b.y, b.width, b.height);

    // Draw the position label only on the first column of a group — per-cell clipping can't paint a
    // true cross-column spanner, so this reads as the group header (the "at ref" detail is on hover).
    const posBandH = 16;
    if (state.firstOfGroup.has(idx)) {
      g.fillStyle = grey6;
      g.font = `600 11px ${GRID_FONT}`;
      g.textAlign = 'left';
      g.textBaseline = 'top';
      // Plain "R": a matrix varies exactly one position, so the number names nothing the reader can
      // act on, and it reads as a claim about which decomposition position this is — arbitrary on a
      // symmetric core, where R1 and R2 are interchangeable.
      g.fillText('R', b.x + 5, b.y + 2, b.width - 10);
    }

    const caption = column.caption;
    const captionBandH = caption ? 14 : 0;
    const depH = Math.max(1, b.height - posBandH - captionBandH);
    const depW = Math.min(b.width, HEADER_W * 2);
    drawDepiction(g, b.x + (b.width - depW) / 2, b.y + posBandH, depW, depH,
      column.substSmiles, null, HEADER_ARGB);

    if (caption) {
      g.fillStyle = grey5;
      g.font = `10px ${GRID_FONT}`;
      g.textAlign = 'center';
      g.textBaseline = 'top';
      g.fillText(caption, b.x + b.width / 2, b.y + posBandH + depH, b.width);
    }
  }

  /**
   * The molblock a row's cells align to: the row's own key in the orientation it is actually drawn
   * in. The key is first aligned to the pane's shared core so rows agree with the header wherever
   * they can, and that result becomes the row's template — deriving the cells' template from a fresh
   * layout of the key instead lets a cell's scaffold sit differently from the core printed beside it.
   */
  private rowTemplate(state: MatrixGridState, rowIdx: number): string | null {
    const cached = state.rowTemplates[rowIdx];
    if (cached !== null)
      return cached;
    const paneRow = state.rows[rowIdx];
    const key = paneRow.matrix.rows[paneRow.rowIndex].keySmiles;
    const aligned = state.paneTemplate !== null ? alignToTemplate(key, state.paneTemplate) : '';
    // A key whose core differs from the shared one cannot align to it, and alignToTemplate hands back
    // the input untouched, so lay that key out on its own instead.
    const template = aligned.includes('V2000') ? aligned : buildAlignmentTemplate(key);
    state.rowTemplates[rowIdx] = template;
    return template;
  }

  /** Paint the pinned-left core cell: the core aligned to its own template + the row label beneath. */
  private paintCore(g: CanvasRenderingContext2D, b: DG.Rect, paneRow: PaneRow, rowIdx: number,
    state: MatrixGridState): void {
    const row = paneRow.matrix.rows[paneRow.rowIndex];
    const template = this.rowTemplate(state, rowIdx);
    const labelH = 16;
    const molH = Math.max(1, b.height - labelH);
    // The default cell paint is suppressed, so the cell fills its own background; it also has to be
    // the colour the depiction bakes in, or the core's anti-aliased bonds carry a fringe of another.
    g.fillStyle = DG.Color.toHtml(WHITE_ARGB);
    g.fillRect(b.x, b.y, b.width, b.height);
    drawDepiction(g, b.x, b.y, b.width, molH, row.keySmiles, template, WHITE_ARGB);
    g.fillStyle = cssColor(this.root, '--grey-6', '#4a4a4a');
    g.font = `600 11px ${GRID_FONT}`;
    g.textAlign = 'center';
    g.textBaseline = 'top';
    g.fillText(paneRow.label, b.x + b.width / 2, b.y + molH + 1, b.width);
  }

  /** Paint one core×substituent cell: potency tint (over white), aligned assembled molecule, value
   *  chip, virtual/deviant markers, and the host-grid selection/current ring. The tint and the
   *  additive-fit check use the row's OWN matrix, so two transfer sides drawn from different matrices
   *  each color against their own activity range. */
  private paintBodyCell(g: CanvasRenderingContext2D, b: DG.Rect, paneRow: PaneRow, rowIdx: number,
    colIndex: number, state: MatrixGridState): void {
    const matrix = paneRow.matrix;
    const ri = paneRow.rowIndex;
    const cell = matrix.cells[ri][paneRow.colIdxs[colIndex]];
    // A cell below the threshold is blanked rather than removed: its row and column are still there
    // because something else in them survived, and leaving the slot empty keeps the grid readable as a
    // table instead of resequencing every neighbour.
    if (cell.kind === 'empty' || cell.value === null || !this.cellVisible(matrix, ri, paneRow.colIdxs[colIndex])) {
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
    // The tint is flattened to an opaque colour and handed to the depiction as its own background, so
    // the structure's anti-aliased bonds blend into the tint rather than into a lighter fringe. The
    // fill still runs first, for cells that have no structure to draw over it.
    const bg = this.flatTint(matrix, value, alpha);
    g.fillStyle = DG.Color.toHtml(bg);
    g.fillRect(b.x, b.y, b.width, b.height);
    if (cell.smiles !== null)
      drawDepiction(g, b.x, b.y, b.width, b.height, cell.smiles, this.rowTemplate(state, rowIdx), bg);

    // '~value' chip, top-left, fainter for a single-observation prediction and green when it is the
    // row's most potent observation (only a prediction-free comparison is worth flagging as "best").
    const isBest = paneRow.markBest === true && !isVirtual && value === this.bestRowValue(paneRow);
    this.paintChip(g, b, `${isVirtual ? '~' : ''}${this.formatActivity(value)}`, {
      corner: 'top-left',
      italic: isVirtual,
      faint: isVirtual && support <= 1,
      color: isBest ? cssColor(this.root, '--green-2', '#1a8a3a') :
        cssColor(this.root, isVirtual ? '--grey-5' : '--grey-6', isVirtual ? '#7d7d7d' : '#4a4a4a'),
    });

    // Id chip, diagonally opposite the potency it belongs to.
    const idText = this.cellIdText(cell);
    if (idText !== null)
      this.paintChip(g, b, idText, {corner: 'bottom-right', color: cssColor(this.root, '--grey-5', '#7d7d7d')});

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

    // Host-grid link: a selected row's cell gets a blue ring.
    if (cell.molIdx !== null && this.dataFrame.selection.get(cell.molIdx)) {
      g.lineWidth = 2;
      g.strokeStyle = cssColor(this.root, '--blue-1', '#2083d5');
      g.strokeRect(b.x + 1, b.y + 1, b.width - 2, b.height - 2);
    }
  }

  /**
   * Draw one corner chip: a translucent white plate carrying a line of text. The potency value and the
   * id share it so the two read as the same kind of annotation over the structure rather than as two
   * unrelated overlays, and the plate is what keeps either legible where a bond runs underneath.
   *
   * The text is ellipsized to the cell instead of being handed to `fillText` as a max width, which
   * condenses the glyphs rather than truncating them — at 10px that turns a long id into a smear.
   */
  private paintChip(g: CanvasRenderingContext2D, b: DG.Rect, text: string, style: ChipStyle): void {
    g.save();
    g.font = `${style.italic === true ? 'italic ' : ''}600 10px ${GRID_FONT}`;
    g.textAlign = 'left';
    g.textBaseline = 'top';
    const trimmed = TextUtils.trimText(text, g, b.width - (CHIP_MARGIN + CHIP_PAD) * 2);
    const chipW = g.measureText(trimmed).width + CHIP_PAD * 2;
    // The far corner is measured from x/y + size rather than through `right`/`bottom`: the grid hands
    // the painter a bare rect whose only own members are x, y, width and height, so those accessors
    // come back undefined and every coordinate derived from them silently becomes NaN.
    const atEnd = style.corner === 'bottom-right';
    const x = atEnd ? b.x + b.width - CHIP_MARGIN - chipW : b.x + CHIP_MARGIN;
    const y = atEnd ? b.y + b.height - CHIP_MARGIN - CHIP_H : b.y + CHIP_MARGIN;
    g.globalAlpha = style.faint === true ? 0.65 : 1;
    g.fillStyle = 'rgba(255, 255, 255, 0.82)';
    g.fillRect(x, y, chipW, CHIP_H);
    g.fillStyle = style.color;
    g.fillText(trimmed, x + CHIP_PAD, y + 2);
    g.restore();
  }

  /** Format an activity value: no decimals at ≥100, one decimal below. */
  private formatActivity(value: number): string {
    return value >= 100 ? value.toFixed(0) : value.toFixed(1);
  }

  /** A context-panel section header with a small provenance badge (mirrors the mockup). */
  private cpSection(title: string, badge?: string): HTMLElement {
    const parts = [ui.divText(title, 'chem-sar-cp-section-title')];
    // Omitted rather than rendered empty: the pill carries its own background and padding, so a blank
    // one still shows as a small coloured box beside the heading.
    if (badge)
      parts.push(ui.divText(badge, 'chem-sar-cp-badge'));
    return ui.divH(parts, 'chem-sar-cp-section');
  }

  /** A label / value row in the context panel. */
  private cpRow(label: string, value: HTMLElement | string): HTMLElement {
    const v = typeof value === 'string' ? ui.divText(value, 'chem-sar-cp-rv') : value;
    return ui.divH([ui.divText(label, 'chem-sar-cp-rl'), v], 'chem-sar-cp-row');
  }

  /**
   * The cell's structure, drawn to the width the panel actually has. A molecule is rasterized at a
   * fixed pixel size, so a hardcoded box leaves the structure stranded at one size while the panel the
   * user drags grows around it. Redrawn on resize, but only past a few pixels of change — each redraw
   * is an RDKit render, and a drag emits a size event per frame.
   */
  private cpStructure(smiles: string): HTMLElement {
    const box = ui.div([], 'chem-sar-cp-structbox');
    let drawnAt = 0;
    const draw = (avail: number): void => {
      const width = Math.min(CP_STRUCT_MAX_W, Math.max(CP_STRUCT_MIN_W, avail));
      if (Math.abs(width - drawnAt) < 8)
        return;
      drawnAt = width;
      ui.empty(box);
      box.appendChild(renderMolecule(smiles,
        {width, height: Math.round(width * CP_STRUCT_ASPECT), popupMenu: false}));
    };
    // A newer panel may claim the single observer slot while this one is still waiting to be laid out.
    const token = ++this.cpStructureToken;
    this.cpStructureSub?.unsubscribe();
    // Wiring up only once the box is actually in the DOM: measured before that, clientWidth is 0 and
    // the first draw lands on the minimum, leaving a too-small structure until something resizes.
    ui.tools.waitForElementInDom(box).then(() => {
      if (token !== this.cpStructureToken)
        return;
      // The panel container is what a drag resizes, and nothing drawn inside can widen it, so it is
      // the authority on how much room there is. The box is measured too and the smaller wins: on its
      // own the box can report the width it was last given rather than the width now available, which
      // is what lets a structure grow with the panel but never shrink back.
      const panel = box.closest('.panel-content') as HTMLElement | null;
      const fit = (): void => draw(panel === null ? Math.floor(box.clientWidth) - BOX_CHROME :
        Math.min(Math.floor(box.clientWidth) - BOX_CHROME, Math.floor(panel.clientWidth) - PANEL_CHROME));
      this.cpStructureSub?.unsubscribe();
      this.cpStructureSub = DG.debounce(ui.onSizeChanged(panel ?? box), 50).subscribe(fit);
      fit();
    });
    return box;
  }

  /** A measured analog that fed a prediction: its structure over its observed value. */

  /**
   * Where a predicted value came from. The additive model is `rowMean + columnMean - grandMean`, so the
   * compounds that determined it are exactly the measured cells sharing this row and this column —
   * listed here with their values, and with the arithmetic spelled out, so the number can be checked
   * rather than trusted.
   */
  private cpPrediction(matrix: SarMatrix, rowIdx: number, colIdx: number): HTMLElement {
    const sameCore: SarMatrixCell[] = [];
    const sameSubstituent: SarMatrixCell[] = [];
    let grandSum = 0;
    let grandN = 0;
    for (let ri = 0; ri < matrix.rows.length; ri++) {
      for (let ci = 0; ci < matrix.columns.length; ci++) {
        const c = matrix.cells[ri][ci];
        if (c.kind !== 'real' || c.value === null)
          continue;
        grandSum += c.value;
        grandN++;
        if (ri === rowIdx)
          sameCore.push(c);
        if (ci === colIdx)
          sameSubstituent.push(c);
      }
    }
    const meanOf = (cells: SarMatrixCell[]): number =>
      cells.length === 0 ? 0 : cells.reduce((s, c) => s + (c.value ?? 0), 0) / cells.length;
    const rowMean = meanOf(sameCore);
    const colMean = meanOf(sameSubstituent);
    const grandMean = grandN === 0 ? 0 : grandSum / grandN;

    // Shown as deviations from the matrix mean, which is what the model actually adds up. Read as raw
    // means it looks wrong whenever both effects are negative — the prediction then falls below every
    // number on screen, because the two deficits compound.
    const coreEffect = rowMean - grandMean;
    const substEffect = colMean - grandMean;
    const signed = (v: number): string => `${v < 0 ? '−' : '+'}${this.formatActivity(Math.abs(v))}`;

    const block = ui.divV([]);
    block.appendChild(this.cpRow(`Matrix mean (n = ${grandN})`, this.formatActivity(grandMean)));
    block.appendChild(this.cpRow(`Core effect (n = ${sameCore.length})`,
      `${signed(coreEffect)}  (mean ${this.formatActivity(rowMean)})`));
    block.appendChild(this.cpRow(`Substituent effect (n = ${sameSubstituent.length})`,
      `${signed(substEffect)}  (mean ${this.formatActivity(colMean)})`));
    block.appendChild(this.cpRow('Sum', `${this.formatActivity(grandMean)} ` +
      `${signed(coreEffect)} ${signed(substEffect)} = ${this.formatActivity(grandMean + coreEffect + substEffect)}`));

    block.appendChild(this.cpReferences(matrix, rowIdx, colIdx));
    return block;
  }

  /**
   * The observed compounds sharing this cell's core and its substituent — the ones a prediction here
   * rests on, and for an existing compound the ones it should be read against. Shown for both kinds of
   * cell: the row and column neighbours are the SAR context whether or not the cell itself was made.
   *
   * The selected cell is left out of its own lists, and the matrix-pane potency threshold applies, so a
   * compound blanked out of the grid does not reappear here. Where the threshold hides some, the badge
   * reads "shown of total" rather than quietly disagreeing with the counts in the arithmetic above.
   */
  private cpReferences(matrix: SarMatrix, rowIdx: number, colIdx: number): HTMLElement {
    const block = ui.divV([]);
    const self = matrix.cells[rowIdx][colIdx];
    const collect = (pick: (ri: number, ci: number) => boolean): {cell: SarMatrixCell, ri: number, ci: number}[] => {
      const found: {cell: SarMatrixCell, ri: number, ci: number}[] = [];
      for (let ri = 0; ri < matrix.rows.length; ri++) {
        for (let ci = 0; ci < matrix.columns.length; ci++) {
          const c = matrix.cells[ri][ci];
          if (c !== self && c.kind === 'real' && c.value !== null && pick(ri, ci))
            found.push({cell: c, ri, ci});
        }
      }
      return found;
    };
    const section = (title: string, all: {cell: SarMatrixCell, ri: number, ci: number}[]): void => {
      const shown = all.filter((e) => this.cellVisible(matrix, e.ri, e.ci));
      if (shown.length === 0)
        return;
      block.appendChild(this.cpSection(title,
        shown.length === all.length ? `${all.length}` : `${shown.length} of ${all.length}`));
      block.appendChild(ui.divH(shown.map((e) => this.cpFragment(e.cell.smiles, e.cell.value ?? 0)),
        'chem-sar-cp-decomp'));
    };
    section('Measured with this core', collect((ri) => ri === rowIdx));
    section('Measured with this substituent', collect((_ri, ci) => ci === colIdx));
    return block;
  }

  /** A small framed fragment (core or substituent) with a caption. */
  /** A framed fragment tile: the structure, and under it the compound's value when one is given. The
   *  structure is omitted for an empty SMILES so a valueless tile is never an empty frame. */
  private cpFragment(smiles: string | null, value?: number): HTMLElement {
    const parts: HTMLElement[] = [];
    if (smiles)
      parts.push(ui.div([renderMolecule(smiles, {width: 78, height: 52, popupMenu: false})], 'chem-sar-cp-frag-box'));
    if (value !== undefined)
      parts.push(ui.divText(this.formatActivity(value), 'chem-sar-cp-rv'));
    return ui.divV(parts, 'chem-sar-cp-frag');
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
    panel.appendChild(header);

    if (cell.smiles)
      panel.appendChild(this.cpStructure(cell.smiles));

    panel.appendChild(this.cpSection(isVirtual ? 'Predicted potency' : 'Potency'));
    panel.appendChild(this.cpRow(isVirtual ? 'Predicted' : 'Observed',
      ui.divText(this.formatActivity(cell.value), 'chem-sar-cp-value')));
    if (isVirtual) {
      panel.appendChild(this.cpRow('Method', FREE_WILSON_METHOD));
      // Both arms, with the split shown: a total of 5 drawn 4 + 1 is a different prediction from one
      // drawn 2 + 3, because the thin arm is what limits it. The weaker arm is called out when it is
      // down to a single compound, which is where the estimate is really an extrapolation.
      const refs = cell.references ?? 0;
      const weakest = cell.support ?? 0;
      panel.appendChild(this.cpRow('Reference points',
        `n = ${refs}${weakest <= 1 ? ' · one arm has a single compound' : ''}`));
      panel.appendChild(this.cpPrediction(matrix, rowIdx, colIdx));
    } else
      panel.appendChild(this.cpReferences(matrix, rowIdx, colIdx));

    panel.appendChild(this.cpSection('Decomposition', 'R-group'));
    // The fragments are drawn, so they caption themselves — the core and each R-group are read from
    // the structures, not from text repeating what the picture already shows.
    const parts = [this.cpFragment(row.keySmiles)];
    matrix.positions.forEach((position) => {
      const v = position === column.position ? column.substSmiles : matrix.refValues[position];
      if (v)
        parts.push(this.cpFragment(v));
    });
    panel.appendChild(ui.divH(parts, 'chem-sar-cp-decomp'));

    if (isVirtual && cell.smiles) {
      panel.appendChild(this.cpSection('Design action'));
      panel.appendChild(ui.divText(
        'This core × substituent is not in the dataset. Add it to the make-list for synthesis triage.',
        'chem-sar-cp-hint'));
      const btn = ui.button('Add to make-list', () => this.addAnalogToMakeList(matrix, rowIdx, colIdx));
      btn.classList.add('chem-sar-cp-generate');
      panel.appendChild(btn);
    }
    return panel;
  }

  /** One displayed row as a transfer-pane grid row: its own matrix, its columns, and its side label. */
  private transferPaneRow(side: TransferSide, colIdxs: number[]): PaneRow {
    return {matrix: this.matrices[side.matrixIndex], rowIndex: side.rowIndex, colIdxs,
      label: this.sideLabel(side), markBest: true};
  }

  /** Most potent OBSERVED value among a row's displayed cells, or null when it has none. Predictions
   *  are excluded: a virtual cell is an estimate, so calling it the row's best would overstate it. */
  private bestRowValue(paneRow: PaneRow): number | null {
    const observed: number[] = [];
    for (const ci of paneRow.colIdxs) {
      const cell = paneRow.matrix.cells[paneRow.rowIndex][ci];
      if (cell.kind === 'real' && cell.value !== null)
        observed.push(cell.value);
    }
    if (!observed.length)
      return null;
    return this.higherIsBetter ? Math.max(...observed) : Math.min(...observed);
  }

  /**
   * The SAR transfer view: two cores from different matrices whose potency trends run in parallel across
   * the R-groups they have both explored, with a trend strip under each so the parallel is visible
   * directly. The two rows go through the same virtualized grid as the matrix pane, each resolving
   * against its own matrix.
   */
  private buildTransferPane(transfer: Transfer): HTMLElement {
    const from = this.sideLabel(transfer.a);
    const to = this.sideLabel(transfer.b);
    const bar = ui.divH([
      ui.divText('Cross-series SAR transfer', 'chem-sar-main-title'),
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

    // Each side keeps its own matrix (a cross-series transfer's two cores live in different ones), so
    // its cells, potency range, core alignment and context panel all resolve against that matrix.
    // Measured pairs first, then the analogs the transfer argues for: one side has the compound, the
    // other only a prediction, so the column shows the evidence and the proposal side by side.
    const rows = [
      this.transferPaneRow(transfer.a, [...transfer.aCols, ...transfer.predictedACols]),
      this.transferPaneRow(transfer.b, [...transfer.bCols, ...transfer.predictedBCols]),
    ];
    // Every column is the same R-position here, and the Label metric is a matrix-pane control, so the
    // headers carry only the position band and the depiction.
    const matchedCount = transfer.substituents.length;
    const columns: PaneColumn[] = [...transfer.substituents, ...transfer.predictedSubstituents]
      .map((substSmiles, i) => ({substSmiles, position: transfer.a.position,
        caption: i < matchedCount ? '' : 'predicted'}));
    const gridHost = this.buildPaneGrid(rows, columns, this.transferSlot);
    // Two rows only: size the grid to its content instead of stretching it, so the trend strips and
    // the statistics stay on screen rather than being pushed below a half-empty grid.
    gridHost.style.flex = '0 0 auto';
    gridHost.style.height = `${COL_HEADER_H + rows.length * CELL_H + GRID_SCROLLBAR_H}px`;

    // A plain div, not ui.box: ui.box pins an explicit pixel width, which stops the content from
    // growing with the pane.
    const scroll = ui.div([this.buildTransferStats(transfer)], 'chem-sar-main-scroll');
    const parts = controlBar ? [bar, controlBar, gridHost, scroll] : [bar, gridHost, scroll];
    return ui.divV(parts, 'chem-sar-main');
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
      row('Shared R-groups', text(`${transfer.substituents.length} · scaffolds ${transfer.similarity.toFixed(2)}`)),
    ];
    if (stats.benefiting) {
      const side = stats.benefiting.side === 'a' ? transfer.a : transfer.b;
      rows.push(row('Benefiting cell', ui.divH([
        text(`${this.sideLabel(side)} ×`),
        renderMoleculeOnColor(stats.benefiting.substSmiles, BENEFIT_MOL_W, BENEFIT_MOL_H, HEADER_ARGB),
      ], 'chem-sar-xfer-ben')));
    }
    return ui.divV([ui.divText('Transfer statistics', 'chem-sar-xfer-title'), ...rows], 'chem-sar-xfer-stats');
  }

  /** Compact summary chips for the current matrix — short labels, full detail on hover: compound
   *  count, cores×substituents, potency range, virtual count, and transfer r. */
  private buildChips(matrix: SarMatrix): HTMLElement {
    const chip = (text: string, tip: string, cls = ''): HTMLElement => {
      const el = ui.divText(text, `chem-sar-chip-badge ${cls}`.trim());
      ui.tooltip.bind(el, () => tip);
      return el;
    };
    const items = [
      chip(`${matrix.realCount} cpd`, `${matrix.realCount} observed compounds`),
      // Reports what is on screen, not what the matrix holds: with a threshold set those differ, and a
      // chip claiming 18×6 over a grid showing 4×2 is the kind of thing nobody notices is wrong.
      ...(() => {
        const {rows, cols} = this.visibleDims(matrix);
        const full = `${matrix.rows.length} cores × ${matrix.columns.length} substituents`;
        const el = chip(`${rows}×${cols}`, rows === matrix.rows.length && cols === matrix.columns.length ?
          full : `${rows} cores × ${cols} substituents shown, filtered from ${full}`);
        this.paneDimsChip = el; // kept so the threshold can retitle it without rebuilding the bar
        return [el];
      })(),
    ];
    if (matrix.realCount) {
      items.push(chip(`${this.formatActivity(matrix.minActivity)}–${this.formatActivity(matrix.maxActivity)}`,
        `Activity range across the matrix, on the ${this.scalingLabel} scale`));
    }
    if (matrix.virtualCount)
      items.push(chip(`${matrix.virtualCount} virtual`, `${matrix.virtualCount} predicted (virtual) analog(s)`));
    const idx = this.matrices.indexOf(matrix);
    // this.transfers is sorted by correlation (computeAllTransfers), so the first match is the strongest.
    const involving = this.transfers.filter((t) => t.a.matrixIndex === idx || t.b.matrixIndex === idx);
    if (involving.length) {
      const best = involving[0];
      items.push(chip(`r ${best.correlation.toFixed(2)}`,
        'Strongest SAR-transfer correlation between this series and another', 'chem-sar-chip-transfer'));
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
          this.renderMatrixPane(); // column choice is a pane concern — the navigator has no stake
        },
      });
      controls.push(varyInput.root);
    }
    const labelInput = ui.input.choice('Caption', {
      value: this.columnCaption,
      items: COLUMN_SORTS,
      onValueChanged: (value) => {
        this.columnCaption = value!;
        this.renderMatrixPane();
      },
    });
    ui.tooltip.bind(labelInput.root, () => 'Annotate each substituent column with a metric — its mean ' +
      'potency (μ) or molecular weight (MW). Columns keep their order; only the caption is added.');
    controls.push(labelInput.root);

    // Every filter lives behind one icon: the two sketchers would dwarf the bar, and splitting the
    // cell cuts away from them left two places to look for "why is this cell blank".
    const filterIcon = ui.icons.filter(() => {
      ui.showPopup(ui.div(this.structureFilterRoot(), 'chem-sar-struct-filters'),
        filterIcon, {vertical: true});
    }, 'Filter cells by potency, reference points, core and R-group');
    filterIcon.classList.add('chem-sar-struct-icon');
    controls.push(filterIcon);
    const controlBar = ui.divH(controls, 'chem-sar-control-bar');

    // Every row of this pane shows the same matrix and the same visible columns; the Vary filter
    // picks the columns and the Label control supplies their captions.
    const visible = this.visibleColIdxs(matrix);
    // Every row is built; the threshold hides them through the grid's row filter afterwards, so moving
    // the slider never rebuilds this pane.
    const rows: PaneRow[] = matrix.rows.map((row, ri) =>
      ({matrix, rowIndex: ri, colIdxs: visible, label: row.label}));
    const columns: PaneColumn[] = visible.map((ci) => ({
      substSmiles: matrix.columns[ci].substSmiles,
      position: matrix.columns[ci].position,
      caption: this.columnSortCaption(matrix, ci),
    }));

    // The grid scrolls and virtualizes internally, so it goes in a flex host that fills the pane (no
    // outer .chem-sar-main-scroll, which would add a second scroll container and its own padding).
    const gridHost = this.buildPaneGrid(rows, columns, this.matrixSlot);
    // Shown in the grid's place when the threshold empties the matrix; an empty grid reads as a broken
    // viewer. Built with the pane and toggled by the filter, so no rebuild is needed to reach it.
    const reach = matrix.realCount ?
      `its compounds run ${this.formatActivity(matrix.minActivity)}–${this.formatActivity(matrix.maxActivity)}` :
      'it has no observed compounds';
    const emptyNote = ui.divText(`No compound in ${matrix.label} is that potent — ${reach}. ` +
      'Lower the threshold, or pick another matrix.', 'chem-sar-empty-note');
    emptyNote.style.display = 'none';
    this.paneGridHost = gridHost;
    this.paneEmptyNote = emptyNote;
    this.applyCellFilter();
    return ui.divV([infoBar, controlBar, gridHost, emptyNote], 'chem-sar-main');
  }

  private render(): void {
    // Keep the navigator scrolled where it was — selecting a card lower down must not jump it to the top.
    const prevNav = this.host.querySelector('.chem-sar-nav-list');
    const navScroll = prevNav instanceof HTMLElement ? prevNav.scrollTop : 0;
    // Before the DOM goes: buildPaneGrid installs a new grid, so the previous one (and everything its
    // handlers close over) has to be let go here.
    this.releaseMatrixGrid();
    ui.empty(this.host);
    // Palette values are memoized by variable name only, so a theme switch would otherwise keep
    // painting the previous theme's colors for the life of the session.
    cssColorCache.clear();
    if (this.matrices.length === 0) {
      this.host.appendChild(ui.divText(this.grouping === SarGrouping.Site ?
        'No SAR matrices found. Try raising the fragment cutoff, or switch grouping to Similarity.' :
        'No SAR matrices found. Try lowering the clustering threshold or raising the fragment cutoff.'));
      return;
    }
    this.host.appendChild(this.buildNavigator());
    const newNav = this.host.querySelector('.chem-sar-nav-list');
    if (newNav instanceof HTMLElement)
      newNav.scrollTop = navScroll;
    // The matrix tab always shows a matrix; SAR transfer has its own tab and never displaces it.
    this.renderMatrixPane();
    // Only when it is the tab on screen: rebuilt hidden, its grid would size against a display:none
    // pane, and rebuilding a pane nobody is looking at wastes the very work the lazy tab defers.
    if (this.transferTabActive)
      this.activateTransferTab();
  }

  /** Rebuild only the matrix pane (right side), leaving the navigator DOM — and its scroll — alone.
   *  This is the whole redraw a card click, a Vary/Caption change or an id-column switch needs. */
  private renderMatrixPane(): void {
    if (this.matrices.length === 0)
      return;
    this.releaseMatrixGrid();
    this.host.querySelector(':scope > .chem-sar-main')?.remove();
    const matrix = this.matrices[Math.min(this.selIndex, this.matrices.length - 1)];
    this.host.appendChild(this.buildMatrixPane(matrix));
    this.syncSelection(); // apply the host grid's current selection/row to the freshly-built cells
  }
}
