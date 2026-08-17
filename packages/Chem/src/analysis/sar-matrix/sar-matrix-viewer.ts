/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../../css/sar-matrix.css';
import {drawMoleculeToCanvas, getRdKitModule, getRdKitService}
  from '../../utils/chem-common-rdkit';
import {getMoleculeRenderer} from '../../package';
import {renderMolecule} from '../../rendering/render-molecule';
import * as TextUtils from '@datagrok-libraries/gridext/src/utils/TextUtils';
import {SCALING_METHODS} from '../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {nestByContainment, rankMatrices, SarRankScheme} from './sar-matrix-ranking';
import {computeAllTransfers, DEFAULT_TRANSFER_SIMILARITY, Transfer, TransferSide, transferStats}
  from './sar-matrix-transfer';
import {MAX_SERIES_LEVELS, runSarMatrix, SarGrouping, SarMatrixParams} from './sar-matrix-run';
import {SarMatrix, SarMatrixCell} from './sar-matrix-types';

/** Transparent (alpha 0) so a drawn core blends with the card/pane instead of showing a white box. */
const CORE_BG_ARGB = 0x00000000;
const HEADER_ARGB = 0xFFF7F7F9;
const WHITE_ARGB = 0xFFFFFFFF;
const ACTIVITY_SCHEME = [DG.Color.red, DG.Color.green];
const CELL_W = 104;
const CELL_H = 76;
const CHIP_H = 13;
const CHIP_PAD = 3;
const CHIP_MARGIN = 3;
const HEADER_W = 78;
const HEADER_H = 46;
const COL_HEADER_H = HEADER_H + 36;
const CORE_W = 132;
/** Room for the grid's own horizontal scrollbar when a grid is sized to its rows, not stretched. */
const GRID_SCROLLBAR_H = 18;
const BENEFIT_MOL_W = 62;
const BENEFIT_MOL_H = 34;
const CARD_CORE_W = 78;
const CARD_CORE_H = 44;
const CELL_W_MAX = 210;
/** Must track the `.chem-sar-nav` width in the stylesheet, or cells are fitted against the wrong pane. */
const NAV_W = 320;
const NAV_COLLAPSED_W = 28;
const TABLE_CHROME = 60;
/** Cell-tint alpha (0-255): solid for observed, fainter for virtual, faintest for thin predictions. */
const REAL_ALPHA = 102;
const VIRTUAL_ALPHA = 46;
const VIRTUAL_ALPHA_MIN = 16;
/** Support at/above which a virtual cell is at full alpha. */
const FULL_SUPPORT = 3;
const FREE_WILSON_METHOD = 'local Free-Wilson (row + column effects)';
const MAKELIST_NAME = 'SAR virtual analogs';

type AnalogPanelBuilder = () => HTMLElement;

const COLSORT_POTENCY = 'Potency';
const COLSORT_MW = 'Molecular weight';
const COLUMN_SORTS = [COLSORT_POTENCY, COLSORT_MW];

/** Properties that only reorder/recolor already-assembled matrices — must NOT re-run fragmentation. */
const RERANK_ONLY_PROPS = ['rankScheme', 'activityDirection'];

/** Properties that only change what is drawn — no re-fragmentation and no re-ranking. */
const RENDER_ONLY_PROPS = ['columnCaption', 'idColumnName'];

/** Properties only the transfer scan reads; matrices are untouched. */
const TRANSFER_ONLY_PROPS = ['transferSimilarity'];

const TAB_MATRIX = 'SAR Matrix';
const TAB_TRANSFER = 'SAR Transfer';

/** Auto derives from scaling (only −lg is higher-is-better); explicit options cover precomputed pIC50 etc. */
const DIR_AUTO = 'Auto (from scaling)';
const DIR_HIGHER = 'Higher is better';
const DIR_LOWER = 'Lower is better';
const ACTIVITY_DIRECTIONS = [DIR_AUTO, DIR_HIGHER, DIR_LOWER];

/** Average MW of a substituent, capping `[*:n]` attachment points with H. Infinity on failure so
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

/** Aligned-molblock cache keyed by template + molecule; cleared wholesale past the cap to avoid leaks. */
const alignCache = new Map<string, string>();
const ALIGN_CACHE_MAX = 4000;
/** Per-core alignment templates by core SMILES. `null` is cached too, so failures aren't retried. */
const templateCache = new Map<string, string | null>();
const TEMPLATE_CACHE_MAX = 4000;

function argbToRgba(argb: number): [number, number, number, number] {
  return [DG.Color.r(argb) / 255, DG.Color.g(argb) / 255, DG.Color.b(argb) / 255, DG.Color.a(argb) / 255];
}

/** Canonical molblock template every cell/core aligns to, so the shared core points the same way.
 *  `[*:n]` dummies are swapped for H so the template is a plain substructure. */
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

/** Draw a molecule with the ARGB background baked into the RDKit draw call — a colour applied
 *  afterwards would leave a pale anti-aliasing fringe around every bond. */
function paintMoleculeOnColor(canvas: HTMLCanvasElement, smiles: string, w: number, h: number,
  argb: number): void {
  try {
    drawMoleculeToCanvas(0, 0, w, h, canvas, smiles, null,
      {normalizeDepiction: true, straightenDepiction: true}, null,
      {clearBackground: true, backgroundColour: argbToRgba(argb)});
  } catch (e) {
  }
}

function renderMoleculeOnColor(smiles: string, w: number, h: number, argb: number): HTMLElement {
  const canvas = ui.canvas(w, h);
  if (smiles)
    paintMoleculeOnColor(canvas as HTMLCanvasElement, smiles, w, h, argb);
  return canvas;
}

/** Drop the per-dataset core-alignment layouts on rebuild; the shared renderer's raster cache is LRU. */
export function clearDepictionCaches(): void {
  alignCache.clear();
}

/** Label every attachment point in a molblock as a plain "R" via an atom alias (RDKit draws a bare
 *  dummy as `*`). Map numbers are meaningless on a one-position matrix. Input unchanged if none. */
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

/** The molecule string a cell/core/header is drawn from: aligned to `template` and R-labelled when
 *  given, else the map-stripped string for the shared renderer. `null` if empty. */
function preparedDepiction(molStr: string, template: string | null): string | null {
  if (!molStr)
    return null;
  const hasAttachment = /\[\*:\d+\]/.test(molStr);
  const plain = molStr.replace(/\[\*:\d+\]/g, '[*]');
  const aligned = template ? alignToTemplate(plain, template) : plain;
  return aligned.includes('V2000') && hasAttachment ? labelAttachmentPoints(aligned) : aligned;
}

/** Scratch canvas the cached `ImageData` is blitted through so the grid's clip is respected
 *  (`putImageData` ignores the clip and would paint over the pinned core column or past the edge). */
let blitCanvas: OffscreenCanvas | null = null;
function ensureBlitCanvas(w: number, h: number): OffscreenCanvas {
  if (!blitCanvas || blitCanvas.width < w || blitCanvas.height < h)
    blitCanvas = new OffscreenCanvas(Math.max(w, blitCanvas?.width ?? 0), Math.max(h, blitCanvas?.height ?? 0));
  return blitCanvas;
}

/** Draw a depiction onto the grid canvas in device pixels, reusing the shared renderer's LRU caches.
 *  The bitmap is produced at device size and blitted 1:1 at a whole-pixel offset so bond lines stay
 *  crisp (drawing through the grid's scaled transform bilinear-filters every bond). The device rect
 *  comes from the context transform, not `devicePixelRatio`, so a grid translate can't displace cells.
 *  Blits via a scratch canvas, not `putImageData`, so the grid's clip is respected. */
function drawDepiction(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number,
  molStr: string, template: string | null, argb: number): void {
  const renderer = getMoleculeRenderer();
  if (!renderer)
    return;
  const molblock = preparedDepiction(molStr, template);
  if (molblock === null)
    return;
  const m = g.getTransform();
  const dw = Math.round(m.a * w);
  const dh = Math.round(m.d * h);
  if (dw < 1 || dh < 1)
    return;
  let image: ImageData;
  try {
    image = renderer.getCachedMolImageData(molblock, dw, dh, argbToRgba(argb));
  } catch {
    return; // malformed structure — leave the cell to its background
  }
  const scratch = ensureBlitCanvas(dw, dh);
  scratch.getContext('2d')!.putImageData(image, 0, 0);
  g.save();
  g.setTransform(1, 0, 0, 1, 0, 0);
  g.drawImage(scratch, 0, 0, dw, dh, Math.round(m.a * x + m.e), Math.round(m.d * y + m.f), dw, dh);
  g.restore();
}

/** Resolve a Datagrok CSS palette variable to a color string, memoized so per-cell paints don't
 *  re-run `getComputedStyle`. */
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

/** Context-panel structure sizing: drawn width is clamped to this range, height follows. */
const CP_STRUCT_MIN_W = 200;
const CP_STRUCT_MAX_W = 520;
const CP_STRUCT_ASPECT = 0.55;
/** Horizontal chrome subtracted from client width to get the width a structure can occupy. */
const BOX_CHROME = 12;
const PANEL_CHROME = 30;

const GRID_FONT = 'Roboto, "Segoe UI", sans-serif';

/** Column captions of the per-cell frame the filter runs over; also the filter widget labels. */
const STRUCT_CORE = 'Core';
const STRUCT_R = 'R';
const STRUCT_POTENCY = 'Potency';
const STRUCT_REFS = 'Reference points';
const STRUCT_MW = 'MW';
const NAV_SERIES = 'Series';
const NAV_CORE = 'Core';
const NAV_BEST = 'Best';
const NAV_MEAN = 'Mean';
const NAV_SPREAD = 'Spread';
const NAV_BEST_R = 'Best R';
const NAV_COMPOUNDS = 'Compounds';
const NAV_CORES = 'Cores';
const NAV_LEVEL = 'Level';

/** How a cell corner chip is drawn. The two corners carry different facts and sit on opposite
 *  diagonals so they never collide however wide either gets. */
interface ChipStyle {
  corner: 'top-left' | 'bottom-right';
  color: string;
  italic?: boolean;
  faint?: boolean;
}

/** Identity of a transfer's source core within its section, keyed on by the nav grouping and the
 *  pane's sibling lookup so they can't disagree about which targets a source reaches. */
function transferSourceKey(t: Transfer): string {
  return `${t.a.matrixIndex}:${t.a.rowIndex}`;
}

/** One displayed grid row. Each row carries its OWN matrix, because a transfer's two sides can come
 *  from different matrices. */
interface PaneRow {
  matrix: SarMatrix;
  rowIndex: number;
  colIdxs: number[];
  label: string;
  /** Highlight the row's most potent observed cell; set for a transfer's two sides. */
  markBest?: boolean;
}

interface PaneColumn {
  substSmiles: string;
  position: string;
  /** Metric caption under the depiction; empty when the Label control is None. */
  caption: string;
}

/** The rendered pane grid plus per-render row/column state, computed once at build time since
 *  `onCellRender` fires per visible cell per repaint. */
interface MatrixGridState {
  grid: DG.Grid;
  /** Scaffold frame backing the grid; held so the grid keeps a live reference for its lifetime. */
  df: DG.DataFrame;
  rows: PaneRow[];
  columns: PaneColumn[];
  /** Grid column NAME (never the grid column index, which pinning shifts) -> index into `columns`. */
  colKeyToIdx: Map<string, number>;
  /** Indices into `columns` that begin an R-position group — only those draw the position label. */
  firstOfGroup: Set<number>;
  /** Structure in the pinned header; null when rows span more than one matrix. */
  headerCore: string | null;
  /** One alignment template for the whole pane so an attachment point sits in the same place
   *  everywhere. Null when rows span several matrices and share no core. */
  paneTemplate: string | null;
  /** Per displayed row, the molblock its cells align to; filled on first paint and kept for the
   *  life of the grid. */
  rowTemplates: (string | null)[];
}

/** One live pane grid plus the subscriptions it owns. Matrix and transfer panels hold one each so a
 *  rebuild of one can't unsubscribe the other's painters or leak its grid. */
interface PaneGridSlot {
  state: MatrixGridState | null;
  subs: {unsubscribe(): void}[];
}

export class SarMatrixViewer extends DG.JsViewer {
  moleculesColumnName: string;
  activityColumnName: string;
  /** Optional column captioning each observed cell. */
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
  /** Correlated core pairs; computed lazily on first SAR-transfer-tab open. */
  private transfers: Transfer[] = [];
  private transfersComputed = false;
  /** Bumped when the transfer list is dropped, so a stale in-flight scan won't publish its result. */
  private transferGeneration = 0;
  /** A transfer scan is already scheduled; a second tab activation must not stack another. */
  private transfersComputing = false;
  private selIndex = 0;
  /** Index into `transfers`; separate from `selIndex` so the two tabs navigate independently. */
  private transferIndex = 0;
  /** "Vary" filter: show only this R-position's column group, or all when empty. */
  private varyPosition = '';
  columnCaption: string;
  private contextCell: {matrix: SarMatrix, ri: number, ci: number} | null = null;
  /** Per-SMILES "SAR analysis" panel builders; cleared on recompute and detach. */
  private readonly analogPanels = new Map<string, AnalogPanelBuilder>();
  private readonly host = ui.divH([], 'chem-sar-matrix');
  private readonly transferHost = ui.divH([], 'chem-sar-xfer-panel');
  private readonly tabs: DG.TabControl;
  private computing = false;
  /** A recompute was requested mid-compute; re-queued when the running one finishes. */
  private dirty = false;
  private computeTimer = 0;
  /** Set in `detach`, so an in-flight compute can't render into a closed viewer. */
  private detached = false;
  private cellW = CELL_W;
  /** Live pane grids, held so selection changes repaint via a cheap `invalidate()`. */
  private readonly matrixSlot: PaneGridSlot = {state: null, subs: []};
  private readonly transferSlot: PaneGridSlot = {state: null, subs: []};
  /** Last pointer event over the grid, so a cell click can honor ctrl/shift (onCellClick carries none). */
  private lastGridMouseEvent: MouseEvent | null = null;
  /** Context-panel structure size observer; replaced per click so they don't stack. */
  private cpStructureSub: {unsubscribe(): void} | null = null;
  /** Guards a slow layout-wait from overwriting a newer panel's observer. */
  private cpStructureToken = 0;
  /** Matrix ids whose children are folded away; kept so re-rank/resize doesn't reopen them. */
  private readonly collapsed = new Set<string>();
  /** Cards per matrix (index-aligned) and their parent chain; built once, then mutated in place. */
  private navCards: HTMLElement[] = [];
  private navParents: number[] = [];
  /** Rasterizes a card's core only when it first scrolls into view. */
  private navCoreObserver: IntersectionObserver | null = null;
  private readonly navPendingCores: Map<Element, () => void> = new Map();
  /** Cell keys passing the potency threshold, or null when unfiltered. */
  private cellPass: Set<string> | null = null;
  /** One row per series, backing the platform filter group over the navigator. */
  private navFrame: DG.DataFrame | null = null;
  private navFilters: DG.FilterGroup | null = null;
  /** Headless view owning `navFilters`; held so it can be closed (dropping the reference alone leaks
   *  the Dart-backed view and frame). */
  private navView: DG.TableView | null = null;
  private navSub: {unsubscribe(): void} | null = null;
  /** Matrix ids the navigator filter admits, or null while it admits everything. */
  private navPass: Set<string> | null = null;
  /** Frame row index -> matrix id. */
  private navKeys: string[] = [];
  /** Refreshes the navigator match count from outside its render. */
  private navMatchCount: (() => void) | null = null;
  /** Frame row index -> cell key. */
  private cellKeys: string[] = [];
  /** One row per (core, substituent) cell, backing the platform filter group over the matrix. */
  private structFrame: DG.DataFrame | null = null;
  private structFilters: DG.FilterGroup | null = null;
  /** Headless view owning `structFilters` (see {@link navView}). */
  private structView: DG.TableView | null = null;
  private structSub: {unsubscribe(): void} | null = null;
  /** Matrix-pane parts the potency threshold updates in place; null while no pane is on screen. */
  private paneGridHost: HTMLElement | null = null;
  private paneEmptyNote: HTMLElement | null = null;
  private paneDimsChip: HTMLElement | null = null;
  /** The matrix pane's filter icon, so the cell filter can flag it as active. */
  private matrixFilterIcon: HTMLElement | null = null;
  /** Whether the navigator is collapsed; on the viewer so it survives a navigator rebuild. */
  private navCollapsed = false;
  /** Guards the one-time collapse-to-roots per analysis. */
  private collapseSeeded = false;
  get helpUrl() {
    return 'https://raw.githubusercontent.com/datagrok-ai/public/refs/heads/master/help/datagrok/solutions/domains/chem/chem.md#sar-matrix';
  }
  constructor() {
    super();
    // COLUMN-typed properties so the picker filters to the right columns; the stored value is the name.
    this.moleculesColumnName = this.addProperty('moleculesColumnName', DG.TYPE.COLUMN, '',
      {semType: DG.SEMTYPE.MOLECULE, category: 'Data', description: 'Structures to analyze'});
    this.activityColumnName = this.addProperty('activityColumnName', DG.TYPE.COLUMN, '',
      {columnTypeFilter: DG.TYPE.NUMERICAL, category: 'Data',
        description: 'Numeric activity or potency the matrices are coloured and ranked by'});
    this.idColumnName = this.addProperty('idColumnName', DG.TYPE.COLUMN, '',
      {nullable: true, category: 'Data',
        description: 'Optional column labelling each measured cell (e.g. compound id)'});
    this.scaling = this.string('scaling', SCALING_METHODS.MINUS_LG, {choices: Object.values(SCALING_METHODS),
      friendlyName: 'Scaling',
      description: 'Activity transform before the additive model: none, log (lg) or −log (-lg)'});
    this.activityDirection = this.string('activityDirection', DIR_AUTO, {choices: ACTIVITY_DIRECTIONS,
      friendlyName: 'Activity direction', description: 'Which end of the scaled activity is more potent'});
    this.fragmentCutoff = this.float('fragmentCutoff', 0.4, {min: 0.1, max: 1, friendlyName: 'Fragment cutoff',
      description: 'Largest substituent kept when fragmenting, as a fraction of the core'});
    this.grouping = this.string('grouping', SarGrouping.Site, {choices: Object.values(SarGrouping),
      friendlyName: 'Grouping',
      description: 'How related cores are gathered into one matrix: Site (shared cut) or Similarity (fingerprint)'});
    this.fragmentationLevels = this.int('fragmentationLevels', 3,
      {min: 1, max: MAX_SERIES_LEVELS, friendlyName: 'Series levels',
        description: 'Nested matrix tiers (L1, L2, …); each level folds matrices one cut broader'});
    this.threshold = this.float('threshold', 0.5, {min: 0, max: 1, friendlyName: 'Similarity threshold',
      description: 'Core-similarity cutoff for Similarity grouping (higher = tighter clusters); ignored for Site'});
    this.transferSimilarity = this.float('transferSimilarity', DEFAULT_TRANSFER_SIMILARITY,
      {min: 0, max: 1, friendlyName: 'Transfer similarity',
        description: 'Optional whole-molecule similarity floor on transfer pairs (0 = off). Transfers ' +
          'are matched by shared R-groups; raise this to restrict them to more alike scaffolds'});
    this.predictVirtual = this.bool('predictVirtual', true, {friendlyName: 'Predict virtual analogs',
      description: 'Fill unmade core × substituent cells with Free-Wilson predictions'});
    this.rankScheme = this.string('rankScheme', SarRankScheme.Potency,
      {choices: [SarRankScheme.Potency, SarRankScheme.Discontinuity, SarRankScheme.Preferred],
        friendlyName: 'Rank by', description: 'How the navigator orders the matrices'});
    this.columnCaption = this.string('columnCaption', COLSORT_POTENCY, {choices: COLUMN_SORTS});
    this.host.style.height = '100%';
    this.transferHost.style.height = '100%';

    // Transfer detection is quadratic in the total row count, so that tab computes on first open.
    this.tabs = ui.tabControl(null, false);
    const matrixPane = this.tabs.addPane(TAB_MATRIX, () => this.host);
    ui.tooltip.bind(matrixPane.header, 'Core × substituent potency matrices, one per series');
    const transferPane = this.tabs.addPane(TAB_TRANSFER, () => this.transferHost);
    ui.tooltip.bind(transferPane.header, 'Pairs of cores whose potency trends run in parallel across the ' +
      'R-groups they have both explored — detected when the tab is first opened');
    // Not on this.subs: tab switching must survive a detach/re-attach cycle.
    this.tabs.onTabChanged.subscribe(() => {
      if (this.tabs.currentPane?.name === TAB_TRANSFER)
        this.activateTransferTab();
      else {
        // A pane rebuilt while hidden sat in a display:none host; repaint now it has real dimensions.
        this.matrixSlot.state?.grid.invalidate();
      }
    });
    this.tabs.root.style.width = '100%';
    this.tabs.root.style.height = '100%';
    this.root.appendChild(this.tabs.root);
  }

  onTableAttached(): void {
    this.detached = false;
    // Re-fit columns of the on-screen grid rather than rebuild the pane, so a splitter drag keeps the
    // user's scroll position (a new grid restarts at the first row).
    this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 200).subscribe(() => {
      if (!this.computing && this.matrices.length)
        this.refitColumns();
    }));
    this.subs.push(this.onContextMenu.subscribe((menu) => this.buildContextMenu(menu)));
    // Two-way link with the host grid. Current row is not watched: nothing paints from it and a
    // repaint re-rasterizes every visible structure.
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50)
      .subscribe(() => this.syncSelection()));
    // Capture-phase reset runs before a cell's bubbling handler, so contextCell reflects only a
    // right-click that landed on a virtual cell.
    this.host.addEventListener('contextmenu', () => this.contextCell = null, true);
    // Surface a clicked analog's SAR context in its panel; scoped to this.subs so it can't inject
    // panes platform-wide after the viewer closes.
    this.subs.push(grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
      const context = acc.context;
      const smiles = context instanceof DG.SemanticValue ? String(context.value) :
        (typeof context === 'string' ? context : null);
      if (!smiles || !this.analogPanels.has(smiles) || acc.getPane('SAR analysis'))
        return;
      const build = this.analogPanels.get(smiles)!;
      acc.addPane('SAR analysis', () => build(), true, acc.panes.length ? acc.panes[0] : null);
    }));
    // Spawn RDKit workers now so their WASM init overlaps the compute debounce; a failure here is
    // re-reported by compute().
    getRdKitService().catch(() => {});
    this.scheduleCompute();
  }

  /** Base `detach` unsubscribes `this.subs`; also stop a pending compute and drop analog-panel
   *  builders so a closed viewer leaves no timer firing or stale entries injecting panes. */
  detach(): void {
    this.detached = true;
    window.clearTimeout(this.computeTimer);
    this.analogPanels.clear();
    this.cpStructureSub?.unsubscribe(); // the panel outlives the viewer; its observer must not
    this.cpStructureSub = null;
    this.releaseMatrixGrid();
    this.releaseSlot(this.transferSlot);
    this.navCoreObserver?.disconnect();
    this.navCoreObserver = null;
    this.navPendingCores.clear();
    this.resetFilters();
    super.detach();
  }

  /** The matrix pane's grid. Named separately from the transfer panel's so a filter can't reach it. */
  private get matrixGrid(): MatrixGridState | null {
    return this.matrixSlot.state;
  }

  /** Release one pane grid: unsubscribing alone isn't enough — the subscriptions keep the grid (and
   *  its DataFrame + SarMatrix) reachable, and a Dart-backed grid is only released when closed. */
  private releaseSlot(slot: PaneGridSlot): void {
    slot.subs.forEach((s) => s.unsubscribe());
    slot.subs = [];
    try {
      slot.state?.grid?.close?.();
    } catch (e) {
      // A view-less grid may not support close; dropping the reference is enough.
    }
    slot.state = null;
  }

  private releaseMatrixGrid(): void {
    this.releaseSlot(this.matrixSlot);
    // Pane is about to be replaced; holding these would let the threshold write into detached elements.
    this.paneGridHost = null;
    this.paneEmptyNote = null;
    this.paneDimsChip = null;
  }

  /** Drop the lazy filter machinery so it rebuilds against the current matrices. Stale keys
   *  (cellKeys/navKeys) partially match reused cluster ids, mis-filtering rather than failing cleanly;
   *  closing the headless views also releases the Dart-backed views + frames. */
  private resetFilters(): void {
    this.resetStructFilter();
    this.navSub?.unsubscribe();
    this.navSub = null;
    try {
      this.navView?.close();
    } catch (e) {
    }
    this.navView = null;
    this.navFrame = null;
    this.navFilters = null;
    this.navMatchCount = null;
    this.navKeys = [];
    this.navPass = null;
  }

  /** Drop only the cell (structure) filter so it rebuilds against the current matrix; its rows and
   *  keys belong to the previous matrix. */
  private resetStructFilter(): void {
    this.structSub?.unsubscribe();
    this.structSub = null;
    try {
      this.structView?.close();
    } catch (e) {
    }
    this.structView = null;
    this.structFrame = null;
    this.structFilters = null;
    this.cellKeys = [];
    this.cellPass = null;
  }

  /** Repaint both panes so a host-grid selection/current-row change shows on the rendered cells. */
  private syncSelection(): void {
    if (!this.dataFrame)
      return;
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

  /** Build a make-list molecule table from virtual cells: predicted activity, support, and full
   *  provenance so each row stands on its own. */
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
    // Some properties don't change the fragmentation, so they must not trigger a full rebuild.
    if (property !== null && RENDER_ONLY_PROPS.includes(property.name)) {
      this.renderMatrixPane();
      return;
    }
    if (property !== null && RERANK_ONLY_PROPS.includes(property.name)) {
      this.reRank();
      return;
    }
    if (property !== null && TRANSFER_ONLY_PROPS.includes(property.name)) {
      this.invalidateTransfers();
      if (this.transferTabActive)
        this.activateTransferTab();
      return;
    }
    this.scheduleCompute();
  }

  /** Re-rank the assembled matrices and redraw without re-fragmenting. */
  private reRank(): void {
    if (!this.matrices.length) {
      this.scheduleCompute();
      return;
    }
    this.matrices = rankMatrices(this.matrices, this.rankScheme as SarRankScheme, this.higherIsBetter);
    // Transfers index into `matrices` by position, so a reorder invalidates them.
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
    // Change arrived mid-compute: re-queue when the running one finishes.
    if (this.computing) {
      this.dirty = true;
      return;
    }

    this.computing = true;
    this.analogPanels.clear();
    clearDepictionCaches();
    // Release the Dart-backed grid so the failure path below doesn't leave one repainting forever.
    this.releaseMatrixGrid();
    // Transfers and filter frames key into the matrices about to be replaced; drop them now so a
    // failed compute leaves nothing pointing at rows that no longer exist. Both rebuild lazily.
    this.invalidateTransfers();
    this.resetFilters();
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
      // Closed while the workers ran; rendering now would build a grid nothing releases.
      if (this.detached)
        return;
      this.matrices = matrices;
      this.collapseSeeded = false;
      this.selIndex = 0;
      // Must clear before render(): render() re-activates the transfer tab, which defers while
      // `computing` is set and would otherwise strand on its "Building..." placeholder.
      this.computing = false;
      this.render();
    } catch (e) {
      if (this.detached)
        return;
      const message = e instanceof Error ? e.message : String(e);
      ui.empty(this.host);
      this.host.appendChild(ui.divText(`SAR Matrix failed: ${message}`));
      // No render() on this path to replace a mid-compute "Building..." placeholder.
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

  /** Whether higher scaled activity is more potent. Explicit Activity direction wins; on Auto only
   *  `-lg` is higher-is-better while raw and `lg` keep the assay's native "lower is more potent". */
  private get higherIsBetter(): boolean {
    if (this.activityDirection === DIR_HIGHER)
      return true;
    if (this.activityDirection === DIR_LOWER)
      return false;
    return this.scaling === SCALING_METHODS.MINUS_LG;
  }

  /** Drop the transfer list and tab content; the matrices it indexed are gone. Recomputes on the
   *  next transfer-tab visit. */
  private invalidateTransfers(): void {
    this.transferGeneration++;
    this.releaseSlot(this.transferSlot);
    ui.empty(this.transferHost);
    this.transfers = [];
    this.transferIndex = 0;
    this.transfersComputed = false;
  }

  /** Potency tint composited over white and returned opaque, for a depiction background: RDKit
   *  anti-aliases bonds against the draw's background, so a translucent tint leaves a pale fringe. */
  private flatTint(matrix: SarMatrix, value: number, alpha: number): number {
    const c = this.tint(matrix, value, alpha);
    const a = alpha / 255;
    const overWhite = (channel: number): number => Math.round(channel * a + 255 * (1 - a));
    return DG.Color.argb(255, overWhite(DG.Color.r(c)), overWhite(DG.Color.g(c)), overWhite(DG.Color.b(c)));
  }

  /** Potency tint at explicit alpha (0-255), green = more potent. `scaleColor`'s own alpha arg is
   *  unusable (transparent mid-range), so RGB is taken and alpha set here; mirrored for lower-is-better. */
  private tint(matrix: SarMatrix, value: number, alpha: number): number {
    // Zero-width range would divide by zero in scaleColor; paint the green end.
    if (matrix.maxActivity === matrix.minActivity)
      return DG.Color.argb(alpha, DG.Color.r(DG.Color.green), DG.Color.g(DG.Color.green), DG.Color.b(DG.Color.green));
    const potency = this.higherIsBetter ? value : matrix.minActivity + matrix.maxActivity - value;
    const base = DG.Color.scaleColor(potency, matrix.minActivity, matrix.maxActivity, undefined, ACTIVITY_SCHEME);
    return DG.Color.argb(alpha, DG.Color.r(base), DG.Color.g(base), DG.Color.b(base));
  }

  // ---- Navigator (left pane) ----------------------------------------------------------------

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

  /** Rows and columns of a matrix that survive the current Vary and cell filters. */
  private visibleDims(matrix: SarMatrix): {rows: number, cols: number} {
    const all = this.visibleColIdxs(matrix);
    if (!this.cellFilterActive)
      return {rows: matrix.rows.length, cols: all.length};
    const cols = all.filter((ci) => matrix.cells.some((_row, ri) => this.cellVisible(matrix, ri, ci)));
    const rows = matrix.rows.filter((_row, ri) =>
      cols.some((ci) => this.cellVisible(matrix, ri, ci))).length;
    return {rows, cols: cols.length};
  }

  private get cellFilterActive(): boolean {
    return this.cellPass !== null;
  }

  private cellVisible(matrix: SarMatrix, ri: number, ci: number): boolean {
    return this.cellPass === null || this.cellPass.has(`${matrix.id}:${ri}:${ci}`);
  }


  /** One row per cell of the SELECTED matrix (core/substituent as molecules, potency/refs/weight as
   *  numbers) so the platform supplies the filters. Scoped to the on-screen matrix to avoid
   *  fingerprinting cells the user can't see; rebuilds on a matrix switch. */
  private buildStructFrame(): DG.DataFrame {
    const cores: string[] = [];
    const subs: string[] = [];
    const potency: number[] = [];
    const refs: number[] = [];
    const weights: number[] = [];
    this.cellKeys = [];
    const matrix = this.matrices[Math.min(this.selIndex, this.matrices.length - 1)];
    if (matrix) {
      for (let ri = 0; ri < matrix.rows.length; ri++) {
        for (let ci = 0; ci < matrix.columns.length; ci++) {
          const cell = matrix.cells[ri][ci];
          if (cell.kind === 'empty' || cell.value === null)
            continue;
          cores.push(matrix.rows[ri].coreSmiles);
          subs.push(matrix.columns[ci].substSmiles);
          potency.push(cell.value);
          // A measured compound has no reference points, so it gets its observed-neighbour count.
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

  /** Read the frame's filter into the passing-cell set, then repaint; null while everything passes. */
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


  /** The cell frame's filter group, built on first use. A filter group belongs to a view, so a
   *  headless one is created and never shown — only its filters are, in a popup. */
  private structureFilterRoot(): HTMLElement {
    if (this.structFilters === null) {
      this.structFrame = this.buildStructFrame();
      this.structView = DG.TableView.create(this.structFrame, false);
      // Default set: the platform already gives molecule columns sketchers and numeric columns histograms.
      this.structFilters = this.structView.getFiltersGroup();
      this.structSub = DG.debounce(this.structFrame.onFilterChanged, 300).subscribe(() => this.syncStructFilter());
    }
    return this.structFilters.root;
  }


  /** One row per series, carrying its core as a molecule and the numbers its card prints, so the
   *  platform supplies sketcher + histogram filters. */
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
      // Direction-adjusted so "more potent" is always the larger number on the histogram.
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
      // Plain string so the platform gives it a searchable category list, for looking up a named series.
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
      this.navView = DG.TableView.create(this.navFrame, false);
      this.navFilters = this.navView.getFiltersGroup();
      // Added explicitly: the default set omits an all-distinct category column, but its search box is
      // how a named series is looked for.
      this.navFilters.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: NAV_SERIES}, false);
      this.navSub = DG.debounce(this.navFrame.onFilterChanged, 300).subscribe(() => this.syncNavFilter());
    }
    return this.navFilters.root;
  }

  private passesFilter(matrix: SarMatrix): boolean {
    return this.navPass === null || this.navPass.has(matrix.id);
  }

  /** The score block on the right of a navigator card. Shared by series and transfer cards so both
   *  lists read as one column of scores. */
  private cardScoreBox(lines: {value: string, label: string}[], tip: () => string): HTMLElement {
    const box = ui.divV(lines.map((ln) => ui.divH([
      ui.divText(ln.value, 'chem-sar-card-score'),
      ui.divText(ln.label, 'chem-sar-card-score-cap'),
    ], 'chem-sar-card-score-line')), 'chem-sar-card-score-box');
    ui.tooltip.bind(box, tip);
    return box;
  }

  /** The navigator card's rank score: best and mean potency on the selected scale, plus a hover tip. */
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
    // Remaining scheme is Discontinuity.
    return {lines: [{value: (matrix.scores[SarRankScheme.Discontinuity] ?? 0).toFixed(2), label: 'spread'}],
      tip: `Largest activity spread within a single core, on the ${this.scalingLabel} scale — ` +
        'how discontinuous the SAR is.'};
  }

  /** The structure a matrix is keyed on, drawable. Isotope-marked inherited sites are renumbered as
   *  attachment points so they draw as R labels rather than numbered stars. */
  private matrixCore(matrix: SarMatrix): string {
    if (!matrix.siteKey)
      return matrix.rows[0]?.coreSmiles ?? '';
    let n = 0;
    return matrix.siteKey.replace(/\[\d+\*\]|\[\*:\d+\]/g, () => `[*:${++n}]`);
  }

  /** Parent index per matrix (`-1` for a root) in ranked order; falls back to containment nesting
   *  when no coarser level was built. */
  private matrixParents(): number[] {
    const byId = new Map<string, number>();
    this.matrices.forEach((matrix, i) => byId.set(matrix.id, i));
    const linked = this.matrices.map((matrix) =>
      matrix.parentId !== undefined ? byId.get(matrix.parentId) ?? -1 : -1);
    return linked.some((p) => p >= 0) ? linked : nestByContainment(this.matrices);
  }

  /** One selectable matrix card: aligned core, a descriptor line, and the rank score. `depth`
   *  indents a folded matrix; `folded` is how many are folded into this one, giving its expander. */
  private buildCard(matrix: SarMatrix, index: number, depth = 0, folded = 0,
    onToggle?: () => void): HTMLElement {
    const desc = `${matrix.rows.length} cores · ${matrix.positions.join('/')} · ${matrix.realCount} cpd` +
      (folded ? ` · ${folded} inside` : '');
    const core = this.lazyCoreDepiction(this.matrixCore(matrix), CARD_CORE_W, CARD_CORE_H);
    core.classList.add('chem-sar-card-core');
    // Level badge, since branches bottom out at different depths so indentation alone doesn't say
    // which level a card sits at.
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
    // The expander is its own click target so collapsing doesn't also select the matrix; leaves get
    // a spacer so every core lines up.
    let twisty: HTMLElement;
    if (folded > 0 && onToggle) {
      twisty = ui.iconFA(!this.collapsed.has(matrix.id) ? 'chevron-down' : 'chevron-right', (e: MouseEvent) => {
        e.stopPropagation();
        onToggle();
        // Card isn't rebuilt on toggle (list mutates in place), so flip the icon here.
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
      // Capped so deep nesting doesn't eat the card.
      card.style.marginLeft = `${Math.min(depth, 3) * 10}px`;
      card.classList.add('chem-sar-card-nested');
    }
    if (index === this.selIndex)
      card.classList.add('selected');
    card.onclick = () => this.selectMatrix(index);
    return card;
  }

  /** Card click: move the highlight and rebuild only the matrix pane, leaving the navigator list
   *  (hundreds of cards) untouched. */
  private selectMatrix(index: number): void {
    if (index === this.selIndex)
      return;
    this.navCards[this.selIndex]?.classList.remove('selected');
    this.selIndex = index;
    this.navCards[index]?.classList.add('selected');
    this.varyPosition = ''; // Vary filter is per-matrix; don't carry it over
    this.resetStructFilter(); // cell filter is scoped to one matrix; rebuild for the new one
    this.renderMatrixPane();
  }

  /** A card-core canvas drawn on first visibility (via `navCoreObserver`); eager fallback when no
   *  observer is active. */
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

  /** A transfer's target label with its R-position; series prefix dropped when the target shares the
   *  source's matrix. */
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

  /** Select a transfer and redraw the trend view; only the transfer tab is rebuilt. */
  private selectTransfer(transfer: Transfer): void {
    this.transferIndex = Math.max(0, this.transfers.indexOf(transfer));
    this.renderTransferPanel();
  }

  /** Whether the SAR Transfer tab is the one on screen. */
  private get transferTabActive(): boolean {
    return this.tabs.currentPane?.name === TAB_TRANSFER;
  }

  /** Bring the SAR Transfer tab up to date, computing the (quadratic) transfer list on first open
   *  behind a spinner. No-op when the list is already current. */
  private activateTransferTab(): void {
    if (this.transfersComputing)
      return;
    if (this.computing) {
      // Matrices are still building; compute() ends in render(), which re-enters here.
      ui.empty(this.transferHost);
      this.transferHost.appendChild(ui.divText('Building SAR matrices...', 'chem-sar-empty-note'));
      return;
    }
    if (this.transfersComputed) {
      // List is current; rebuild the pane only if it's not on screen, else repaint the live grid.
      if (this.transferHost.childElementCount === 0)
        this.renderTransferPanel();
      else
        this.transferSlot.state?.grid.invalidate();
      return;
    }
    this.transfersComputing = true;
    ui.empty(this.transferHost);
    ui.setUpdateIndicator(this.transferHost, true, 'Detecting SAR transfers...');
    // Timeout so the tab switch and indicator paint before the synchronous quadratic pass.
    setTimeout(async () => {
      // A matrix rebuild started meanwhile would make these transfers stale; its render() re-enters.
      if (this.detached || this.computing) {
        this.transfersComputing = false;
        ui.setUpdateIndicator(this.transferHost, false);
        return;
      }
      // Captured before the await so a tab bounce / threshold change in the window can't publish a stale scan.
      const generation = this.transferGeneration;
      let detected: Transfer[] | null = null;
      try {
        detected = await computeAllTransfers(this.matrices, this.transferSimilarity, this.higherIsBetter);
      } catch (e) {
        console.error('SAR Matrix | transfer detection failed', e);
      }
      this.transfersComputing = false;
      ui.setUpdateIndicator(this.transferHost, false);
      // A bumped generation means these indices no longer address the on-screen matrices.
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

  /** Rebuild the SAR-transfer tab: the transfer list plus the trend view of the selected one. */
  private renderTransferPanel(): void {
    // Release the tab's own grid (not the matrix's) so switching transfers doesn't stack up grids.
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

  /** Every transfer sharing this one's source core — the pane's target-dropdown alternatives. */
  private transferSiblings(transfer: Transfer): Transfer[] {
    const key = transferSourceKey(transfer);
    return this.transfers.filter((t) => transferSourceKey(t) === key);
  }

  /** One selectable card for a group of transfers sharing a source core. */
  private buildTransferCard(group: Transfer[]): HTMLElement {
    const selected = this.transfers[this.transferIndex] ?? null;
    // Selected transfer when it belongs to this group, otherwise the strongest (first).
    const shown = (selected && group.includes(selected)) ? selected : group[0];

    const body = ui.divV([
      ui.divText(this.sideLabel(shown.a), 'chem-sar-card-title'),
    ], 'chem-sar-card-body');
    // Painted up front, not via the navigator's deferred observer, which belongs to the matrix list
    // and would leave these blank on refill.
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
      // A field assignment raises no onPropertyChanged, so re-rank explicitly here.
      onValueChanged: (value) => {
        this.rankScheme = value!;
        this.reRank();
      },
    });
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
    const navIdleTip = 'Filter series by core structure, potency, SAR spread, size';
    const navIcon = ui.icons.filter(() => {
      ui.showPopup(ui.div(this.navFilterRoot(), 'chem-sar-struct-filters'), navIcon, {vertical: true});
    }, navIdleTip);
    navIcon.classList.add('chem-sar-struct-icon');
    const refill = (): void => {
      this.updateNavVisibility();
      const filtered = this.navPass !== null;
      const hits = this.matrices.filter((matrix) => this.passesFilter(matrix)).length;
      matchCount.innerText = filtered ? `${hits} of ${this.matrices.length} match` : '';
      this.markFilterIcon(navIcon, filtered,
        `${hits} of ${this.matrices.length} series shown — filter active. Click to edit or clear.`, navIdleTip);
    };
    this.navMatchCount = refill;

    const header = ui.divV([sub, ui.divH([ui.form([rankInput]), navIcon], 'chem-sar-nav-controls'), matchCount],
      'chem-sar-nav-header');
    this.fillNavList(list, parents);
    refill(); // initial visibility + match count
    const nav = ui.divV([header, list], 'chem-sar-nav');
    const collapseBtn = ui.iconFA('chevron-left', () => {
      this.navCollapsed = !this.navCollapsed;
      nav.classList.toggle('chem-sar-collapsed', this.navCollapsed);
      this.refitColumns();
    }, 'Collapse / expand the series panel');
    collapseBtn.classList.add('chem-sar-nav-collapse');
    nav.appendChild(collapseBtn);
    nav.classList.toggle('chem-sar-collapsed', this.navCollapsed);
    return nav;
  }

  /** The navigator's child lists, from the parent chain the cards were built against. */
  private navChildren(parents: number[]): Map<number, number[]> {
    const children = new Map<number, number[]>();
    parents.forEach((p, i) => {
      if (p >= 0)
        children.set(p, [...(children.get(p) ?? []), i]);
    });
    return children;
  }

  /** Build every navigator card once (parent immediately followed by its subtree); visibility is
   *  handled separately by {@link updateNavVisibility}. Structures rasterize lazily, so this costs DOM only. */
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
    }, {root: list, rootMargin: '300px'}); // pre-draw past the viewport for smooth scrolling

    const children = this.navChildren(parents);
    const descendants = (i: number): number =>
      (children.get(i) ?? []).reduce((n, child) => n + 1 + descendants(child), 0);

    // A fresh analysis opens collapsed to its roots. Seeded once, so a redraw after a toggle keeps
    // whatever the user has opened since.
    if (!this.collapseSeeded) {
      this.collapsed.clear();
      for (const parent of children.keys())
        this.collapsed.add(this.matrices[parent].id);
      // Raise the selection to its root rather than expanding that root, so the pane opens on the
      // broadest matrix of the family.
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

  /** Show/hide the already-built cards per collapse state and filter (pure display toggles). A parent
   *  that fails the filter still shows when a descendant passes; collapse is ignored while filtering. */
  private updateNavVisibility(): void {
    const parents = this.navParents;
    const filtering = this.navPass !== null;
    const children = this.navChildren(parents);
    const hits = new Map<number, boolean>();
    const anyHit = (i: number): boolean => {
      const cached = hits.get(i);
      if (cached !== undefined)
        return cached;
      hits.set(i, true); // guard against a malformed parent chain looping back on itself
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

  /** Grow cells to fill the pane, never below CELL_W (the table scrolls horizontally instead). */
  private fitCellWidth(nCols: number): number {
    if (nCols <= 0)
      return CELL_W;
    const navW = this.navCollapsed ? NAV_COLLAPSED_W : NAV_W;
    const avail = (this.root.clientWidth || 900) - navW - CORE_W - TABLE_CHROME - nCols * 6;
    return Math.max(CELL_W, Math.min(CELL_W_MAX, Math.floor(avail / nCols)));
  }

  /** Flag a filter icon as active: a dot plus the active tooltip, or the idle tooltip otherwise. */
  private markFilterIcon(icon: HTMLElement, active: boolean, activeTip: string, idleTip: string): void {
    icon.classList.toggle('chem-sar-filter-on', active);
    ui.tooltip.bind(icon, active ? activeTip : idleTip);
  }

  /** Apply the cell filter to the on-screen grid instead of rebuilding: rows via the grid's row
   *  filter, columns via `GridColumn.visible`, failing cells blanked by a repaint. */
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
    const idleTip = 'Filter cells by potency, reference points, core and R-group';
    if (matrix !== null) {
      const {rows, cols} = this.visibleDims(matrix);
      const full = `${matrix.rows.length} cores × ${matrix.columns.length} substituents`;
      const filtered = rows !== matrix.rows.length || cols !== matrix.columns.length;
      if (this.paneDimsChip !== null) {
        this.paneDimsChip.innerText = `${rows}×${cols}`;
        ui.tooltip.bind(this.paneDimsChip, () => filtered ?
          `${rows} cores × ${cols} substituents shown, filtered from ${full}` : full);
      }
      if (this.matrixFilterIcon !== null) {
        this.markFilterIcon(this.matrixFilterIcon, this.cellFilterActive,
          `${rows}×${cols} of ${full} shown — filter active. Click to edit or clear.`, idleTip);
      }
    } else if (this.matrixFilterIcon !== null)
      this.markFilterIcon(this.matrixFilterIcon, this.cellFilterActive, 'Filter active. Click to edit or clear.', idleTip);
    const emptied = state.df.filter.trueCount === 0 || this.visibleGridCols(state) === 0;
    if (this.paneGridHost !== null)
      this.paneGridHost.style.display = emptied ? 'none' : '';
    if (this.paneEmptyNote !== null)
      this.paneEmptyNote.style.display = emptied ? '' : 'none';
  }

  private visibleGridCols(state: MatrixGridState): number {
    let n = 0;
    state.colKeyToIdx.forEach((_idx, key) => {
      if (state.grid.col(key)?.visible !== false)
        n++;
    });
    return n;
  }

  /** Re-fit the live grid's columns to the pane width without a rebuild, so scroll position survives
   *  a resize. Skips sub-pixel changes to avoid repainting on every debounce tick. */
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

  /** Column indices to show: all, or only the "Vary" position's group when one is chosen. */
  private visibleColIdxs(matrix: SarMatrix): number[] {
    const all = matrix.columns.map((_, ci) => ci);
    // Threshold is NOT applied here; it hides columns via `GridColumn.visible` on the live grid.
    return this.varyPosition && matrix.positions.includes(this.varyPosition) ?
      all.filter((ci) => matrix.columns[ci].position === this.varyPosition) : all;
  }

  /** Mean observed potency of a column (real cells only), or null when none. */
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

  /** Caption under the substituent header: mean potency or substituent MW; empty when labelling is off. */
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

  /** Virtualized render of a row/column spec: a scaffold DataFrame backs a DG.Grid whose every cell
   *  is hand-painted in `onCellRender`, so only viewport cells draw. Each displayed row resolves
   *  against its own matrix, which lets a cross-series transfer show rows from two matrices. */
  private buildPaneGrid(rows: PaneRow[], columns: PaneColumn[], slot: PaneGridSlot): HTMLElement {
    this.cellW = this.fitCellWidth(columns.length);
    // Group boundaries computed once here; rescanning in `onCellRender` would be quadratic per repaint.
    const firstOfGroup = new Set<number>();
    columns.forEach((column, i) => {
      if (i === 0 || columns[i - 1].position !== column.position)
        firstOfGroup.add(i);
    });

    const df = DG.DataFrame.create(rows.length);
    df.columns.addNewString('Core');
    const colKeyToIdx = new Map<string, number>();
    // Stable string keys (never the grid column idx, which pinning shifts).
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
    // Shared template only when every row genuinely shares one concrete core; otherwise null and
    // those panes align per row (a single template would silently fail to align dissimilar cores).
    const firstCore = rows.length > 0 ? rows[0].matrix.rows[rows[0].rowIndex]?.coreSmiles ?? null : null;
    const sharedCore = firstCore !== null &&
      rows.every((row) => row.matrix.rows[row.rowIndex]?.coreSmiles === firstCore) ? firstCore : null;
    // Cap attachment points off the header so it names the shared scaffold without labelling filled sites.
    const headerCore = sharedCore === null ? null : sharedCore.replace(/\[\*:\d+\]/g, '[H]');
    const paneTemplate = headerCore !== null ? buildAlignmentTemplate(headerCore) : null;
    const state: MatrixGridState = {grid, df, rows, columns, colKeyToIdx, firstOfGroup, headerCore, paneTemplate,
      rowTemplates: new Array(rows.length).fill(null)};

    // Owned by this grid instance; unsubscribed when replaced or on detach so renders don't leak handlers.
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
      // The context already carries the DPR scale, so paint in CSS coords; save/restore isolates state.
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

    // Track the pointer for ctrl/shift extend, and decode right-clicks to the virtual cell under them.
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

    // Plain flex host that fills the pane, not ui.box (which would pin a fixed width).
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    return ui.div([grid.root], 'chem-sar-grid-host');
  }

  /** Resolve a grid body cell to its row descriptor and matrix column index; null for 'Core' and
   *  anything outside the built set. */
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

  /** Click on a grid body cell: open the Molecule context, set the host grid's current row/selection,
   *  and register the cell's SAR panel against whichever matrix the clicked side belongs to. */
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
    // A threshold-blanked cell draws empty and isn't selectable.
    if (cell.kind === 'empty' || cell.value === null || !this.cellVisible(matrix, ri, ci))
      return;
    grok.shell.windows.showContextPanel = true;
    // A real compound becomes the current row and extends selection; a virtual analog has no row.
    if (cell.molIdx !== null) {
      this.dataFrame.currentRowIdx = cell.molIdx;
      const event = this.lastGridMouseEvent ?? new MouseEvent('click');
      this.dataFrame.selection.handleClick((i) => i === cell.molIdx, event, true);
    } else
      this.dataFrame.currentRowIdx = -1;
    // Assembled cell opens the Molecule context (SAR context gated into the "SAR analysis" pane);
    // an unassembled virtual falls back to the standalone SAR panel.
    if (cell.smiles) {
      this.analogPanels.set(cell.smiles, () => this.buildCellPanel(matrix, ri, ci));
      grok.shell.o = DG.SemanticValue.fromValueType(cell.smiles, DG.SEMTYPE.MOLECULE);
    } else
      grok.shell.o = this.buildCellPanel(matrix, ri, ci);
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
    // A threshold-blanked cell must hover as empty too, or the hidden value returns under the cursor.
    if (cell.kind === 'empty' || cell.value === null || !this.cellVisible(paneRow.matrix, paneRow.rowIndex, ci))
      return false;
    ui.tooltip.show(ui.divText(this.cellTooltipText(cell)), x, y);
    return true;
  }

  /** The id chip's text for an observed cell; null for virtual analogs (no source row to identify). */
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

  /** Paint an R-group column header (position band + substituent depiction + sort caption) or the
   *  'Core' column's "Aligned core" label. */
  private paintHeader(g: CanvasRenderingContext2D, b: DG.Rect, name: string, state: MatrixGridState): void {
    const grey6 = cssColor(this.root, '--grey-6', '#4a4a4a');
    const grey5 = cssColor(this.root, '--grey-5', '#7d7d7d');
    if (name === 'Core') {
      // Each header paints its own background (default paint suppressed), else it keeps stale pixels.
      g.fillStyle = DG.Color.toHtml(HEADER_ARGB);
      g.fillRect(b.x, b.y, b.width, b.height);
      const captionH = 14;
      if (state.headerCore !== null) {
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

    // Position label only on the first column of a group (per-cell clipping can't span columns).
    const posBandH = 16;
    if (state.firstOfGroup.has(idx)) {
      g.fillStyle = grey6;
      g.font = `600 11px ${GRID_FONT}`;
      g.textAlign = 'left';
      g.textBaseline = 'top';
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

  /** The molblock a row's cells align to: the row key aligned to the pane's shared core (so cells
   *  match the core printed beside them), else laid out on its own. */
  private rowTemplate(state: MatrixGridState, rowIdx: number): string | null {
    const cached = state.rowTemplates[rowIdx];
    if (cached !== null)
      return cached;
    const paneRow = state.rows[rowIdx];
    const key = paneRow.matrix.rows[paneRow.rowIndex].keySmiles;
    const aligned = state.paneTemplate !== null ? alignToTemplate(key, state.paneTemplate) : '';
    // A key that can't align to the shared core comes back untouched, so lay it out on its own.
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
    // Fill in the depiction's baked-in colour so anti-aliased bonds don't fringe.
    g.fillStyle = DG.Color.toHtml(WHITE_ARGB);
    g.fillRect(b.x, b.y, b.width, b.height);
    drawDepiction(g, b.x, b.y, b.width, molH, row.keySmiles, template, WHITE_ARGB);
    g.fillStyle = cssColor(this.root, '--grey-6', '#4a4a4a');
    g.font = `600 11px ${GRID_FONT}`;
    g.textAlign = 'center';
    g.textBaseline = 'top';
    g.fillText(paneRow.label, b.x + b.width / 2, b.y + molH + 1, b.width);
  }

  /** Paint one core×substituent cell: potency tint, aligned molecule, value/id chips, markers, and
   *  selection ring. Tint uses the row's OWN matrix, so transfer sides colour against their own range. */
  private paintBodyCell(g: CanvasRenderingContext2D, b: DG.Rect, paneRow: PaneRow, rowIdx: number,
    colIndex: number, state: MatrixGridState): void {
    const matrix = paneRow.matrix;
    const ri = paneRow.rowIndex;
    const cell = matrix.cells[ri][paneRow.colIdxs[colIndex]];
    // A cell below the threshold is blanked, not removed, so the grid stays a readable table.
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
    // Flattened opaque and used as the depiction background so bonds blend into the tint; fill first
    // for cells with no structure.
    const bg = this.flatTint(matrix, value, alpha);
    g.fillStyle = DG.Color.toHtml(bg);
    g.fillRect(b.x, b.y, b.width, b.height);
    if (cell.smiles !== null)
      drawDepiction(g, b.x, b.y, b.width, b.height, cell.smiles, this.rowTemplate(state, rowIdx), bg);

    // '~value' chip, top-left; green when it's the row's most potent observation.
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

  /** Draw one corner chip: a translucent white plate carrying a line of text (kept legible over a
   *  bond). Text is ellipsized to the cell rather than condensed by `fillText`'s max width. */
  private paintChip(g: CanvasRenderingContext2D, b: DG.Rect, text: string, style: ChipStyle): void {
    g.save();
    g.font = `${style.italic === true ? 'italic ' : ''}600 10px ${GRID_FONT}`;
    g.textAlign = 'left';
    g.textBaseline = 'top';
    const trimmed = TextUtils.trimText(text, g, b.width - (CHIP_MARGIN + CHIP_PAD) * 2);
    const chipW = g.measureText(trimmed).width + CHIP_PAD * 2;
    // Measure the far corner from x/y + size: the rect has no right/bottom accessors (would be NaN).
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

  /** A context-panel section header with an optional provenance badge. */
  private cpSection(title: string, badge?: string): HTMLElement {
    const parts = [ui.divText(title, 'chem-sar-cp-section-title')];
    if (badge)
      parts.push(ui.divText(badge, 'chem-sar-cp-badge'));
    return ui.divH(parts, 'chem-sar-cp-section');
  }

  /** A label / value row in the context panel. */
  private cpRow(label: string, value: HTMLElement | string): HTMLElement {
    const v = typeof value === 'string' ? ui.divText(value, 'chem-sar-cp-rv') : value;
    return ui.divH([ui.divText(label, 'chem-sar-cp-rl'), v], 'chem-sar-cp-row');
  }

  /** The cell's structure, drawn to the panel's actual width and redrawn on resize past a few pixels
   *  (each redraw is an RDKit render). */
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
    // A newer panel may claim the single observer slot while this one is still awaiting layout.
    const token = ++this.cpStructureToken;
    this.cpStructureSub?.unsubscribe();
    // Wired up only once the box is in the DOM, else clientWidth is 0 and the first draw is too small.
    ui.tools.waitForElementInDom(box).then(() => {
      if (token !== this.cpStructureToken)
        return;
      // Measure the panel container (what a drag resizes) and the box, smaller wins, so the structure
      // grows with the panel but never gets stuck at a stale box width.
      const panel = box.closest('.panel-content') as HTMLElement | null;
      const fit = (): void => draw(panel === null ? Math.floor(box.clientWidth) - BOX_CHROME :
        Math.min(Math.floor(box.clientWidth) - BOX_CHROME, Math.floor(panel.clientWidth) - PANEL_CHROME));
      this.cpStructureSub?.unsubscribe();
      this.cpStructureSub = DG.debounce(ui.onSizeChanged(panel ?? box), 50).subscribe(fit);
      fit();
    });
    return box;
  }

  /** Where a predicted value came from: the additive model `rowMean + colMean - grandMean` with the
   *  contributing cells and the arithmetic spelled out. */
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

    // Shown as deviations from the matrix mean, which is what the model actually adds up.
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

  /** Observed compounds sharing this cell's core and substituent — its SAR context. The cell itself
   *  is excluded and the cell filter applies; the badge reads "shown of total" when some are hidden. */
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

  /** A framed fragment tile: the structure with the compound's value beneath; structure omitted for
   *  an empty SMILES. */
  private cpFragment(smiles: string | null, value?: number): HTMLElement {
    const parts: HTMLElement[] = [];
    if (smiles)
      parts.push(ui.div([renderMolecule(smiles, {width: 78, height: 52, popupMenu: false})], 'chem-sar-cp-frag-box'));
    if (value !== undefined)
      parts.push(ui.divText(this.formatActivity(value), 'chem-sar-cp-rv'));
    return ui.divV(parts, 'chem-sar-cp-frag');
  }

  /** SAR context for the selected cell: header, structure, potency, the prediction/references, and —
   *  for a virtual analog — a make-list action. */
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
      // Call out a single-compound arm: that's where the estimate is really an extrapolation.
      const refs = cell.references ?? 0;
      const weakest = cell.support ?? 0;
      panel.appendChild(this.cpRow('Reference points',
        `n = ${refs}${weakest <= 1 ? ' · one arm has a single compound' : ''}`));
      panel.appendChild(this.cpPrediction(matrix, rowIdx, colIdx));
    } else
      panel.appendChild(this.cpReferences(matrix, rowIdx, colIdx));

    panel.appendChild(this.cpSection('Decomposition', 'R-group'));
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

  private transferPaneRow(side: TransferSide, colIdxs: number[]): PaneRow {
    return {matrix: this.matrices[side.matrixIndex], rowIndex: side.rowIndex, colIdxs,
      label: this.sideLabel(side), markBest: true};
  }

  /** Most potent OBSERVED value among a row's displayed cells, or null; predictions excluded. */
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

  /** The SAR transfer view: two cores whose potency trends run in parallel, as two rows of the same
   *  virtualized grid, each resolving against its own matrix. */
  private buildTransferPane(transfer: Transfer): HTMLElement {
    const from = this.sideLabel(transfer.a);
    const to = this.sideLabel(transfer.b);
    const bar = ui.divH([
      ui.divText('Cross-series SAR transfer', 'chem-sar-main-title'),
      ui.divText(`${from} → ${to} · across ${transfer.a.position} · r = ${transfer.correlation.toFixed(2)}`),
    ], 'chem-sar-main-bar');

    // Dropdown to switch targets when this source transfers to more than one.
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

    // Each side keeps its own matrix. Measured pairs first, then the predicted analogs the transfer argues for.
    const rows = [
      this.transferPaneRow(transfer.a, [...transfer.aCols, ...transfer.predictedACols]),
      this.transferPaneRow(transfer.b, [...transfer.bCols, ...transfer.predictedBCols]),
    ];
    const matchedCount = transfer.substituents.length;
    const columns: PaneColumn[] = [...transfer.substituents, ...transfer.predictedSubstituents]
      .map((substSmiles, i) => ({substSmiles, position: transfer.a.position,
        caption: i < matchedCount ? '' : 'predicted'}));
    const gridHost = this.buildPaneGrid(rows, columns, this.transferSlot);
    // Two rows only: size to content, not stretch, so the statistics stay on screen.
    gridHost.style.flex = '0 0 auto';
    gridHost.style.height = `${COL_HEADER_H + rows.length * CELL_H + GRID_SCROLLBAR_H}px`;

    // Plain div, not ui.box (which would pin a fixed width).
    const scroll = ui.div([this.buildTransferStats(transfer)], 'chem-sar-main-scroll');
    const parts = controlBar ? [bar, controlBar, gridHost, scroll] : [bar, gridHost, scroll];
    return ui.divV(parts, 'chem-sar-main');
  }

  /** The "Transfer statistics" block: correlation, fold-change match, and the benefiting cell. */
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

  /** Compact summary chips for the current matrix, full detail on hover. */
  private buildChips(matrix: SarMatrix): HTMLElement {
    const chip = (text: string, tip: string, cls = ''): HTMLElement => {
      const el = ui.divText(text, `chem-sar-chip-badge ${cls}`.trim());
      ui.tooltip.bind(el, () => tip);
      return el;
    };
    const items = [
      chip(`${matrix.realCount} cpd`, `${matrix.realCount} observed compounds`),
      // Reports what is on screen, not what the matrix holds (a threshold makes them differ).
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
    // this.transfers is sorted by correlation, so the first match is the strongest.
    const involving = this.transfers.filter((t) => t.a.matrixIndex === idx || t.b.matrixIndex === idx);
    if (involving.length) {
      const best = involving[0];
      items.push(chip(`r ${best.correlation.toFixed(2)}`,
        'Strongest SAR-transfer correlation between this series and another', 'chem-sar-chip-transfer'));
    }
    return ui.divH(items, 'chem-sar-chips');
  }

  private buildMatrixPane(matrix: SarMatrix): HTMLElement {
    const infoBar = ui.divH([ui.divText(matrix.label, 'chem-sar-main-title'), this.buildChips(matrix)],
      'chem-sar-main-bar');

    const controls: HTMLElement[] = [];
    if (matrix.positions.length >= 2) {
      const varyValue = this.varyPosition && matrix.positions.includes(this.varyPosition) ?
        this.varyPosition : 'All';
      const varyInput = ui.input.choice('Vary', {
        value: varyValue,
        items: ['All', ...matrix.positions],
        onValueChanged: (value) => {
          this.varyPosition = value === 'All' ? '' : value!;
          this.renderMatrixPane();
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

    // All filters behind one icon; the sketchers would dwarf the bar inline.
    const filterIcon = ui.icons.filter(() => {
      ui.showPopup(ui.div(this.structureFilterRoot(), 'chem-sar-struct-filters'),
        filterIcon, {vertical: true});
    }, 'Filter cells by potency, reference points, core and R-group');
    filterIcon.classList.add('chem-sar-struct-icon');
    this.matrixFilterIcon = filterIcon;
    controls.push(filterIcon);
    const controlBar = ui.divH(controls, 'chem-sar-control-bar');

    const visible = this.visibleColIdxs(matrix);
    // Every row is built; the filter hides them via the grid's row filter, so no rebuild on filter change.
    const rows: PaneRow[] = matrix.rows.map((row, ri) =>
      ({matrix, rowIndex: ri, colIdxs: visible, label: row.label}));
    const columns: PaneColumn[] = visible.map((ci) => ({
      substSmiles: matrix.columns[ci].substSmiles,
      position: matrix.columns[ci].position,
      caption: this.columnSortCaption(matrix, ci),
    }));

    const gridHost = this.buildPaneGrid(rows, columns, this.matrixSlot);
    // Shown in the grid's place when the filter empties the matrix; toggled by the filter, no rebuild.
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
    // Preserve the navigator scroll position across the rebuild.
    const prevNav = this.host.querySelector('.chem-sar-nav-list');
    const navScroll = prevNav instanceof HTMLElement ? prevNav.scrollTop : 0;
    // Release the old grid before its DOM goes.
    this.releaseMatrixGrid();
    ui.empty(this.host);
    // Palette cache is keyed by name only, so clear it on a possible theme switch.
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
    this.renderMatrixPane();
    // Only when it's the tab on screen: rebuilt hidden, its grid would size against a display:none pane.
    if (this.transferTabActive)
      this.activateTransferTab();
  }

  /** Rebuild only the matrix pane, leaving the navigator DOM and its scroll alone. */
  private renderMatrixPane(): void {
    if (this.matrices.length === 0)
      return;
    this.releaseMatrixGrid();
    this.host.querySelector(':scope > .chem-sar-main')?.remove();
    const matrix = this.matrices[Math.min(this.selIndex, this.matrices.length - 1)];
    this.host.appendChild(this.buildMatrixPane(matrix));
    this.syncSelection();
  }
}
