/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import '../../../css/sar-matrix.css';
import {_package} from '../../package';
import {getRdKitModule, getRdKitService} from '../../utils/chem-common-rdkit';
import {SCALING_METHODS} from '../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {scaleActivity} from '../molecular-matched-pairs/mmp-viewer/mmpa-utils';
import {nestByContainment, rankMatrices, SarRankScheme} from './build/sar-matrix-ranking';
import {DEFAULT_TRANSFER_SIMILARITY} from './transfer/sar-matrix-transfer';
import {SarFragmentColumns, standardizeFragment} from './build/sar-matrix-columns';
import {MAX_SERIES_LEVELS, runSarMatrix, SarGrouping, SarMatrixParams} from './build/sar-matrix-run';
import {closeGridQuietly, finiteOrNaN, observedMolecules, SarMatrix, SarMatrixCell, SarMatrixCellKind}
  from './sar-matrix-types';
import {CARD_CORE_H, CARD_CORE_W, CELL_H, CELL_W, CELL_W_MAX, chipBadge,
  clearCssColorCache, COL_HEADER_H, CORE_BG_ARGB, CORE_W, MatrixCellRef, MatrixGridState,
  NAV_COLLAPSED_W, NAV_W,
  FrameFilter, paintMoleculeOnColor, PaneColumn, PaneGridSlot, PaneRow, renderMoleculeOnColor,
  scrollWithin, TAB_MAKELIST, TAB_MATRIX, TAB_SUMMARY, TAB_TRANSFER, TABLE_CHROME} from './sar-matrix-ui-common';
import {buildAlignmentTemplate, clearDepictionCaches, coreDepictionBlock, matrixCore} from './render/sar-matrix-depict';
import {MakeListPanel} from './sar-matrix-make-list';
import {MatrixPainter} from './render/sar-matrix-paint';
import {SummaryPanel} from './summary/sar-matrix-summary-panel';
import {TransferPanel} from './transfer/sar-matrix-transfer-panel';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';

// Above the filter popup that opened the sketcher, which the platform stacks at 10000.
const SKETCHER_DIALOG_Z = 10001;

type AnalogPanelBuilder = () => HTMLElement;

/** How long the current object is defended after a cell click. Long enough to outlast the platform's
 *  restore, short enough that the next thing the user selects is theirs. */
const CONTEXT_GUARD_MS = 400;

/** Scaffold keys are derived per matrix and reused across renders; one entry per matrix in play. */
const SCAFFOLD_KEY_CACHE_MAX = 4000;

const COLSORT_POTENCY = 'Potency';
const COLSORT_MW = 'Molecular weight';
const COLUMN_SORTS = [COLSORT_POTENCY, COLSORT_MW];

/** Properties that only reorder/recolor already-assembled matrices — must NOT re-run fragmentation. */
const RERANK_ONLY_PROPS = ['rankScheme', 'activityDirection'];

/** Properties that only change what is drawn — no re-fragmentation and no re-ranking. */
const RENDER_ONLY_PROPS = ['columnCaption', 'idColumnName'];

/** Properties only the transfer scan reads; matrices are untouched. */
const TRANSFER_ONLY_PROPS = ['transferSimilarity'];
/** Read only when fragmenting, so they cannot change matrices built from existing R-groups. */
const FRAGMENT_ONLY_PROPS = ['fragmentCutoff', 'grouping', 'fragmentationLevels', 'threshold',
  'useMcsAnchors'];


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

/** Put the grid and the matrix in one tab group rather than side by side. Docking both FILL at the root
 *  gives each the view's full width; a split leaves the matrix showing two of its columns, and the
 *  navigator alone claims a fixed 320px of it. The table stays a tab away instead of being closed. */
export function dockSarMatrixTabs(view: DG.TableView, viewer: DG.Viewer): void {
  // Anchored to the grid's own node: a FILL into the container root re-forms the tab group from
  // whatever the container holds at that moment, so a container-level change can leave one panel and
  // the tab-handle strip with it.
  const gridNode = view.dockManager.findNode(view.grid.root);
  // A miss falls back to the container root, which is the shape above — so it is said out loud rather
  // than left looking like the anchored call.
  if (gridNode === undefined) {
    _package.logger.warning('SAR Matrix | the grid has no dock node; docking into the container root, ' +
      'which re-forms its tab group');
  }
  view.dockManager.dock(viewer, DG.DOCK_TYPE.FILL, gridNode ?? null, 'SAR Matrix');
}


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


export class SarMatrixViewer extends DG.JsViewer {
  moleculesColumnName: string;
  /** Empty unless the user grouped the series themselves; then it replaces the automatic grouping. */
  seriesColumnName: string;
  activityColumnName: string;
  /** Optional column captioning each observed cell. */
  idColumnName: string;
  /** Empty unless the matrices come from R-group columns already in the table. The rows are the
   *  R-groups other than the axis. */
  coreColumnName: string;
  rGroupColumnNames: string[];
  axisColumnName: string;
  scaling: string;
  activityDirection: string;
  fragmentCutoff: number;
  grouping: string;
  fragmentationLevels: number;
  threshold: number;
  transferSimilarity: number;
  predictVirtual: boolean;
  predictUnmeasured: boolean;
  minCompounds: number;
  useMcsAnchors: boolean;
  rankScheme: string;

  matrices: SarMatrix[] = [];
  private selIndex = 0;
  /** "Vary" filter: show only this R-position's column group, or all when empty. */
  private varyPosition = '';
  columnCaption: string;
  /** The assembled set, serialized so a saved project restores it instead of spending minutes
   *  rebuilding. Kept out of layouts: it is megabytes and only lines up with the table it was built
   *  from, while a layout is meant to stay small and portable across tables. */
  matricesData: string;
  private contextCell: MatrixCellRef | null = null;
  /** Last cell clicked in either pane's grid, so the make-list tab can act on it without the Context
   *  Panel. Dropped on recompute: it holds a matrix that the new set has replaced. */
  selectedCell: MatrixCellRef | null = null;
  /** Per-SMILES "SAR analysis" panel builders; cleared on recompute and detach. */
  private readonly analogPanels = new Map<string, AnalogPanelBuilder>();
  /** Armed only for the moment after a cell click, to take the context back if the platform restores
   *  the table view's previous current object over it. */
  private contextGuard: {unsubscribe(): void} | null = null;
  private readonly host = ui.divH([], 'chem-sar-matrix');
  private readonly tabs: DG.TabControl;
  /** The three panels split out of this class. They reach back through the members below, which are
   *  public for that reason rather than for outside callers. */
  private readonly painter = new MatrixPainter(this);
  private readonly makeListPanel = new MakeListPanel(this);
  private readonly transferPanel = new TransferPanel(this);
  private readonly summaryPanel = new SummaryPanel(this);
  computing = false;
  /** A recompute was requested mid-compute; re-queued when the running one finishes. */
  private dirty = false;
  private computeTimer = 0;
  private reported = '';
  /** Set in `detach`, so an in-flight compute can't render into a closed viewer. */
  detached = false;
  private cellW = CELL_W;
  /** Live pane grids, held so selection changes repaint via a cheap `invalidate()`. */
  private readonly matrixSlot: PaneGridSlot = {state: null, subs: []};
  /** Last pointer event over the grid, so a cell click can honor ctrl/shift (onCellClick carries none). */
  private lastGridMouseEvent: MouseEvent | null = null;
  /** Matrix ids whose children are folded away; kept so re-rank/resize doesn't reopen them. */
  private readonly collapsed = new Set<string>();
  /** Cards per matrix (index-aligned) and their parent chain; built once, then mutated in place. */
  private navCards: HTMLElement[] = [];
  private navParents: number[] = [];
  /** Rasterizes a card's core only when it first scrolls into view. */
  private navCoreObserver: IntersectionObserver | null = null;
  private readonly navPendingCores: Map<Element, () => void> = new Map();
  /** The matrix the cell filter was built from; cells of any other are outside it. */
  private cellPassMatrixId: string | null = null;
  /** One row per series. */
  private readonly navFilter = new FrameFilter<string>(
    () => this.buildNavFrame(),
    (group) => {
      // Explicit: the default set omits an all-distinct category, whose search box finds a named
      // series, and the sketcher matches the core rather than the compounds under it.
      group.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: NAV_SERIES}, false);
      group.updateOrAdd({type: 'Chem:substructureFilter', column: NAV_CORE}, false);
    },
    () => {
      this.updateNavVisibility();
      this.navMatchCount?.();
    });
  /** Refreshes the navigator match count from outside its render. */
  private navMatchCount: (() => void) | null = null;
  /** One row per (core, substituent) cell of the matrix on screen. */
  private readonly cellFilter = new FrameFilter<string>(
    () => this.buildStructFrame(),
    null,
    () => {
      this.cellPassMatrixId = this.cellFilter.active ?
        this.matrices[Math.min(this.selIndex, this.matrices.length - 1)]?.id ?? null : null;
      this.applyCellFilter();
    });
  /** Matrix-pane parts the potency threshold updates in place; null while no pane is on screen. */
  private paneGridHost: HTMLElement | null = null;
  private paneEmptyNote: HTMLElement | null = null;
  private paneDimsChip: HTMLElement | null = null;
  /** The matrix pane's filter icon, so the cell filter can flag it as active. */
  private matrixFilterIcon: HTMLElement | null = null;
  /** Whether the navigator is collapsed; on the viewer so it survives a navigator rebuild. */
  private navCollapsed = false;
  /** Header rows standing above the roots that key on one scaffold; see {@link scaffoldKeyOf}. */
  private navGroups: {header: HTMLElement, members: number[]}[] = [];
  private navGroupOfRoot = new Map<number, string>();
  private collapsedScaffolds = new Set<string>();
  private scaffoldKeys = new DG.LruCache<string, string>(SCAFFOLD_KEY_CACHE_MAX);
  private scaffoldReps = new Map<string, string>();
  private scaffoldRepsFor = '';
  /** Guards the one-time collapse-to-roots per analysis. */
  private collapseSeeded = false;
  /** See {@link hierarchy}; null whenever `matrices` has been replaced or reordered. */
  private hierarchyCache: {parents: number[], roots: string[], tiers: number[]} | null = null;
  /** Values the chosen scaling cannot represent. Counted here rather than read back from the run: a
   *  set carried in by a layout skips `runSarMatrix` entirely and would report none. */
  unscalableCount = 0;
  /** Rows the assay reached at all, whatever the scaling can then do with them. A compound nobody
   *  measured and one the decomposition could not place are different problems with different fixes,
   *  and the coverage bar can only separate them if it knows this. */
  assayedCount = 0;
  get helpUrl() {
    return 'https://raw.githubusercontent.com/datagrok-ai/public/refs/heads/master/help/datagrok/solutions/domains/chem/chem.md#sar-matrix';
  }

  /** The analysis and the series on screen: where the matrices came from and what their cells hold. */
  getWidgetStatus(): DG.U2.IWidgetStatus {
    const base = super.getWidgetStatus();
    const matrix = this.matrices[Math.min(this.selIndex, this.matrices.length - 1)];
    const cells = matrix?.cells.flat() ?? [];
    const count = (kind: SarMatrixCellKind): number => cells.filter((cell) => cell.kind === kind).length;
    return {...base, values: {...base.values,
      'computing': this.computing,
      'problem': this.rGroupFault() ?? '',
      'source': this.rGroupColumns() !== null ? 'R-group columns' : 'fragments',
      'matrices': this.matrices.length,
      'series': matrix?.label ?? '',
      'rows': matrix?.rows.length ?? 0,
      'columns': matrix?.columns.length ?? 0,
      'measured': count('real'),
      'predicted': count('virtual'),
      'predicted with structure': cells.filter((cell) => cell.kind === 'virtual' && cell.smiles !== null).length,
      'impossible': count('impossible'),
      'row labels': matrix?.rows.slice(0, 20).map((row) => row.label).join(', ') ?? '',
    }};
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
    this.seriesColumnName = this.addProperty('seriesColumnName', DG.TYPE.COLUMN, '',
      {nullable: true, category: 'Data', friendlyName: 'Series column',
        description: 'Optional: series you have assigned yourself. Compounds sharing a value make one ' +
          'matrix, named with that value. Leave empty to group by structure instead'});
    this.coreColumnName = this.addProperty('coreColumnName', DG.TYPE.COLUMN, '',
      {nullable: true, category: 'R-groups', friendlyName: 'Core',
        description: 'Column with the core, its attachment points marked [*:1], [*:2], …. Set it with ' +
          'R-groups and Matrix columns to use an existing R-group decomposition instead of fragmenting'});
    this.rGroupColumnNames = this.columnList('rGroupColumnNames', [],
      {nullable: true, columnTypeFilter: DG.TYPE.CATEGORICAL, category: 'R-groups', friendlyName: 'R-groups',
        description: 'Columns with the substituent at each attachment point of the core'});
    this.axisColumnName = this.addProperty('axisColumnName', DG.TYPE.COLUMN, '',
      {nullable: true, category: 'R-groups', friendlyName: 'Matrix columns',
        description: 'The R-group whose substituents become the matrix columns. The core and the other ' +
          'R-groups make up the rows'});
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
    this.predictUnmeasured = this.bool('predictUnmeasured', true, {friendlyName: 'Predict untested compounds',
      description: 'Also predict compounds the dataset already holds but has no activity for. They ' +
        'stay marked as made, so they are offered for testing rather than synthesis'});
    this.minCompounds = this.int('minCompounds', 3, {min: 1, max: 10, friendlyName: 'Min compounds',
      description: 'Fewest DISTINCT measured compounds a matrix must hold. Two compounds spread over ' +
        'a large grid predict every cell from the same pair of numbers, which reads as a dense matrix ' +
        'built from no SAR at all'});
    this.useMcsAnchors = this.bool('useMcsAnchors', false, {friendlyName: 'Group leftovers by MCS',
      description: 'Also cover the compounds no shared core could group, by searching those for a ' +
        'common core. No matrix is lost, though a series whose core is too small to anchor is pooled ' +
        'with the leftovers and can end up varying several positions. Costly on large sets'});
    this.rankScheme = this.string('rankScheme', SarRankScheme.Potency,
      {choices: [SarRankScheme.Potency, SarRankScheme.Discontinuity, SarRankScheme.Preferred],
        friendlyName: 'Rank by', description: 'How the navigator orders the matrices'});
    this.columnCaption = this.string('columnCaption', COLSORT_POTENCY, {choices: COLUMN_SORTS});
    this.matricesData = this.string('matricesData', '', {userEditable: false, includeInLayout: false});
    this.host.style.height = '100%';
    this.transferPanel.root.style.height = '100%';
    this.makeListPanel.root.style.height = '100%';
    this.summaryPanel.root.style.height = '100%';
    // Capture-phase reset runs before a cell's bubbling handler, so contextCell reflects only a
    // right-click on a virtual cell. Bound once: `host` outlives every attach, so re-binding stacks copies.
    this.host.addEventListener('contextmenu', () => this.contextCell = null, true);

    // Transfer detection is quadratic in the total row count, so that tab computes on first open.
    this.tabs = ui.tabControl(null, {vertical: false});
    // First, so the analysis opens on what it found rather than on one arbitrary series.
    const summaryPane = this.tabs.addPane(TAB_SUMMARY, () => this.summaryPanel.root);
    ui.tooltip.bind(summaryPane.header, 'What the whole analysis found, on one screen: the scale it is ' +
      'read on, the three answers, and which series to open — with the fitted effects, the analogs ' +
      'worth making and the methodology each a segment away');
    const matrixPane = this.tabs.addPane(TAB_MATRIX, () => this.host);
    ui.tooltip.bind(matrixPane.header, 'Core × substituent potency matrices, one per series');
    const transferPane = this.tabs.addPane(TAB_TRANSFER, () => this.transferPanel.root);
    ui.tooltip.bind(transferPane.header, 'Pairs of cores whose potency trends run in parallel across the ' +
      'R-groups they have both explored — detected when the tab is first opened');
    const makeListPane = this.tabs.addPane(TAB_MAKELIST, () => this.makeListPanel.root);
    ui.tooltip.bind(makeListPane.header, 'Compounds you have collected — predicted analogs and made ones alike');
    this.makeListPanel.renderMakeList();
    // Not on this.subs: tab switching must survive a detach/re-attach cycle.
    this.tabs.onTabChanged.subscribe(() => {
      const name = this.tabs.currentPane?.name;
      if (name === TAB_TRANSFER)
        this.transferPanel.activateTransferTab();
      else if (name === TAB_SUMMARY)
        this.summaryPanel.activateSummaryTab();
      else if (name === TAB_MATRIX) {
        // The matrix pane is usually built while hidden behind the Summary tab, so its grid took both
        // its column widths and its canvas metrics against a display:none host.
        this.refitColumns();
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
    // A sketcher opened from the filter popup lands under it. Doing this in CSS would restyle every
    // sketcher dialog in the platform, since dialogs are body-level and the classes are the platform's.
    this.subs.push(grok.events.onDialogShown.subscribe((dialog: DG.Dialog) => {
      // Matched by class, not a tracked element: both the navigator's and the matrix's filter popups
      // host a sketcher, and only one of them is opened through showFilterPopup.
      if (document.querySelector('.d4-popup-host .chem-sar-struct-filters') === null)
        return;
      const sketcher = dialog.root.querySelector('.chem-sketcher-host');
      if (!(sketcher instanceof HTMLElement))
        return;
      dialog.root.style.zIndex = `${SKETCHER_DIALOG_Z}`;
    }));
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
    this.contextGuard?.unsubscribe();
    this.contextGuard = null;
    this.analogPanels.clear();
    this.selectedCell = null;
    this.releaseMatrixGrid();
    this.transferPanel.release();
    this.makeListPanel.release();
    this.summaryPanel.release();
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
  releaseSlot(slot: PaneGridSlot): void {
    slot.subs.forEach((s) => s.unsubscribe());
    slot.subs = [];
    closeGridQuietly(slot.state?.grid);
    slot.state = null;
  }

  private releaseMatrixGrid(): void {
    this.releaseSlot(this.matrixSlot);
    // Pane is about to be replaced; holding these would let the threshold write into detached elements.
    this.paneGridHost = null;
    this.paneEmptyNote = null;
    this.paneDimsChip = null;
  }

  /** Rebuild the filters against the current matrices: stale keys partially match reused cluster ids,
   *  which mis-filters rather than failing cleanly. */
  private resetFilters(): void {
    this.resetStructFilter();
    this.navFilter.reset();
    this.navMatchCount = null;
  }

  /** Drop only the cell filter so it rebuilds against the current matrix; its rows key the previous one. */
  private resetStructFilter(): void {
    this.cellFilter.reset();
    this.cellPassMatrixId = null;
  }

  /** Repaint both panes so a host-grid selection/current-row change shows on the rendered cells. */
  private syncSelection(): void {
    if (!this.dataFrame)
      return;
    this.matrixSlot.state?.grid.invalidate();
    this.transferPanel.invalidateGrid();
  }

  // ---- Export virtual analogs to a make-list --------------------------------------------------

  /** Right-click menu: collect this matrix's (or every matrix's) predicted analogs into the make-list. */
  private buildContextMenu(menu: DG.Menu): void {
    if (!this.matrices.length)
      return;
    const group = menu.group('SAR Matrix');
    if (this.contextCell) {
      const {matrix, ri, ci} = this.contextCell;
      group.item('Add this compound to make list', () => this.makeListPanel.addAnalogToMakeList(matrix, ri, ci));
    }
    const current = this.matrices[Math.min(this.selIndex, this.matrices.length - 1)];
    if (current?.virtualCount)
      group.item(`Add virtual analogs to make-list (${current.label})`, () => this.makeListPanel.addMatrixAnalogsToMakeList([current]));
    if (this.matrices.reduce((n, m) => n + m.virtualCount, 0) > (current?.virtualCount ?? 0)) {
      group.item(`Add virtual analogs to make-list (all ${this.matrices.length} matrices)`,
        () => this.makeListPanel.addMatrixAnalogsToMakeList(this.matrices));
    }
  }


  onPropertyChanged(property: DG.Property | null): void {
    super.onPropertyChanged(property);
    // Some properties don't change the fragmentation, so they must not trigger a full rebuild.
    if (property !== null && RENDER_ONLY_PROPS.includes(property.name)) {
      this.renderMatrixPane();
      // Cards name rows by the Id column, and this branch never reaches render().
      this.refreshSummary();
      return;
    }
    if (property !== null && RERANK_ONLY_PROPS.includes(property.name)) {
      this.reRank();
      return;
    }
    if (property !== null && TRANSFER_ONLY_PROPS.includes(property.name)) {
      this.transferPanel.invalidateTransfers();
      if (this.transferTabActive)
        this.transferPanel.activateTransferTab();
      // The summary's transfer card just went back to "not scanned".
      this.refreshSummary();
      return;
    }
    if (property !== null && FRAGMENT_ONLY_PROPS.includes(property.name) && this.rGroupColumns() !== null)
      return;
    this.scheduleCompute();
  }

  private refreshSummary(): void {
    this.summaryPanel.invalidate();
    if (this.summaryTabActive)
      this.summaryPanel.activateSummaryTab();
  }

  /** Re-rank the assembled matrices and redraw without re-fragmenting. */
  private reRank(): void {
    // A compute in flight already re-ranks against the live properties when it lands, so re-sorting
    // the set it is about to replace would only repaint a stale analysis over its own loader.
    if (this.computing)
      return;
    if (!this.matrices.length) {
      this.scheduleCompute();
      return;
    }
    this.matrices = rankMatrices(this.matrices, this.rankScheme as SarRankScheme, this.higherIsBetter);
    this.hierarchyCache = null;
    // Transfers index into `matrices` by position, so a reorder invalidates them.
    this.transferPanel.invalidateTransfers();
    // No re-activation here: render() below does it once, and doing it twice repaints every depiction.
    this.summaryPanel.invalidate();
    // Both filter frames are stale, not just the cell one: the navigator's Best/Mean fold the activity
    // direction in and Best R is a rank score, so a kept frame would filter the series on numbers the
    // new scheme or direction has already replaced. The cell filter is keyed by matrix id, which
    // re-ranking invalidates as well.
    this.resetFilters();
    this.selIndex = 0;
    this.render();
  }

  /** Why the R-group settings cannot build a matrix, or null when they can or are all empty. */
  private rGroupFault(): string | null {
    const core = this.coreColumnName;
    const axis = this.axisColumnName;
    if (!core && !axis)
      return null;
    if (!core || !axis)
      return core ? 'pick the R-group for the matrix columns' : 'pick the core column';
    if (core === axis)
      return 'the core and the matrix columns are the same column';
    const missing = [core, axis, ...(this.rGroupColumnNames ?? [])].filter((name) => !this.dataFrame?.col(name));
    return missing.length > 0 ? `${missing.map((name) => `"${name}"`).join(', ')} not found in the table` : null;
  }

  /** The R-group columns the matrices are built from, or null to fragment the molecules. */
  private rGroupColumns(): SarFragmentColumns | null {
    if (!this.dataFrame || this.rGroupFault() !== null || !this.coreColumnName)
      return null;
    const core = this.dataFrame.col(this.coreColumnName)!;
    const column = this.dataFrame.col(this.axisColumnName)!;
    const rows = (this.rGroupColumnNames ?? []).filter((name) => name !== column.name && name !== core.name)
      .map((name) => this.dataFrame.col(name)!);
    return {core, rows, column};
  }

  /** Warns once per message: the property panel fires on every keystroke. */
  private report(message: string): void {
    if (message === this.reported)
      return;
    this.reported = message;
    grok.shell.warning(`SAR Matrix: ${message}.`);
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
    // Incomplete R-group settings run nothing rather than fragmenting the whole table; a saved set still restores.
    const fault = this.rGroupFault();
    const restorable = !this.matrices.length && this.matricesData !== '' && this.matricesData !== '[]';
    if (fault !== null && !restorable) {
      this.report(fault);
      if (!this.matrices.length) {
        ui.empty(this.host);
        this.host.appendChild(ui.divText(`SAR Matrix: ${fault}.`));
        this.summaryPanel.showMessage(`SAR Matrix: ${fault}.`);
      }
      return;
    }

    this.reported = '';
    this.computing = true;
    this.analogPanels.clear();
    this.selectedCell = null;
    this.unscalableCount = this.countUnscalable(activity as DG.Column<number>);
    this.assayedCount = (activity as DG.Column<number>).stats.valueCount;
    clearDepictionCaches();
    // Release the Dart-backed grid so the failure path below doesn't leave one repainting forever.
    this.releaseMatrixGrid();
    // Transfers and filter frames key into the matrices about to be replaced; drop them now so a
    // failed compute leaves nothing pointing at rows that no longer exist. Both rebuild lazily.
    this.transferPanel.invalidateTransfers();
    // Invalidation empties the panel, so a user sitting on the Summary tab would watch a blank pane
    // for the whole build unless the loader is put back now.
    this.refreshSummary();
    this.resetFilters();
    ui.empty(this.host);
    this.host.appendChild(ui.loader());
    const progress = DG.TaskBarProgressIndicator.create('Building SAR matrices...');
    try {
      const params: SarMatrixParams = {
        scaling: this.scaling as SCALING_METHODS,
        fragmentCutoff: this.fragmentCutoff,
        predictVirtual: this.predictVirtual,
        predictUnmeasured: this.predictUnmeasured,
        minCompounds: this.minCompounds,
        useMcsAnchors: this.useMcsAnchors,
        grouping: this.grouping as SarGrouping,
        fragmentationLevels: this.fragmentationLevels,
        higherIsBetter: this.higherIsBetter,
        threshold: this.threshold,
        rankScheme: this.rankScheme as SarRankScheme,
        seriesColumn: this.seriesColumnName ? this.dataFrame.col(this.seriesColumnName) : null,
        fragmentColumns: this.rGroupColumns(),
      };
      // A set carried in by a layout is the analysis already; rebuilding it would cost minutes. An
      // empty one is not a set: restoring it would leave the settings that produced it unable to run.
      const carried = this.matricesData && !this.matrices.length ?
        JSON.parse(this.matricesData) as SarMatrix[] : null;
      // Re-ranked on the way in: reRank() reorders without re-serializing, so a carried order goes stale.
      const matrices = carried?.length ? rankMatrices(carried, this.rankScheme as SarRankScheme, this.higherIsBetter) :
        await runSarMatrix(molecules, activity as DG.Column<number>, params);
      // Closed while the workers ran; rendering now would build a grid nothing releases.
      if (this.detached)
        return;
      // Ranked against the LIVE properties, not the ones snapshotted into `params` before the await:
      // a scheme or direction changed while the workers ran would otherwise leave the order and the
      // scores disagreeing with the property panel until the next rank change. rankMatrices recomputes
      // the scores, so this reconciles both.
      this.matrices = rankMatrices(matrices, this.rankScheme as SarRankScheme, this.higherIsBetter);
      this.hierarchyCache = null;
      this.matricesData = JSON.stringify(this.matrices);
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
      this.transferPanel.showMessage(`SAR Matrix failed: ${message}`);
      this.summaryPanel.showMessage(`SAR Matrix failed: ${message}`);
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
  get higherIsBetter(): boolean {
    if (this.activityDirection === DIR_HIGHER)
      return true;
    if (this.activityDirection === DIR_LOWER)
      return false;
    return this.scaling === SCALING_METHODS.MINUS_LG;
  }

  /** Whether a difference in activity units is already a log ratio. A transformed scale is; so is a
   *  raw column the user declared higher-is-better, which is a precomputed pIC50 rather than a
   *  concentration, and where a fold change would understate a four-log swing as "2x". */
  get activityIsLog(): boolean {
    return this.scalingLabel !== 'raw' || this.activityDirection === DIR_HIGHER;
  }


  // ---- Navigator (left pane) ----------------------------------------------------------------

  get scalingLabel(): string {
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
    return this.cellFilter.active;
  }

  /** Whether a cell survives the cell filter. The filter is built from one matrix and the transfer
   *  view draws rows from others, which are outside it and so pass — a keyed lookup alone would miss
   *  every one and blank the grid. */
  cellVisible(matrix: SarMatrix, ri: number, ci: number): boolean {
    if (!this.cellFilter.active || matrix.id !== this.cellPassMatrixId)
      return true;
    return this.cellFilter.passes(`${matrix.id}:${ri}:${ci}`);
  }

  /** One row per cell of the SELECTED matrix (core/substituent as molecules, potency/refs/weight as
   *  numbers) so the platform supplies the filters. Scoped to the on-screen matrix to avoid
   *  fingerprinting cells the user can't see; rebuilds on a matrix switch. */
  private buildStructFrame(): {frame: DG.DataFrame, keys: string[]} {
    const cores: string[] = [];
    const subs: string[] = [];
    const potency: number[] = [];
    const refs: number[] = [];
    const weights: number[] = [];
    const cellKeys: string[] = [];
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
          cellKeys.push(`${matrix.id}:${ri}:${ci}`);
        }
      }
    }
    const coreCol = DG.Column.fromStrings(STRUCT_CORE, cores);
    const subCol = DG.Column.fromStrings(STRUCT_R, subs);
    coreCol.semType = DG.SEMTYPE.MOLECULE;
    subCol.semType = DG.SEMTYPE.MOLECULE;
    return {frame: DG.DataFrame.fromColumns([
      coreCol,
      subCol,
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, STRUCT_POTENCY, potency),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, STRUCT_REFS, refs),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, STRUCT_MW, weights),
    ]), keys: cellKeys};
  }

  /** Measured compounds sharing this cell's row or column — what a prediction here would rest on. */
  observedNeighbours(matrix: SarMatrix, ri: number, ci: number): number {
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

  /** One row per series, carrying its core as a molecule and the numbers its card prints, so the
   *  platform supplies sketcher + histogram filters. */
  private buildNavFrame(): {frame: DG.DataFrame, keys: string[]} {
    const series: string[] = [];
    const cores: string[] = [];
    const best: number[] = [];
    const mean: number[] = [];
    const spread: number[] = [];
    const bestR: number[] = [];
    const compounds: number[] = [];
    const coreCount: number[] = [];
    const level: number[] = [];
    const navKeys: string[] = [];
    for (const matrix of this.matrices) {
      series.push(matrix.label);
      cores.push(matrixCore(matrix));
      // Direction-adjusted so "more potent" is always the larger number on the histogram.
      const top = matrix.realCount ? (this.higherIsBetter ? matrix.maxActivity : -matrix.minActivity) : NaN;
      best.push(top);
      const m = this.meanRealActivity(matrix);
      mean.push(m === null ? NaN : (this.higherIsBetter ? m : -m));
      spread.push(finiteOrNaN(matrix.scores[SarRankScheme.Discontinuity]));
      // Already higher-is-better: preferredScore folds the direction in, so re-applying it would invert it.
      bestR.push(finiteOrNaN(matrix.scores[SarRankScheme.Preferred]));
      compounds.push(matrix.realCount);
      coreCount.push(matrix.rows.length);
      level.push(matrix.level - 1);
      navKeys.push(matrix.id);
    }
    const coreCol = DG.Column.fromStrings(NAV_CORE, cores);
    coreCol.semType = DG.SEMTYPE.MOLECULE;
    return {frame: DG.DataFrame.fromColumns([
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
    ]), keys: navKeys};
  }

  private passesFilter(matrix: SarMatrix): boolean {
    return this.navFilter.passes(matrix.id);
  }

  /** The score block on the right of a navigator card. Shared by series and transfer cards so both
   *  lists read as one column of scores. */
  cardScoreBox(lines: {value: string, label: string}[], tip: () => string): HTMLElement {
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
      const raw = matrix.scores[scheme];
      // No column here has the two observations a mean needs, so there is no score to show.
      if (raw === undefined || !Number.isFinite(raw)) {
        return {lines: [{value: '—', label: 'best R'}],
          tip: 'No substituent in this series is measured on two cores, so it has no preferred-R score.'};
      }
      return {lines: [{value: this.formatActivity(this.higherIsBetter ? raw : -raw), label: 'best R'}],
        tip: `Best mean potency of any single substituent across the cores (${unit}, ${this.scalingLabel} scale).`};
    }
    // Remaining scheme is Discontinuity.
    return {lines: [{value: (matrix.scores[SarRankScheme.Discontinuity] ?? 0).toFixed(2), label: 'spread'}],
      tip: `Largest activity spread within a single core, on the ${this.scalingLabel} scale — ` +
        'how discontinuous the SAR is.'};
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

  /**
   * The one hierarchy every part of the viewer reads: parent index, lineage root id and displayed tier
   * per matrix, in ranked order.
   *
   * Cached because `matrixParents` falls back to a quadratic containment scan, and because deriving a
   * second hierarchy elsewhere would disagree with the navigator's on exactly the runs that take that
   * fallback. Dropped whenever `matrices` is replaced or reordered — the indices are positional.
   */
  private get hierarchy(): {parents: number[], roots: string[], tiers: number[]} {
    if (this.hierarchyCache !== null)
      return this.hierarchyCache;
    const parents = this.matrixParents();
    const children = this.navChildren(parents);
    // Tiers as the list actually shows them, not as the fragmentation built them: a matrix folded from
    // groups that all assembled into nothing has no tier under it, and calling it L2 sends the reader
    // looking for one.
    const tierOf = (i: number): number =>
      1 + (children.get(i) ?? []).reduce((deepest, child) => Math.max(deepest, tierOf(child)), 0);
    const roots = this.matrices.map((_matrix, i) => {
      let root = i;
      for (let guard = 0; parents[root] >= 0 && guard <= parents.length; guard++)
        root = parents[root];
      return this.matrices[root].id;
    });
    return this.hierarchyCache = {parents, roots, tiers: this.matrices.map((_m, i) => tierOf(i))};
  }

  /** Lineage root id per matrix — what a finding is counted under, so one compound set re-cut at three
   *  tiers cannot read as three confirmations. */
  get matrixRoots(): string[] {
    return this.hierarchy.roots;
  }

  /** Displayed tier per matrix: the L1/L2/L3 the navigator cards print. */
  get matrixTiers(): number[] {
    return this.hierarchy.tiers;
  }

  /**
   * The fragment column every matrix varies, or null when these came from fragmentation. Only when
   * this is set does one substituent label mean the same thing in two series — a fragmented matrix
   * names its positions from its own decomposition.
   *
   * Read off the matrices rather than the properties: `readFragmentColumns` warns as a side effect,
   * and a set carried in by a layout never ran the decomposition the properties describe.
   */
  get axisRole(): string | null {
    const axis = this.axisColumnName;
    return this.coreColumnName && axis && this.matrices.length > 0 &&
      this.matrices.every((m) => m.positions[0] === axis) ? axis : null;
  }

  /** Select the compounds carrying one value of one component column. An offset is a statement about
   *  the subset of the table that carries a value, and the only way to check one is to look at it. */
  selectRoleValue(role: string, value: string): void {
    const df = this.dataFrame;
    const column = df?.col(role);
    if (df == null || column == null)
      return;
    // The level came off the decomposition, which standardizes every cell it reads, so the raw text
    // of a molblock or a non-canonically written SMILES never equals it. Standardized per distinct
    // category rather than per row: a degrader set has hundreds of warheads over thousands of rows.
    const wanted = standardizeFragment(value).value;
    const matches = new Set(column.categories
      .filter((category) => standardizeFragment(category).value === wanted));
    df.selection.init((i) => matches.has(column.getString(i)));
  }

  get transferSummary(): {scanned: boolean, count: number} {
    return this.transferPanel.scanSummary;
  }

  /** The column the cores came from, under the same condition as {@link axisRole}: a decomposition
   *  the user handed over names its own parts, and a fragmented one has no name to give. */
  get coreRole(): string | null {
    return this.axisRole === null ? null : this.coreColumnName;
  }

  /** Whether a matrix IS one core, so its fitted mean compares cores rather than groupings. A series
   *  column or an empty row axis groups by something else, and then the cores are not the series. */
  get coresAreSeries(): boolean {
    return this.matrices.length > 0 &&
      this.matrices.every((m) => new Set(m.rows.map((r) => r.coreSmiles)).size === 1);
  }

  /** Fragment columns that could run across the top, the current axis included. */
  get roleColumns(): string[] {
    const names = [this.axisColumnName, ...(this.rGroupColumnNames ?? [])];
    return [...new Set(names)].filter((name) =>
      name && name !== this.coreColumnName && this.dataFrame?.col(name) != null);
  }

  /** Transpose the analysis onto another fragment column. Through `setOptions` rather than the field:
   *  a plain assignment raises no property change, so nothing would recompute. */
  setColumnAxis(name: string): void {
    this.setOptions({axisColumnName: name});
  }

  private countUnscalable(activity: DG.Column<number>): number {
    if (this.scaling === SCALING_METHODS.NONE)
      return 0;
    const scaled = scaleActivity(activity, this.scaling as SCALING_METHODS).getRawData();
    let n = 0;
    for (let i = 0; i < activity.length; i++) {
      if (!activity.isNone(i) && !Number.isFinite(scaled[i]))
        n++;
    }
    return n;
  }

  /** The matrix's scaffold with every attachment made alike, canonicalized. Matrices varying
   *  different positions of one scaffold key on the same structure and differ only in which
   *  attachment is inherited, so normalizing the marks away collapses them onto a single key. */
  private scaffoldKeyOf(matrix: SarMatrix): string {
    const source = matrix.siteKey || matrix.rows[0]?.coreSmiles || '';
    if (!source)
      return '';
    return this.scaffoldKeys.getOrCreate(source, () => {
      const bare = source.replace(/\[\d+\*\]|\[\*:\d+\]/g, '[*]');
      let mol = null;
      let key = bare;
      try {
        mol = getRdKitModule().get_mol(bare);
        if (mol?.is_valid())
          key = mol.get_smiles();
      } catch {
        key = bare;
      } finally {
        mol?.delete();
      }
      return key;
    });
  }

  /**
   * Folds every scaffold onto the most general one present that contains it. One ring system appears
   * under several keys when different positions are pinned rather than left open, so keying on
   * structure alone would split it into separate groups. "More general" means strictly more open
   * positions, which prunes the pairwise test and keeps the fold acyclic.
   *
   * Each key is parsed once as both query and molecule: parsing per pair costs tens of thousands of
   * main-thread RDKit calls on a few hundred scaffolds.
   */
  private scaffoldRepresentatives(keys: string[]): Map<string, string> {
    const distinct = [...new Set(keys.filter((k) => k))];
    const signature = [...distinct].sort().join('\n');
    if (this.scaffoldRepsFor === signature)
      return this.scaffoldReps;
    const rdkit = getRdKitModule();
    const openness = new Map(distinct.map((key) => [key, (key.match(/\*/g) ?? []).length]));
    const ordered = [...distinct].sort((a, b) =>
      openness.get(b)! - openness.get(a)! || (a < b ? -1 : a > b ? 1 : 0));

    const queries = new Map<string, RDMol | null>();
    const targets = new Map<string, RDMol | null>();
    const parse = (cache: Map<string, RDMol | null>, key: string, asQuery: boolean): RDMol | null => {
      if (!cache.has(key)) {
        let mol: RDMol | null = null;
        try {
          mol = asQuery ? rdkit.get_qmol(key) : rdkit.get_mol(key);
          if (!asQuery && !mol?.is_valid()) {
            mol?.delete();
            mol = null;
          }
        } catch {
          mol = null;
        }
        cache.set(key, mol);
      }
      return cache.get(key) ?? null;
    };
    const contains = (query: string, target: string): boolean => {
      const q = parse(queries, query, true);
      const t = parse(targets, target, false);
      if (!q || !t)
        return false;
      try {
        const match = JSON.parse(t.get_substruct_match(q));
        return Array.isArray(match?.atoms) && match.atoms.length > 0;
      } catch {
        return false;
      }
    };

    const reps = new Map<string, string>();
    try {
      for (let i = 0; i < ordered.length; i++) {
        const key = ordered[i];
        let rep = key;
        for (let j = 0; j < i; j++) {
          const candidate = ordered[j];
          if (openness.get(candidate)! <= openness.get(key)!)
            continue;
          if (contains(candidate, key)) {
            rep = reps.get(candidate) ?? candidate;
            break;
          }
        }
        reps.set(key, rep);
      }
    } finally {
      for (const mol of [...queries.values(), ...targets.values()])
        mol?.delete();
    }
    this.scaffoldReps = reps;
    this.scaffoldRepsFor = signature;
    return reps;
  }

  /** Row standing above the matrices that share one scaffold. Not a matrix — it has no rows of its
   *  own — so it selects nothing and only folds the set away. Its every attachment draws unlabelled,
   *  since at this level no one position is the axis. */
  private buildScaffoldCard(scaffold: string, count: number): HTMLElement {
    const core = this.lazyCoreDepiction(coreDepictionBlock(scaffold, -1), CARD_CORE_W, CARD_CORE_H);
    core.classList.add('chem-sar-card-core');
    const body = ui.divV([
      ui.divText('Scaffold', 'chem-sar-card-name'),
      ui.divText(`${count} series`, 'chem-sar-card-desc'),
    ], 'chem-sar-card-body');
    const twisty = ui.iconFA(this.collapsedScaffolds.has(scaffold) ? 'chevron-right' : 'chevron-down',
      (e: MouseEvent) => {
        e.stopPropagation();
        if (!this.collapsedScaffolds.delete(scaffold))
          this.collapsedScaffolds.add(scaffold);
        this.updateNavVisibility();
        const open = !this.collapsedScaffolds.has(scaffold);
        twisty.classList.toggle('fa-chevron-down', open);
        twisty.classList.toggle('fa-chevron-right', !open);
      });
    twisty.classList.add('chem-sar-card-twisty');
    ui.tooltip.bind(twisty, () => this.collapsedScaffolds.has(scaffold) ?
      `Show the ${count} series on this scaffold` : 'Hide the series on this scaffold');
    const card = ui.divH([twisty, core, body], 'chem-sar-card chem-sar-scaffold-card');
    ui.tooltip.bind(card, () => 'One scaffold, one series per position it varies. Each keeps its own ' +
      'column axis — the series are listed rather than merged, since substituents at different ' +
      'positions do not belong in one column.');
    return card;
  }

  /** One selectable matrix card: aligned core, a descriptor line, and the rank score. `depth`
   *  indents a folded matrix; `folded` is how many are folded into this one, giving its expander. */
  private buildCard(matrix: SarMatrix, index: number, depth = 0, folded = 0, shown = 1,
    onToggle?: () => void): HTMLElement {
    // The folded count rides on the level badge instead of the descriptor, which ran out of room for it.
    const desc = `${matrix.rows.length} cores · ${matrix.realCount} cpd`;
    const core = this.lazyCoreDepiction(matrixCore(matrix), CARD_CORE_W, CARD_CORE_H);
    core.classList.add('chem-sar-card-core');
    // Level badge, since branches bottom out at different depths so indentation alone doesn't say
    // which level a card sits at.
    const level = ui.divText(folded ? `L${shown}·${folded}` : `L${shown}`, 'chem-sar-card-level');
    // Carried by the badge alone. The descriptor has one line to hold cores and compounds, and a third
    // term there is what pushed it into an ellipsis.
    const gap = this.subSeriesGap(matrix);
    if (gap !== null)
      level.classList.add('chem-sar-card-level-partial');
    ui.tooltip.bind(level, () => (shown === 1 ?
      'Level 1 — no finer series sit under this one. Each of its rows is a core with its substituents.' :
      `Level ${shown} — a coarser matrix, holding the compounds of the level-${shown - 1} matrices ` +
      'below it, whose cores agree one further cut deeper.') +
      (folded ? ` The number after the level is how many it folds in: ${folded}.` : '') +
      (gap !== null ? ` ${gap.tip}` : ''));
    const body = ui.divV([
      ui.divH([ui.divText(matrix.label, 'chem-sar-card-name'), level], 'chem-sar-card-title'),
      ui.divText(desc, 'chem-sar-card-desc'),
    ], 'chem-sar-card-body');
    const sc = this.cardScore(matrix);
    const scoreBox = this.cardScoreBox(sc.lines, () => sc.tip);
    // Held to the body's two rows so a single-line score (Best R, Spread) still starts level with the
    // series name instead of floating between the two text rows.
    scoreBox.classList.add('chem-sar-card-score-2row');
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
    this.applyMatrixSelection(index, '');
  }

  /** Move the navigator highlight and rebuild the matrix pane, varying `position` (empty = all
   *  columns). Runs even on the already-selected index: a landing may need a different vary axis,
   *  and it needs a freshly built grid to scroll into. */
  private applyMatrixSelection(index: number, position: string): void {
    this.navCards[this.selIndex]?.classList.remove('selected');
    this.selIndex = index;
    this.navCards[index]?.classList.add('selected');
    this.varyPosition = position;
    // Cell filter is scoped to one matrix; rebuild for the new one.
    this.resetStructFilter();
    this.renderMatrixPane();
  }

  /** Whether the Summary tab is the one on screen. */
  get summaryTabActive(): boolean {
    return this.tabs.currentPane?.name === TAB_SUMMARY;
  }

  /** Bring a matrix on screen: the SAR Matrix tab, its navigator card unfolded and scrolled to, and
   *  the pane rebuilt for it. Resolved by id, not by a held index — reRank reorders `matrices`.
   *  `position` pins the Vary axis when the caller needs two specific columns on screen together.
   *  Returns the matrix's index, or -1 when it is not in the current analysis. */
  revealMatrix(matrix: SarMatrix, position?: string): number {
    const index = this.matrices.findIndex((m) => m.id === matrix.id);
    if (index < 0) {
      grok.shell.info('That series is not in the current analysis.');
      return -1;
    }
    this.tabs.currentPane = this.tabs.getPane(TAB_MATRIX);
    this.unfoldNavTo(index);
    this.applyMatrixSelection(index,
      position !== undefined && this.matrices[index].positions.includes(position) ? position : '');
    this.scrollNavTo(index);
    return index;
  }

  /** Land on one cell. The scroll and the open are deferred: the grid is created in this tick and has
   *  no layout yet, so scrollToCell would measure against zero. */
  revealCell(matrix: SarMatrix, ri: number, ci: number, position?: string): void {
    const index = this.revealMatrix(matrix, position);
    if (index < 0)
      return;
    const target = this.matrices[index];
    window.setTimeout(() => {
      const state = this.matrixGrid;
      if (state === null || this.detached)
        return;
      this.scrollPaneTo(state, ri, ci);
      this.openCell(target, ri, ci, state.grid, null);
    }, 0);
  }

  /** Bring one navigator card into view by moving the navigator's own scroller. */
  private scrollNavTo(index: number): void {
    const card = this.navCards[index];
    const list = card?.closest('.chem-sar-nav-list');
    if (card === undefined || !(list instanceof HTMLElement))
      return;
    scrollWithin(list, card);
  }

  /** Open every collapsed ancestor and scaffold header hiding a matrix's card. A navigator filter is
   *  the user's, so it is reported rather than cleared. */
  private unfoldNavTo(index: number): void {
    for (let p = this.navParents[index]; p >= 0; p = this.navParents[p])
      this.collapsed.delete(this.matrices[p].id);
    let root = index;
    for (let guard = 0; this.navParents[root] >= 0 && guard <= this.navParents.length; guard++)
      root = this.navParents[root];
    const key = this.navGroupOfRoot.get(root);
    if (key !== undefined)
      this.collapsedScaffolds.delete(key);
    this.updateNavVisibility();
    if (!this.passesFilter(this.matrices[index])) {
      grok.shell.info(`${this.matrices[index].label} is hidden by the series filter — clear it to see ` +
        'the card in the list.');
    }
  }

  /** Scroll the matrix pane to one cell. The displayed column carries a synthetic key, so the matrix
   *  column index has to be mapped back through the pane's own column list. */
  private scrollPaneTo(state: MatrixGridState, ri: number, ci: number): void {
    const paneRowIdx = state.rows.findIndex((row) => row.rowIndex === ri);
    if (paneRowIdx < 0)
      return;
    const colIdxs = state.rows[paneRowIdx].colIdxs;
    for (const [key, idx] of state.colKeyToIdx) {
      if (colIdxs[idx] === ci) {
        state.grid.scrollToCell(key, paneRowIdx);
        return;
      }
    }
  }

  /** Bring a tab to the front, for a panel handing work to another one. */
  showTab(name: string): void {
    this.tabs.currentPane = this.tabs.getPane(name);
  }

  addCellsToMakeList(cells: MatrixCellRef[], emptyMessage: string): void {
    this.makeListPanel.addCellsToMakeList(cells, emptyMessage);
  }

  get hostRowCount(): number {
    return this.dataFrame?.rowCount ?? 0;
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


  /** Whether the SAR Transfer tab is the one on screen. */
  get transferTabActive(): boolean {
    return this.tabs.currentPane?.name === TAB_TRANSFER;
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
    // Named by its tooltip rather than a label, keeping the header one row of aligned controls. Rebuilt
    // per hover to mark the active scheme, so not a setTooltip string; bound to the editor because the
    // select covers the root and swallows its hover.
    rankInput.caption = '';
    ui.tooltip.bind(rankInput.input, () => {
      const unit = this.activityColumnName || 'activity';
      const mark = (scheme: string): string => scheme === this.rankScheme ? '▸ ' : '   ';
      return ui.divV([
        ui.divText('Rank by — how the navigator orders the series:'),
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

    const parents = this.hierarchy.parents;
    const roots = parents.filter((p) => p < 0).length;
    const deepest = this.matrices.reduce((m, matrix) => Math.max(m, matrix.level), 2) - 1;
    // The breakdown is read once, so it hangs off a hover rather than costing a row above the list.
    const countBadge = ui.divText(
      `${this.matrices.length} matri${this.matrices.length === 1 ? 'x' : 'ces'}`, 'chem-sar-nav-badge');
    ui.tooltip.bind(countBadge, () =>
      `${this.matrices.length} matrices in ${roots} famil${roots === 1 ? 'y' : 'ies'} · ` +
      `${deepest} level${deepest === 1 ? '' : 's'}`);
    const list = ui.div([], 'chem-sar-nav-list');
    const matchCount = ui.divText('', 'chem-sar-nav-matches');
    const navIdleTip = 'Filter series by core structure, potency, SAR spread, size';
    const navIcon = ui.icons.filter(() => {
      ui.showPopup(ui.div(this.navFilter.root(), 'chem-sar-struct-filters'), navIcon, {vertical: true});
    }, navIdleTip);
    navIcon.classList.add('chem-sar-struct-icon');
    const refill = (): void => {
      this.updateNavVisibility();
      const filtered = this.navFilter.active;
      const hits = this.matrices.filter((matrix) => this.passesFilter(matrix)).length;
      matchCount.innerText = filtered ? `${hits} of ${this.matrices.length} match` : '';
      this.markFilterIcon(navIcon, filtered,
        `${hits} of ${this.matrices.length} series shown — filter active. Click to edit or clear.`, navIdleTip);
    };
    this.navMatchCount = refill;

    const header = ui.divV([
      ui.divH([rankInput.root, countBadge, navIcon], 'chem-sar-nav-controls'),
      matchCount,
    ], 'chem-sar-nav-header');
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
    const tiers = this.matrixTiers;

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
      const card = this.buildCard(this.matrices[i], i, depth, descendants(i), tiers[i], () => {
        if (!this.collapsed.delete(id))
          this.collapsed.add(id);
        this.updateNavVisibility();
      });
      this.navCards[i] = card;
      list.appendChild(card);
      for (const child of children.get(i) ?? [])
        emit(child, depth + 1);
    };
    // Roots keyed on one scaffold differ only in which of its positions they vary, so they are listed
    // under a shared header instead of scattered through the ranking. Bucket order follows the first
    // root of each bucket, keeping the ranking's order of families intact.
    this.navGroups = [];
    this.navGroupOfRoot.clear();
    const rootKeys = new Map<number, string>();
    parents.forEach((p, root) => {
      if (p < 0)
        rootKeys.set(root, this.scaffoldKeyOf(this.matrices[root]));
    });
    // Folded first, so a ring system whose positions are pinned differently reads as one scaffold.
    const reps = this.scaffoldRepresentatives([...rootKeys.values()]);
    const byScaffold = new Map<string, number[]>();
    for (const [root, rawKey] of rootKeys) {
      const key = reps.get(rawKey) ?? rawKey;
      const bucket = byScaffold.get(key);
      if (bucket === undefined)
        byScaffold.set(key, [root]);
      else
        bucket.push(root);
    }
    for (const [key, members] of byScaffold) {
      // A lone root needs no header: there is nothing to group it with.
      if (members.length < 2 || !key) {
        for (const root of members)
          emit(root, 0);
        continue;
      }
      const header = this.buildScaffoldCard(key, members.length);
      list.appendChild(header);
      this.navGroups.push({header, members});
      for (const root of members) {
        this.navGroupOfRoot.set(root, key);
        emit(root, 1);
      }
    }
  }

  /** Show/hide the already-built cards per collapse state and filter (pure display toggles). A parent
   *  that fails the filter still shows when a descendant passes; collapse is ignored while filtering. */
  private updateNavVisibility(): void {
    const parents = this.navParents;
    const filtering = this.navFilter.active;
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
    // A folded scaffold header hides the families under it, whatever their own collapse state.
    const underFoldedScaffold = (i: number): boolean => {
      let root = i;
      for (let guard = 0; parents[root] >= 0 && guard <= parents.length; guard++)
        root = parents[root];
      const key = this.navGroupOfRoot.get(root);
      return key !== undefined && this.collapsedScaffolds.has(key);
    };
    this.navCards.forEach((card, i) => {
      if (!card)
        return;
      const visible = filtering ? anyHit(i) : (!hiddenByCollapse(i) && !underFoldedScaffold(i));
      card.style.display = visible ? '' : 'none';
    });
    for (const group of this.navGroups)
      group.header.style.display = !filtering || group.members.some((m) => anyHit(m)) ? '' : 'none';
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
  buildPaneGrid(rows: PaneRow[], columns: PaneColumn[], slot: PaneGridSlot): HTMLElement {
    this.cellW = this.fitCellWidth(columns.length);
    // Group boundaries computed once here; rescanning in `onCellRender` would be quadratic per repaint.
    const firstOfGroup = new Set<number>();
    columns.forEach((column, i) => {
      if (i === 0 || columns[i - 1].position !== column.position)
        firstOfGroup.add(i);
    });
    const axisIsChemical = columns.some((column) => column.substSmiles.includes('[*:'));

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
    const onlyMatrix = rows.length > 0 && rows.every((row) => row.matrix === rows[0].matrix) ?
      rows[0].matrix : null;
    // The H-capped shared core draws no attachment point, so it is only the fallback across matrices.
    const headerDepiction = onlyMatrix !== null ? matrixCore(onlyMatrix) : headerCore;
    const state: MatrixGridState = {grid, df, rows, columns, colKeyToIdx, firstOfGroup, axisIsChemical,
      headerDepiction, paneTemplate, rowTemplates: new Array(rows.length).fill(null)};

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
          this.painter.paintHeader(g, b, name, state);
        else {
          const ri = grid.gridRowToTable(c.gridRow);
          if (ri >= 0 && ri < rows.length) {
            if (name === 'Core')
              this.painter.paintCore(g, b, rows[ri], ri, state);
            else {
              const idx = colKeyToIdx.get(name);
              if (idx !== undefined)
                this.painter.paintBodyCell(g, b, rows[ri], ri, idx, state);
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
    if (!this.cellVisible(paneRow.matrix, paneRow.rowIndex, ci))
      return;
    this.openCell(paneRow.matrix, paneRow.rowIndex, ci, state.grid,
      this.lastGridMouseEvent ?? new MouseEvent('click'));
  }

  /** Open a cell: context panel, host current row, and the selected cell the cart acts on. `event` is
   *  the click that reached the grid; null arrives from a landing, which moves the current row without
   *  rewriting a selection the user made elsewhere. The grid is passed because both panes' grids come
   *  from buildPaneGrid and each must repaint its own. */
  private openCell(matrix: SarMatrix, ri: number, ci: number, grid: DG.Grid,
    event: MouseEvent | null): void {
    const cell = matrix.cells[ri][ci];
    // A threshold-blanked cell draws empty and isn't selectable.
    if (cell.kind === 'empty' || cell.value === null)
      return;
    this.selectedCell = {matrix, ri, ci};
    grok.shell.windows.showContextPanel = true;
    // A real compound becomes the current row and extends selection; a virtual analog has no row.
    if (cell.molIdx !== null) {
      this.dataFrame.currentRowIdx = cell.molIdx;
      if (event !== null)
        this.dataFrame.selection.handleClick((i) => i === cell.molIdx, event, true);
    } else
      this.dataFrame.currentRowIdx = -1;
    // Assembled cell opens the Molecule context (SAR context gated into the "SAR analysis" pane);
    // an unassembled virtual falls back to the standalone SAR panel.
    if (cell.smiles)
      this.showMoleculeContext(cell.smiles, () => this.makeListPanel.buildCellPanel(matrix, ri, ci));
    else
      grok.shell.o = this.makeListPanel.buildCellPanel(matrix, ri, ci);
    grid.invalidate();
  }

  /** Collect an explicit set of cells, for the transfer pane's cart. */
  /** The cart both panes carry. Adding sits with the grid rather than only in the make-list tab: the
   *  cell to add is picked here, and a control that lives one tab away asks the user to leave what
   *  they are looking at. Both grids come from `buildPaneGrid`, so a click in either sets the same
   *  selected cell and one icon serves both. Left enabled with no selection so the reason lands as a
   *  message instead of a dead icon. */
  cartIcon(): HTMLElement {
    const icon = ui.iconFA('cart-plus', () => this.makeListPanel.addSelectedToMakeList());
    // Second class so the cart and the filter funnel, which share the icon styling, stay separable.
    icon.classList.add('chem-sar-struct-icon', 'chem-sar-cart-icon');
    ui.tooltip.bind(icon, () => {
      const cell = this.selectedCell;
      if (cell === null)
        return 'Add to make list — click a cell first.';
      const target = cell.matrix.cells[cell.ri][cell.ci];
      if (target.smiles === null)
        return 'Add to make list — the selected cell has no structure.';
      const kind = target.kind === 'virtual' ? 'virtual analog' :
        target.kind === 'unmeasured' ? 'untested compound' : 'synthesized compound';
      return `Add ${cell.matrix.label} · ${cell.matrix.rows[cell.ri].label} × ` +
        `${cell.matrix.columns[cell.ci].position} (${kind}) to the make list`;
    });
    return icon;
  }

  /** Open the platform's molecule context for a structure. The builder, when given, is gated into it
   *  as the "SAR analysis" pane — the accordion hook above resolves it by SMILES.
   *
   *  A click that moves the source table's current row is followed by the platform restoring the
   *  current object to the one the table view held before, which would drop the molecule just
   *  opened. Take it back the first time that happens, then stop listening: staying armed would let
   *  a later, deliberate selection elsewhere get yanked back here. */
  showMoleculeContext(smiles: string, build?: AnalogPanelBuilder): void {
    grok.shell.windows.showContextPanel = true;
    if (build)
      this.analogPanels.set(smiles, build);
    const value = DG.SemanticValue.fromValueType(smiles, DG.SEMTYPE.MOLECULE);
    this.contextGuard?.unsubscribe();
    // Compared by SMILES, not by identity: reading the current object hands back a fresh wrapper
    // around the same value, so an identity test would call our own assignment a replacement.
    const sub = grok.events.onCurrentObjectChanged.subscribe(() => {
      const current = grok.shell.o;
      if (current instanceof DG.SemanticValue && String(current.value) === smiles)
        return;
      this.dropContextGuard(sub);
      grok.shell.o = value;
    });
    this.contextGuard = sub;
    window.setTimeout(() => this.dropContextGuard(sub), CONTEXT_GUARD_MS);
    grok.shell.o = value;
  }

  /** Retire one guard, leaving a newer one that has since replaced it alone. */
  private dropContextGuard(sub: {unsubscribe(): void}): void {
    sub.unsubscribe();
    if (this.contextGuard === sub)
      this.contextGuard = null;
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
    if (cell.kind === 'impossible') {
      ui.tooltip.show(ui.divText('This row has no attachment point for this substituent.'), x, y);
      return true;
    }
    // A threshold-blanked cell must hover as empty too, or the hidden value returns under the cursor.
    if (cell.kind === 'empty' || cell.value === null || !this.cellVisible(paneRow.matrix, paneRow.rowIndex, ci))
      return false;
    ui.tooltip.show(ui.divText(this.cellTooltipText(cell)), x, y);
    return true;
  }

  /** The id chip's text for an observed cell; null for virtual analogs (no source row to identify). */
  cellIdText(cell: SarMatrixCell): string | null {
    // An untested compound has a registration id like any other; showing it is what says "you have
    // this one already".
    if (!this.idColumnName || (cell.kind !== 'real' && cell.kind !== 'unmeasured') || cell.molIdx === null)
      return null;
    const column = this.dataFrame.col(this.idColumnName);
    if (column === null || column.isNone(cell.molIdx))
      return null;
    return String(column.get(cell.molIdx));
  }

  /** Cell hover text: predicted (+ support) for a virtual analog, observed for a real one, and for an
   *  untested compound both — the compound is in hand, only the number is estimated. */
  private cellTooltipText(cell: SarMatrixCell): string {
    const value = cell.value!;
    const support = cell.support ?? 0;
    const predicted = `Predicted ${this.formatActivity(value)} · support n=${support}` +
      `${support <= 1 ? ' (low)' : ''}`;
    if (cell.kind === 'virtual')
      return predicted;
    if (cell.kind === 'unmeasured')
      return `In the dataset, never tested — ${predicted.toLowerCase()}`;
    return `Observed ${this.formatActivity(value)}`;
  }


  /** Format an activity value: no decimals at ≥100, one decimal below. */
  formatActivity(value: number): string {
    return value >= 100 ? value.toFixed(0) : value.toFixed(1);
  }


  /** Most potent OBSERVED value among a row's displayed cells, or null; predictions excluded. */
  bestRowValue(paneRow: PaneRow): number | null {
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
   * How many of a matrix's compounds its sub-series reach, when they do not reach all of them.
   *
   * Sub-series are finer views derived from this matrix, not a division of it: a folded group that
   * assembled into nothing still contributed its compounds here, so expanding shows fewer than the
   * parent holds. Null when there are no sub-series or they cover it fully.
   */
  private subSeriesGap(matrix: SarMatrix): {covered: number, tip: string} | null {
    const children = this.matrices.filter((m) => m.parentId === matrix.id);
    // Nothing to reconcile with no sub-series: the card reads L1 and promises none.
    if (!children.length)
      return null;
    const covered = new Set<number>(children.flatMap((child) => [...observedMolecules(child)]));
    if (covered.size >= matrix.realCount)
      return null;
    return {
      covered: covered.size,
      tip: `${children.length} sub-series cover ${covered.size} of this matrix's ${matrix.realCount} ` +
        `compounds. The other ${matrix.realCount - covered.size} sit in groups too thin to form a ` +
        'matrix of their own, so they are counted here but appear in none of the series below.',
    };
  }

  /** Compact summary chips for the current matrix, full detail on hover. */
  private buildChips(matrix: SarMatrix): HTMLElement {
    // Reports what is on screen, not what the matrix holds (a threshold makes them differ), and is
    // kept so the threshold can retitle it without rebuilding the bar.
    const {rows, cols} = this.visibleDims(matrix);
    const full = `${matrix.rows.length} cores × ${matrix.columns.length} substituents`;
    this.paneDimsChip = chipBadge(`${rows}×${cols}`, rows === matrix.rows.length && cols === matrix.columns.length ?
      full : `${rows} cores × ${cols} substituents shown, filtered from ${full}`);
    const items = [
      chipBadge(`${matrix.realCount} cpd`, `${matrix.realCount} observed compounds`),
      this.paneDimsChip,
    ];
    if (matrix.realCount) {
      items.push(chipBadge(`${this.formatActivity(matrix.minActivity)}–${this.formatActivity(matrix.maxActivity)}`,
        `Activity range across the matrix, on the ${this.scalingLabel} scale`));
    }
    if (matrix.virtualCount)
      items.push(chipBadge(`${matrix.virtualCount} virtual`, `${matrix.virtualCount} predicted (virtual) analog(s)`));
    // Compounds in hand with no assay value are the cheapest thing a matrix can point at, so they are
    // called out rather than left to be inferred from the cells that carry a '~' on a solid frame.
    const untested = matrix.cells.reduce((n, row) =>
      n + row.filter((c) => c.kind === 'unmeasured').length, 0);
    if (untested)
      items.push(chipBadge(`${untested} untested`, `${untested} compound(s) the dataset holds with no activity value`));
    const gap = this.subSeriesGap(matrix);
    if (gap !== null)
      items.push(chipBadge(`${gap.covered}/${matrix.realCount} in sub-series`, gap.tip, 'chem-sar-chip-partial'));
    // Whether the additive model this matrix's predictions rest on actually holds here. Cross-validated
    // out of sample, so it is the one number that can contradict the predictions rather than restate
    // them, and a low value marks a matrix whose virtual cells should not be trusted.
    const conf = matrix.confidence;
    if (conf) {
      const tip = `Leave-one-out R² ${conf.r2.toFixed(2)}, RMSE ${this.formatActivity(conf.rmse)} over ` +
        `${conf.n} of ${conf.total} measured cells. Near 1 the substituent effects add up and the ` +
        'predictions are trustworthy; near 0 or below they do not, and this matrix\'s analogs are guesses.';
      items.push(chipBadge(`R² ${conf.r2.toFixed(2)}`, tip, conf.r2 < 0.5 ? 'chem-sar-chip-partial' : ''));
    }
    const idx = this.matrices.indexOf(matrix);
    // The panel keeps them sorted by correlation, so the first match is the strongest.
    const involving = this.transferPanel.transfersInvolving(idx);
    if (involving.length) {
      const best = involving[0];
      items.push(chipBadge(`r ${best.correlation.toFixed(2)}`,
        'Strongest SAR-transfer correlation between this series and another', 'chem-sar-chip-transfer'));
    }
    return ui.divH(items, 'chem-sar-chips');
  }

  /** A filter popup, pinned to the viewer's left edge so it lands over the navigator. Where the anchor
   *  puts it, it covers the cells the filter acts on — the one thing that has to stay visible while
   *  the filter is being adjusted. */
  private showFilterPopup(content: HTMLElement, icon: HTMLElement): void {
    ui.showPopup(content, icon, {vertical: true});
    const popup = content.closest('.d4-popup-host');
    if (popup instanceof HTMLElement)
      popup.style.left = `${Math.round(this.host.getBoundingClientRect().left)}px`;
  }

  private buildMatrixPane(matrix: SarMatrix): HTMLElement {
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
    // Named by its tooltip rather than a label, which keeps the header one row of aligned controls.
    labelInput.caption = '';
    labelInput.setTooltip('Caption — annotate each substituent column with a metric: its mean potency ' +
      '(μ) or molecular weight (MW). Columns keep their order; only the caption is added.');
    controls.push(labelInput.root);

    controls.push(this.cartIcon());

    // All filters behind one icon; the sketchers would dwarf the bar inline.
    const filterIcon = ui.icons.filter(() => {
      this.showFilterPopup(ui.div(this.cellFilter.root(), 'chem-sar-struct-filters'), filterIcon);
    }, 'Filter cells by potency, reference points, core and R-group');
    filterIcon.classList.add('chem-sar-struct-icon', 'chem-sar-filter-icon');
    this.matrixFilterIcon = filterIcon;
    controls.push(filterIcon);
    // The controls ride in the title row instead of a row of their own: a second bar costs vertical
    // space the matrix itself needs, and the pane is often only a few rows tall.
    const controlBar = ui.divH(controls, 'chem-sar-control-bar');
    controlBar.classList.add('chem-sar-control-inline');
    const infoBar = ui.divH([ui.divText(matrix.label, 'chem-sar-main-title'), this.buildChips(matrix), controlBar],
      'chem-sar-main-bar');

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
    return ui.divV([infoBar, gridHost, emptyNote], 'chem-sar-main');
  }

  /** Why the analysis produced nothing, named by the grouping that produced it. */
  noMatricesMessage(): string {
    return this.grouping === SarGrouping.Site ?
      'No SAR matrices found. Try raising the fragment cutoff, or switch grouping to Similarity.' :
      'No SAR matrices found. Try lowering the clustering threshold or raising the fragment cutoff.';
  }

  private render(): void {
    // Preserve the navigator scroll position across the rebuild.
    const prevNav = this.host.querySelector('.chem-sar-nav-list');
    const navScroll = prevNav instanceof HTMLElement ? prevNav.scrollTop : 0;
    // Release the old grid before its DOM goes.
    this.releaseMatrixGrid();
    ui.empty(this.host);
    // Palette cache is keyed by name only, so clear it on a possible theme switch.
    clearCssColorCache();
    if (this.summaryTabActive)
      this.summaryPanel.activateSummaryTab();
    if (this.matrices.length === 0) {
      this.host.appendChild(ui.divText(this.noMatricesMessage(), 'chem-sar-empty-note'));
      return;
    }
    this.host.appendChild(this.buildNavigator());
    const newNav = this.host.querySelector('.chem-sar-nav-list');
    if (newNav instanceof HTMLElement)
      newNav.scrollTop = navScroll;
    this.renderMatrixPane();
    // Only when it's the tab on screen: rebuilt hidden, its grid would size against a display:none pane.
    if (this.transferTabActive)
      this.transferPanel.activateTransferTab();
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
