import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

import {_package} from '../../package';
import {getRdKitModule} from '../../utils/chem-common-rdkit';
import {AdditiveFit, fitAdditiveFromTriples} from './sar-matrix-assemble';
import {cachedAtomCount, median} from './sar-matrix-decompose';
import {fitRoleEffects, ROLE_FIT_MAX_SWEEPS, RoleFit, RoleLevel,
  RoleSummary} from './sar-matrix-role-fit';
import {closeGridQuietly, logSarTime, SarMatrix, SarMatrixCell} from './sar-matrix-types';
import {ANALOG_W, BENEFIT_MOL_H, BENEFIT_MOL_W, CARD_CORE_H, CARD_CORE_W, CELL_H, CORE_BG_ARGB,
  CORE_W, MatrixCellRef, paintMoleculeOnColor, STRIP_MOL_H, STRIP_MOL_W,
  TAB_TRANSFER} from './sar-matrix-ui-common';

const SUM_ROWS = 3;
const SUM_POOL = 24;
/** Below this the leave-one-out fit does not support acting on a prediction. */
const TRUST_R2 = 0.5;
/** A prediction resting on fewer measured compounds than this on either axis is extrapolation. */
const MIN_SUPPORT = 3;
/** Measured pairs a pooled swap needs, and lineages it must span, before it is worth a row. */
const SWAP_MIN_PAIRS = 3;
const SWAP_MIN_SERIES = 2;
/** Measured cells kept per row when enumerating swaps; both extremes survive the trim. */
const SWAP_ROW_CAP = 32;
/** Lineages an R-group must win in before it is ranked; two is the mean of two numbers. */
const RGROUP_MIN_SERIES = 3;
/** Cross-validatable cells a fit needs before its R² is quotable as "best validated". */
const BEST_FIT_MIN_N = 8;
const TRUST_LIST_MAX = 10;
/** Dock widths and height below which the tab drops marks for their text, in px. */
const NARROW_PX = 720;
const XNARROW_PX = 560;
const SHORT_PX = 340;
/** Component rows the landing band shows before folding the rest away, so the swap row under them and
 *  the band under that stay on the screen the reader lands on. */
const FINDING_ROWS = 4;
/** Matrices a per-core outcome strip can carry before a row of squares stops being a mark. */
const STRIP_SLOTS = 12;
/** Reliable losers shown under the winners; a loss over two lineages is not a rule, so the tier the
 *  winners get is not offered here. */
const LOSER_ROWS = 2;
/** Shortest a strip square's directional fill may draw and still read as a value rather than the
 *  mid-rule itself. */
const STRIP_MIN_FILL = 3;
/** Rows the analog list holds. A chemist does not browse ten thousand; two hundred is more than a
 *  quarter's synthesis and fits one frame of rendered structures. */
const ANALOG_LIST_MAX = 200;
/** Multiples of a series' own leave-one-out error a tick run shows before a count stops reading at a
 *  glance. */
const GAIN_TICKS = 4;
/** Difference in noise-corrected spread below which two component columns are one answer with two
 *  names: offsets print at two decimals, so anything under this is not on screen at all. */
const ROLE_SPREAD_TIE = 0.01;

/** Said on the Overview and again on the core card, which a reader can arrive at either way round. */
const CORES_NOT_COMPARABLE = 'Cores are not comparable across series here — a core is one row of one ' +
  'matrix and recurs only inside its own fold lineage. Each series\' best core is in its Start-here ' +
  'expand.';

const PANE_OVERVIEW = 'Overview';
const PANE_EFFECTS = 'Effects';
const PANE_MAKING = 'Worth making';
const PANE_METHOD = 'Method';
/** The Effects segment's last tab: what is read off measured compounds inside each series — which
 *  substituent came first where, which core scored best, which swap was actually made — as against the
 *  component tabs, which are one fitted model over the whole table. */
const PANE_SERIES = 'Measured in series';
const PANES = [PANE_OVERVIEW, PANE_EFFECTS, PANE_MAKING, PANE_METHOD];

/** Fixed slots of the reason lane, in this order on every row — the lane is only comparable down a
 *  column if slot 3 means the same thing everywhere. */
const REASON_GLYPHS = ['n', '↔', '★', '✚', '✓'];
const REASON_WORDS = ['most measured compounds', 'widest measured range',
  'holds the best measured compound', 'most predictions worth making', 'best-validated fit'];

const ANALOG_COLS = {
  analog: 'Analog', predicted: 'Predicted', gain: 'Gain', interest: 'Gain / error', support: 'n',
  r2: 'R²', rmse: '± error', neighbours: 'Tried around it', series: 'Series', tier: 'Tier',
  evidence: 'Evidence', core: 'Core', fixed: 'Fixed R-groups', rgroup: 'R-group',
};

/** What the summary reads back from the viewer, and the things it asks the viewer to do. */
export interface SummaryHost {
  readonly matrices: SarMatrix[];
  /** Lineage root id per matrix, index-aligned with `matrices` — the one hierarchy the navigator uses. */
  readonly matrixRoots: string[];
  /** Displayed tier per matrix (the L1/L2/L3 on the navigator cards), index-aligned with `matrices`. */
  readonly matrixTiers: number[];
  /** Rows carrying an activity value at all. A compound the assay never reached and one no core could
   *  group are both outside every matrix, and only this separates them. */
  readonly assayedCount: number;
  /** Whether SAR transfer detection has run, and what it found. Detection is lazy and lives on its own
   *  tab, so the landing screen reports its state rather than triggering it. */
  readonly transferSummary: {scanned: boolean, count: number};
  readonly higherIsBetter: boolean;
  readonly scalingLabel: string;
  /** Whether a difference in activity units is already a log ratio. A raw column declared
   *  higher-is-better is a precomputed pIC50: no transform applied, but the numbers are logs. */
  readonly activityIsLog: boolean;
  readonly activityColumnName: string;
  /** Rows of the host table — the denominator the coverage line is a fraction of. */
  readonly hostRowCount: number;
  /** Assayed values the chosen scaling cannot represent; they are in no matrix and are not untested. */
  readonly unscalableCount: number;
  /** The fragment column every matrix varies, or null when these came from fragmentation. Only when
   *  it is set does one substituent label mean the same thing in two series. */
  readonly axisRole: string | null;
  /** The column the cores came from, or null when they came from fragmentation. A degrader set names
   *  it Linker, and "Best core" over an unnamed scaffold is the same finding nobody can act on. */
  readonly coreRole: string | null;
  /** Whether a matrix IS one core, so its fitted mean compares cores rather than groupings. */
  readonly coresAreSeries: boolean;
  /** Fragment columns that could run across the top, the current axis included. */
  readonly roleColumns: string[];
  setColumnAxis(name: string): void;
  /** Select the compounds carrying one value of one component column, in the table the analysis ran on. */
  selectRoleValue(role: string, value: string): void;
  readonly predictVirtual: boolean;
  readonly predictUnmeasured: boolean;
  readonly computing: boolean;
  noMatricesMessage(): string;
  formatActivity(value: number): string;
  cellIdText(cell: SarMatrixCell): string | null;
  observedNeighbours(matrix: SarMatrix, ri: number, ci: number): number;
  revealCell(matrix: SarMatrix, ri: number, ci: number, position?: string): void;
  revealMatrix(matrix: SarMatrix): void;
  showTab(name: string): void;
  addCellsToMakeList(cells: MatrixCellRef[], emptyMessage: string): void;
}

interface SummaryRow {
  /** What the pool deduplicates on; held rather than re-derived, since only the pool knows it. */
  key: string;
  matrix: SarMatrix;
  ri: number;
  ci: number;
  /** Always higher-is-better, so one comparator orders every pool. */
  score: number;
}

/** Both extremes by value once a group offers more pairs than the cap allows. A cap taken in the
 *  group's own order would keep the commonest substituents and truncate away the rare one that jumped
 *  two logs; this loses only mid-range pairs. */
function keepExtremes<T>(items: T[], valueOf: (item: T) => number): T[] {
  if (items.length <= SWAP_ROW_CAP)
    return items;
  const half = SWAP_ROW_CAP >> 1;
  const sorted = [...items].sort((a, b) => valueOf(a) - valueOf(b));
  return [...sorted.slice(0, half), ...sorted.slice(sorted.length - half)];
}

function supportOf(row: SummaryRow): number {
  return row.matrix.cells[row.ri][row.ci].support ?? 0;
}

/** Ties end at `matrix.id`: fragments arrive in worker-completion order, so anything resolved by
 *  index or Map order would make the cards depend on scheduling. */
function finerSeries(a: SummaryRow, b: SummaryRow): boolean {
  return a.matrix.level !== b.matrix.level ? a.matrix.level < b.matrix.level : a.matrix.id < b.matrix.id;
}

function betterSupported(a: SummaryRow, b: SummaryRow): boolean {
  const sa = supportOf(a);
  const sb = supportOf(b);
  return sa !== sb ? sa > sb : finerSeries(a, b);
}

function betterEvidenced(a: SummaryRow, b: SummaryRow): boolean {
  const sa = supportOf(a);
  const sb = supportOf(b);
  if (sa !== sb)
    return sa > sb;
  const ra = a.matrix.confidence?.r2 ?? Number.NEGATIVE_INFINITY;
  const rb = b.matrix.confidence?.r2 ?? Number.NEGATIVE_INFINITY;
  return ra !== rb ? ra > rb : a.matrix.id < b.matrix.id;
}

/** A bounded best-of pool that deduplicates on insert. Deduplicating a pool of cells afterwards can
 *  leave it shorter than the card has room for, and the same compound occupies a cell in every tier
 *  that folded it. */
class TopList {
  private readonly rows: SummaryRow[] = [];
  private readonly byKey = new Map<string, SummaryRow>();
  /** Keys that cleared every gate and are not in the pool for want of room. A set, not a counter: one
   *  structure is offered by every tier that folded its core, and a count would report it once per
   *  tier under a label the reader hears as structures. */
  private readonly missed = new Set<string>();

  /** `better` picks which occurrence of one key to keep; the score alone cannot, since a compound
   *  carries the same activity wherever it appears. */
  constructor(private readonly better: (a: SummaryRow, b: SummaryRow) => boolean,
    private readonly cap = SUM_POOL) {}

  /** The whole pool, best first — what a browsable list shows, as against the few rows a card has
   *  room for. */
  get all(): readonly SummaryRow[] {
    return this.rows;
  }

  /** How many further structures qualified. A capped list presented as a total states the cap as the
   *  finding, so every screen printing `all.length` has to be able to name this too. */
  get overflow(): number {
    return this.missed.size;
  }

  /** Drop the keys a late-arriving fact disqualifies. The structures the dataset already holds are
   *  only fully known once every matrix has been walked, and a pool offered to during that walk cannot
   *  have tested against the complete set. */
  prune(drop: (key: string) => boolean): void {
    for (let i = this.rows.length - 1; i >= 0; i--) {
      if (!drop(this.rows[i].key))
        continue;
      this.byKey.delete(this.rows[i].key);
      this.rows.splice(i, 1);
    }
    // The squeezed-out keys are disqualified by the same fact, and counting one of them as a structure
    // the list has no room for would offer to make something the dataset already holds.
    for (const key of this.missed) {
      if (drop(key))
        this.missed.delete(key);
    }
  }

  /** Primitive parameters: a full pool rejects a weak cell before anything is allocated. */
  offer(key: string, score: number, matrix: SarMatrix, ri: number, ci: number): void {
    const held = this.byKey.get(key);
    if (held === undefined && this.rows.length === this.cap && score <= this.rows[this.cap - 1].score) {
      this.missed.add(key);
      return;
    }
    const row: SummaryRow = {key, matrix, ri, ci, score};
    if (held !== undefined) {
      if (!this.better(row, held))
        return;
      this.rows.splice(this.rows.indexOf(held), 1);
    }
    let at = 0;
    while (at < this.rows.length && this.rows[at].score >= score)
      at++;
    this.rows.splice(at, 0, row);
    this.byKey.set(key, row);
    this.missed.delete(key);
    if (this.rows.length > this.cap) {
      const out = this.rows.pop()!;
      this.byKey.delete(out.key);
      this.missed.add(out.key);
    }
  }

  /** One row per series first — three analogs of one chemotype is one finding, not three — then top
   *  up from what the cap held back, so a run with fewer series than rows still fills the card.
   *  `skip` runs inside both passes: a card that filters afterwards renders short while the pool
   *  still holds rows it could have shown. */
  take(skip?: (row: SummaryRow) => boolean): SummaryRow[] {
    const out: SummaryRow[] = [];
    const seen = new Set<string>();
    for (const row of this.rows) {
      if (out.length === SUM_ROWS)
        break;
      if (seen.has(row.matrix.id) || skip?.(row))
        continue;
      seen.add(row.matrix.id);
      out.push(row);
    }
    for (const row of this.rows) {
      if (out.length === SUM_ROWS)
        break;
      if (!out.includes(row) && !skip?.(row))
        out.push(row);
    }
    return out;
  }
}

/** One measured cell of the row being walked, kept only long enough to enumerate that row's pairs. */
interface RowCell {
  ci: number;
  value: number;
  molIdx: number;
}

/** A measured swap pooled across series: one R-group exchanged inside one row, everything else
 *  identical. `from`/`to` are the lexically ordered fragment pair, so a swap cannot split into two
 *  mirror pools oriented by an arbitrary column index; the display direction is chosen at render. */
interface SwapPool {
  from: string;
  to: string;
  n: number;
  sum: number;
  min: number;
  max: number;
  /** Instances where `from → to` improved potency; the reverse count is `n - nUp`. */
  nUp: number;
  roots: Set<string>;
  /** Unordered molIdx pairs already counted — one compound pair recurs in every tier that folded it. */
  seen: Set<string>;
  best: {matrix: SarMatrix, ri: number, ciFrom: number, ciTo: number, delta: number} | null;
  sampled: boolean;
}

/** A pool read in the direction that improves potency. The pool is keyed on the lexically ordered
 *  fragment pair, so the improving direction is whichever end of its range survives; the reverse swap
 *  is the same measurements with every sign flipped. */
function orient(pool: SwapPool): {forward: boolean, from: string, to: string, worst: number,
  widest: number, mean: number, up: number} {
  const forward = pool.min >= -pool.max;
  return {
    forward,
    from: forward ? pool.from : pool.to,
    to: forward ? pool.to : pool.from,
    worst: forward ? pool.min : -pool.max,
    widest: forward ? pool.max : -pool.min,
    mean: (forward ? pool.sum : -pool.sum) / pool.n,
    up: forward ? pool.nUp : pool.n - pool.nUp,
  };
}

/** One matrix where a substituent came first at its position, under the fitted model. */
interface RGroupWin {
  matrix: SarMatrix;
  ri: number;
  ci: number;
  root: string;
  /** Measured cells behind the column — the tiebreak where no reference is recorded to compare on. */
  n: number;
  /** `dir * (colEffect[c] − colEffect[reference])` — gauge-free, so this one pools. */
  refDelta: number | null;
  fitHolds: boolean;
}

interface RGroupRow {
  subst: string;
  win: RGroupWin;
  /** Lineages won after the trust gate, and before it. */
  k: number;
  m: number;
  /** Lineages that carried this substituent as a comparable column at all — the denominator a reader
   *  hears in "best in 5 of 6". Wins alone would read as attempts and flatter every row. */
  tried: number;
  magnitude: number | null;
  lo: number;
  hi: number;
  magnitudeSeries: number;
}

/** One matrix's columns keyed by substituent, for the per-core strip. Only the columns the pooled
 *  comparison also counts are recorded, so the strip and the coverage sentence beneath it describe one
 *  set. `refDelta` is null only where the series records no reference substituent, and 0 where the
 *  column IS that reference. */
type StripCols = Map<string, {ci: number, n: number, refDelta: number | null}>;

/** What the R-group card accumulates over the one fit per matrix. */
interface RGroupAcc {
  wins: Map<string, RGroupWin[]>;
  losses: Map<string, RGroupWin[]>;
  tried: Map<string, Set<string>>;
}

/** Everything one matrix contributes, extracted from its fit and its cells; the fit arrays are
 *  discarded with it so 345 of them are never live at once. */
interface SeriesStat {
  matrix: SarMatrix;
  root: string;
  tier: number;
  /** Distinct measured molecules from the walk — never `realCount`, which is captured pre-prune. */
  cpd: number;
  realCells: number;
  totalCells: number;
  /** Cells the decomposition cannot express, so the fill denominator is not `rows × columns`. Zero in
   *  a matrix whose assembly never marked them, which is NOT the same as a complete grid. */
  impossibleCells: number;
  lo: number;
  hi: number;
  best: {ri: number, ci: number, value: number} | null;
  /** The fit's own count-weighted mean. Unlike the centred effects it is one activity column under one
   *  scaling, so it is the one fitted quantity comparable between matrices. */
  typical: number;
  /** Most potent prediction that passes the trust gate — the model's own ceiling for this series. */
  bestVirtual: {ri: number, ci: number, value: number} | null;
  /** Unfilled cells this series carries a prediction for. With the matrix's own R2 it says WHY there
   *  is nothing to aim at. */
  virtualCells: number;
  trusted: number;
  converged: boolean;
  /** Filled only for the matrices a strip can draw; null everywhere else. */
  stripCols: StripCols | null;
  /** Range of the fitted substituent and core effects — centred, same units, same fit, so the two
   *  are comparable within this matrix and the verdict pools as a count. */
  colRange: number | null;
  rowRange: number | null;
  bestRow: {ri: number, effect: number, n: number} | null;
}

interface StartRow {
  primary: SeriesStat;
  /** Indices into the fixed reason-lane slots, not phrases: the lane's shape is what a reader
   *  compares down the column. */
  reasons: number[];
}

interface SummaryData {
  compounds: number;
  untested: number;
  measuredCells: number;
  minObserved: number | null;
  maxObserved: number | null;
  trustedCells: number;
  trustedStructures: number;
  /** Predictions that pass the same gate but carry no structure: the core holds an attachment point
   *  none of the picked fragment columns fills, so nothing can be completed over it. */
  trustedNoStructure: number;
  /** Typical leave-one-out prediction error across series — model error, never an assay σ. Null when
   *  no series has a cross-validated fit, so it is never printed as zero. */
  modelError: number | null;
  fitHolds: number;
  unchecked: number;
  lowR2: SarMatrix[];
  lowR2Virtual: number;
  families: number;
  axisRole: string | null;
  coreRole: string | null;
  coresAreSeries: boolean;
  /** One additive fit over every role column at once; null outside fragment-columns mode, where a
   *  substituent label is local to its own series and one pooled offset would average unrelated
   *  quantities. */
  roleFit: RoleFit | null;
  /** Best measured cell per value of each role column, so a leaderboard row lands on a compound.
   *  Keyed by role name and in the order the fit was given the roles, which the fit's own ranking
   *  discards. */
  roleBest: Map<string, Map<string, MatrixCellRef & {value: number}>>;
  tierCounts: {tier: number, n: number}[];
  /** Series whose additive fit stopped short of its tolerance. They are left out of every pooled
   *  R-group comparison, so an empty leaderboard has to be able to name this as the cause. */
  nonConverged: number;
  series: SeriesStat[];
  startHere: StartRow[];
  swaps: SwapPool[];
  /** The same, per component column, grouped straight off the component values rather than off one
   *  matrix's columns — so every component's measured pairs are in hand at once and none of them costs
   *  a rebuild. Empty outside fragment-columns mode, where a substituent label is local to its series. */
  swapsByRole: Map<string, SwapPool[]>;
  /** Whether any swap candidate existed at all — "none clears the gate" is a different sentence from
   *  "no row carries three measured substituents". */
  swapCandidates: number;
  /** The largest single measured move anywhere, so "nothing pooled clears the gate" can still name
   *  what the data does hold. */
  swapBest: SwapPool | null;
  rgroups: RGroupRow[];
  rgroupsThin: RGroupRow[];
  rgroupLosers: RGroupRow[];
  test: TopList;
  /** Virtual analogs ranked on predicted gain in their own series' leave-one-out error, deduplicated
   *  on structure — tiers overlap, so one molecule is proposed by several matrices. */
  analogs: TopList;
  /** Same gates bar the cross-validatable-cell floor: a tiny fit's error is not a scale to divide by,
   *  so these are listed without a rank rather than mixed into one. */
  analogsThin: TopList;
  /** The best of what the trust gate turned down, by raw gain over its own series' best measured
   *  compound. Shown only when the two lists above are empty, so the screen always names a next
   *  structure rather than reporting that nothing qualifies. */
  analogsAny: TopList;

  /** Why a prediction that passes support and fit is still not offered. Never one total: they call for
   *  different things from the reader, and only `withheldFitFails` has a list to open — a fit that was
   *  never checked and a fit that stopped short of its tolerance appear nowhere. */
  withheldBelowError: number;
  withheldThinSupport: number;
  withheldFitFails: number;
  withheldUnchecked: number;
  withheldNotConverged: number;
  /** Predicted structures the dataset already holds — an assay plate, not a synthesis. Structures, so
   *  it stands beside the analog count rather than multiplying it by the tiers that folded the core. */
  alreadyHeld: number;
  /** Structures that cleared every gate and are past the list's row cap, so "N worth making" can say
   *  what it is the top N of. */
  analogOverflow: number;
}

/** One row of the Worth-making grid and which list it came off: `ranked` cleared every gate, `thin` rests
 *  on a fit with too few cross-validatable cells to rank by, `ungated` cleared no gate at all and is
 *  shown only when the other two lists are empty. */
interface AnalogRow {
  row: SummaryRow;
  evidence: string;
}

/** One answer tile's content: the structures the conclusion is about, the conclusion, and the one clause
 *  that stops the conclusion being over-read. */
interface Answer {
  answer: string;
  negative: string;
  /** Drawn above the conclusion. A substituent written out is truncated to fit the tile, and three
   *  analogs of one scaffold truncate to the same string. */
  art?: HTMLElement;
}

/** Digits grouped with a space rather than the locale's separator: a locale that groups with a period
 *  renders 2749 compounds as "2.749" directly beside an activity of "-5.9". */
function count(n: number): string {
  return n.toLocaleString('en-US').replace(/,/g, ' ');
}

/** A fragment SMILES cut to fit a one-line descriptor; the full string stays in the tooltip. */
function shortSmiles(smiles: string): string {
  return smiles.length > 20 ? `${smiles.slice(0, 19)}…` : smiles;
}

/** A blank role value is the unsubstituted parent, which is a level like any other — rendered as a
 *  name rather than as the empty string it is stored as.
 *
 *  Not shortened: a fixed character cap renders three analogs of one scaffold as the same string, since
 *  what distinguishes them is past the cap. A row ellipsizes in CSS instead, to whatever width it has. */
/** One half of a measured pair: the component value it carries, what it measured, and which compound it
 *  is — the three things a pool needs from each side, whether the side came off a matrix cell or off a
 *  row of component columns. */
interface SwapSide {
  value: string;
  activity: number;
  mol: number;
}

/** Where a component's swap card sits on the Effects segment, so the band row for that component scrolls
 *  to the card about that component rather than to whichever one was built first. */
function swapAnchor(role: string): string {
  return `swaps-${role}`;
}

function roleValueName(value: string): string {
  return value === '' ? 'nothing at this position' : value;
}

/** Whether this role's own values may be put in an order at all: two of them have to be readable, and
 *  a split of credit that does not survive refitting each half of the table is not this role's. */
function roleRanks(role: RoleSummary): boolean {
  return role.levels.length >= 2 && (role.repeat === null || role.repeat >= TRUST_R2);
}

/**
 * The Summary tab: the landing screen of the analysis. The scale and direction are pinned above
 * everything, because they invalidate every ranking under them; under that a segmented control opens
 * one of four panes, of which the first — the answers, which series to open, and what the analysis
 * covers — fits without scrolling at a normal dock size. Each row lands on the cell it describes.
 *
 * Holds one Dart-backed object, the analog grid, and it is built only when its own segment is opened;
 * every teardown path closes it. It holds no view and no `TableView`, so it can reach no dock node.
 */
export class SummaryPanel {
  readonly root = ui.divV([], 'chem-sar-main');
  private data: SummaryData | null = null;
  /** Depictions are painted one tick after the DOM lands, so a tab switch is not held up by RDKit. */
  private paintTimer = 0;
  /** Collecting over every cell of every matrix holds the main thread, so it runs off the activation
   *  stack — which also lets the loader paint. */
  private collectTimer = 0;
  private pendingPaints: (() => void)[] = [];

  /** Set by the answer tile so the R-group leaderboard opens its top row when the Effects pane builds. */
  private expandTopRGroup = false;
  /** Which of the Effects segment's own tabs is open: an index into its role cards, or past the last
   *  of them the per-series evidence. Held across pane rebuilds so leaving the segment and coming back
   *  does not lose the component the reader was on. */
  private effectsTab = 0;
  /** Whether the poorly-fitting list starts open: a reader who followed a "the fit does not hold"
   *  link came for that list, not for the explanation above it. */
  private openTrustList = false;
  /**
   * The fold tier every ranking on this tab is read at, or null for all of them together.
   *
   * Not a cut of the data: the matrices are built once and all of them stay. This only decides which
   * of them the tab's rankings walk — a tier holds the same compounds as the tier below it over cores
   * cut one bond broader, so reading at one tier answers "what does the SAR look like at this breadth"
   * without the other breadths mixed in.
   */
  private tierFilter: number | null = null;

  /** Which halves of the landing band the reader has shut, by key. Held across pane rebuilds, so a
   *  collapsed half stays collapsed when the segment is left and come back to. */
  private readonly foldedBands = new Set<string>();
  /** The one line that reports the transfer scan, kept so its state can be refreshed without
   *  rebuilding the pane around it. */
  private transferLine: HTMLElement | null = null;
  /** A fault line stands until something invalidates it: collecting over it would replace the reason
   *  the analysis failed with a generic empty-state note. */
  private message = false;
  private paneHost: HTMLElement | null = null;
  private currentPane = PANE_OVERVIEW;
  private readonly segButtons = new Map<string, HTMLElement>();
  /** The one scrolling region of the pane on screen, or null where the pane does not scroll. Moving
   *  this rather than calling `scrollIntoView` is what keeps a scroll inside the panel. */
  private scroller: HTMLElement | null = null;
  private analogGrid: DG.Grid | null = null;
  private analogSub: Subscription | null = null;
  private analogSources: MatrixCellRef[] = [];
  private sizeSub: Subscription | null = null;
  /** The width breakpoints the pane on screen was built against; empty before the first build. */
  private widthKey = '';

  constructor(private readonly host: SummaryHost) {}

  activateSummaryTab(): void {
    if (this.host.computing) {
      this.reset();
      this.root.appendChild(ui.divV([ui.loader(), ui.divText('Building SAR matrices...')],
        'chem-sar-empty-note'));
      return;
    }
    // Occupancy alone cannot say the tab is current: the loader above is a child too, and the
    // re-activation that follows a compute would then leave it standing forever.
    if (this.message || this.collectTimer !== 0 ||
      (this.data !== null && this.root.childElementCount > 0)) {
      // Detection runs on the Transfer tab and this pane is re-shown rather than rebuilt, so the one
      // line that reports it would otherwise keep inviting a scan that has already run.
      this.syncTransferLine();
      return;
    }
    if (this.data !== null) {
      this.render();
      return;
    }
    // Letting the loader paint first also lets the dock finish laying out its tab strip before this
    // pane's content lands.
    ui.empty(this.root);
    this.root.appendChild(ui.divV([ui.loader(), ui.divText('Reading the analysis...')],
      'chem-sar-empty-note'));
    this.collectTimer = window.setTimeout(() => {
      this.collectTimer = 0;
      if (this.host.computing || this.message)
        return;
      this.data = this.collect();
      this.render();
    }, 0);
  }

  invalidate(): void {
    this.reset();
    this.message = false;
  }

  /** Re-read the analysis at another fold tier and rebuild the segment on screen. Nothing is recomputed
   *  in the viewer: the matrices are already built, and this only changes which of them are walked. */
  private setTierFilter(tier: number | null): void {
    if (this.tierFilter === tier)
      return;
    this.tierFilter = tier;
    this.closeAnalogGrid();
    this.data = this.collect();
    this.render(true);
  }

  showMessage(text: string): void {
    this.reset();
    this.message = true;
    this.root.appendChild(ui.divText(text, 'chem-sar-empty-note'));
  }

  release(): void {
    this.reset();
    this.message = false;
    this.sizeSub?.unsubscribe();
    this.sizeSub = null;
  }

  /** Every timer, every Dart-backed object and every rendered node this panel owns. The grid is the
   *  one thing a dropped reference does not release, so it is closed on all three teardown paths. */
  private reset(): void {
    window.clearTimeout(this.paintTimer);
    window.clearTimeout(this.collectTimer);
    this.collectTimer = 0;
    this.closeAnalogGrid();
    this.data = null;
    this.paneHost = null;
    this.scroller = null;
    this.currentPane = PANE_OVERVIEW;
    // A new analysis has its own tiers. Kept, a tier the new run does not hold filters every matrix
    // out, which reads as an empty dataset — and the chip bar that would clear it is built from the
    // tiers that survived the filter, so there is none.
    this.tierFilter = null;
    this.widthKey = '';
    this.segButtons.clear();
    this.pendingPaints = [];
    ui.empty(this.root);
  }

  private closeAnalogGrid(): void {
    this.analogSub?.unsubscribe();
    this.analogSub = null;
    closeGridQuietly(this.analogGrid);
    this.analogGrid = null;
    this.analogSources = [];
  }

  /** `scrollIntoView` scrolls every scrollable ancestor, the dock container included, so it can move
   *  the layout around the panel; only this pane's own scroller may move. */
  private scrollTo(el: HTMLElement): void {
    const s = this.scroller;
    if (s === null)
      return;
    s.scrollTop += el.getBoundingClientRect().top - s.getBoundingClientRect().top;
  }

  // ---- Collection -----------------------------------------------------------------------------

  /**
   * One pass over every cell of every matrix, feeding the totals, the two pools and the per-series
   * statistics.
   *
   * The observed cells are collected as triples on the way past and the additive fit is run from
   * those, so the fit does not walk the grid a second time. The buffers are reused across matrices:
   * 345 fresh sets would trade the scan cost for GC cost, and nothing caps the SUM of cells across
   * matrices — only each matrix.
   */
  private collect(): SummaryData {
    const t0 = performance.now();
    const host = this.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const log = host.activityIsLog;
    const roots = host.matrixRoots;
    const tiers = host.matrixTiers;

    const compounds = new Set<number>();
    const untested = new Set<number>();
    const trustedStructures = new Set<string>();
    // Assembled from the source table's molecule column while a virtual SMILES comes from RDKit's
    // linker, so the two canonicalisations need not agree: this catches matches, and is not a proof of
    // novelty — which is why a listed analog says "not in this dataset" and never "never made".
    const measuredStructures = new Set<string>();
    const test = new TopList(betterSupported);
    const analogs = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    const analogsThin = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    // Everything the gate turned down, ranked on raw gain. A screen that answers "nothing qualifies"
    // and stops leaves the reader with no next structure at all, when what the gate actually withheld
    // is the confidence, not the ranking.
    const analogsAny = new TopList(betterEvidenced, ANALOG_LIST_MAX);
    const swaps = new Map<string, SwapPool>();
    const acc: RGroupAcc = {wins: new Map(), losses: new Map(), tried: new Map()};
    const series: SeriesStat[] = [];
    const lowR2: SarMatrix[] = [];
    // A strip needs a fixed partner axis and few enough slots to read as a mark; fragmented matrices
    // overlap through tiers and have neither.
    const wantStrip = host.axisRole !== null && host.matrices.length <= STRIP_SLOTS;

    let measuredCells = 0;
    let trustedNoStructure = 0;
    let trustedCells = 0;
    let withheldBelowError = 0;
    let withheldThinSupport = 0;
    let withheldFitFails = 0;
    let withheldUnchecked = 0;
    let withheldNotConverged = 0;
    const alreadyHeld = new Set<string>();
    let lowR2Virtual = 0;
    let unchecked = 0;
    let fitHolds = 0;
    let nonConverged = 0;
    let sampledRows = 0;
    let swapCandidates = 0;
    let minObserved = Infinity;
    let maxObserved = -Infinity;

    const obsRow: number[] = [];
    const obsCol: number[] = [];
    const obsVal: number[] = [];
    const mols = new Set<number>();
    const rowCells: RowCell[] = [];
    // Admitting an analog on its gain needs this matrix's best measured cell, which is final only when
    // its row loop ends — so the candidates wait here rather than being offered inside it. Reused
    // across matrices: one triple of arrays, never 345 of them.
    const riBuf: number[] = [];
    const ciBuf: number[] = [];
    const valBuf: number[] = [];
    /** The same, for the candidates the trust gate turned down — they still need this matrix's best
     *  measured cell to carry a gain. */
    const riAlt: number[] = [];
    const ciAlt: number[] = [];
    const valAlt: number[] = [];

    // `roleColumns` excludes the core column, so the core role is prepended: left out, the report omits
    // the one role the request names first. Empty outside fragment-columns mode, which is what keeps the
    // fit off a decomposition whose labels are local to one series.
    const roleNames = host.axisRole === null ? [] : [host.coreRole!, ...host.roleColumns];
    // Where each role's value lives, decided once: the core on the row, the axis on the column, and
    // every other role folded into the row identity under its own column name.
    const roleKind = roleNames.map((name) => name === host.coreRole ? 0 : name === host.axisRole ? 1 : 2);
    const roleValues: string[][] = roleNames.map(() => []);
    const roleActivity: number[] = [];
    const roleMol: number[] = [];
    const roleBest = roleNames.map(() => new Map<string, MatrixCellRef & {value: number}>());

    for (let mi = 0; mi < host.matrices.length; mi++) {
      // Reading at one tier: skipping the others here is what keeps every ranking below — components,
      // swaps, cores, Start here, Worth making — on the same set of series, since all of them are built
      // from this one walk.
      if (this.tierFilter !== null && tiers[mi] !== this.tierFilter)
        continue;
      const matrix = host.matrices[mi];
      const nRows = matrix.rows.length;
      const nCols = matrix.columns.length;
      const r2 = matrix.confidence?.r2 ?? null;
      const holds = r2 !== null && r2 >= TRUST_R2;

      obsRow.length = 0;
      obsCol.length = 0;
      obsVal.length = 0;
      riBuf.length = 0;
      ciBuf.length = 0;
      valBuf.length = 0;
      riAlt.length = 0;
      ciAlt.length = 0;
      valAlt.length = 0;
      mols.clear();
      let realCells = 0;
      let impossibleCells = 0;
      let trusted = 0;
      let lo = Infinity;
      let hi = -Infinity;
      let best: {ri: number, ci: number, value: number} | null = null;
      let bestVirtual: {ri: number, ci: number, value: number} | null = null;
      let virtualCells = 0;

      for (let ri = 0; ri < nRows; ri++) {
        rowCells.length = 0;
        for (let ci = 0; ci < nCols; ci++) {
          const cell = matrix.cells[ri][ci];
          const value = cell.value;
          if (cell.kind === 'real' && value !== null && cell.molIdx !== null) {
            compounds.add(cell.molIdx);
            mols.add(cell.molIdx);
            measuredCells++;
            realCells++;
            obsRow.push(ri);
            obsCol.push(ci);
            obsVal.push(value);
            rowCells.push({ci, value, molIdx: cell.molIdx});
            if (roleNames.length > 0) {
              for (let r = 0; r < roleNames.length; r++) {
                const level = roleKind[r] === 0 ? matrix.rows[ri].coreSmiles :
                  roleKind[r] === 1 ? matrix.columns[ci].substSmiles :
                    matrix.rows[ri].foldedValues[roleNames[r]] ?? '';
                roleValues[r].push(level);
                const held = roleBest[r].get(level);
                if (held === undefined || dir * value > dir * held.value)
                  roleBest[r].set(level, {matrix, ri, ci, value});
              }
              roleActivity.push(value);
              roleMol.push(cell.molIdx);
            }
            if (value < lo)
              lo = value;
            if (value > hi)
              hi = value;
            if (value < minObserved)
              minObserved = value;
            if (value > maxObserved)
              maxObserved = value;
            if (best === null || dir * value > dir * best.value)
              best = {ri, ci, value};
            if (cell.smiles !== null)
              measuredStructures.add(cell.smiles);
          } else if (cell.kind === 'unmeasured' && cell.molIdx !== null) {
            untested.add(cell.molIdx);
            if (cell.smiles !== null)
              measuredStructures.add(cell.smiles);
            if (value !== null)
              test.offer(String(cell.molIdx), dir * value, matrix, ri, ci);
          } else if (cell.kind === 'virtual' && value !== null) {
            const support = cell.support ?? 0;
            virtualCells++;
            const passes = support >= MIN_SUPPORT && holds;
            if (passes && (bestVirtual === null || dir * value > dir * bestVirtual.value))
              bestVirtual = {ri, ci, value};
            // A structureless prediction is neither actionable nor withheld for want of evidence: the
            // core carries an attachment point no picked column fills, which is a settings fix.
            if (cell.smiles === null) {
              if (passes)
                trustedNoStructure++;
            } else if (passes) {
              trustedCells++;
              trusted++;
              trustedStructures.add(cell.smiles);
              if (measuredStructures.has(cell.smiles))
                alreadyHeld.add(cell.smiles);
              else {
                riBuf.push(ri);
                ciBuf.push(ci);
                valBuf.push(value);
              }
            } else {
              if (!measuredStructures.has(cell.smiles)) {
                riAlt.push(ri);
                ciAlt.push(ci);
                valAlt.push(value);
              }
              if (support < MIN_SUPPORT)
                withheldThinSupport++;
              else if (r2 === null)
                withheldUnchecked++;
              else
                withheldFitFails++;
            }
          } else if (cell.kind === 'impossible')
            impossibleCells++;
        }
        if (rowCells.length >= 2) {
          swapCandidates++;
          if (rowCells.length > SWAP_ROW_CAP)
            sampledRows++;
          this.poolSwaps(swaps, matrix, ri, rowCells, dir, log, roots[mi]);
        }
      }

      const fit = fitAdditiveFromTriples(obsRow, obsCol, obsVal, nRows, nCols);
      if (!fit.converged)
        nonConverged++;
      // Shrinkage toward the zero initialisation differs per design, so a gain read off a fit that
      // stopped short is not comparable with one read off a converged fit.
      const conf = matrix.confidence ?? null;
      if (!fit.converged || best === null || conf === null)
        withheldNotConverged += riBuf.length;
      else {
        // A fit whose error rests on a handful of cross-validatable cells is not a scale to divide by:
        // a plane through four points has a tiny RMSE and would top every ranking.
        const thin = conf.n < BEST_FIT_MIN_N || conf.rmse <= 0;
        for (let k = 0; k < riBuf.length; k++) {
          const gain = dir * (valBuf[k] - best.value);
          if (thin) {
            if (gain > 0) {
              analogsThin.offer(matrix.cells[riBuf[k]][ciBuf[k]].smiles!, gain,
                matrix, riBuf[k], ciBuf[k]);
            }
          } else if (gain < conf.rmse)
            withheldBelowError++;
          else {
            analogs.offer(matrix.cells[riBuf[k]][ciBuf[k]].smiles!, gain / conf.rmse,
              matrix, riBuf[k], ciBuf[k]);
          }
        }
      }
      // Ranked on raw gain, never on gain over error: the fits behind these are the ones that were not
      // trusted, so their error is not a scale to divide by.
      if (best !== null) {
        for (let k = 0; k < riAlt.length; k++) {
          const gain = dir * (valAlt[k] - best.value);
          if (gain > 0) {
            analogsAny.offer(matrix.cells[riAlt[k]][ciAlt[k]].smiles!, gain,
              matrix, riAlt[k], ciAlt[k]);
          }
        }
      }
      const stat: SeriesStat = {
        matrix, root: roots[mi], tier: tiers[mi], cpd: mols.size, realCells, impossibleCells,
        totalCells: nRows * nCols, lo: realCells ? lo : 0, hi: realCells ? hi : 0, best, bestVirtual,
        typical: fit.grandMean, trusted, converged: fit.converged, virtualCells,
        stripCols: this.recordRGroupExtremes(acc, matrix, fit, dir, roots[mi], holds, wantStrip),
        ...this.extractRanges(matrix, fit, dir),
      };
      series.push(stat);

      if (r2 === null)
        unchecked++;
      else if (r2 < TRUST_R2) {
        lowR2.push(matrix);
        lowR2Virtual += matrix.virtualCount;
      } else
        fitHolds++;
    }

    if (nonConverged > 0) {
      _package.logger.warning(`SAR Matrix | the additive fit of ${nonConverged} series did not reach ` +
        'its tolerance; they are left out of the pooled R-group comparison');
    }
    if (sampledRows > 0) {
      _package.logger.warning(`SAR Matrix | ${sampledRows} rows carry more than ${SWAP_ROW_CAP} measured ` +
        'cells; the swap pool keeps their most and least potent halves, so mid-range pairs are dropped');
    }

    // The walk offers a candidate against the structures seen so far, so a molecule first measured in a
    // later matrix is only caught here.
    const dropHeld = (key: string): boolean => {
      if (!measuredStructures.has(key))
        return false;
      alreadyHeld.add(key);
      return true;
    };
    analogs.prune(dropHeld);
    analogsThin.prune(dropHeld);
    analogsAny.prune(dropHeld);

    // Over the matrices this walk read, not over every matrix: at a tier filter none of the others
    // reaches any other number on the tab, and a model error pooled across tiers is the error of no
    // ranking shown — it sets the band every effect is read against and the Start-here threshold.
    const walked = host.matrices.filter((_, mi) =>
      this.tierFilter === null || tiers[mi] === this.tierFilter);
    const confN = walked.map((m) => m.confidence?.n ?? null)
      .filter((n): n is number => n !== null);
    const rmses = walked.map((m) => m.confidence?.rmse ?? null)
      .filter((r): r is number => r !== null);
    const families = new Set(roots).size;
    const roleFit = fitRoleEffects({names: roleNames, values: roleValues, activity: roleActivity,
      molIdx: roleMol, minSupport: MIN_SUPPORT, higherIsBetter: host.higherIsBetter});
    const data: SummaryData = {
      compounds: compounds.size,
      untested: untested.size,
      measuredCells,
      minObserved: measuredCells ? minObserved : null,
      maxObserved: measuredCells ? maxObserved : null,
      trustedCells,
      trustedStructures: trustedStructures.size,
      trustedNoStructure,
      modelError: rmses.length ? median(rmses) : null,
      fitHolds, unchecked, lowR2, lowR2Virtual,
      families,
      axisRole: host.axisRole,
      coreRole: host.coreRole,
      coresAreSeries: host.coresAreSeries,
      roleFit,
      roleBest: new Map(roleNames.map((name, r) => [name, roleBest[r]])),
      tierCounts: this.countTiers(tiers),
      nonConverged,
      series,
      startHere: this.mergeStartHere(series, confN, dir),
      swaps: this.rankSwaps(swaps),
      swapsByRole: this.poolRoleSwaps(roleNames, roleValues, roleActivity, roleMol, dir, log),
      swapCandidates,
      swapBest: this.largestSwap(swaps),
      ...this.rankRGroups(acc),
      test, analogs, analogsThin, analogsAny,
      withheldBelowError, withheldThinSupport, withheldFitFails, withheldUnchecked, withheldNotConverged,
      alreadyHeld: alreadyHeld.size,
      analogOverflow: analogs.overflow + analogsThin.overflow,
    };
    logSarTime('summary', t0);
    return data;
  }

  private countTiers(tiers: number[]): {tier: number, n: number}[] {
    const byTier = new Map<number, number>();
    for (const tier of tiers)
      byTier.set(tier, (byTier.get(tier) ?? 0) + 1);
    return [...byTier.entries()].sort((a, b) => a[0] - b[0]).map(([tier, n]) => ({tier, n}));
  }

  /** The scalars the screen needs out of one fit, so the Float64Arrays can be dropped with it. Both
   *  ranges are of centred effects, so they are comparable to each other inside this matrix — and to
   *  nothing outside it. One support floor on both axes: an effect estimated from two observations is
   *  noisier than one estimated from three, so admitting columns at a lower floor than rows would
   *  widen the column range before any chemistry entered. */
  private extractRanges(matrix: SarMatrix, fit: ReturnType<typeof fitAdditiveFromTriples>, dir: number):
    {colRange: number | null, rowRange: number | null, bestRow: {ri: number, effect: number, n: number} | null} {
    let colLo = Infinity;
    let colHi = -Infinity;
    let colN = 0;
    for (let ci = 0; ci < matrix.columns.length; ci++) {
      if (fit.colN[ci] < MIN_SUPPORT)
        continue;
      colN++;
      colLo = Math.min(colLo, fit.colEffect[ci]);
      colHi = Math.max(colHi, fit.colEffect[ci]);
    }
    let rowLo = Infinity;
    let rowHi = -Infinity;
    let rowN = 0;
    let bestRow: {ri: number, effect: number, n: number} | null = null;
    for (let ri = 0; ri < matrix.rows.length; ri++) {
      if (fit.rowN[ri] < MIN_SUPPORT)
        continue;
      rowN++;
      rowLo = Math.min(rowLo, fit.rowEffect[ri]);
      rowHi = Math.max(rowHi, fit.rowEffect[ri]);
      if (bestRow === null || dir * fit.rowEffect[ri] > dir * bestRow.effect)
        bestRow = {ri, effect: fit.rowEffect[ri], n: fit.rowN[ri]};
    }
    return {
      colRange: colN >= 2 ? colHi - colLo : null,
      rowRange: rowN >= 2 ? rowHi - rowLo : null,
      bestRow,
    };
  }

  /**
   * Measured pairs inside ONE row index, pooled on the unordered fragment pair.
   *
   * The row index, not the core: a row is keyed on `[coreSmiles, ...foldedValues]`, so two rows can
   * share a core and differ at a folded position — pooling per core would break "everything else
   * identical", which is the whole claim of the card.
   */
  private poolSwaps(pools: Map<string, SwapPool>, matrix: SarMatrix, ri: number, cells: RowCell[],
    dir: number, log: boolean, root: string): void {
    const sampled = cells.length > SWAP_ROW_CAP;
    const kept = keepExtremes(cells, (cell) => cell.value);
    for (let a = 0; a < kept.length; a++) {
      for (let b = a + 1; b < kept.length; b++) {
        const added = this.addSwapPair(pools,
          {value: matrix.columns[kept[a].ci].substSmiles, activity: kept[a].value, mol: kept[a].molIdx},
          {value: matrix.columns[kept[b].ci].substSmiles, activity: kept[b].value, mol: kept[b].molIdx},
          root, sampled, dir, log);
        if (added === null)
          continue;
        const lower = added.flip ? kept[b] : kept[a];
        const upper = added.flip ? kept[a] : kept[b];
        if (added.pool.best === null || Math.abs(added.delta) > Math.abs(added.pool.best.delta))
          added.pool.best = {matrix, ri, ciFrom: lower.ci, ciTo: upper.ci, delta: added.delta};
      }
    }
  }

  /**
   * Fold one measured pair into its pool, or decline it. Null when the two sides carry the same value,
   * when a fold cannot be taken, or when this compound pair is already counted.
   *
   * The caller records `best` rather than this: a pair found through a matrix can name the two cells it
   * came from, and a pair found by grouping the component columns cannot.
   */
  private addSwapPair(pools: Map<string, SwapPool>, a: SwapSide, b: SwapSide,
    root: string, sampled: boolean, dir: number, log: boolean):
    {pool: SwapPool, delta: number, flip: boolean} | null {
    const {value: sa, activity: va, mol: molA} = a;
    const {value: sb, activity: vb, mol: molB} = b;
    if (sa === sb)
      return null;
    const flip = sa > sb;
    const lowerVal = flip ? vb : va;
    const upperVal = flip ? va : vb;
    // A fold multiplies, so its statistics are geometric: pooling the log of the ratio makes the mean,
    // min and max the fold statistics a chemist would quote, and lets one accumulator serve both
    // scales. A non-positive raw value has no fold at all.
    if (!log && (lowerVal <= 0 || upperVal <= 0))
      return null;
    const delta = log ? dir * (upperVal - lowerVal) :
      dir * (Math.log10(upperVal) - Math.log10(lowerVal));
    const key = `${flip ? sb : sa}\0${flip ? sa : sb}`;
    let pool = pools.get(key);
    if (pool === undefined) {
      pool = {from: flip ? sb : sa, to: flip ? sa : sb, n: 0, sum: 0, min: Infinity, max: -Infinity,
        nUp: 0, roots: new Set(), seen: new Set(), best: null, sampled: false};
      pools.set(key, pool);
    }
    const pair = `${Math.min(molA, molB)}:${Math.max(molA, molB)}`;
    if (pool.seen.has(pair))
      return null;
    pool.seen.add(pair);
    pool.n++;
    pool.sum += delta;
    pool.min = Math.min(pool.min, delta);
    pool.max = Math.max(pool.max, delta);
    if (delta > 0)
      pool.nUp++;
    pool.roots.add(root);
    pool.sampled = pool.sampled || sampled;
    return {pool, delta, flip};
  }

  /**
   * Measured swaps for every component at once, read off the component values of each measured cell
   * rather than off one matrix's columns.
   *
   * Nothing is rebuilt to reach another component. A swap is two compounds alike in every component but
   * one; where the components are given as columns that is a grouping, not a decomposition — group the
   * measured cells on every component except the one being swapped, and each pair inside a group that
   * differs in it is a matched pair. The matrix columns were only ever one route to the same thing, and
   * being one route is what made the other components cost a rebuild.
   */
  private poolRoleSwaps(roleNames: string[], roleValues: string[][], roleActivity: number[],
    roleMol: number[], dir: number, log: boolean): Map<string, SwapPool[]> {
    const out = new Map<string, SwapPool[]>();
    for (let r = 0; r < roleNames.length; r++) {
      const buckets = new Map<string, number[]>();
      // One compound occupies a cell in every tier that folded it, so the same partner tuple arrives
      // several times; the first occurrence can form every pair the later ones could.
      const placed = new Map<string, Set<number>>();
      for (let k = 0; k < roleActivity.length; k++) {
        let key = '';
        for (let j = 0; j < roleNames.length; j++) {
          if (j !== r)
            key += `${roleValues[j][k]}\u0001`;
        }
        let bucket = buckets.get(key);
        if (bucket === undefined) {
          bucket = [];
          buckets.set(key, bucket);
          placed.set(key, new Set());
        }
        const seen = placed.get(key)!;
        if (seen.has(roleMol[k]))
          continue;
        seen.add(roleMol[k]);
        bucket.push(k);
      }
      const pools = new Map<string, SwapPool>();
      for (const [key, bucket] of buckets) {
        const sampled = bucket.length > SWAP_ROW_CAP;
        const kept = keepExtremes(bucket, (i) => roleActivity[i]);
        for (let a = 0; a < kept.length; a++) {
          for (let b = a + 1; b < kept.length; b++) {
            this.addSwapPair(pools,
              {value: roleValues[r][kept[a]], activity: roleActivity[kept[a]], mol: roleMol[kept[a]]},
              {value: roleValues[r][kept[b]], activity: roleActivity[kept[b]], mol: roleMol[kept[b]]},
              key, sampled, dir, log);
          }
        }
      }
      out.set(roleNames[r], this.rankSwaps(pools));
    }
    return out;
  }

  /** The swap's worth in its better direction: what it bought in EVERY pair we have. Ranking a mean
   *  over many small pools selects the noisiest pool instead. */
  private swapScore(pool: SwapPool): number {
    return Math.max(pool.min, -pool.max);
  }

  /** The single widest measured move anywhere, gate or no gate — so a card with no qualifying pool can
   *  still say what the data does hold rather than only what it does not. */
  private largestSwap(pools: Map<string, SwapPool>): SwapPool | null {
    let best: SwapPool | null = null;
    for (const pool of pools.values()) {
      const reach = Math.max(pool.max, -pool.min);
      if (best === null || reach > Math.max(best.max, -best.min))
        best = pool;
    }
    return best;
  }

  private rankSwaps(pools: Map<string, SwapPool>): SwapPool[] {
    return [...pools.values()]
      .filter((pool) => pool.n >= SWAP_MIN_PAIRS && pool.roots.size >= SWAP_MIN_SERIES)
      .sort((a, b) => this.swapScore(b) - this.swapScore(a) ||
        (a.from < b.from ? -1 : a.from > b.from ? 1 : a.to < b.to ? -1 : 1))
      .slice(0, SUM_ROWS);
  }

  /**
   * The substituents the fitted model puts first and last at this matrix's position, recorded once per
   * lineage, plus every substituent this series tried and — where a strip will draw — its columns.
   *
   * The fit, not the raw column mean: a column mean confounds a substituent with whichever cores it
   * happened to be made on, which is the normal case in a modular library. The magnitude is a
   * within-series difference against that series' own reference substituent, so it survives the fact
   * that two matrices centre their effects over different substituent menus. An ordinal survives it
   * too: adding a constant to every `colEffect` of a matrix does not reorder them, which is what lets
   * the last place pool the same way the first does.
   */
  private recordRGroupExtremes(acc: RGroupAcc, matrix: SarMatrix, fit: AdditiveFit, dir: number,
    root: string, fitHolds: boolean, strip: boolean): StripCols | null {
    const reference = matrix.refValues[matrix.positions[0] ?? ''];
    const found = reference ? matrix.columns.findIndex((c) => c.substSmiles === reference) : -1;
    // A reference measured once has a fitted effect made of one residual, so a difference against it
    // is noise wearing a comparator's name.
    const refCi = found >= 0 && fit.colN[found] >= 2 ? found : -1;
    // A series whose fit stopped short contributes to no pool, so it takes no slot either: an empty
    // strip entry would draw the "never tried" mark over a series that tried the group and was dropped.
    const cols: StripCols | null = strip && fit.converged ? new Map() : null;
    let bestCi = -1;
    let worstCi = -1;
    for (let c = 0; c < matrix.columns.length; c++) {
      const subst = matrix.columns[c].substSmiles;
      // Only a converged fit may be compared with another matrix's, and a column measured once has a
      // fitted effect made of one residual. The strip is filled under the same test, so a square and
      // the coverage sentence under it can never describe different sets of columns.
      if (fit.colN[c] < 2 || !fit.converged)
        continue;
      if (cols !== null) {
        // 0 where the column IS the reference — the honest difference against itself, which is not the
        // "no comparator at all" the null stands for.
        cols.set(subst, {ci: c, n: fit.colN[c],
          refDelta: refCi < 0 ? null : refCi === c ? 0 :
            dir * (fit.colEffect[c] - fit.colEffect[refCi])});
      }
      let attempts = acc.tried.get(subst);
      if (attempts === undefined)
        acc.tried.set(subst, attempts = new Set());
      attempts.add(root);
      if (bestCi < 0 || dir * fit.colEffect[c] > dir * fit.colEffect[bestCi] ||
        (fit.colEffect[c] === fit.colEffect[bestCi] && subst < matrix.columns[bestCi].substSmiles))
        bestCi = c;
      if (worstCi < 0 || dir * fit.colEffect[c] < dir * fit.colEffect[worstCi] ||
        (fit.colEffect[c] === fit.colEffect[worstCi] && subst < matrix.columns[worstCi].substSmiles))
        worstCi = c;
    }
    if (bestCi >= 0)
      this.pushExtreme(acc.wins, matrix, fit, dir, root, fitHolds, bestCi, refCi, 1);
    // One comparable column is first and last at once, which is not a loss.
    if (worstCi >= 0 && worstCi !== bestCi)
      this.pushExtreme(acc.losses, matrix, fit, dir, root, fitHolds, worstCi, refCi, -1);
    return cols;
  }

  /** Record one column as this lineage's extreme, keeping the occurrence with the stronger
   *  within-series margin. `keep` is +1 for the winner pool and −1 for the loser pool, so one dedup
   *  serves both. */
  private pushExtreme(pool: Map<string, RGroupWin[]>, matrix: SarMatrix, fit: AdditiveFit, dir: number,
    root: string, fitHolds: boolean, ci: number, refCi: number, keep: number): void {
    // The most potent measured cell of the column, so the row lands on a compound rather than on a
    // hole the reader has to hunt through.
    const ri = this.bestMeasuredRow(matrix, ci);
    if (ri < 0)
      return;
    const win: RGroupWin = {
      matrix, ri, ci, root, n: fit.colN[ci],
      refDelta: refCi >= 0 && refCi !== ci ? dir * (fit.colEffect[ci] - fit.colEffect[refCi]) : null,
      fitHolds,
    };
    const subst = matrix.columns[ci].substSmiles;
    const held = pool.get(subst);
    if (held === undefined) {
      pool.set(subst, [win]);
      return;
    }
    // One entry per lineage: a substituent extreme in a matrix and again in the tier that folded it has
    // been confirmed once, over one set of compounds, and must not read as two series.
    const at = held.findIndex((other) => other.root === root);
    if (at < 0)
      held.push(win);
    else if (this.strongerWin(win, held[at], keep))
      held[at] = win;
  }

  /**
   * Orders two occurrences of one substituent, best of `keep`'s pool first.
   *
   * Never on the fitted level itself: `colEffect` is centred over whichever substituents its own
   * matrix happened to make, so two matrices' levels have different zeros and picking a maximum over
   * them is arithmetic on two different gauges. The difference against each matrix's own reference has
   * no gauge in it; where a matrix records no reference, the measured cells behind the column decide.
   */
  private strongerWin(a: RGroupWin, b: RGroupWin, keep: number): boolean {
    if (a.refDelta !== null && b.refDelta !== null && a.refDelta !== b.refDelta)
      return keep * a.refDelta > keep * b.refDelta;
    if ((a.refDelta === null) !== (b.refDelta === null))
      return b.refDelta === null;
    return a.n !== b.n ? a.n > b.n : a.matrix.id < b.matrix.id;
  }

  private rankRGroups(acc: RGroupAcc):
    {rgroups: RGroupRow[], rgroupsThin: RGroupRow[], rgroupLosers: RGroupRow[]} {
    const winners = this.rankPool(acc.wins, acc.tried, 1);
    const losers = this.rankPool(acc.losses, acc.tried, -1);
    return {
      rgroups: winners.filter((row) => row.k >= RGROUP_MIN_SERIES).slice(0, SUM_ROWS + 2),
      rgroupsThin: winners.filter((row) => row.k < RGROUP_MIN_SERIES).slice(0, SUM_ROWS),
      // No thin tier for losers: a loss over two lineages is not a rule, and striking a building block
      // from the next enumeration is acted on without being checked.
      rgroupLosers: losers.filter((row) => row.k >= RGROUP_MIN_SERIES).slice(0, LOSER_ROWS),
    };
  }

  /** One pool ordered best-first, `keep` being +1 for wins and −1 for losses so the strongest of each
   *  pool comes out on top of its own list. */
  private rankPool(pool: Map<string, RGroupWin[]>, tried: Map<string, Set<string>>,
    keep: number): RGroupRow[] {
    const rows: RGroupRow[] = [];
    for (const [subst, entries] of pool) {
      const trusted = entries.filter((win) => win.fitHolds);
      if (trusted.length < SWAP_MIN_SERIES)
        continue;
      const deltas = trusted.map((win) => win.refDelta).filter((d): d is number => d !== null);
      const ranked = [...trusted].sort((a, b) => this.strongerWin(a, b, keep) ? -1 : 1);
      rows.push({
        subst, win: ranked[0], k: trusted.length, m: entries.length,
        tried: tried.get(subst)?.size ?? entries.length,
        magnitude: deltas.length ? median(deltas) : null,
        lo: deltas.length ? Math.min(...deltas) : 0,
        hi: deltas.length ? Math.max(...deltas) : 0,
        magnitudeSeries: deltas.length,
      });
    }
    // Rows with no comparator sort last in either pool: they still rank on how often the group placed.
    const strength = (row: RGroupRow): number =>
      row.magnitude === null ? Number.NEGATIVE_INFINITY : keep * row.magnitude;
    rows.sort((a, b) => b.k - a.k || strength(b) - strength(a) || (a.subst < b.subst ? -1 : 1));
    return rows;
  }

  /**
   * Five independent argmaxes, then merged by lineage root.
   *
   * Independent because a series holding the most compounds, the widest range AND the best compound is
   * the answer, and pushing it out of a taken slot would put a weaker series under a label it does not
   * earn. Merged because those are one place to look, not three.
   */
  private mergeStartHere(series: SeriesStat[], confN: number[], dir: number): StartRow[] {
    if (series.length === 0)
      return [];
    // Floored well above MIN_CV_POINTS: an R² ranked at the cross-validation floor is an extreme-value
    // draw from the noisiest end, on a row captioned "best validated".
    const nFloor = Math.max(BEST_FIT_MIN_N, median(confN));
    const pick = (score: (s: SeriesStat) => number | null): SeriesStat | null => {
      let winner: SeriesStat | null = null;
      let bestScore = 0;
      for (const stat of series) {
        const value = score(stat);
        if (value === null)
          continue;
        if (winner === null || value > bestScore ||
          (value === bestScore && stat.matrix.id < winner.matrix.id)) {
          winner = stat;
          bestScore = value;
        }
      }
      return winner;
    };
    // Index-aligned with REASON_GLYPHS / REASON_WORDS: the lane is a fixed set of slots, so the order
    // here is the order on screen.
    const slots: (SeriesStat | null)[] = [
      pick((s) => s.cpd || null),
      pick((s) => s.realCells ? s.hi - s.lo : null),
      pick((s) => s.best === null ? null : dir * s.best.value),
      pick((s) => s.trusted || null),
      pick((s) => s.matrix.confidence && s.matrix.confidence.n >= nFloor ? s.matrix.confidence.r2 : null),
    ];
    const byRoot = new Map<string, StartRow>();
    for (let slot = 0; slot < slots.length; slot++) {
      const stat = slots[slot];
      if (stat === null)
        continue;
      const held = byRoot.get(stat.root);
      if (held === undefined)
        byRoot.set(stat.root, {primary: stat, reasons: [slot]});
      else
        held.reasons.push(slot);
    }
    return [...byRoot.values()];
  }

  // ---- Layout ---------------------------------------------------------------------------------

  /**
   * @param keepPlace Keep the segment and sub-tab the reader is on. Set when the same analysis is
   * re-read — changing the fold tier — and left alone for a new one.
   */
  private render(keepPlace = false): void {
    const data = this.data!;
    const host = this.host;
    window.clearTimeout(this.paintTimer);
    this.pendingPaints = [];
    this.transferLine = null;
    this.closeAnalogGrid();
    this.segButtons.clear();
    this.paneHost = null;
    this.scroller = null;
    // A new analysis is a new landing: a reader returning to a different run must not arrive mid-way
    // into someone else's orientation. Re-reading the same one at another tier is the opposite case —
    // the reader asked a question about the segment they are looking at and has to stay on it.
    const pane = keepPlace ? this.currentPane : PANE_OVERVIEW;
    if (!keepPlace) {
      this.expandTopRGroup = false;
      this.effectsTab = 0;
      this.openTrustList = false;
    }
    this.currentPane = PANE_OVERVIEW;
    ui.empty(this.root);

    if (host.matrices.length === 0) {
      // No segment bar: four panes over nothing is four ways to read one message.
      this.root.appendChild(ui.divText(host.noMatricesMessage(), 'chem-sar-empty-note'));
      return;
    }

    this.root.appendChild(this.buildScaleBand(data));
    this.root.appendChild(this.buildSegBar(data));
    this.paneHost = ui.div([], 'chem-sar-sum-pane');
    this.root.appendChild(this.paneHost);
    // Before the first build, not after: a pane decides here whether it draws a mark or its text
    // fallback, and it can only read a class that is already on the root.
    this.widthKey = '';
    this.applySize();
    this.showPane(pane);

    if (this.sizeSub === null) {
      // The pane is docked, so its width is not the viewport's: a media query would key on the wrong
      // number entirely.
      this.sizeSub = DG.debounce(ui.onSizeChanged(this.root), 150).subscribe(() => this.applySize());
    }
  }

  /**
   * The breakpoint classes, and a rebuild of the pane on screen whenever a width one crosses.
   *
   * Height is read as well as width: the Overview's no-scroll look is an `overflow: hidden`, and a
   * dock that is wide but short would otherwise clip its bottom row with no scrollbar to say so.
   * The rebuild is what keeps a mark and its text fallback from both vanishing — CSS can hide the
   * reason lane, only a rebuild can put the phrases back.
   */
  private applySize(): void {
    const w = this.root.clientWidth;
    if (w === 0)
      return;
    const narrow = this.root.classList.toggle('chem-sar-sum-narrow', w < NARROW_PX);
    const xnarrow = this.root.classList.toggle('chem-sar-sum-xnarrow', w < XNARROW_PX);
    this.root.classList.toggle('chem-sar-sum-short', this.root.clientHeight < SHORT_PX);
    const key = `${narrow}|${xnarrow}`;
    const crossed = this.widthKey !== '' && this.widthKey !== key;
    this.widthKey = key;
    if (crossed && this.paneHost !== null)
      this.showPane(this.currentPane);
  }

  /**
   * The second level of tabs, built from plain divs rather than `ui.tabControl`.
   *
   * A nested tab control renders the platform's own tab chrome 26px under the viewer's, and two
   * identical strips is exactly the confusion a second level has to avoid. This differs in shape, in
   * active state, in size and position, and in carrying a readout on its right — which a tab strip
   * never does and a toolbar always does.
   */
  private tierChip(tier: number | null, label: string, n: number): HTMLElement {
    // The tier name goes in its own element so a reader or a test can address it without the count
    // beside it, which changes with the data.
    const chip = ui.divH([ui.divText(label, 'chem-sar-sum-tier-name'),
      ui.divText(`· ${count(n)}`, 'chem-sar-sum-tier-n')], 'chem-sar-chip-badge chem-sar-sum-role');
    chip.classList.toggle('chem-sar-sum-role-on', this.tierFilter === tier);
    ui.tooltip.bind(chip, () => (tier === null ?
      `Read every ranking on this tab over all ${count(n)} series at once. A compound cut at several ` +
      'breadths is then counted at each of them.' :
      `Read every ranking on this tab over the ${count(n)} series cut at level ${tier} only — the ` +
      'components, the swaps, the cores, Start here and Worth making.') +
      ' L1 is a leaf series; each level above it holds the same compounds over cores that agree one ' +
      'further cut deeper. Nothing is rebuilt: the series already exist either way.');
    chip.onclick = () => this.setTierFilter(tier);
    return chip;
  }

  private buildSegBar(data: SummaryData): HTMLElement {
    const segs = PANES.map((name) => {
      const el = ui.divText(name, 'chem-sar-sum-seg');
      el.onclick = () => this.showPane(name);
      this.segButtons.set(name, el);
      return el;
    });
    // On the segment bar rather than on one pane: the tier decides what every segment is reading, so it
    // has to be visible and changeable from all of them. A decomposition handed over as columns is all
    // one tier, and there the chips would be a control with one setting.
    const readout: HTMLElement[] = [];
    if (data.tierCounts.length > 1) {
      const total = data.tierCounts.reduce((sum, {n}) => sum + n, 0);
      readout.push(this.tierChip(null, 'All', total));
      for (const {tier, n} of data.tierCounts)
        readout.push(this.tierChip(tier, `L${tier}`, n));
    } else {
      const line = ui.divText(`${count(this.host.matrices.length)} series · ` +
        `${count(data.families)} families`, 'chem-sar-sum-seg-readout');
      ui.tooltip.bind(line, 'Matrices over fold lineages. The same compounds re-cut, not this ' +
        'many findings.');
      readout.push(line);
    }
    return ui.divH([ui.divH(segs, 'chem-sar-sum-seg-group'),
      ui.divH(readout, 'chem-sar-sum-seg-tiers')], 'chem-sar-sum-seg-bar');
  }

  /** Build one segment and discard the last. Panes are never held: an unvisited Worth-making segment
   *  then never creates a DataFrame, and a revisited one is rebuilt from data that has not changed. */
  private showPane(name: string, anchor?: string): void {
    const paneHost = this.paneHost;
    const data = this.data;
    if (paneHost === null || data === null)
      return;
    // Unconditional: re-opening the segment that holds it builds a second grid, and only closing the
    // first releases the Dart-backed object behind it.
    this.closeAnalogGrid();
    this.currentPane = name;
    for (const [key, el] of this.segButtons)
      el.classList.toggle('chem-sar-sum-seg-on', key === name);
    window.clearTimeout(this.paintTimer);
    this.pendingPaints = [];
    this.scroller = null;
    ui.empty(paneHost);
    // An anchor on this segment names a card, and a card now lives on a tab rather than at a scroll
    // offset, so it selects the tab before the pane is built.
    if (name === PANE_EFFECTS && anchor !== undefined)
      this.effectsTab = this.effectsTabOf(data, anchor);
    paneHost.classList.toggle('chem-sar-sum-pane-fixed', name === PANE_OVERVIEW);
    paneHost.appendChild(name === PANE_EFFECTS ? this.buildEffectsPane(data) :
      name === PANE_MAKING ? this.buildMakingPane(data) :
        name === PANE_METHOD ? this.buildMethodPane(data) : this.buildOverviewPane(data));

    this.paintTimer = window.setTimeout(() => {
      this.flushPaints();
      if (anchor === undefined)
        return;
      const target = paneHost.querySelector(`[data-anchor="${anchor}"]`);
      if (target instanceof HTMLElement)
        this.scrollTo(target);
    }, 0);
  }

  /** What the transfer scan has found, said in the one line that offers it. Detection is lazy and runs
   *  on its own tab, so this is read at activation rather than when the pane was built. */
  private syncTransferLine(): void {
    const line = this.transferLine;
    if (line === null)
      return;
    const scan = this.host.transferSummary;
    line.textContent = !scan.scanned ?
      'SAR transfers — pairs of cores whose trends run in parallel · open the tab to detect →' :
      scan.count === 0 ?
        'No SAR transfers: no two cores share enough R-groups, on scaffolds alike enough · open the ' +
        'tab to change the similarity →' :
        `${count(scan.count)} SAR transfers — pairs of cores whose trends run in parallel · open the ` +
        'tab →';
  }

  /** Switch to the methodology segment with the trust list on screen — where every "the fit does not
   *  hold here" on the tab resolves to. */
  private revealTrust(): void {
    this.openTrustList = true;
    this.showPane(PANE_METHOD, 'trust');
  }

  private anchored(el: HTMLElement, anchor: string): HTMLElement {
    el.dataset.anchor = anchor;
    return el;
  }

  // ---- The four panes ---------------------------------------------------------------------------

  private buildOverviewPane(data: SummaryData): HTMLElement {
    const cores = this.topCores(data);
    const list = ui.divV(this.buildStartHere(data), 'chem-sar-sum-list');
    this.scroller = list;
    return ui.divV([
      this.coverageBar(data),
      this.setupLine(data),
      this.buildAnswers(data, cores),
      this.buildListHeader(data),
      list,
      this.buildTotals(data),
    ], 'chem-sar-sum-overview');
  }

  /**
   * One tab per component column, then the per-series evidence.
   *
   * Stacked, this segment is as tall as the decomposition is wide — five component columns is five
   * cards before the two that were already here — and a landing screen that has to be scrolled to be
   * read is not one. The ordering sentence stays above the strip because it is the only finding here
   * that is about every tab at once, and each tab carries its own spread so the comparison between
   * components does not cost a click each.
   */
  private buildEffectsPane(data: SummaryData): HTMLElement {
    // The R-group card is dropped wherever the fit already ranks the axis role, for the reason the core
    // card is: two rankings of one column, one a count of within-series wins and one an adjusted effect,
    // read as a disagreement rather than as two kinds of evidence.
    const series: HTMLElement[] = this.roleAnswer(data, data.axisRole) === null ?
      [this.anchored(this.buildRGroupCard(data), 'rgroup')] : [];
    // Dropped only where the core column is already ranked by the fit, adjusted for the components it
    // was paired with: two rankings of one column, one confounded and one adjusted, in different units,
    // is a worse screen than either alone. Where the fit declines to rank it, this is the only one left.
    let coreNote: HTMLElement | null = null;
    if (this.roleAnswer(data, data.coreRole) === null) {
      if (data.coresAreSeries)
        series.push(this.anchored(this.buildCoreCard(data, this.topCores(data)), 'cores'));
      else {
        // Where cores are not comparable the card has no ranking it could ever carry, and a card slot
        // spent on one sentence is a slot. It keeps the anchor so the tile that points here still lands.
        coreNote = this.anchored(ui.divText(CORES_NOT_COMPARABLE,
          'chem-sar-cp-hint chem-sar-sum-sub-lead'), 'cores');
      }
    }
    // One card per component, matching the band that links here. A single card carried whichever
    // component the matrix columns enumerated, so clicking the Linker row landed on a card about
    // Warhead — the band's own answer contradicted on arrival.
    const roles = data.roleFit === null ? [] : data.roleFit.roles.map((role) => role.name);
    if (roles.length === 0)
      series.push(this.anchored(this.buildSwapCard(data), 'swaps'));
    else {
      for (const name of roles)
        series.push(this.anchored(this.buildSwapCard(data, name), swapAnchor(name)));
    }

    const roleCards = this.buildRoleCards(data);
    if (roleCards.length === 0) {
      return coreNote === null ? this.effectsScroll(series) :
        ui.divV([coreNote, this.effectsScroll(series)], 'chem-sar-sum-effects');
    }

    const tabs = [...roleCards.map((card) => ({label: card.label, note: card.note, cards: [card.el]})),
      {label: PANE_SERIES, note: 'counted, not fitted' as string | null, cards: series}];
    const open = Math.min(Math.max(this.effectsTab, 0), tabs.length - 1);
    const bar = ui.divH(tabs.map((tab, i) => {
      const parts = [ui.divText(tab.label, 'chem-sar-sum-sub-name')];
      if (tab.note !== null)
        parts.push(ui.divText(tab.note, 'chem-sar-sum-sub-note'));
      const el = ui.divH(parts, 'chem-sar-sum-sub');
      el.classList.toggle('chem-sar-sum-sub-on', i === open);
      el.onclick = () => {
        this.effectsTab = i;
        this.showPane(PANE_EFFECTS);
      };
      return el;
    }), 'chem-sar-sum-sub-bar');

    const lead = this.roleOrdering(data);
    return ui.divV([...(lead === null ? [] : [lead]), bar, this.effectsScroll(tabs[open].cards)],
      'chem-sar-sum-effects');
  }

  private effectsScroll(cards: HTMLElement[]): HTMLElement {
    const grid = ui.div(cards, 'chem-sar-sum-grid');
    // A lone card has no neighbour to sit beside, so left at one track's width it strands the rest of
    // the pane. The width is not decoration here: these row names are SMILES, and at one track three
    // analogs of one scaffold truncate to the same prefix and read as the same structure.
    grid.classList.toggle('chem-sar-sum-grid-one', cards.length === 1);
    const scroll = ui.div([grid], 'chem-sar-sum-scroll');
    this.scroller = scroll;
    return scroll;
  }

  /** The Effects tab an anchor lives on: a role anchor carries its own index, and everything else is
   *  on the per-series tab, which sits after the role cards. */
  private effectsTabOf(data: SummaryData, anchor: string): number {
    const role = /^role-(\d+)$/.exec(anchor);
    if (role !== null)
      return Number(role[1]);
    if (data.axisRole === null)
      return 0;
    return this.roleFitRefusal(data) === null ? data.roleFit!.roles.length : 1;
  }

  /** Not a scroll pane: the ranked list is a DG.Grid, which needs a definite height to virtualize
   *  against, and a grid inside a scrolling parent resolves to its own minimum however tall the dock
   *  is. The shelf above it keeps a bounded share and scrolls inside it. */
  private buildMakingPane(data: SummaryData): HTMLElement {
    return ui.divV([this.buildShelfSegment(data), this.buildAnalogBlock(data)], 'chem-sar-sum-making');
  }

  private buildMethodPane(data: SummaryData): HTMLElement {
    const parts: HTMLElement[] = [this.buildChips(data),
      this.anchored(this.buildTrustSection(data), 'trust')];
    const scroll = ui.divV(parts, 'chem-sar-sum-scroll');
    this.scroller = scroll;
    return scroll;
  }

  /** The five reason-lane glyphs printed once, directly above the lane they label and always on
   *  screen — a mark whose legend is a scroll away is a legend hunt. */
  private buildListHeader(data: SummaryData): HTMLElement {
    const label = ui.divText('Start here', 'chem-sar-sum-card-title');
    if (!this.lanesFit(data))
      return ui.divH([label], 'chem-sar-sum-list-head');
    const keys = REASON_GLYPHS.map((glyph, i) =>
      ui.divH([ui.divText(glyph, 'chem-sar-sum-slot chem-sar-sum-slot-on'),
        ui.divText(REASON_WORDS[i], 'chem-sar-sum-legend-word')], 'chem-sar-sum-legend-key'));
    return ui.divH([label, ui.divH(keys, 'chem-sar-sum-legend')], 'chem-sar-sum-list-head');
  }

  /** Three totals that each open the segment answering them. Magnitudes that only describe the run
   *  stay on Method: a landing screen carries decisions, not census figures. */
  private buildTotals(data: SummaryData): HTMLElement {
    const host = this.host;
    const ranked = data.analogs.all.length;
    const thin = data.analogsThin.all.length;
    const gated = data.analogs.all.concat(data.analogsThin.all);
    const shown = gated.length > 0 ? gated : data.analogsAny.all;
    const structures = new Set(shown.map((r) => r.key)).size;
    // A capped list printed as a total states the cap as the finding, so the cap is named wherever the
    // length is.
    const capped = gated.length > 0 && data.analogOverflow > 0 ? `top ${count(ranked + thin)} of ` +
      `${count(ranked + thin + data.analogOverflow)}` : count(shown.length);
    const making = ui.divText(shown.length === 0 ? 'Nothing is predicted above what is already made →' :
      gated.length === 0 ? `${capped} best candidates, none past the trust gate →` :
        `${capped} worth making · ${count(structures)} distinct structures →`, 'chem-sar-sum-total');
    making.onclick = () => this.showPane(PANE_MAKING);
    ui.tooltip.bind(making, () => 'Predicted analogs this dataset has no row for, ranked on the gain ' +
      'they buy over the best compound their own series has already made.' + (data.analogOverflow === 0 ?
      '' : ` The ranked and the thin list hold ${count(ANALOG_LIST_MAX)} rows each: ` +
      `${count(data.analogOverflow)} further structures cleared the same gate and are in neither.`));

    const trust = ui.divText(`${count(data.fitHolds)} fits hold · ${count(data.unchecked)} unchecked · ` +
      `${count(data.lowR2.length)} do not →`, 'chem-sar-sum-total');
    trust.onclick = () => this.revealTrust();
    ui.tooltip.bind(trust, 'A verdict on each series\' fit, not a partition of your library — one ' +
      'compound sits in several series and is routinely on both sides. Unchecked means unverified, ' +
      'not wrong.');

    const parts = [making, trust];
    if (host.unscalableCount > 0) {
      const bad = ui.divText(`${count(host.unscalableCount)} values cannot be scaled by ` +
        `${host.scalingLabel}`, 'chem-sar-sum-total chem-sar-chip-partial');
      ui.tooltip.bind(bad, () => `${host.unscalableCount} of the "${host.activityColumnName}" values ` +
        `were assayed but cannot be scaled by ${host.scalingLabel}, which needs positive numbers. ` +
        'They are excluded from every matrix — set Scaling to "none" to use them as they are.');
      parts.push(bad);
    }
    return ui.divH(parts, 'chem-sar-sum-totals');
  }

  // ---- Band A: scale, then method ---------------------------------------------------------------

  /**
   * The two things that invalidate every ranking below, and nothing else — this band does not scroll,
   * so anything pinned here is width the answers never get back.
   */
  private buildScaleBand(data: SummaryData): HTMLElement {
    const host = this.host;
    const lines: HTMLElement[] = [];
    const line = (text: string): HTMLElement => ui.divText(text, 'chem-sar-sum-orient-line');

    const transform = host.scalingLabel === 'raw' ? 'untransformed' : `${host.scalingLabel} applied`;
    const direction = host.higherIsBetter ? 'higher is better' : 'lower is better';
    // Never printed as 0 when no series has a fit: a floor of zero reads as "everything is resolved".
    const floor = data.modelError === null ? '' :
      ` · ± ${host.formatActivity(data.modelError)} model error`;
    // Log-ness of a raw column is inferred from the declared direction, not read off the data, and
    // every fold and log-unit claim on the tab rests on it. Percent inhibition, ΔTm and ΔG are raw and
    // higher-is-better too, so the inference has to be on screen rather than in a getter.
    // Stated on the hover rather than on the band: it qualifies the fold figures further down, and on
    // the band it read as one more property of the column.
    const assumedLog = host.scalingLabel === 'raw' && host.activityIsLog;
    const scaleLine = line(`${host.activityColumnName || 'Activity'} · ${transform} · ${direction}`);
    ui.tooltip.bind(scaleLine, () => 'The scale and direction every ranking here depends on.' +
      (!assumedLog ? '' : ' Fold changes on this tab assume the column is already a log, inferred from ' +
      'its being untransformed and higher-is-better — a percent inhibition, ΔTm or ΔG column is the ' +
      'same shape and is not a log.') +
      (data.modelError === null ? '' : ' Model error, not assay error: no replicates here, so it is an ' +
      'upper bound on the assay.'));
    const head: HTMLElement[] = [scaleLine];
    if (data.minObserved !== null)
      head.push(this.rangeRule(data));
    else
      head.push(line('nothing observed'));
    if (floor !== '')
      head.push(line(floor.replace(/^ · /, '')));
    lines.push(ui.divH(head, 'chem-sar-sum-orient-row'));
    // Only where this analysis did transform: a column it left alone is on the scale it was measured
    // on, and log permeability or ΔG is negative throughout without anything being wrong.
    if (data.maxObserved !== null && data.maxObserved < 0 && host.scalingLabel !== 'raw') {
      // Redundant with the rule's own zero mark, and deliberately so: a double transform inverts every
      // ranking on the tab, which is the one consequence worth saying twice.
      lines.push(ui.divText('Every observed value is negative — check Scaling.',
        'chem-sar-sum-orient-line chem-sar-sum-alarm'));
    }
    return ui.divV(lines, 'chem-sar-sum-orient');
  }

  // ---- Marks ------------------------------------------------------------------------------------

  /**
   * The observed range as a rule, spanning exactly lo..hi. It must not extend to zero: on a column that
   * is negative throughout — log solubility, ΔG — the right-hand label would sit at the track's end
   * where the value is zero and not the label's number.
   */
  private rangeRule(data: SummaryData): HTMLElement {
    const host = this.host;
    const lo = data.minObserved!;
    const hi = data.maxObserved!;
    const span = hi - lo || 1;
    const track = ui.div([], 'chem-sar-sum-rule');
    const bar = ui.div([], 'chem-sar-sum-rule-bar');
    bar.style.left = '0%';
    bar.style.width = '100%';
    track.appendChild(bar);
    // Only where the measurements actually straddle it, which is the only place the mark can be read
    // off the rule at all.
    const showZero = lo < 0 && hi > 0;
    if (showZero) {
      const zero = ui.div([], 'chem-sar-sum-rule-zero');
      zero.style.left = `${((0 - lo) / span) * 100}%`;
      track.appendChild(zero);
    }
    const box = ui.divH([
      ui.divText(host.formatActivity(lo), 'chem-sar-sum-rule-end'),
      track,
      ui.divText(host.formatActivity(hi), 'chem-sar-sum-rule-end'),
    ], 'chem-sar-sum-rule-box');
    ui.tooltip.bind(box, () => `Observed ${host.formatActivity(lo)} to ${host.formatActivity(hi)} over ` +
      `${count(data.measuredCells)} measured cells. The rule spans what was measured, end to end.` +
      (showZero ? ' The tick is zero.' : ''));
    return box;
  }

  /**
   * MK-cov: what fraction of the table the analysis sees, in four segments in a fixed order, over a
   * label line naming them in the same order — so the label is the legend.
   */
  private coverageBar(data: SummaryData): HTMLElement {
    const host = this.host;
    const total = Math.max(1, host.hostRowCount);
    // The two ways a row misses out are split because they call for different things: one is a plate
    // to run, the other a library too diverse for any core to recur. Lumped, a sparse assay and a
    // diverse library read the same.
    const unplaced = Math.max(0, host.assayedCount - data.compounds - host.unscalableCount);
    const never = Math.max(0, total - host.assayedCount - data.untested);
    const segs = [
      {n: data.compounds, word: 'measured', cls: 'chem-sar-sum-cov-real'},
      {n: data.untested, word: 'held, never assayed', cls: 'chem-sar-sum-cov-held'},
      {n: host.unscalableCount, word: 'assayed, cannot be scaled', cls: 'chem-sar-sum-cov-bad'},
      {n: unplaced, word: 'assayed, no analog to pair with', cls: 'chem-sar-sum-cov-out'},
      {n: never, word: 'never assayed', cls: 'chem-sar-sum-cov-none'},
      // A zero-count segment is omitted rather than drawn at the minimum width, which would state a
      // quantity that is not there.
    ].filter((s) => s.n > 0);
    const bar = ui.divH(segs.map((s) => {
      const el = ui.div([], `chem-sar-sum-cov-seg ${s.cls}`);
      el.style.flexGrow = `${s.n}`;
      return el;
    }), 'chem-sar-sum-cov');
    const label = ui.divText(segs.map((s) => `${count(s.n)} ${s.word}`).join(' · '),
      'chem-sar-sum-cov-label');
    const box = ui.divV([bar, label], 'chem-sar-sum-cov-box');
    ui.tooltip.bind(box, () => `${count(total)} table rows. Rows, not compounds: the denominator is the ` +
      'table\'s raw row count, so this ignores any filter on the source grid. A segment narrower than ' +
      'a couple of pixels is drawn at that floor and stops being proportional, which is why every ' +
      'count is printed.');
    return box;
  }

  /**
   * MK-D: a signed magnitude against the screen's own model error.
   *
   * `scale` is the largest magnitude in the set this bar is compared against — the card's own rows,
   * except where one fit produced the rows of several cards and they are genuinely on one scale. A bar
   * ending inside the grey band is one the analysis cannot resolve.
   */
  private effectBar(value: number, scale: number, modelError: number | null): HTMLElement {
    const track = ui.div([], 'chem-sar-sum-effbar');
    const half = scale > 0 ? scale : 1;
    // No band at all rather than a zero-width one: a zero-width band reads as "everything here is
    // resolved", which is the opposite of what a missing model error means.
    if (modelError !== null && modelError > 0) {
      const band = ui.div([], 'chem-sar-sum-effband');
      const half2 = Math.min(50, (modelError / half) * 50);
      band.style.left = `${50 - half2}%`;
      band.style.width = `${2 * half2}%`;
      track.appendChild(band);
    }
    track.appendChild(ui.div([], 'chem-sar-sum-effzero'));
    const w = Math.min(50, (Math.abs(value) / half) * 50);
    const bar = ui.div([], `chem-sar-sum-effval chem-sar-sum-eff-${value >= 0 ? 'pos' : 'neg'}`);
    bar.style.left = `${value >= 0 ? 50 : 50 - w}%`;
    bar.style.width = `${w}%`;
    track.appendChild(bar);
    ui.tooltip.bind(track, () => `${this.formatEffect(value)}, against the widest bar it is compared ` +
      `with (${this.formatEffect(half)}).` + (modelError === null ? ' No series has a cross-validated fit, ' +
      'so there is no error band to read it against.' :
      ` The grey band is ± ${this.host.formatActivity(modelError)}, this analysis\' own typical ` +
      'prediction error: a bar ending inside it is a difference this analysis cannot resolve.'));
    return track;
  }

  /** MK-E: whether this series' fit was checked, and whether it held. No colour — `confidence` is null
   *  where too few cells could be cross-validated or every measured value was identical, and a traffic
   *  light would make "unchecked" read as a failing grade. */
  private trustDot(matrix: SarMatrix): HTMLElement {
    const conf = matrix.confidence;
    const dot = ui.div([], 'chem-sar-sum-dot' + (!conf ? ' chem-sar-sum-dot-none' :
      conf.r2 < TRUST_R2 ? ' chem-sar-sum-dot-open' : ''));
    ui.tooltip.bind(dot, () => !conf ?
      'Unchecked: this series has fewer than four measured cells that could be held out, or every ' +
      'measured value in it is identical. Its predictions are unverified, not wrong.' :
      `R² ${conf.r2.toFixed(2)} ± ${this.host.formatActivity(conf.rmse)} predicting cells held out of ` +
      `the fit, over ${conf.n} of ${conf.total} measured — the fit ` +
      `${conf.r2 >= TRUST_R2 ? 'holds' : 'does not hold'} at R² ${TRUST_R2}. A property of the model, ` +
      'not of the data.');
    return dot;
  }

  /** MK-F: how many of a series' own typical prediction errors a predicted gain is worth. Unfilled
   *  slots are drawn, not omitted, so "one of four" never looks like "one of one". */
  private tickRun(gain: number, rmse: number | null): HTMLElement {
    const filled = rmse !== null && rmse > 0 ?
      Math.max(0, Math.min(GAIN_TICKS, Math.floor(gain / rmse))) : 0;
    const ticks: HTMLElement[] = [];
    for (let i = 0; i < GAIN_TICKS; i++)
      ticks.push(ui.div([], `chem-sar-sum-tick${i < filled ? ' chem-sar-sum-tick-on' : ''}`));
    const run = ui.divH(ticks, 'chem-sar-sum-ticks');
    ui.tooltip.bind(run, () => rmse === null || rmse <= 0 ?
      'This series has no cross-validated error to measure the gain against.' :
      `The predicted gain is ${(gain / rmse).toFixed(1)} times this series' own leave-one-out error. ` +
      'Under one error the model cannot tell this analog from the best compound already made.');
    return run;
  }

  /** MK-B: the five fixed reasons a series is worth opening, won or not won, always in the same order
   *  — the lane's shape is what a reader compares down the column, which five variable-length phrases
   *  can never be. */
  private reasonLane(reasons: number[]): HTMLElement {
    const won = new Set(reasons);
    const slots = REASON_GLYPHS.map((glyph, i) => ui.divText(won.has(i) ? glyph : '',
      `chem-sar-sum-slot${won.has(i) ? ' chem-sar-sum-slot-on' : ''}`));
    const lane = ui.divH(slots, 'chem-sar-sum-lane');
    ui.tooltip.bind(lane, () => reasons.map((i) => REASON_WORDS[i]).join(' · '));
    return lane;
  }

  /** MK-C: one slot per series that tried this group — solid where it won and that series' fit holds,
   *  hatched where it won and the fit does not, outlined where it was tried and did not win. */
  private winTally(row: RGroupRow, placed: string): HTMLElement {
    const shown = Math.min(row.tried, 10);
    const slots: HTMLElement[] = [];
    for (let i = 0; i < shown; i++) {
      const state = i < row.k ? ' chem-sar-sum-slot-on' :
        i < row.m ? ' chem-sar-sum-slot-hatch' : '';
      const slot = ui.div([], `chem-sar-sum-wslot${state}`);
      if (i >= row.k && i < row.m) {
        slot.onclick = (e: MouseEvent) => {
          e.stopPropagation();
          this.revealTrust();
        };
      }
      slots.push(slot);
    }
    const parts: HTMLElement[] = [ui.divH(slots, 'chem-sar-sum-tally')];
    if (row.tried > shown)
      parts.push(ui.divText(`+${row.tried - shown}`, 'chem-sar-sum-tally-more'));
    const box = ui.divH(parts, 'chem-sar-sum-tally-box');
    ui.tooltip.bind(box, () => `${placed} in ${row.k} of the ${row.tried} series that tried it` +
      (row.m > row.k ? `, and in ${row.m - row.k} more whose additive fit does not hold — a group that ` +
      `is ${placed} in additive series and not in non-additive ones is core-dependent. Click a hatched ` +
      'slot to see those series.' : '. Every fit behind those placings holds.'));
    return box;
  }

  private buildChips(data: SummaryData): HTMLElement {
    const chip = (text: string, tip: string, cls = ''): HTMLElement => {
      const el = ui.divText(text, `chem-sar-chip-badge ${cls}`.trim());
      ui.tooltip.bind(el, () => tip);
      return el;
    };
    const items = [
      chip(`${count(data.trustedCells)} predictions pass the trust gate · ` +
        `${count(data.trustedStructures)} distinct structures`,
      // The second constant is concatenated rather than interpolated: the bundler folds two adjacent
      // template operands that each carry a compile-time constant and drops the first one's trailing
      // text, which loses a clause without failing the build.
      `Gate: ${MIN_SUPPORT} supporting compounds on both axes, in a series whose own R² ≥ ` +
      TRUST_R2 + '. The distinct-structure count is the actionable one — it deduplicates the same ' +
      'analog predicted in several tiers.'),
    ];
    if (data.trustedNoStructure > 0) {
      items.push(chip(`${count(data.trustedNoStructure)} predictions have a value and no structure`,
        'These pass the same trust gate as the analogs in Do next, but the core carries an attachment ' +
        'point none of the picked fragment columns fills, so no structure can be completed over it. ' +
        'Add that column to the R-group columns and they appear in Do next.', 'chem-sar-chip-partial'));
    }
    return ui.divH(items, 'chem-sar-sum-chips');
  }

  // ---- Band A′: the three answers -------------------------------------------------------------

  /**
   * What the analysis answers, in words: which R-group, which core, and which swap buys the most.
   *
   * Every number here is already on a card below, so a reader who stops at this band has the finding
   * and one who does not finds the evidence under it. Each tile carries its own negative, because a
   * leaderboard with no stated limit is read as a ranking of the chemistry rather than of what was made.
   */
  private buildAnswers(data: SummaryData, cores: SeriesStat[]): HTMLElement {
    const findings = this.buildFindings(data);
    if (findings !== null)
      return findings;
    return ui.divH([
      this.answerTile(`Best ${data.axisRole ?? 'R-group'}`, this.rgroupAnswer(data), 'rgroup',
        () => this.expandTopRGroup = true),
      this.answerTile(`Best ${data.coreRole ?? 'core'}`, this.coreAnswer(data, cores), 'cores'),
      this.answerTile('What improves potency', this.potencyAnswer(data), 'swaps'),
    ], 'chem-sar-sum-answers');
  }

  /**
   * One row per component column and one for the best measured swap, in a single full-width band.
   *
   * Every component, not two of them: tiles for the axis and the core left the rest to be hunted for,
   * which sent readers to the axis switch on Method — a rebuild of every matrix — to see a ranking that
   * was already computed. One band rather than a strip beside a tile, because the tile's column was too
   * narrow for a structure and fell back to a truncated SMILES, which reads identically for every analog
   * of one scaffold.
   *
   * Null outside fragment-columns mode and wherever the fit is refused — a fragmented substituent label
   * is local to its series, so there is no component to name and the three tiles stay.
   */
  private buildFindings(data: SummaryData): HTMLElement | null {
    const fit = data.roleFit;
    if (fit === null || this.roleFitRefusal(data) !== null)
      return null;
    const rows = fit.roles.map((role, index) => this.componentRow(role, index));
    const body: HTMLElement[] = [
      ui.divText(`spans = the range of ${this.host.activityColumnName} across that component's ` +
        `values · offsets are against the library mean of ${this.host.formatActivity(fit.mean)}`,
      'chem-sar-cp-hint'),
      ...rows.slice(0, FINDING_ROWS),
    ];
    // A six-component table would otherwise push the swap row, and the band under it, off the screen the
    // reader lands on. The ranking is by spread, so what folds away is what moves the endpoint least.
    if (rows.length > FINDING_ROWS)
      body.push(this.moreComponents(rows.slice(FINDING_ROWS)));
    // One row per component, not one for whichever component the matrix columns happened to enumerate:
    // a matched pair is a grouping on the other components, so every component's pairs are in hand.
    const swapRows = fit.roles.map((role) => this.swapRow(data, role.name));
    const withPairs = fit.roles.filter((role) => (data.swapsByRole.get(role.name) ?? []).length > 0);
    // Each half folds to its own headline, so collapsing one to read the other never costs the answer.
    const lead = fit.roles[0];
    return ui.divV([
      this.bandGroup('comp', `What to change · ${count(fit.roles.length)} components`,
        `${lead.name} leads, spans ${lead.spread.toFixed(2)}`, ui.divV(body)),
      this.bandGroup('swap', `Best measured swap · ${count(fit.roles.length)} components`,
        withPairs.length === 0 ? 'none pooled' :
          `${withPairs[0].name} ${this.swapHeadline(data, withPairs[0].name)}`, ui.divV(swapRows)),
    ], 'chem-sar-sum-comp');
  }

  /**
   * Which column is the core, which one the matrix columns enumerate, and which fold into the row.
   *
   * The band below names components without saying what each one is to the analysis, and the three are
   * not interchangeable: the core is the scaffold every row is drawn from, the one across the columns is
   * the only one whose values sit side by side in the grid, and the rest are part of a row's identity.
   */
  private setupLine(data: SummaryData): HTMLElement {
    const host = this.host;
    if (data.axisRole === null) {
      return this.hint('Cores and substituents found by cutting the molecules — no component columns ' +
        'were given.', 'Every core here is a fragment the cutting produced, so its substituent labels ' +
        'are local to the series they came from and mean nothing across series.');
    }
    const folded = host.roleColumns.filter((name) => name !== data.axisRole);
    const parts = [`core: ${data.coreRole}`, `across the matrix columns: ${data.axisRole}`];
    if (folded.length > 0)
      parts.push(`folded into the row: ${folded.join(', ')}`);
    return this.hint(parts.join(' · '),
      'The core is the scaffold every row is drawn from. The component across the columns is the one ' +
      'whose values sit side by side in the SAR Matrix grid; the rest are part of what identifies a ' +
      'row. Every component is ranked and has its swaps pooled either way — the choice only decides ' +
      'the grid\'s layout.');
  }

  /** One foldable half of the landing band: a heading that states its own answer when shut, so the band
   *  can be narrowed to the half the reader wants without hiding what the other half concluded. */
  private bandGroup(key: string, title: string, headline: string, body: HTMLElement): HTMLElement {
    const open = !this.foldedBands.has(key);
    const lead = ui.divText(headline, 'chem-sar-cp-hint');
    lead.style.display = open ? 'none' : '';
    const head = this.foldHead(title, body, open, 'chem-sar-sum-band-head', (show) => {
      lead.style.display = show ? 'none' : '';
      if (show)
        this.foldedBands.delete(key);
      else
        this.foldedBands.add(key);
    });
    head.appendChild(lead);
    return ui.divV([head, body]);
  }

  /** The best pooled swap of one component, or of whichever component the matrix columns enumerated
   *  when the components are not columns and there is only the one pool. */
  private topSwap(data: SummaryData, role?: string): SwapPool | undefined {
    return role === undefined ? data.swaps[0] : (data.swapsByRole.get(role) ?? [])[0];
  }

  /** The swap's conclusion in one phrase, for the group heading when it is shut. */
  private swapHeadline(data: SummaryData, role?: string): string {
    const top = this.topSwap(data, role);
    if (top === undefined)
      return 'none pooled';
    const {worst, widest} = orient(top);
    return worst > 0 ? `≥ ${this.formatDelta(worst)} in ${top.n} pairs` :
      `${this.formatDelta(worst)} to ${this.formatDelta(widest)}`;
  }

  /** One component: how much it moves the endpoint, and the best value of it.
   *
   * The range is written out rather than drawn against the widest component. A bar whose scale is the
   * other rows reads as a fraction of something, and with three components within 0.03 of one another it
   * draws a ranking the numbers do not support. */
  private componentRow(role: RoleSummary, index: number): HTMLElement {
    const top = roleRanks(role) ? role.levels[0] : undefined;
    const parts: HTMLElement[] = [
      ui.divText(role.name, 'chem-sar-sum-comp-name'),
      ui.divText(`spans ${role.spread.toFixed(2)}`, 'chem-sar-sum-comp-spread'),
    ];
    if (top === undefined)
      parts.push(ui.divH([ui.divText('not ranked', 'chem-sar-cp-hint')], 'chem-sar-sum-comp-slot'));
    else {
      parts.push(ui.divH([top.value === '' ? this.hint('nothing at this position',
        `The compounds that leave ${role.name} empty score best. The column is blank for them, so there ` +
        'is no group to draw — in a decomposition that is hydrogen, in a list of named components it ' +
        'means nobody recorded one.') :
        this.depiction(top.value, STRIP_MOL_W, STRIP_MOL_H)], 'chem-sar-sum-comp-slot'));
      parts.push(ui.divText(this.formatEffect(top.coef), 'chem-sar-sum-comp-best'));
      parts.push(ui.divText(`over ${count(top.n)} compounds`, 'chem-sar-cp-hint'));
    }
    const row = ui.divH(parts, 'chem-sar-sum-comp-row');
    ui.tooltip.bind(row, () => `${role.name} moves ${this.host.activityColumnName} by ` +
      `${role.spread.toFixed(2)} across its values. Click for its full ranking.`);
    row.onclick = () => this.showPane(PANE_EFFECTS, `role-${index}`);
    return row;
  }

  /** The components that move the endpoint least, folded away behind their own count. */
  private moreComponents(rows: HTMLElement[]): HTMLElement {
    const body = ui.divV(rows);
    const head = this.foldHead(`${rows.length} more`, body, false, 'chem-sar-sum-comp-more');
    return ui.divV([head, body]);
  }

  /**
   * The best measured swap, on the same row grammar as the components above it. Structures and two
   * numbers rather than a sentence, which has no natural length and clips wherever the band is narrow.
   */
  private swapRow(data: SummaryData, role?: string): HTMLElement {
    // Named like a component row, because that is what it reports on: the pairs are grouped on the other
    // components, so each component has its own and the row has to say which.
    const parts: HTMLElement[] = [
      ui.divText(role ?? 'Best swap', 'chem-sar-sum-comp-name'),
      ui.divText(role === undefined ? '' : 'swapped', 'chem-sar-sum-comp-spread'),
    ];
    const top = this.topSwap(data, role);
    const only = role === undefined ? '' :
      ` Pairs alike in every component but ${role}, pooled across the rest of the table.`;
    let tip: string;
    if (top !== undefined) {
      const {from, to, worst, widest, mean} = orient(top);
      parts.push(this.pairArt(from, to, 'chem-sar-sum-comp-slot chem-sar-sum-comp-pair'));
      parts.push(ui.divText(worst > 0 ? `≥ ${this.formatDelta(worst)}` :
        `mean ${this.formatDelta(mean)}`, 'chem-sar-sum-comp-best'));
      parts.push(ui.divText(`${this.formatDelta(worst)} to ${this.formatDelta(widest)} · ` +
        `${count(top.n)} pairs · ${count(top.roots.size)} ${role === undefined ? 'series' : 'contexts'}`,
      'chem-sar-cp-hint'));
      tip = (worst > 0 ? 'This swap improved potency in every pair it was measured in.' :
        'No swap improved potency in every pair it was measured in; this one has the best floor.') +
        ' Both compounds were made and measured; nothing here is fitted or predicted.' + only +
        ' Click for the full pool.';
    } else {
      const best = data.swapBest;
      parts.push(ui.divH([ui.divText(data.swapCandidates === 0 && role === undefined ?
        'no row carries two measured values' :
        `nothing clears ${SWAP_MIN_PAIRS} pairs in ${SWAP_MIN_SERIES} ` +
        `${role === undefined ? 'series' : 'contexts'}`, 'chem-sar-cp-hint')],
      'chem-sar-sum-comp-slot'));
      if (role === undefined && best?.best != null) {
        parts.push(ui.divText(`largest single move ${this.formatDelta(Math.abs(best.best.delta))} ` +
          `in ${best.best.matrix.label}`, 'chem-sar-cp-hint'));
      }
      tip = 'A swap is pooled only where the same substitution was measured against several different ' +
        'backgrounds, so one pair cannot carry it.' + only + ' Click for what was rejected.';
    }
    const row = ui.divH(parts, 'chem-sar-sum-comp-row');
    ui.tooltip.bind(row, () => tip);
    row.onclick = () => this.showPane(PANE_EFFECTS, role === undefined ? 'swaps' : swapAnchor(role));
    return row;
  }

  /** One conclusion: the name, the number, its bar against the model error, and one clause. The
   *  evidence is a segment away, so nothing here has to carry it. */
  private answerTile(title: string, text: Answer, anchor: string, onOpen?: () => void): HTMLElement {
    const go = ui.iconFA('chevron-right');
    go.classList.add('chem-sar-sum-go');
    const body: HTMLElement[] = [
      ui.divH([ui.divText(title, 'chem-sar-sum-card-title'), go], 'chem-sar-sum-answer-head'),
    ];
    if (text.art !== undefined)
      body.push(text.art);
    body.push(ui.divH([ui.divText(text.answer, 'chem-sar-sum-answer-value')], 'chem-sar-sum-answer-line'));
    // Not `chem-sar-card-desc`: this is a sentence, and that class clips to one line.
    body.push(ui.divText(text.negative, 'chem-sar-cp-hint'));
    const tile = ui.divV(body, 'chem-sar-sum-answer');
    // The answer itself, not only the routing: the value line clips, and a conclusion that fits nowhere
    // on the tile has to be readable somewhere.
    ui.tooltip.bind(tile, () => `${text.answer} — ${text.negative}. The evidence behind this answer is ` +
      'on the Effects segment; click to go to it.');
    tile.onclick = () => {
      onOpen?.();
      this.showPane(PANE_EFFECTS, anchor);
    };
    return tile;
  }

  /** Whether the front-runner may be named alone: its margin over the runner-up has to clear the
   *  screen's own model error, or the two are one answer with two names. */
  private leads(first: number | null, second: number | null, modelError: number | null): boolean {
    if (second === null)
      return true;
    if (first === null)
      return false;
    return first - second > (modelError ?? 0);
  }

  /** The substituent a tile's conclusion is about, drawn. Undefined where there is no one structure —
   *  a tile that names two groups or declines to rank has nothing to draw. */
  private answerArt(smiles: string | null): HTMLElement | undefined {
    return smiles ? this.depiction(smiles, STRIP_MOL_W, STRIP_MOL_H) : undefined;
  }

  /** One role's tile, read off the global fit — the same answer for the axis role and the core role,
   *  since in fragment-columns mode one fit ranks both. Null where that role's card is not ranking. */
  private roleAnswer(data: SummaryData, name: string | null): Answer | null {
    const fit = data.roleFit;
    if (fit === null || this.roleFitRefusal(data) !== null)
      return null;
    const role = fit.roles.find((other) => other.name === name);
    if (role === undefined || !roleRanks(role))
      return null;
    const [top, second] = role.levels;
    // Signed by direction: `levels` is best-first, so on a lower-is-better column the leader's
    // coefficient is the smaller number.
    const dir = this.host.higherIsBetter ? 1 : -1;
    if (!this.leads(dir * top.coef, dir * second.coef, fit.cvRmse)) {
      // A bimodal component has an answer no leaderboard can state: which group, not which value. The
      // top two being inseparable is exactly the case where the split is the finding.
      const split = role.split;
      if (split !== null) {
        return {answer: `two groups: ${count(split.hiCount)} at ${this.formatEffect(split.hiMean)}, ` +
          `${count(split.loCount)} at ${this.formatEffect(split.loMean)}`,
        negative: `against the library mean of ${this.host.formatActivity(fit.mean)}; no single ` +
          `${role.name} leads, but the groups are ${split.gap.toFixed(2)} apart — which group a ` +
          'compound carries matters more than the choice inside it'};
      }
      return {answer: `No single ${role.name} leads`,
        negative: `top two ${this.formatEffect(dir * (top.coef - second.coef))} apart, within the ` +
          `± ${this.host.formatActivity(fit.cvRmse!)} this fit resolves`};
    }
    return {answer: `${this.formatEffect(top.coef)} · over ${count(top.n)} compounds`,
      negative: `against the library mean of ${this.host.formatActivity(fit.mean)}, adjusted for the ` +
        'other components — observational, not a potency', art: this.answerArt(top.value)};
  }

  private rgroupAnswer(data: SummaryData): Answer {
    if (!this.host.activityIsLog) {
      // There is no mark for "not answerable", so this branch stays a sentence.
      return {answer: 'Not answerable on a raw scale',
        negative: 'a fitted effect is a difference in assay units and says nothing about how large ' +
          'the change is — set Scaling to lg or −lg, or declare the column higher-is-better if it is ' +
          'already a pIC50'};
    }
    const loser = data.rgroupLosers[0];
    const negative = loser === undefined ?
      'over each series\' own most-common substituent, not over your whole library' :
      `stop making ${shortSmiles(loser.subst)} — last in ${loser.k} of ${loser.tried} tried`;
    const rows = data.rgroups;
    if (rows.length === 0) {
      return {answer: `None in ${RGROUP_MIN_SERIES} series`,
        negative: 'no substituent comes first at its position in that many series whose additive fit holds'};
    }
    const top = rows[0];
    const partner = data.axisRole === null ? 'series' : 'core';
    if (!this.leads(top.magnitude, rows[1]?.magnitude ?? null, data.modelError)) {
      return {answer: `Two lead equally`,
        negative: `the ${partner} decides — ${top.k} of ${top.tried} against ${rows[1].k} of ` +
          `${rows[1].tried}`,
        art: ui.divH([this.depiction(top.subst, STRIP_MOL_W, STRIP_MOL_H),
          ui.divText('or', 'chem-sar-sum-comp-arrow'),
          this.depiction(rows[1].subst, STRIP_MOL_W, STRIP_MOL_H)], 'chem-sar-sum-comp-pair')};
    }
    return {answer: (top.magnitude === null ? '' : `${this.formatEffect(top.magnitude)} · `) +
      `first on ${top.k} of ${top.tried}`, negative, art: this.answerArt(top.subst)};
  }

  private coreAnswer(data: SummaryData, cores: SeriesStat[]): Answer {
    const host = this.host;
    if (!data.coresAreSeries) {
      // A core recurs only inside its own lineage, so there is no ranking a mark could carry.
      return {answer: 'Not comparable across series',
        // Two lines at tile width; a third is clipped, and a clipped sentence is what the tile is for.
        negative: 'a core recurs only inside its own lineage — nothing to pool across series'};
    }
    if (cores.length === 0)
      return {answer: `None with ${MIN_SUPPORT} compounds`, negative: 'mean of what was made'};
    // Which cut depth was ranked, because only cores cut alike are comparable and the card lets the
    // reader change it.
    const depth = this.coreTiers(data).length > 1 ? `among the L${cores[0].tier} cores; ` : '';
    const negative = `${depth}mean of what was made — a core only ever paired with good caps scores high`;
    const dir = host.higherIsBetter ? 1 : -1;
    const top = cores[0];
    if (!this.leads(dir * top.typical, cores.length > 1 ? dir * cores[1].typical : null, data.modelError)) {
      return {answer: `No core runs clear`,
        negative: `${depth}the top ${cores.length} sit inside ` +
          `${host.formatActivity(data.modelError ?? 0)} model error of each other`};
    }
    return {answer: `${top.matrix.label} · typical ${host.formatActivity(top.typical)} over ` +
      `${count(top.cpd)} cpd`, negative,
    art: this.answerArt(top.matrix.rows[0]?.keySmiles ?? null)};
  }

  private potencyAnswer(data: SummaryData): Answer {
    const negative = 'both compounds measured, nothing fitted — one R-group swapped, the rest identical';
    const top = data.swaps[0];
    if (top !== undefined) {
      const {from, to, worst, widest, mean} = orient(top);
      const art = this.pairArt(from, to, 'chem-sar-sum-comp-pair');
      if (worst > 0)
        return {answer: `≥ ${this.formatDelta(worst)} in ${top.n} pairs`, negative, art};
      // The pool is ranked on its floor, and a negative floor under this heading would state the
      // reverse of the question: the swap lost in at least one pair it was measured in.
      return {answer: `${this.formatDelta(worst)} to ${this.formatDelta(widest)} over ${top.n} pairs`,
        negative: `no swap improves potency in every pair it was measured in — this one is the best ` +
          `floor; mean ${this.formatDelta(mean)}`, art};
    }
    // Nothing pooled: what the gate rejected is the finding, and there is nothing to draw.
    const gate = data.swapCandidates === 0 ? 'No row carries two measured substituents' :
      `Nothing clears ${SWAP_MIN_PAIRS} pairs in ${SWAP_MIN_SERIES} series`;
    const best = data.swapBest;
    if (best === null || best.best === null)
      return {answer: gate, negative};
    return {answer: gate, negative: `the largest single measured move is ` +
      `${this.formatDelta(Math.abs(best.best.delta))} in ${best.best.matrix.label}`};
  }

  // ---- Rows -----------------------------------------------------------------------------------

  /** Canvas now, RDKit later: a synchronous pass over every structure on the tab would stall the tab
   *  switch that asked for them. */
  private depiction(smiles: string | null, w: number, h: number): HTMLElement {
    const canvas = ui.canvas(w, h);
    canvas.classList.add('chem-sar-card-core');
    if (smiles)
      this.pendingPaints.push(() => paintMoleculeOnColor(canvas, smiles, w, h, CORE_BG_ARGB));
    return canvas;
  }

  /** Two structures and the arrow between them: the shape every swap is drawn in. */
  private pairArt(from: string, to: string, cls: string): HTMLElement {
    return ui.divH([this.depiction(from, STRIP_MOL_W, STRIP_MOL_H),
      ui.divText('→', 'chem-sar-sum-comp-arrow'),
      this.depiction(to, STRIP_MOL_W, STRIP_MOL_H)], cls);
  }

  private flushPaints(): void {
    const paints = this.pendingPaints;
    this.pendingPaints = [];
    for (const paint of paints)
      paint();
  }

  private badge(text: string, tip: string, partial = false): HTMLElement {
    const el = ui.divText(text, `chem-sar-chip-badge${partial ? ' chem-sar-chip-partial' : ''}`);
    ui.tooltip.bind(el, tip);
    return el;
  }

  /** Signed, and finer than `formatActivity`: fitted effects and pooled differences live in tenths,
   *  where one decimal renders three differently-ranked R-groups as three rows all reading "+0.1". */
  private formatEffect(value: number): string {
    const abs = Math.abs(value);
    return `${value < 0 ? '−' : '+'}${abs < 10 ? abs.toFixed(2) : abs.toFixed(1)}`;
  }

  /** A measured difference in the reader's own units: a signed log difference where the scale is a
   *  log, and the fold it corresponds to where it is not. A loss is written as its reciprocal fold:
   *  10^value rounds every loss past one log to "0.0×", which destroys the number it is reporting. */
  private formatDelta(value: number): string {
    if (this.host.activityIsLog)
      return this.formatEffect(value);
    const fold = Math.pow(10, Math.abs(value));
    const shown = fold < 10 ? fold.toFixed(1) : fold.toFixed(0);
    return value < 0 ? `1/${shown}×` : `${shown}×`;
  }

  /** Every clickable row on the tab. `.chem-sar-card` carries the hover, the left border and the row
   *  metric the two text lines and the value column are all sized against. */
  private summaryRow(opts: {
    depiction: HTMLElement | null,
    name: string,
    badges: HTMLElement[],
    desc: string,
    /** A mark lane between the text column and the value, at a fixed width so it reads down a list. */
    mark?: HTMLElement,
    value: string,
    caption: string,
    valueTip?: string,
    faint?: boolean,
    cart?: MatrixCellRef,
    tip: string,
    onClick: () => void,
    /** Makes the chevron a toggle of its own rather than a repeat of the row's click. */
    onChevron?: () => void,
  }): HTMLElement {
    const value = ui.divText(opts.value,
      `chem-sar-sum-value${opts.faint ? ' chem-sar-sum-faint' : ''}`);
    const caption = ui.divText(opts.caption, 'chem-sar-sum-cap');
    const valueBox = ui.divV([value, caption], 'chem-sar-sum-valuebox');
    if (opts.valueTip)
      ui.tooltip.bind(valueBox, () => opts.valueTip!);
    const parts: HTMLElement[] = [];
    if (opts.depiction !== null)
      parts.push(opts.depiction);
    parts.push(ui.divV([
      ui.divH([ui.divText(opts.name, 'chem-sar-card-name'), ...opts.badges], 'chem-sar-card-title'),
      ui.divText(opts.desc, 'chem-sar-card-desc'),
    ], 'chem-sar-card-body'));
    if (opts.mark !== undefined)
      parts.push(opts.mark);
    parts.push(valueBox);
    if (opts.cart !== undefined) {
      const ref = opts.cart;
      const cart = ui.iconFA('cart-plus', (e: MouseEvent) => {
        e.stopPropagation();
        this.host.addCellsToMakeList([ref], 'This cell has no structure to add.');
      }, 'Add this compound to the Make list');
      cart.classList.add('chem-sar-struct-icon', 'chem-sar-cart-icon');
      parts.push(cart);
    }
    const toggle = opts.onChevron;
    const go = toggle === undefined ? ui.iconFA('chevron-right') :
      ui.iconFA('chevron-down', (e: MouseEvent) => {
        e.stopPropagation();
        toggle();
      }, 'What is worth opening this series for');
    go.classList.add('chem-sar-sum-go');
    parts.push(go);
    const row = ui.divH(parts, 'chem-sar-card chem-sar-sum-row');
    ui.tooltip.bind(row, () => opts.tip);
    row.onclick = opts.onClick;
    return row;
  }

  private cellLocation(matrix: SarMatrix, ri: number, ci: number): string {
    return `${matrix.label} · ${matrix.rows[ri].label} × ${matrix.columns[ci].position}`;
  }

  private openTip(matrix: SarMatrix, ri: number, ci: number): string {
    return `Open ${this.cellLocation(matrix, ri, ci)} in the SAR Matrix`;
  }

  private card(title: string, subtitle: string, body: HTMLElement[], footer?: HTMLElement): HTMLElement {
    const parts = [ui.divText(title, 'chem-sar-sum-card-title'),
      ui.divText(subtitle, 'chem-sar-cp-hint'), ...body];
    if (footer !== undefined)
      parts.push(footer);
    return ui.divV(parts, 'chem-sar-sum-card');
  }

  /** A card with no qualifying rows still renders: a vanishing card reflows the grid and hides the
   *  fact that the analysis found nothing of that kind. */
  private reason(text: string): HTMLElement[] {
    return [ui.divText(text, 'chem-sar-cp-hint')];
  }

  // ---- Band B: start here ---------------------------------------------------------------------

  /**
   * Whether the reason lane is drawn on this render — the header's key, the per-row lane and the text
   * fallback all read this one answer, so the mark can never lose its legend or the phrases with it.
   *
   * A lane that wins every slot carries no information, so with one series the phrases come back as
   * text: a constant mark is decoration. Below the narrow breakpoint five 14px slots and a readable row
   * cannot both fit, and `applySize` rebuilds the pane when that breakpoint is crossed.
   */
  private lanesFit(data: SummaryData): boolean {
    return data.startHere.length > 1 && !this.root.classList.contains('chem-sar-sum-xnarrow');
  }

  private buildStartHere(data: SummaryData): HTMLElement[] {
    const host = this.host;
    const parts: HTMLElement[] = [];
    if (data.startHere.length === 0)
      parts.push(...this.reason('No series holds a measured compound.'));
    const lanes = this.lanesFit(data);
    for (const row of data.startHere) {
      const stat = row.primary;
      const {matrix} = stat;
      const rowKey = stat.best !== null ? matrix.rows[stat.best.ri].keySmiles : matrix.rows[0]?.keySmiles;
      const neighbours = stat.best === null ? 0 :
        host.observedNeighbours(matrix, stat.best.ri, stat.best.ci);
      const expand = ui.divV([]);
      expand.style.display = 'none';
      let filled = false;
      parts.push(this.summaryRow({
        depiction: this.depiction(rowKey ?? null, CARD_CORE_W, CARD_CORE_H),
        name: matrix.label,
        badges: [
          this.badge(`L${stat.tier}`, stat.tier === 1 ?
            'A leaf series: no finer series sit under this one.' :
            `A coarser series, holding the compounds of the L${stat.tier - 1} series below it, whose ` +
            'cores agree one further cut deeper.'),
          this.trustDot(matrix),
          this.gainBadge(stat),
        ],
        mark: lanes ? this.reasonLane(row.reasons) : undefined,
        // Cells the decomposition cannot express are never predicted and never offered, so they are no
        // part of the denominator — but only a matrix whose assembly marked them can say so.
        desc: `${count(stat.cpd)} cpd · ${count(stat.realCells)} of ` +
          `${count(stat.totalCells - stat.impossibleCells)} ` +
          `${stat.impossibleCells > 0 ? 'makeable cells' : 'cells'} measured · ` +
          `${host.formatActivity(stat.lo)}–${host.formatActivity(stat.hi)}` +
          (stat.realCells ? ` · typical ${host.formatActivity(stat.typical)}` : ''),
        value: stat.best === null ? '—' : host.formatActivity(stat.best.value),
        caption: 'best',
        valueTip: stat.best === null ? 'This series holds no measured compound.' :
          `${host.formatActivity(stat.best.value)} measured. ${neighbours} measured cells share this ` +
          'core or this substituent.',
        tip: `Open ${matrix.label} in the SAR Matrix. "Typical" is this series' fitted mean — the mean ` +
          'of what was made, count-weighted, so a series whose chemists only made their good ideas ' +
          'scores high and one heavily resampled core dominates it. The gap between typical and best ' +
          'is the lucky-singleton signal.' + (stat.impossibleCells > 0 ? '' :
          ' No cell of this series was marked unmakeable, which can mean the decomposition never ' +
          'assessed them rather than that the grid is complete.'),
        onClick: () => stat.best === null ? host.revealMatrix(matrix) :
          host.revealCell(matrix, stat.best.ri, stat.best.ci),
        onChevron: () => {
          if (!filled) {
            filled = true;
            for (const part of this.startExpandParts(stat))
              expand.appendChild(part);
            // Rows added after the paint tick has run, so these depictions need one of their own.
            this.flushPaints();
          }
          expand.style.display = expand.style.display === 'none' ? '' : 'none';
        },
      }));
      if (!lanes)
        parts.push(ui.divText(row.reasons.map((i) => REASON_WORDS[i]).join(' · '), 'chem-sar-cp-hint'));
      parts.push(expand);
    }
    const transfers = ui.divText('', 'chem-sar-cp-hint');
    this.transferLine = transfers;
    this.syncTransferLine();
    transfers.classList.add('chem-sar-sum-click');
    transfers.onclick = () => host.showTab(TAB_TRANSFER);
    parts.push(transfers);
    return parts;
  }

  /**
   * The best unmade compound the model names in this series, against the best one already measured.
   *
   * An additive FLOOR, not a ceiling on the chemistry — a cliff is exactly what this model cannot see,
   * so a gain near zero means "no further additive gain", never "no further potency".
   */
  /**
   * Which gate left this series with nothing to aim at. The gate is one condition — a predicted cell
   * needs MIN_SUPPORT measured compounds on both axes AND a series whose leave-one-out R² holds — but
   * it fails in ways that call for different things: a setting to change, more compounds, or nothing
   * at all because the grid is already full.
   */
  private noGainReason(stat: SeriesStat): {label: string, tip: string} {
    const r2 = stat.matrix.confidence?.r2 ?? null;
    if (!this.host.predictVirtual) {
      return {label: 'predictions off', tip: '"Predict virtual analogs" is off, so this series has no ' +
        'predicted cells to compare against. Turn it on in the viewer\'s properties.'};
    }
    if (stat.virtualCells === 0) {
      return {label: 'nothing unfilled', tip: 'Every cell this series can express is already measured ' +
        'or cannot be made, so there is no unfilled cell to predict.'};
    }
    if (r2 === null) {
      return {label: 'fit unchecked', tip: `This series carries ${count(stat.virtualCells)} predicted ` +
        'cells, but its additive fit was never cross-validated — too few measured cells to leave one ' +
        'out — so none of them can be trusted. Unchecked, not wrong.'};
    }
    if (r2 < TRUST_R2) {
      return {label: 'fit fails', tip: `This series carries ${count(stat.virtualCells)} predicted ` +
        `cells, but its own R² is ${r2.toFixed(2)}, below ${TRUST_R2}: substituent effects ` +
        'do not add up here, so its predictions are arithmetic rather than estimates.'};
    }
    return {label: 'thin support', tip: `All ${count(stat.virtualCells)} predicted cells here rest on ` +
      `fewer than ${MIN_SUPPORT} measured compounds on one of their two axes. The fit holds; the ` +
      'evidence behind each individual cell does not.'};
  }

  private gainBadge(stat: SeriesStat): HTMLElement {
    const host = this.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const reach = stat.bestVirtual !== null && stat.best !== null ?
      dir * (stat.bestVirtual.value - stat.best.value) : null;
    let tip = 'The best unfilled cell the model names, against the best compound already measured. ' +
      'Negative means the best compound already beats anything the additive model can build from this ' +
      'series\' parts — what a cliff looks like. Near zero means no further ADDITIVE gain, not that ' +
      'the series is exhausted.';
    let label = `${this.formatEffect(reach ?? 0)} vs best made`;
    if (reach === null)
      ({label, tip} = this.noGainReason(stat));
    const badge = this.badge(label, tip, reach === null);
    if (stat.bestVirtual !== null) {
      badge.classList.add('chem-sar-sum-click');
      badge.onclick = (e: MouseEvent) => {
        e.stopPropagation();
        host.revealCell(stat.matrix, stat.bestVirtual!.ri, stat.bestVirtual!.ci);
      };
    }
    // The one thing the number itself cannot say: whether +0.4 is four of this series' own prediction
    // errors or a quarter of one.
    const rmse = stat.matrix.confidence?.rmse ?? null;
    if (reach === null || rmse === null || !host.activityIsLog)
      return badge;
    return ui.divH([badge, this.tickRun(reach, rmse)], 'chem-sar-sum-gain');
  }

  private startExpandParts(stat: SeriesStat): HTMLElement[] {
    const host = this.host;
    const {matrix} = stat;
    const parts: HTMLElement[] = [];

    if (stat.colRange !== null && stat.rowRange !== null) {
      const story = stat.colRange > stat.rowRange ? 'this is a column story' : 'this is a row story';
      const spread = ui.divText(`Spread: substituent choice spans ${stat.colRange.toFixed(2)} · core ` +
        `choice spans ${stat.rowRange.toFixed(2)} — ${story}`, 'chem-sar-cp-hint');
      ui.tooltip.bind(spread, () => 'Ranges of this series\' own fitted effects, over substituents ' +
        `measured on at least ${MIN_SUPPORT} cores and cores measured at at least ${MIN_SUPPORT} ` +
        'substituents — one floor on both sides, or the looser side of the comparison would come out ' +
        'wider on noise alone. Both are centred on the same fit, so they compare to each other — and ' +
        'to no other series\' numbers.');
      parts.push(spread);
    }

    if (stat.bestRow !== null) {
      const row = matrix.rows[stat.bestRow.ri];
      const folded = Object.values(row.foldedValues);
      const oneStructure = matrix.rows.length > 1 &&
        new Set(matrix.rows.map((r) => r.keySmiles)).size === 1;
      const heading = folded.length === 0 ? 'Best core' :
        oneStructure ? 'Best R-group combination' : 'Best core + fixed R-groups';
      // With R-group columns every row can be labelled "Row N" and draw the same structure, so what
      // distinguishes the winning row is its folded substituents, not its picture.
      const art = folded.length > 0 && oneStructure ?
        ui.divH(folded.slice(0, 3).map((subst) =>
          this.depiction(subst, BENEFIT_MOL_W, BENEFIT_MOL_H)), 'chem-sar-sum-pair') :
        this.depiction(row.keySmiles, CARD_CORE_W, CARD_CORE_H);
      const text = ui.divV([
        ui.divText(`${heading}: ${row.label} · ${this.formatEffect(stat.bestRow.effect)} ` +
          `(n = ${stat.bestRow.n})`, 'chem-sar-sum-prose'),
        ui.divText('Best on average across the substituents it was actually measured at — not the ' +
          'single best compound, and not comparable to another series.', 'chem-sar-cp-hint'),
      ]);
      const block = ui.divH([art, text], 'chem-sar-sum-detail');
      block.classList.add('chem-sar-sum-click');
      const ci = this.bestMeasuredCol(matrix, stat.bestRow.ri);
      if (ci >= 0)
        block.onclick = () => host.revealCell(matrix, stat.bestRow!.ri, ci);
      parts.push(block);
    }

    // The one place on this screen where non-additivity is a finding rather than a trust problem, and
    // where the instruction is the opposite of the leaderboards': hold the pair, vary elsewhere.
    const conf = matrix.confidence;
    const dir = host.higherIsBetter ? 1 : -1;
    // `?? null`, not a strict read: a matrix carried in by a project or layout was serialized by
    // whichever Chem wrote it, and one without these keys yields undefined, which `!== null` passes.
    const outlier = (conf ? (dir > 0 ? conf.hi : conf.lo) : null) ?? null;
    if (conf && outlier !== null && Math.abs(outlier.residual) > conf.rmse &&
      outlier.ri < matrix.rows.length && outlier.ci < matrix.columns.length) {
      const beat = ui.divText(`${matrix.rows[outlier.ri].label} × ` +
        `${shortSmiles(matrix.columns[outlier.ci].substSmiles)} beats its additive expectation by ` +
        `${this.formatEffect(dir * outlier.residual)}, against this series' own ` +
        `±${host.formatActivity(conf.rmse)} — hold that pair and vary elsewhere.`, 'chem-sar-cp-hint');
      beat.classList.add('chem-sar-sum-click');
      ui.tooltip.bind(beat, 'Out-of-sample: the cell was held out and predicted from the rest, so ' +
        'the model could not pull itself toward it. It beats the additive sum in this series\' own ' +
        'units; it does not say the mechanism is understood.');
      beat.onclick = () => host.revealCell(matrix, outlier.ri, outlier.ci);
      parts.push(beat);
    }

    const position = matrix.positions[0] ?? '';
    const reference = matrix.refValues[position];
    const refLine = ui.divText(reference ? `Varies ${position} · reference substituent:` :
      `Varies ${position} · no reference substituent recorded`, 'chem-sar-cp-hint');
    ui.tooltip.bind(refLine, 'The position this series explores, and the most frequently observed ' +
      'substituent at it — the comparator its fitted effects are quoted against. Position labels come ' +
      'from each series\' own decomposition, so R1 here and R1 elsewhere are unrelated.');
    parts.push(reference ?
      ui.divH([refLine, this.depiction(reference, BENEFIT_MOL_W, BENEFIT_MOL_H)], 'chem-sar-sum-detail') :
      refLine);
    return parts;
  }

  /** The most potent measured cell of one column, so a leaderboard row lands on a compound. */
  private bestMeasuredRow(matrix: SarMatrix, ci: number): number {
    const dir = this.host.higherIsBetter ? 1 : -1;
    let ri = -1;
    for (let r = 0; r < matrix.rows.length; r++) {
      const cell = matrix.cells[r][ci];
      if (cell.kind === 'real' && cell.value !== null &&
        (ri < 0 || dir * cell.value > dir * matrix.cells[ri][ci].value!))
        ri = r;
    }
    return ri;
  }

  /** The most potent measured cell of one row, so an expand lands on a compound. */
  private bestMeasuredCol(matrix: SarMatrix, ri: number): number {
    const dir = this.host.higherIsBetter ? 1 : -1;
    let ci = -1;
    for (let c = 0; c < matrix.columns.length; c++) {
      const cell = matrix.cells[ri][c];
      if (cell.kind === 'real' && cell.value !== null &&
        (ci < 0 || dir * cell.value > dir * matrix.cells[ri][ci].value!))
        ci = c;
    }
    return ci;
  }

  // ---- Band C: the three cards ----------------------------------------------------------------

  private buildSwapCard(data: SummaryData, role?: string): HTMLElement {
    const host = this.host;
    const pools = role === undefined ? data.swaps : data.swapsByRole.get(role) ?? [];
    const title = role === undefined ? 'What buys potency inside these series' :
      `Swapping ${role} — what it bought`;
    const subtitle = 'Differences between two compounds that were both made and both measured, alike in ' +
      'everything but one R-group. Nothing here is fitted or predicted: if the pair was never made, it ' +
      'is not on this card.';
    const unit = role === undefined ? 'series' : 'contexts';
    const rows = pools.map((pool) => {
      const {forward, from, to, worst, widest: bestOf, mean, up} = orient(pool);
      const allUp = up === pool.n;
      // Null for a pool grouped off the component columns: the pair was found by matching component
      // values, so there is no one cell of one matrix it can be opened at.
      const best = pool.best;
      const ci = best === null ? -1 : forward ? best.ciTo : best.ciFrom;
      const position = best === null ? '' : best.matrix.columns[ci].position;
      const show = (value: number): string => this.formatDelta(value);
      // Badges, not descriptor tail: the descriptor clips to one line, and these three are exactly the
      // qualifiers that stop the row being over-read — appended last they are the first characters lost.
      const badges: HTMLElement[] = [];
      if (allUp) {
        badges.push(this.badge('all up', 'Every measured pair moved the same way. At 3 pairs that is ' +
          'one-in-four under a coin-flip null.'));
      }
      if (pool.sampled) {
        badges.push(this.badge('sampled row', 'One contributing row carried more than ' +
          `${SWAP_ROW_CAP} measured cells; its most and least potent halves were kept and mid-range ` +
          'pairs dropped.', true));
      }
      // Only on a log scale: the pooled delta is a log ratio there and the model error is in the same
      // units, while on a raw scale the two are different quantities.
      if (host.activityIsLog && data.modelError !== null && Math.abs(worst) < data.modelError) {
        badges.push(this.badge('under resolution', 'The worst case is smaller than this analysis\' own ' +
          'typical prediction error, so the whole row may be noise.', true));
      }
      const pair = ui.divH([
        this.depiction(from, BENEFIT_MOL_W, BENEFIT_MOL_H),
        ui.divText('→', 'chem-sar-sum-arrow'),
        this.depiction(to, BENEFIT_MOL_W, BENEFIT_MOL_H),
      ], 'chem-sar-sum-pair');
      return this.summaryRow({
        depiction: pair,
        name: best === null ? `${role} swap` : `${best.matrix.label} · ${position}`,
        badges,
        desc: `mean ${show(mean)} (${show(worst)} → ${show(bestOf)}) · ${pool.n} pairs · ` +
          `${pool.roots.size} ${unit}`,
        value: show(worst),
        caption: 'worst case',
        valueTip: `This swap was worth at least ${show(worst)} in every one of the ${pool.n} measured ` +
          `pairs, across ${pool.roots.size} ${unit}. Each pair is two measured compounds alike in ` +
          'everything else, so what they do not share cancels exactly and no fit is involved. At 3 ' +
          'pairs "all up" is one-in-four under a coin-flip null.' +
          (pool.sampled ? ` One group carried more than ${SWAP_ROW_CAP} measured compounds; its most ` +
            'and least potent halves were kept and mid-range pairs dropped.' : ''),
        // A pool grouped off the component columns has no one cell to open, so the row selects the
        // compounds that carry the better side instead of going nowhere.
        tip: best === null ? `Select every compound whose ${role} is the one on the right` :
          this.openTip(best.matrix, best.ri, ci),
        onClick: best === null ? () => host.selectRoleValue(role!, to) :
          () => host.revealCell(best.matrix, best.ri, ci, position),
      });
    });
    if (rows.length > 0)
      return this.card(title, subtitle, rows);
    return this.card(title, subtitle, this.reason(data.swapCandidates === 0 && role === undefined ?
      'No row of any series carries two measured substituents, so there is no swap to measure.' :
      `No ${role === undefined ? '' : `${role} `}swap clears ${SWAP_MIN_PAIRS} measured pairs in ` +
      `${SWAP_MIN_SERIES} ${unit}.`));
  }

  private buildRGroupCard(data: SummaryData): HTMLElement {
    const host = this.host;
    const title = data.axisRole === null ? 'R-groups — within-series ranking' :
      `${data.axisRole} — within-series ranking`;
    const subtitle = 'Ranked by how often the fitted model places this group first at its position. The ' +
      'number is its margin over that series\' own most-common substituent — a within-series ' +
      'comparison, so two series\' numbers are only loosely comparable.';
    // An additive effect in raw assay units cannot be re-expressed as a fold, which is the one thing
    // that would make it mean something; the swap card refuses the same arithmetic.
    if (!host.activityIsLog) {
      return this.card(title, subtitle, this.reason('The activity is on a raw scale, where a fitted ' +
        'effect is a difference in assay units and says nothing about how large the change is. Set ' +
        'Scaling to lg or −lg, or declare the column higher-is-better if it is already a pIC50.'));
    }
    const slots = data.series.filter((stat) => stat.stripCols !== null);
    // One scale for every bar on this card, so the bars compare down the card — and, deliberately, to
    // nothing on any other card.
    const shown = [...data.rgroups.slice(0, SUM_ROWS), ...data.rgroupsThin, ...data.rgroupLosers];
    const scale = shown.reduce((m, r) => Math.max(m, Math.abs(r.magnitude ?? 0)), 0);
    // With one row the bar is full-length and says nothing.
    const bars = shown.filter((r) => r.magnitude !== null).length > 1 ? scale : 0;
    // One strip scale for the whole card, like the bars: normalised per row, every row's own best
    // difference draws at full height, so a +1.8 row and a +0.05 row terminate alike.
    let stripScale = 0;
    for (const stat of slots) {
      for (const r of shown) {
        const delta = stat.stripCols!.get(r.subst)?.refDelta;
        if (delta != null)
          stripScale = Math.max(stripScale, Math.abs(delta));
      }
    }
    const rows = data.rgroups.slice(0, SUM_ROWS).map((row, i) =>
      this.rgroupBlock(data, slots, row, bars, stripScale, {primary: i === 0}));
    const thin = data.rgroupsThin.map((row) =>
      this.rgroupBlock(data, slots, row, bars, stripScale, {thin: true}));
    if (rows.length === 0 && thin.length === 0) {
      // Convergence is named because it is the one cause the reader cannot see anywhere else: a fit
      // that stopped short is silently excluded from every pool, and R² says nothing about it.
      return this.card(title, subtitle, this.reason(`No substituent comes first at its position in ` +
        `${SWAP_MIN_SERIES} series whose additive fit holds.` + (data.nonConverged === 0 ? '' :
        ` The additive fit of ${count(data.nonConverged)} of the ${count(this.host.matrices.length)} ` +
        'series did not reach its tolerance, and those are left out of this comparison.')));
    }
    const body: HTMLElement[] = [];
    // Two numbers for one substituent are a tab apart, and nothing else says which question each
    // answers. Only where that second number exists: where the fit declined to rank this column there
    // is no other card to reconcile with.
    if (this.roleAnswer(data, data.axisRole) !== null) {
      body.push(ui.divText('A value can place first in most series and still sit at the library mean: ' +
        'this card counts wins within a series against that series\' own reference, while the ' +
        `${data.axisRole} tab reports one offset pooled over the table.`,
      'chem-sar-sum-prose'));
    }
    if (slots.length > 0 && data.nonConverged > 0) {
      body.push(ui.divText(`${count(data.nonConverged)} further series carry no square and count in ` +
        'no row: their additive fit did not reach its tolerance, so their effects cannot be ordered ' +
        'against another series\'.', 'chem-sar-cp-hint'));
    }
    body.push(...rows);
    if (thin.length > 0) {
      body.push(ui.divText(`Won in only two series — a median over two lineages is the mean of two ` +
        'numbers:', 'chem-sar-cp-hint'));
      body.push(...thin);
    }
    if (data.rgroupLosers.length > 0) {
      const head = data.rgroupLosers[0];
      body.push(ui.divText(`Stop making: ${shortSmiles(head.subst)} came last at its position in ` +
        `${head.k} of the ${head.tried} series it was tried in` +
        (head.magnitude === null ? '.' : `, costing ${this.formatEffect(head.magnitude)} against each ` +
        'series\' own reference.'), 'chem-sar-cp-hint'));
      body.push(...data.rgroupLosers.map((row) =>
        this.rgroupBlock(data, slots, row, bars, stripScale, {loser: true})));
    }
    return this.card(title, subtitle, body, this.sizeVerdict(data));
  }

  /** What one strip slot is. A strip is drawn for any fragment-columns run, but only where a matrix IS
   *  one core is a slot a core — with a series column, or with nothing on the row axis, one slot holds
   *  several cores. */
  private slotNoun(data: SummaryData): string {
    return data.coresAreSeries ? 'core' : 'series';
  }

  /** One leaderboard entry: the row, its per-slot strip where one can be drawn, and an expand that
   *  says what the row does and does not claim. */
  private rgroupBlock(data: SummaryData, slots: SeriesStat[], row: RGroupRow, scale: number,
    stripScale: number, opts: {thin?: boolean, loser?: boolean, primary?: boolean}): HTMLElement {
    const parts: HTMLElement[] = [this.rgroupRow(data, slots, row, scale, opts)];
    if (slots.length > 0)
      parts.push(this.buildStrip(data, slots, row, stripScale));
    const body = ui.divV([]);
    body.style.display = 'none';
    const toggle = ui.divText('What this row claims →', 'chem-sar-cp-hint');
    toggle.classList.add('chem-sar-sum-click');
    let filled = false;
    const open = (): void => {
      if (!filled) {
        filled = true;
        for (const part of this.rgroupExpandParts(data, row))
          body.appendChild(part);
        this.flushPaints();
      }
      body.style.display = '';
    };
    toggle.onclick = () => {
      if (body.style.display === 'none')
        open();
      else
        body.style.display = 'none';
    };
    // The answer tile above asked for its evidence, and the pane it asked for is only built now.
    if (opts.primary && this.expandTopRGroup) {
      this.expandTopRGroup = false;
      open();
    }
    parts.push(toggle, body);
    return ui.divV(parts);
  }

  private rgroupExpandParts(data: SummaryData, row: RGroupRow): HTMLElement[] {
    const slot = this.slotNoun(data);
    // A reading instruction, not a restatement of the marks above it.
    const read = ui.divText(`How to read this: first on most of the ${slot}s it was tried on and never ` +
      `worse than the reference means the group travels — carry it forward. First on one ${slot} while ` +
      `costing potency on the rest means hold it to that ${slot} and vary elsewhere.`,
    'chem-sar-sum-prose');
    const art = ui.divH([this.depiction(row.subst, CARD_CORE_W * 2, CARD_CORE_H * 2), read],
      'chem-sar-sum-detail');
    return [art];
  }

  /**
   * MK-A. One square per slot, in a fixed order that is never re-sorted per row — the row is only
   * readable against its neighbours if the same slot means the same series everywhere.
   *
   * Sign is carried by the direction the fill grows from the mid-rule, and only redundantly by hue:
   * a tint alone renders +0.9 and −0.9 as the same grey. The mid-rule is the legend.
   *
   * `widest` is the card's largest difference, not this row's, so fill heights compare down the card
   * the way the slot order lets the squares compare across it.
   */
  private buildStrip(data: SummaryData, slots: SeriesStat[], row: RGroupRow,
    widest: number): HTMLElement {
    const host = this.host;
    const slot = this.slotNoun(data);
    const squares = slots.map((stat) => {
      const col = stat.stripCols!.get(row.subst);
      if (col === undefined) {
        const none = ui.divText('–', 'chem-sar-sum-sq chem-sar-sum-sq-none');
        ui.tooltip.bind(none, () => `${stat.matrix.label} — never tried on this ${slot}`);
        return none;
      }
      const square = ui.div([], 'chem-sar-sum-sq');
      const delta = col.refDelta;
      if (delta === null) {
        square.classList.add('chem-sar-sum-sq-hatch');
        ui.tooltip.bind(square, () => `${stat.matrix.label} — tried here over ${col.n} measured cells, ` +
          'but this series records no reference substituent to compare against');
      } else if (delta === 0 && stat.matrix.refValues[stat.matrix.columns[col.ci].position] === row.subst) {
        square.classList.add('chem-sar-sum-sq-ref');
        ui.tooltip.bind(square, () => `${stat.matrix.label} — this group IS this series' reference ` +
          `substituent, over ${col.n} measured cells, so there is nothing to compare it against here`);
      } else {
        square.classList.add('chem-sar-sum-sq-val');
        // No difference draws no fill: the floor exists so a small value still reads as a direction,
        // and applied to an exact zero it would draw a win that is not there.
        if (delta !== 0) {
          const fill = ui.div([], `chem-sar-sum-fill chem-sar-sum-${delta > 0 ? 'up' : 'down'}`);
          fill.style.height = `${Math.max(STRIP_MIN_FILL,
            Math.round((widest > 0 ? Math.abs(delta) / widest : 0) * 9))}px`;
          square.appendChild(fill);
        }
        ui.tooltip.bind(square, () => `${stat.matrix.label} · ${this.formatEffect(delta)} against this ` +
          `series' own reference, over ${col.n} measured cells` + (delta === 0 ?
          ' — no difference at all, so nothing grows from the mid-rule.' :
          ` — ${delta > 0 ? 'better than' : 'worse than'} it, so the fill grows ` +
          `${delta > 0 ? 'up' : 'down'} from the mid-rule. Heights compare down this card, against its ` +
          `widest difference (${this.formatEffect(widest)}).`));
      }
      square.classList.add('chem-sar-sum-click');
      square.onclick = () => {
        const ri = this.bestMeasuredRow(stat.matrix, col.ci);
        if (ri >= 0)
          host.revealCell(stat.matrix, ri, col.ci, stat.matrix.columns[col.ci].position);
      };
      return square;
    });
    return ui.divH(squares, 'chem-sar-sum-strip');
  }

  /** One verdict for the whole card, not a number per row: five heavy-atom counts do not add
   *  themselves up in a reader's head, and it is the sum that ends a series at MW 620. */
  private sizeVerdict(data: SummaryData): HTMLElement | undefined {
    const rows = data.rgroups;
    if (rows.length === 0)
      return undefined;
    const line = ui.divText('', 'chem-sar-cp-hint');
    this.pendingPaints.push(() => {
      const rdkit = getRdKitModule();
      let groups = 0;
      let references = 0;
      let n = 0;
      for (const row of rows) {
        const {matrix, ci} = row.win;
        const reference = matrix.refValues[matrix.columns[ci].position];
        const atoms = cachedAtomCount(row.subst, rdkit);
        const refAtoms = reference ? cachedAtomCount(reference, rdkit) : 0;
        if (atoms > 0 && refAtoms > 0) {
          groups += atoms;
          references += refAtoms;
          n++;
        }
      }
      line.innerText = n === 0 ? '' : `Your top ${n} groups average ${Math.round(groups / n)} heavy ` +
        `atoms; the series references they beat average ${Math.round(references / n)}.`;
    });
    ui.tooltip.bind(line, 'Heavy atoms of the fragment, the attachment point included — not ' +
      'molecular weight, and not lipophilicity, which is what actually kills the compound. A smoke ' +
      'alarm, not a property model.');
    return line;
  }

  private rgroupRow(data: SummaryData, slots: SeriesStat[], row: RGroupRow, scale: number,
    opts: {thin?: boolean, loser?: boolean}): HTMLElement {
    const host = this.host;
    const {win} = row;
    const position = win.matrix.columns[win.ci].position;
    const reference = win.matrix.refValues[position] ?? '';
    const atoms = this.badge('', 'Heavy atoms of the fragment, the attachment point included. The ' +
      'model scores potency alone and carries no efficiency term, so an unlabelled leaderboard is a ' +
      'heavy-atom leaderboard wearing a potency label.');
    this.pendingPaints.push(() => {
      const n = cachedAtomCount(row.subst, getRdKitModule());
      atoms.innerText = n ? `${n} atoms` : '—';
    });
    // Attempts, not wins: "best in 5 of 6" is heard as "tried in 6", and a denominator of wins makes
    // every row on the card read as near-perfect.
    const placed = opts.loser ? 'last' : 'best';
    const badges = [this.winTally(row, placed), atoms];
    if (opts.thin)
      badges.push(this.badge('2 series', 'Two lineages only — the number below is the mean of two.', true));
    // Only a single contributing series may have its reference named: each series is measured against
    // its own most-common substituent, and those genuinely differ, so one SMILES cannot stand for the
    // comparator of a median pooled over several.
    const comparator = row.magnitudeSeries > 1 ? 'each series\' own reference' : shortSmiles(reference);
    return this.summaryRow({
      depiction: this.depiction(row.subst, CARD_CORE_W, CARD_CORE_H),
      name: `${win.matrix.label} · ${position}`,
      badges,
      // The strip below carries the per-slot picture where one can be drawn, and a range line beside it
      // would say the same thing twice and worse.
      desc: row.magnitude === null ? 'no reference substituent to compare against' :
        slots.length > 0 ? `vs ${comparator} over ${row.magnitudeSeries} series` :
          `vs ${comparator} · ${this.formatEffect(row.lo)} → ${this.formatEffect(row.hi)} ` +
          `over ${row.magnitudeSeries} series`,
      mark: row.magnitude === null || scale <= 0 ? undefined :
        this.effectBar(row.magnitude, scale, data.modelError),
      value: row.magnitude === null ? '—' : this.formatEffect(row.magnitude),
      caption: 'vs reference',
      valueTip: row.magnitude === null ?
        'No series contributing this placing also carries its reference substituent as a column, so ' +
        'there is no within-series comparison to quote. The row still ranks on how often it placed.' :
        'Median difference against each series\' own most-common substituent' +
        (row.magnitudeSeries > 1 ? ', which differs from series to series' : ` (${reference})`) +
        '. A within-series difference, so it is free of the fact that two series centre their effects ' +
        `over different substituent menus — but only ${row.magnitudeSeries} of the ${row.k} ` +
        'contributing series carry one, and two series\' numbers are still only loosely comparable.',
      tip: this.openTip(win.matrix, win.ri, win.ci),
      onClick: () => host.revealCell(win.matrix, win.ri, win.ci, position),
    });
  }

  // ---- Band B′: one fit over every component column ---------------------------------------------

  /**
   * One card per role column, ordered by that role's noise-corrected spread, or one card naming why
   * nothing is ranked.
   *
   * Absent rather than empty outside fragment-columns mode: a substituent label discovered by
   * fragmentation is local to its own series, so there is no scale on which one pooled offset could be
   * read and the question does not apply at all.
   */
  private buildRoleCards(data: SummaryData): {label: string, note: string | null, el: HTMLElement}[] {
    if (data.axisRole === null)
      return [];
    const refusal = this.roleFitRefusal(data);
    if (refusal !== null) {
      return [{label: 'Components', note: null,
        el: this.anchored(this.card('Component contributions', refusal, this.roleDropped(data)),
          'role-0')}];
    }
    const roles = data.roleFit!.roles;
    // One scale across every role card, against the tab's usual rule: each of these offsets comes from
    // one fit against one library average, so they are on one scale and a shared bar is the honest
    // rendering. Per-matrix effects elsewhere are centred inside their own matrix and are not. Over
    // the rows that will be drawn only — a coefficient on no card would shorten every bar on every
    // other one for a reason the reader cannot see.
    let scale = 0;
    for (const role of roles.filter(roleRanks)) {
      for (const {level} of this.roleRows(role))
        scale = Math.max(scale, Math.abs(level.coef));
    }
    return roles.map((role, i) => ({
      label: role.name,
      note: roleRanks(role) ? role.spread.toFixed(2) : null,
      el: this.anchored(this.buildRoleCard(data, role, scale), `role-${i}`),
    }));
  }

  /** The one finding that is about every component at once, so it sits above the tabs rather than on
   *  any one of them. The order of the roles is stable; the size of the gap between two of them is not,
   *  so it is named as a rank and never quoted as a ratio. Omitted where a role is refusing to rank its
   *  own values, since that refusal is what invalidates the spread ordering. */
  private roleOrdering(data: SummaryData): HTMLElement | null {
    const roles = data.roleFit?.roles;
    if (roles === undefined || this.roleFitRefusal(data) !== null || !roles.every(roleRanks))
      return null;
    // Consecutive roles closer than the tie threshold are joined by an approximation sign rather than
    // an inequality: the order is stable, a gap that small is not, and one chain says so however many
    // component columns there are.
    const chain = roles.map((role, i) => i === 0 ? role.name :
      `${roles[i - 1].spread - role.spread < ROLE_SPREAD_TIE ? '≈' : '>'} ${role.name}`).join(' ');
    const el = this.hint(`Changing ${roles[0].name} moves ${this.host.activityColumnName} most: ${chain}.`,
      'Each component\'s range between its best and worst value, corrected for the fact that a column ' +
      'with more values gets a wider range by chance alone. That correction is what makes columns with ' +
      'different numbers of values comparable. "≈" marks a gap too small to order.');
    el.classList.add('chem-sar-sum-sub-lead');
    return el;
  }

  /** Why the whole fit cannot be read, or null when it can. */
  private roleFitRefusal(data: SummaryData): string | null {
    const fit = data.roleFit;
    if (fit === null) {
      return `No component column carries two values with ${MIN_SUPPORT} or more compounds each, so ` +
        'there is nothing to compare.';
    }
    if (!fit.converged) {
      return `The fit did not settle within ${ROLE_FIT_MAX_SWEEPS} passes. Each component's offsets are ` +
        'still shrunk toward zero by a different amount and cannot be compared, so nothing is ranked.';
    }
    if (fit.cvR2 === null) {
      return 'Too few compounds to cross-validate a fit over ' +
        `${count(fit.roles.reduce((sum, role) => sum + role.fitted, 0))} component values. Unchecked, ` +
        'not wrong — nothing is ranked.';
    }
    if (fit.cvR2 < TRUST_R2) {
      return `Cross-validated R² ${fit.cvR2.toFixed(2)} — below ${TRUST_R2}, the additive reading does ` +
        'not hold on this table: knowing every component does not predict a compound the fit had not ' +
        'seen. The offsets are not ranked.';
    }
    return null;
  }

  /** Stated wherever the fit is reported: the connected-block prune drops measured compounds by design,
   *  so the compound count the card prints is not the table. */
  private roleDropped(data: SummaryData): HTMLElement[] {
    const dropped = data.roleFit?.dropped ?? 0;
    if (dropped === 0)
      return [];
    return [ui.divText(`${count(dropped)} compounds share no component value with the rest of the table; ` +
      'only the largest connected block is fitted, because offsets from separate blocks rest on ' +
      'unrelated baselines.', 'chem-sar-sum-prose')];
  }

  /** The rows one role card shows: its best three, then its worst two. Two losers rather than one —
   *  where a role falls into two groups the single lowest value is routinely its least-supported
   *  member, and one row would then stand for a whole cluster by its smallest. */
  private roleRows(role: RoleSummary): {level: RoleLevel, faint: boolean}[] {
    const top = role.levels.slice(0, SUM_ROWS);
    return [...top.map((level) => ({level, faint: false})),
      ...role.levels.slice(Math.max(top.length, role.levels.length - LOSER_ROWS))
        .map((level) => ({level, faint: true}))];
  }

  /**
   * Which component the matrix columns actually enumerate, and the one control that changes it.
   *
   * A card ranking R4 sits beside a grid whose columns are R1, which reads as the two disagreeing. They
   * are answering different questions: the card is the fit over the whole table, the grid enumerates one
   * component at a time.
   */
  private axisSwitchNote(name: string, data: SummaryData): HTMLElement {
    const host = this.host;
    const note = this.hint('Ranked from the fit over the whole table. The SAR Matrix columns enumerate ' +
      `${data.axisRole}, so these values are not the ones across the top there.`,
    `Putting ${name} across the columns means rebuilding: every matrix is reassembled and every fit ` +
      'recomputed, which takes as long as the first run did.');
    const pill = ui.divText(`Put ${name} across the columns — rebuilds`,
      'chem-sar-chip-badge chem-sar-sum-role');
    // What will and will not move, on the control: this ranking is the same fit either way, so a reader
    // who clicks expecting these rows to change waits out a rebuild for a screen that looks identical.
    ui.tooltip.bind(pill, () => `Reassembles every matrix with ${name} across the columns, and pools ` +
      `its measured pairs instead of ${data.axisRole}'s. The offsets on this card do not change — they ` +
      'come from one fit over the whole table, whichever component the columns enumerate.');
    pill.onclick = () => {
      if (!host.computing)
        host.setColumnAxis(name);
    };
    return ui.divH([note, pill], 'chem-sar-sum-chips');
  }

  private buildRoleCard(data: SummaryData, role: RoleSummary, scale: number): HTMLElement {
    const host = this.host;
    const fit = data.roleFit!;
    const title = `${role.name} — offsets from the additive fit`;
    const subtitle = `Additive fit over ${fit.roles.length} component columns. Each offset is against ` +
      `the library mean of ${host.formatActivity(fit.mean)}, adjusted for the other components. ` +
      `R² ${fit.cvR2!.toFixed(2)} ± ${host.formatActivity(fit.cvRmse!)} predicting unseen compounds, ` +
      `n = ${count(fit.compounds)}.`;
    const footer = this.roleFooter();
    // On every card, because each is now the only one its reader can see: a subtitle counting fewer
    // compounds than the tab does, on a card that never explains the difference, is the one omission
    // this note exists to prevent.
    const dropped = this.roleDropped(data);
    if (role.levels.length < 2) {
      return this.card(title, subtitle, [...this.reason(`${role.levels.length === 1 ? 'Only one' : 'No'} ` +
        `${role.name} value carries ${MIN_SUPPORT} or more compounds, so there is nothing to compare.`),
      ...dropped], footer);
    }
    // Prediction is nearly unique even where the split of credit between the roles is not, so this is
    // the only check that this role's share of it means anything.
    if (role.repeat !== null && role.repeat < TRUST_R2) {
      return this.card(title, subtitle, [...this.reason(`${role.name} offsets do not repeat: fitted on ` +
        `each half of the table separately they correlate at r ${role.repeat.toFixed(2)}, below ` +
        `${TRUST_R2}. The other components in this table predict which ${role.name} a compound carries ` +
        'closely enough that the fit cannot separate their contributions on a subset, so these values ' +
        'are not ranked.'), ...dropped], footer);
    }

    const chips: HTMLElement[] = [
      this.badge(`spread ${role.spread.toFixed(2)}`, 'Count-weighted sd of this component\'s offsets, ' +
        'with estimation noise subtracted so columns with different numbers of values compare. ' +
        'Not a range.'),
      this.badge(`${role.levels.length} of ${role.fitted} values ranked`,
        `${role.thin} values seen fewer than ${MIN_SUPPORT} times stay in the fit but carry no readable ` +
        'offset of their own.'),
    ];
    if (role.repeat !== null) {
      chips.push(this.badge(`repeats at r ${role.repeat.toFixed(2)}`,
        'Fitted on each half of the table separately, these offsets correlate this closely. Below ' +
        TRUST_R2 + ' they would not be ranked.'));
    }
    const body: HTMLElement[] = [ui.divH(chips, 'chem-sar-sum-chips')];
    if (role.name !== data.axisRole && host.roleColumns.includes(role.name))
      body.push(this.axisSwitchNote(role.name, data));
    const split = role.split;
    if (split !== null) {
      body.push(ui.divText(`Bimodal: ${split.hiCount} values at ${this.formatEffect(split.hiMean)} ` +
        `(n = ${count(split.hiN)}) and ${split.loCount} at ${this.formatEffect(split.loMean)} ` +
        `(n = ${count(split.loN)}), a gap of ${split.gap.toFixed(2)} against the ` +
        `± ${host.formatActivity(fit.cvRmse!)} this fit resolves.`, 'chem-sar-sum-prose'));
    }

    const best = data.roleBest.get(role.name)!;
    const others = fit.roles.filter((other) => other !== role).map((other) => other.name).join(', ');
    for (const {level, faint} of this.roleRows(role)) {
      const ref = best.get(level.value)!;
      // The axis is the one role pinned in the Vary filter, so its row can land on the position too.
      const position = role.name === data.axisRole ? ref.matrix.columns[ref.ci].position : undefined;
      body.push(this.summaryRow({
        depiction: level.value === '' ? null : this.depiction(level.value, CARD_CORE_W, CARD_CORE_H),
        name: roleValueName(level.value),
        badges: [this.roleCountBadge(role.name, level)],
        desc: `best measured ${host.formatActivity(ref.value)} · ${ref.matrix.label}`,
        mark: this.effectBar(level.coef, scale, fit.cvRmse),
        value: this.formatEffect(level.coef),
        caption: 'offset',
        valueTip: `${this.formatEffect(level.coef)} against the library mean of ` +
          `${host.formatActivity(fit.mean)}, over ${count(level.n)} compounds, adjusted for ${others}. ` +
          'An offset, not a potency.',
        faint,
        tip: role.name === data.coreRole ? `Open ${ref.matrix.label} in the SAR Matrix` :
          this.openTip(ref.matrix, ref.ri, ref.ci),
        onClick: () => role.name === data.coreRole ? host.revealMatrix(ref.matrix) :
          host.revealCell(ref.matrix, ref.ri, ref.ci, position),
      }));
    }
    body.push(...dropped);
    return this.card(title, subtitle, body, footer);
  }

  /** The compound count, which also selects those compounds. An offset is a statement about a subset
   *  of the table, and the only way to check one is to look at the subset it is about. */
  private roleCountBadge(role: string, level: RoleLevel): HTMLElement {
    const badge = this.badge(`${count(level.n)} cpd`, 'Measured compounds carrying this value — click ' +
      'to select them in the table. No ± is printed: one from this count alone would ignore the ' +
      'degrees of freedom the other components cost.');
    badge.classList.add('chem-sar-sum-click');
    badge.onclick = (e: MouseEvent) => {
      e.stopPropagation();
      this.host.selectRoleValue(role, level.value);
    };
    return badge;
  }

  /** The negative clause and the reconciliation, on every role card: each is now the only card its
   *  reader can see, so neither can be carried by a neighbour. */
  private roleFooter(): HTMLElement {
    return ui.divV([
      this.prose('Offset = the mean difference of the compounds carrying this value, adjusted for the ' +
        'other components.',
      'Observational, not causal: it is what the compounds that happen to carry this value did, so it ' +
        'is neither a potency nor a prediction.'),
      this.prose(`Fitted over the whole table. The "${PANE_SERIES}" tab counts instead.`,
        `A count on "${PANE_SERIES}" only sees comparisons that were actually made inside one row — a ` +
        'substitution nobody tried side by side is invisible there, and visible here.'),
    ]);
  }

  // ---- Cores ----------------------------------------------------------------------------------

  /** Series a core leaderboard may rank at all: a matrix has to BE one core, since elsewhere a core is
   *  one row of one matrix and its fitted mean compares groupings, not scaffolds. */
  private rankableCores(data: SummaryData): SeriesStat[] {
    return data.coresAreSeries ?
      data.series.filter((stat) => stat.realCells > 0 && stat.cpd >= MIN_SUPPORT) : [];
  }

  /** The fold tiers that carry rankable cores, shallowest first. */
  private coreTiers(data: SummaryData): number[] {
    return [...new Set(this.rankableCores(data).map((stat) => stat.tier))].sort((a, b) => a - b);
  }

  /** The tier ranked when the reader has not picked one: the one holding the most cores, which is the
   *  cut depth this library actually recurs at. */
  private defaultCoreTier(data: SummaryData): number | null {
    const tiers = this.coreTiers(data);
    if (tiers.length === 0)
      return null;
    const counts = new Map(tiers.map((tier) =>
      [tier, this.rankableCores(data).filter((stat) => stat.tier === tier).length]));
    return tiers.reduce((best, tier) => counts.get(tier)! > counts.get(best)! ? tier : best, tiers[0]);
  }

  /**
   * One tier's cores, best first.
   *
   * Never two tiers at once: a folded tier holds the same compounds over a core cut one bond broader,
   * so a mixed ranking enters one compound twice under two cores and reads as two scaffolds agreeing.
   */
  private topCores(data: SummaryData): SeriesStat[] {
    const tier = this.tierFilter ?? this.defaultCoreTier(data);
    if (tier === null)
      return [];
    const dir = this.host.higherIsBetter ? 1 : -1;
    return this.rankableCores(data).filter((stat) => stat.tier === tier)
      .sort((a, b) => dir * (b.typical - a.typical) || (a.matrix.id < b.matrix.id ? -1 : 1))
      .slice(0, SUM_ROWS);
  }

  private buildCoreCard(data: SummaryData, cores: SeriesStat[]): HTMLElement {
    const host = this.host;
    const title = `${data.coreRole ?? 'Cores'} — the average of what was made on each`;
    // Not matrix.scores: the preferred score is a raw column mean, which is exactly the statistic the
    // fitted cards were built to replace, and the potency score reads counts captured before the prune.
    const subtitle = 'The mean activity of each core\'s measured compounds, in the activity column\'s ' +
      'own units. Not corrected for which substituents each core was paired with, so a core only ever ' +
      'tried with good groups scores high — but unlike a fitted offset, these do compare between cores.';
    if (!data.coresAreSeries)
      return this.card(title, subtitle, this.reason(CORES_NOT_COMPARABLE));
    if (cores.length === 0)
      return this.card(title, subtitle, this.reason(`No core holds ${MIN_SUPPORT} measured compounds.`));
    const body: HTMLElement[] = [];
    body.push(...cores.map((stat) => this.summaryRow({
      depiction: this.depiction(stat.matrix.rows[0]?.coreSmiles ?? null, CARD_CORE_W, CARD_CORE_H),
      name: stat.matrix.label,
      badges: [this.trustDot(stat.matrix)],
      desc: `${count(stat.cpd)} cpd · best ${host.formatActivity(stat.best!.value)} · ` +
        `${host.formatActivity(stat.lo)}–${host.formatActivity(stat.hi)}`,
      value: host.formatActivity(stat.typical),
      caption: 'typical',
      valueTip: 'The fitted mean of this core\'s measured compounds, count-weighted. It is the mean of ' +
        'what was made: a core only ever paired with good caps scores high, and a heavily resampled ' +
        'substituent dominates it.',
      tip: `Open ${stat.matrix.label} in the SAR Matrix`,
      onClick: () => host.revealMatrix(stat.matrix),
    })));
    const dir = host.higherIsBetter ? 1 : -1;
    let holder: SeriesStat | null = null;
    // Inside the ranked tier: a compound on a core cut at another depth is not an alternative to the
    // cores listed here, so naming it would set up a comparison the card refuses to make.
    for (const stat of this.rankableCores(data).filter((s) => s.tier === cores[0].tier)) {
      if (stat.best !== null && (holder === null || dir * stat.best.value > dir * holder.best!.value))
        holder = stat;
    }
    // The carry-forward decision nothing else on the screen states: the best single compound and the
    // best scaffold on average are routinely not the same core.
    if (holder !== null && holder.matrix.id !== cores[0].matrix.id) {
      body.push(ui.divText(`${cores[0].matrix.label} runs highest on average. Your single best compound ` +
        `is on ${holder.matrix.label}, whose typical is ${host.formatActivity(holder.typical)}: one ` +
        'lucky substituent, not a better scaffold.', 'chem-sar-cp-hint'));
    }
    return this.card(title, subtitle, body);
  }

  /** Compounds the dataset already holds with no activity value. Deliberately not gated on R²: a plate
   *  is cheap and a synthesis is not. */
  private buildShelfSegment(data: SummaryData): HTMLElement {
    const host = this.host;
    const body: HTMLElement[] = [];
    body.push(ui.divText('Already in the table, never assayed — test these first',
      'chem-sar-sum-card-title'));
    body.push(this.prose('No synthesis needed: these compounds are rows of your table that carry no ' +
      'value in the activity column, ranked by what the model predicts for them.',
    'They pass no fit gate, only support ≥ 2 measured compounds on both axes — running an assay on a ' +
      'compound you already have is cheap, so the bar it has to clear is lower than for a synthesis.'));
    // Filtered inside the pool's own passes, not afterwards: a pure extrapolation at support 1 would
    // otherwise top the first card a chemist sees, with nothing but a faint value to say so.
    const shelf = data.test.take((row) => supportOf(row) < 2);
    if (shelf.length === 0) {
      body.push(...this.reason(host.predictUnmeasured ?
        'No untested compound rests on two measured compounds on both axes.' :
        '"Predict untested compounds" is off, so untested compounds carry no predicted value to rank by.'));
    }
    for (const row of shelf) {
      const {matrix, ri, ci} = row;
      const cell = matrix.cells[ri][ci];
      const support = cell.support ?? 0;
      body.push(this.summaryRow({
        depiction: this.depiction(cell.smiles, CARD_CORE_W, CARD_CORE_H),
        name: host.cellIdText(cell) ?? '(no id)',
        badges: [
          this.badge(`n=${support}`, 'Measured compounds backing the prediction on the weaker of this ' +
            'core and this substituent.'),
          this.trustDot(matrix),
        ],
        desc: this.cellLocation(matrix, ri, ci),
        value: `~${host.formatActivity(cell.value!)}`,
        caption: 'predicted',
        cart: {matrix, ri, ci},
        valueTip: 'Predicted, not measured — the compound exists, only the number is estimated. This ' +
          'segment is deliberately not gated on R²: a plate is cheap and a synthesis is not.',
        tip: this.openTip(matrix, ri, ci),
        onClick: () => host.revealCell(matrix, ri, ci),
      }));
    }
    if (shelf.length > 0) {
      const refs: MatrixCellRef[] = shelf.map((row) => ({matrix: row.matrix, ri: row.ri, ci: row.ci}));
      body.push(ui.divH([ui.button(refs.length === 1 ? 'Queue this compound for testing' :
        `Queue these ${refs.length} for testing`,
      () => host.addCellsToMakeList(refs, 'Nothing to add.'))], 'chem-sar-sum-foot'));
    }

    return ui.divV(body, 'chem-sar-sum-shelf');
  }

  // ---- Worth making ------------------------------------------------------------------------------

  /**
   * The analogs the model says are worth making, ranked on the gain each buys over the best
   * compound its own series has already made, in that series' own leave-one-out error.
   *
   * Not on the predicted value. Under an additive model the global maximum is the best row crossed
   * with the best column — one obvious corner per series, which a chemist reads straight off the grid;
   * the intercept it carries is the series' own observed mean, so ranking on it ranks series by what
   * was already made; and a high prediction is a reason to synthesise only if it beats the best
   * compound already measured. The division by the series' own error is a reliability weighting:
   * a +1.2 where predictions are typically off by 1.0 is a worse offer than a +0.7 where they are off
   * by 0.15.
   */
  private buildAnalogBlock(data: SummaryData): HTMLElement {
    const host = this.host;
    const body: HTMLElement[] = [
      ui.divText('Worth making — ranked by gain over what this series has already made',
        'chem-sar-sum-card-title'),
      ui.divText('Every number here is a model output: this series\' fitted core effect plus its ' +
        'fitted substituent effect. Nothing here has been measured. The additive model cannot see a ' +
        'cliff, so a large gain is a hypothesis to test and a small one is not a refutation.',
      'chem-sar-cp-hint'),
    ];
    if (data.trustedNoStructure > 0) {
      body.push(ui.divText(`${count(data.trustedNoStructure)} further predictions pass the same gate but ` +
        'have no structure to make: the core carries an attachment point none of the picked fragment ' +
        'columns fills. Add that column and they appear here.', 'chem-sar-cp-hint'));
    }

    const ranked = data.analogs.all;
    const thin = data.analogsThin.all;
    const fallback = ranked.length === 0 && thin.length === 0;
    const spare = fallback ? data.analogsAny.all : [];
    if (fallback && spare.length === 0) {
      body.push(...this.reason(!host.predictVirtual ?
        'Predicted analogs are off — turn on "Predict virtual analogs" to see what to make.' :
        'No series predicts anything above the best compound it has already measured.'));
      return ui.divV(body, 'chem-sar-sum-making-block');
    }
    if (fallback) {
      // The gate withholds confidence, not the ranking. Naming the best candidates and saying what they
      // are short of leaves a next structure on the screen; refusing to name any leaves none.
      body.push(ui.divText(`Nothing clears the trust gate, so these are the largest predicted gains ` +
        `without it — ranked on gain alone, not on gain over error. ` +
        `${this.gateShortfall(data)}`, 'chem-sar-sum-warn'));
    }

    const rows: AnalogRow[] = [...ranked.map((row) => ({row, evidence: 'ranked'})),
      ...thin.map((row) => ({row, evidence: 'thin'})),
      ...spare.map((row) => ({row, evidence: 'ungated'}))];
    const structures = new Set(rows.map((r) => r.row.key)).size;
    const series = new Set(rows.map((r) => r.row.matrix.id)).size;
    const dropped = data.withheldBelowError + data.withheldThinSupport + data.withheldFitFails +
      data.withheldUnchecked + data.withheldNotConverged;
    const chips: HTMLElement[] = [
      this.badge(`${count(structures)} analogs · ${count(series)} series`,
        'Molecules, not cells: one structure is proposed by every tier that folded its core, and only ' +
        'the best-evidenced occurrence is listed.'),
    ];
    if (dropped > 0) {
      chips.push(this.badge(`${count(dropped)} dropped: ${count(data.withheldBelowError)} gain under ` +
        `this series' own error · ${count(data.withheldThinSupport)} thin support`,
      'A prediction whose gain is smaller than its series\' own typical prediction error is one the ' +
      'model cannot distinguish from the best compound that series has already measured. Thin support ' +
      `means fewer than ${MIN_SUPPORT} measured compounds on one of the two axes.`, true));
    }
    // Only the checked-and-failed cause has a list: a fit that was never cross-validated and a fit that
    // stopped short of its tolerance appear on no screen, so they say so rather than pointing at one.
    if (data.withheldFitFails > 0) {
      const drop = this.badge(`${count(data.withheldFitFails)} in a fit checked and failed →`,
        `These series' own R² is below ${TRUST_R2}. Click for the list.`, true);
      drop.onclick = () => this.revealTrust();
      chips.push(drop);
    }
    if (data.withheldUnchecked + data.withheldNotConverged > 0) {
      chips.push(this.badge(`${count(data.withheldUnchecked)} fit never checked · ` +
        `${count(data.withheldNotConverged)} fit did not converge`,
      'Unchecked means the series has too few cross-validatable cells to verify — unverified, not ' +
      'wrong. A fit that did not reach its tolerance is excluded from every pooled comparison, and R² ' +
      'says nothing about it. Neither has a list to open.', true));
    }
    if (data.alreadyHeld > 0) {
      chips.push(this.badge(`${count(data.alreadyHeld)} already in the table — see the block above`,
        'The structure is one the table already carries, so this is an assay to run rather than a ' +
        'synthesis.'));
    }
    body.push(ui.divH(chips, 'chem-sar-sum-chips'));

    if (data.analogOverflow > 0) {
      body.push(ui.divText(`Ranked and thin hold the top ${count(ANALOG_LIST_MAX)} each by gain in ` +
        `their own series' error: ${count(data.analogOverflow)} further structures cleared the same ` +
        'gate and are in neither. Narrow the analysis to reach them.', 'chem-sar-cp-hint'));
    }

    const refs: MatrixCellRef[] = rows.map((r) => ({matrix: r.row.matrix, ri: r.row.ri, ci: r.row.ci}));
    body.push(ui.divH([
      // "These", not "all": the list is capped, and a button promising all of them would be the cap
      // stated as a total for the third time.
      ui.button(`Add these ${refs.length} to Make list`,
        () => host.addCellsToMakeList(refs, 'Nothing to add.')),
      ui.button('Add selected', () => host.addCellsToMakeList(this.selectedAnalogs(),
        'Select rows in the list first.')),
    ], 'chem-sar-sum-foot'));
    if (thin.length > 0) {
      body.push(ui.divText(ranked.length === 0 ?
        'Thin evidence — fewer than ' + BEST_FIT_MIN_N + ' cross-validatable cells behind every one of ' +
        'these series\' fits, so none of them is ranked.' :
        `The ${count(thin.length)} rows marked thin rest on fewer than ${BEST_FIT_MIN_N} ` +
        'cross-validatable cells and are not ranked against the rows above them.', 'chem-sar-cp-hint'));
    }
    body.push(this.buildAnalogGrid(data, rows));
    return ui.divV(body, 'chem-sar-sum-making-block');
  }

  /** What the best ungated candidates are short of, in the order a reader can act on: a fit nobody
   *  checked is a different request from a fit that was checked and failed. */
  private gateShortfall(data: SummaryData): string {
    if (this.host.matrices.every((matrix) => !matrix.confidence))
      return 'No series has a cross-validated fit, so none of them could be checked.';
    const parts: string[] = [];
    if (data.withheldThinSupport > 0) {
      parts.push(`${count(data.withheldThinSupport)} rest on fewer than ${MIN_SUPPORT} measured ` +
        'compounds on one axis');
    }
    if (data.withheldFitFails > 0)
      parts.push(`${count(data.withheldFitFails)} sit in a series whose fit was checked and failed`);
    if (data.withheldUnchecked > 0)
      parts.push(`${count(data.withheldUnchecked)} sit in a series with too few cells to check`);
    return parts.length === 0 ? '' : `Of what was withheld, ${parts.join('; ')}.`;
  }

  private selectedAnalogs(): MatrixCellRef[] {
    const grid = this.analogGrid;
    if (grid === null)
      return [];
    const selection = grid.dataFrame.selection;
    return this.analogSources.filter((_ref, i) => selection.get(i));
  }

  /** A DG.Grid, because browsing is the point: the chemist sorts by support, filters to one series,
   *  selects fifteen and sends them to the Make list. */
  private buildAnalogGrid(data: SummaryData, rows: AnalogRow[]): HTMLElement {
    const host = this.host;
    const dir = host.higherIsBetter ? 1 : -1;
    const byMatrix = new Map(host.matrices.map((m, i) => [m.id, host.matrixTiers[i]]));
    const best = new Map(data.series.map((s) => [s.matrix.id, s.best]));
    // With fragment columns every row of a series can carry the same core structure, and drawing the
    // same picture on every grid row is worse than drawing none.
    const degenerate = rows.every(({row}) =>
      new Set(row.matrix.rows.map((r) => r.keySmiles)).size === 1);

    const analog: string[] = [];
    const predicted: number[] = [];
    const gains: number[] = [];
    const interest: number[] = [];
    const support: number[] = [];
    const r2: number[] = [];
    const rmse: number[] = [];
    const neighbours: number[] = [];
    const seriesName: string[] = [];
    const tier: string[] = [];
    const evidence: string[] = [];
    const core: string[] = [];
    const rgroup: string[] = [];
    const sources: MatrixCellRef[] = [];

    for (const {row, evidence: kind} of rows) {
      const {matrix, ri, ci} = row;
      const cell = matrix.cells[ri][ci];
      // Nullable: an ungated row can come off a series whose fit was never cross-validated at all.
      const conf = matrix.confidence ?? null;
      const anchor = best.get(matrix.id) ?? null;
      const gain = anchor === null ? 0 : dir * (cell.value! - anchor.value);
      analog.push(cell.smiles ?? '');
      predicted.push(cell.value!);
      gains.push(gain);
      interest.push(conf !== null && conf.rmse > 0 ? gain / conf.rmse : NaN);
      support.push(cell.support ?? 0);
      r2.push(conf?.r2 ?? NaN);
      rmse.push(conf?.rmse ?? NaN);
      // Only for the rows that survived the ranking: this is O(rows + columns) per call.
      neighbours.push(host.observedNeighbours(matrix, ri, ci));
      seriesName.push(matrix.label);
      tier.push(`L${byMatrix.get(matrix.id) ?? 1}`);
      evidence.push(kind);
      core.push(degenerate ? Object.values(matrix.rows[ri].foldedValues).join(' · ') :
        matrix.rows[ri].keySmiles);
      rgroup.push(matrix.columns[ci].substSmiles);
      sources.push({matrix, ri, ci});
    }

    const columns = [
      DG.Column.fromStrings(ANALOG_COLS.analog, analog),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.predicted, predicted),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.gain, gains),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.interest, interest),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, ANALOG_COLS.support, support),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.r2, r2),
      DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, ANALOG_COLS.rmse, rmse),
      DG.Column.fromList(DG.COLUMN_TYPE.INT, ANALOG_COLS.neighbours, neighbours),
      DG.Column.fromStrings(ANALOG_COLS.series, seriesName),
      DG.Column.fromStrings(ANALOG_COLS.tier, tier),
      DG.Column.fromStrings(ANALOG_COLS.evidence, evidence),
      DG.Column.fromStrings(degenerate ? ANALOG_COLS.fixed : ANALOG_COLS.core, core),
      DG.Column.fromStrings(ANALOG_COLS.rgroup, rgroup),
    ];
    const descriptions: {[name: string]: string} = {
      [ANALOG_COLS.analog]: 'The assembled structure — a combination this dataset has no row for. It ' +
        'knows nothing of your collection, the catalogue or the literature, so this is "not in this ' +
        'dataset", never "never made".',
      [ANALOG_COLS.predicted]: 'A model output in the host column\'s units: this series\' fitted core ' +
        'effect plus its fitted substituent effect. Nothing here was measured.',
      [ANALOG_COLS.gain]: 'Against the best MEASURED compound of this same series — never against ' +
        'another series, and never against the dataset as a whole.',
      [ANALOG_COLS.interest]: 'The ranking quantity: the gain in multiples of this series\' own ' +
        'leave-one-out error. Under one error the model cannot tell this analog from the compound you ' +
        'already have.',
      [ANALOG_COLS.support]: 'Measured compounds on the weaker of this core and this substituent.',
      [ANALOG_COLS.r2]: 'How well this series\' model predicts a measured compound when that compound ' +
        'is left out of the fit. A property of the model, not of the data.',
      [ANALOG_COLS.rmse]: 'Typical size of this series\' prediction error.',
      [ANALOG_COLS.neighbours]: 'Measured cells sharing this core OR this substituent — how much has ' +
        'already been tried around this combination. Not the size of the series.',
      [ANALOG_COLS.evidence]: `"thin" means fewer than ${BEST_FIT_MIN_N} cross-validatable cells behind ` +
        'this series\' fit, so its error is not a scale to rank by; those rows are listed, not ranked. ' +
        '"ungated" means the row cleared no trust gate and is shown because nothing else did.',
      [degenerate ? ANALOG_COLS.fixed : ANALOG_COLS.core]: degenerate ?
        'The substituents this row holds fixed — every row of these series draws the same core.' :
        'The core this analog is built on.',
    };
    for (const column of columns) {
      const text = descriptions[column.name];
      if (text !== undefined)
        column.setTag(DG.TAGS.DESCRIPTION, text);
    }
    for (const name of [ANALOG_COLS.analog, ANALOG_COLS.rgroup, degenerate ? '' : ANALOG_COLS.core]) {
      const column = columns.find((c) => c.name === name);
      if (column !== undefined)
        column.semType = DG.SEMTYPE.MOLECULE;
    }

    const frame = DG.DataFrame.fromColumns(columns);
    this.analogSources = sources;
    const grid = DG.Viewer.grid(frame);
    this.analogGrid = grid;
    // The grid root is a ui-box, which pins itself to a fixed size and leaves the rest blank.
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    grid.setOptions({rowHeight: CELL_H});
    for (const [name, width] of [[ANALOG_COLS.analog, ANALOG_W], [ANALOG_COLS.core, CORE_W],
      [ANALOG_COLS.rgroup, CARD_CORE_W]] as [string, number][]) {
      const gridCol = grid.col(name);
      if (gridCol)
        gridCol.width = width;
    }
    // A click, not the current row: setting a current row is something the grid does to itself while it
    // is built, and landing on a cell switches the outer tab.
    this.analogSub = grid.onCellClick.subscribe((cell: DG.GridCell) => {
      const at = cell.tableRowIndex ?? -1;
      if (at >= 0 && at < sources.length) {
        const ref = sources[at];
        host.revealCell(ref.matrix, ref.ri, ref.ci);
      }
    });
    return ui.div([grid.root], 'chem-sar-sum-analog-grid');
  }

  // ---- Trust section --------------------------------------------------------------------------

  /**
   * What the two R² on this tab mean and what the gate does with them, then the series at each end of it.
   * The poorly-fitting list opens itself when something elsewhere on the tab links here, since that is
   * the one arrival where the reader came for the list rather than the explanation.
   */
  private buildTrustSection(data: SummaryData): HTMLElement {
    const host = this.host;
    const parts: HTMLElement[] = [
      ui.divText('Fit quality', 'chem-sar-sum-card-title'),
      // Two quantities with one name on one screen, and nothing else says they are not the same number.
      this.prose('Two different R² appear on this tab: one for the whole table, one for each series.',
        'Both are scored by predicting compounds that were left out of the fit, so they say how well ' +
        'the model predicts rather than how well it describes what it was shown. In a series, 1 means ' +
        'its substituent effects add perfectly, 0 means the model does no better than simply using that ' +
        'series\' average, and below 0 it does worse than that average.'),
      this.prose(`Trust gate: ${MIN_SUPPORT} measured compounds on both axes, in a series whose own R² ` +
        'is at least ' + TRUST_R2 + '.',
      'Predictions that fail it are still computed and still drawn in the matrix. They are left out of ' +
        'the counts above, and out of Worth making unless nothing at all clears the gate.'),
    ];
    const worst = [...data.lowR2].sort((a, b) =>
      (a.confidence!.r2 - b.confidence!.r2) || (a.id < b.id ? -1 : 1));
    if (worst.length === 0) {
      parts.push(this.prose('Every cross-validated fit holds up.',
        'Substituent effects add in each of them, so their predictions can be read as estimates.'));
    } else {
      parts.push(this.foldable(`Where the additive model does not hold · ${count(worst.length)}`,
        `${count(data.lowR2Virtual)} predicted cells inherit a fit that reproduces its own measured ` +
        'cells poorly.', this.trustList(worst), this.openTrustList));
    }
    const held = host.matrices.filter((matrix) => matrix.confidence != null &&
      matrix.confidence.r2 >= TRUST_R2)
      .sort((a, b) => (b.confidence!.r2 - a.confidence!.r2) || (a.id < b.id ? -1 : 1));
    if (held.length > 0) {
      parts.push(this.foldable(`Where it holds best · ${count(held.length)}`,
        'Best first. A prediction in one of these rests on a fit that reproduced its own measured cells.',
        this.trustList(held), false));
    }
    return ui.divV(parts, 'chem-sar-sum-trust');
  }

  /** A statement on the pane and its qualification on the hover. Method is read for the one fact it is
   *  opened for, and a paragraph per fact buries that fact in the others. */
  private prose(text: string, tip: string): HTMLElement {
    const el = ui.divText(text, 'chem-sar-sum-prose');
    ui.tooltip.bind(el, () => tip);
    return el;
  }

  /** The same, in the faint style a row uses for its trailing clause. */
  private hint(text: string, tip: string): HTMLElement {
    const el = ui.divText(text, 'chem-sar-cp-hint');
    ui.tooltip.bind(el, () => tip);
    return el;
  }

  /** A chevron heading that shows and hides `body`; the caller owns whatever else sits on the row. */
  private foldHead(title: string, body: HTMLElement, open: boolean, cls = '',
    onToggle?: (open: boolean) => void): HTMLElement {
    const icon = ui.iconFA('chevron-right');
    icon.classList.add('chem-sar-sum-fold-icon');
    icon.classList.toggle('chem-sar-sum-fold-open', open);
    body.style.display = open ? '' : 'none';
    const head = ui.divH([icon, ui.divText(title, 'chem-sar-sum-fold-title')],
      `chem-sar-sum-fold ${cls}`.trim());
    head.onclick = () => {
      const show = body.style.display === 'none';
      body.style.display = show ? '' : 'none';
      icon.classList.toggle('chem-sar-sum-fold-open', show);
      onToggle?.(show);
    };
    return head;
  }

  /** The lead line above the fold stays, since a count is worth reading without opening anything. */
  private foldable(title: string, lead: string, body: HTMLElement, open: boolean): HTMLElement {
    return ui.divV([this.foldHead(title, body, open), ui.divText(lead, 'chem-sar-cp-hint'), body]);
  }

  private trustList(matrices: SarMatrix[]): HTMLElement {
    const host = this.host;
    const list = ui.div(matrices.slice(0, TRUST_LIST_MAX).map((matrix) => {
      const conf = matrix.confidence!;
      return this.summaryRow({
        // "Series 7" names a series without showing one. Every other list on the tab draws the
        // chemistry it is talking about, and a reader deciding whether to trust a series' predictions
        // wants to see which scaffold they are about.
        depiction: this.depiction(matrix.rows[0]?.coreSmiles ?? null, CARD_CORE_W, CARD_CORE_H),
        name: matrix.label,
        badges: [this.badge(`R² ${conf.r2.toFixed(2)} ± ${host.formatActivity(conf.rmse)}`,
          'How well this series\' model predicts a measured cell held out of its own fit, and the ' +
          'typical size of that error. Near 0 or below, substituent effects do not add up here.',
          conf.r2 < TRUST_R2)],
        desc: `${conf.n} of ${conf.total} measured cells checked`,
        value: `${count(matrix.virtualCount)}`,
        caption: 'predicted cells',
        tip: `Open ${matrix.label} in the SAR Matrix`,
        onClick: () => host.revealMatrix(matrix),
      });
    }));
    if (matrices.length > TRUST_LIST_MAX)
      list.appendChild(ui.divText(`+${matrices.length - TRUST_LIST_MAX} more`, 'chem-sar-cp-hint'));
    return list;
  }
}
