/* The Summary tab's data model: what one walk over the matrices produces, the contract it reads
   from the viewer, the thresholds that gate it and the small pure helpers it is built from. A leaf
   module, so the collector, the swap and R-group algorithms and the panel can all depend on it. */
import {RoleFit} from './sar-matrix-role-fit';
import {SarMatrix, SarMatrixCell} from './sar-matrix-types';
import {MatrixCellRef} from './sar-matrix-ui-common';

export const SUM_ROWS = 3;

const SUM_POOL = 24;

/** Below this the leave-one-out fit does not support acting on a prediction. */
export const TRUST_R2 = 0.5;

/** A prediction resting on fewer measured compounds than this on either axis is extrapolation. */
export const MIN_SUPPORT = 3;

/** Measured pairs a pooled swap needs, and lineages it must span, before it is worth a row. */
export const SWAP_MIN_PAIRS = 3;

export const SWAP_MIN_SERIES = 2;

/** Measured cells kept per row when enumerating swaps; both extremes survive the trim. */
export const SWAP_ROW_CAP = 32;

/** Lineages an R-group must win in before it is ranked; two is the mean of two numbers. */
export const RGROUP_MIN_SERIES = 3;

/** Cross-validatable cells a fit needs before its R² is quotable as "best validated". */
export const BEST_FIT_MIN_N = 8;

/** Matrices a per-core outcome strip can carry before a row of squares stops being a mark. */
export const STRIP_SLOTS = 12;

/** Reliable losers shown under the winners. */
export const LOSER_ROWS = 2;

/** Rows the analog list holds. A chemist does not browse ten thousand; two hundred is more than a
 *  quarter's synthesis and fits one frame of rendered structures. */
export const ANALOG_LIST_MAX = 200;

export interface SummaryRow {
  /** What the pool deduplicates on; held rather than re-derived, since only the pool knows it. */
  key: string;
  matrix: SarMatrix;
  ri: number;
  ci: number;
  /** Always higher-is-better, so one comparator orders every pool. */
  score: number;
}

/** What the summary reads back from the viewer, and the things it asks the viewer to do. */
export interface SummaryHost {
  readonly matrices: SarMatrix[];
  /** Lineage root id per matrix, index-aligned with `matrices` — the one hierarchy the navigator uses. */
  readonly matrixRoots: string[];
  readonly matrixTiers: number[];
  readonly assayedCount: number;
  /** Whether SAR transfer detection has run, and what it found. Detection is lazy and lives on its own
   *  tab, so the landing screen reports its state rather than triggering it. */
  readonly transferSummary: {scanned: boolean, count: number};
  readonly higherIsBetter: boolean;
  readonly scalingLabel: string;
  readonly activityIsLog: boolean;
  readonly activityColumnName: string;
  /** Rows of the host table — the denominator the coverage line is a fraction of. */
  readonly hostRowCount: number;
  /** Assayed values the chosen scaling cannot represent; they are in no matrix and are not untested. */
  readonly unscalableCount: number;
  readonly axisRole: string | null;
  /** The column the cores came from, or null when they came from fragmentation. A degrader set names
   *  it Linker, and "Best core" over an unnamed scaffold is the same finding nobody can act on. */
  readonly coreRole: string | null;
  readonly coresAreSeries: boolean;
  readonly roleColumns: string[];
  setColumnAxis(name: string): void;
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

export function supportOf(row: SummaryRow): number {
  return row.matrix.cells[row.ri][row.ci].support ?? 0;
}

/** A bounded best-of pool that deduplicates on insert. Deduplicating a pool of cells afterwards can
 *  leave it shorter than the card has room for, and the same compound occupies a cell in every tier
 *  that folded it. */
export class TopList {
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
export interface RowCell {
  ci: number;
  value: number;
  molIdx: number;
}

/** A measured swap pooled across series: one R-group exchanged inside one row, everything else
 *  identical. `from`/`to` are the lexically ordered fragment pair, so a swap cannot split into two
 *  mirror pools oriented by an arbitrary column index; the display direction is chosen at render. */
export interface SwapPool {
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

/** One matrix where a substituent came first at its position, under the fitted model. */
export interface RGroupWin {
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

export interface RGroupRow {
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
export type StripCols = Map<string, {ci: number, n: number, refDelta: number | null}>;

/** What the R-group card accumulates over the one fit per matrix. */
export interface RGroupAcc {
  wins: Map<string, RGroupWin[]>;
  losses: Map<string, RGroupWin[]>;
  tried: Map<string, Set<string>>;
}

/** Everything one matrix contributes, extracted from its fit and its cells; the fit arrays are
 *  discarded with it, so the walk never holds one per matrix at once. */
export interface SeriesStat {
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
  /** Filled only for the matrices a strip can draw; null everywhere else. */
  stripCols: StripCols | null;
  /** Range of the fitted substituent and core effects — centred, same units, same fit, so the two
   *  are comparable within this matrix and the verdict pools as a count. */
  colRange: number | null;
  rowRange: number | null;
  bestRow: {ri: number, effect: number, n: number} | null;
}

export interface StartRow {
  primary: SeriesStat;
  /** Indices into the fixed reason-lane slots, not phrases: the lane's shape is what a reader
   *  compares down the column. */
  reasons: number[];
}

export interface SummaryData {
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

/** One half of a measured pair: the component value it carries, what it measured, and which compound it
 *  is — the three things a pool needs from each side, whether the side came off a matrix cell or off a
 *  row of component columns. */
export interface SwapSide {
  value: string;
  activity: number;
  mol: number;
}

/** Index of the most potent measured cell of one matrix row or column, so a row of the tab lands on a
 *  compound; -1 where none of them is measured. */
export function bestMeasured(cells: SarMatrixCell[], dir: number): number {
  let best = -1;
  cells.forEach((cell, i) => {
    if (cell.kind === 'real' && cell.value !== null && (best < 0 || dir * cell.value > dir * cells[best].value!))
      best = i;
  });
  return best;
}
