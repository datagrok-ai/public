/* The Summary's data model: what the walk produces, what it reads from the viewer, and its thresholds.
   A leaf module, so everything else in the Summary can depend on it. */
import {RoleFit} from './sar-matrix-role-fit';
import {SarMatrix, SarMatrixCell} from '../sar-matrix-types';
import {MatrixCellRef} from '../sar-matrix-ui-common';

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

/** Reliable losers shown under the winners. */
export const LOSER_ROWS = 2;

/** Rows the analog list holds: more than a quarter's synthesis. */
export const ANALOG_LIST_MAX = 200;

/** A sparse series needs hundreds of sweeps to reach the fit tolerance; under the matrix fill's cap of 50
 *  many would read as unconverged and drop out of every pooled comparison. */
export const SERIES_FIT_SWEEPS = 2000;

/** Measured compounds a prediction must rest on, on both axes, before it is listed at all. */
export const MIN_LISTED_SUPPORT = 2;

export interface SummaryRow {
  /** What the pool deduplicates on. */
  key: string;
  matrix: SarMatrix;
  ri: number;
  ci: number;
  /** Always higher-is-better, so one comparator orders every pool. */
  score: number;
}

/** What the Summary reads from the viewer, and asks it to do. */
export interface SummaryHost {
  readonly matrices: SarMatrix[];
  /** Lineage root id per matrix, index-aligned with `matrices`. */
  readonly matrixRoots: string[];
  readonly matrixTiers: number[];
  readonly assayedCount: number;
  readonly higherIsBetter: boolean;
  readonly scalingLabel: string;
  readonly activityIsLog: boolean;
  readonly activityColumnName: string;
  /** The denominator of the coverage line. */
  readonly hostRowCount: number;
  /** Assayed values the scaling cannot represent: in no matrix, yet not untested either. */
  readonly unscalableCount: number;
  readonly axisRole: string | null;
  /** The column the cores came from, or null when they came from fragmentation; it names the core tile. */
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
  addCellsToMakeList(cells: MatrixCellRef[], emptyMessage: string): void;
}

export function supportOf(row: SummaryRow): number {
  return row.matrix.cells[row.ri][row.ci].support ?? 0;
}

/** A bounded best-of pool that deduplicates on insert: one compound has a cell in every tier that folded
 *  it, and deduplicating afterwards would leave the list short. */
export class TopList {
  private readonly rows: SummaryRow[] = [];
  private readonly byKey = new Map<string, SummaryRow>();
  /** Keys that qualified but did not fit. A set, not a counter, so one structure counts once across tiers. */
  private readonly missed = new Set<string>();

  /** `better` picks which occurrence of a key to keep; the score cannot, being the same everywhere. */
  constructor(private readonly better: (a: SummaryRow, b: SummaryRow) => boolean,
    private readonly cap = SUM_POOL) {}

  /** The whole pool, best first. */
  get all(): readonly SummaryRow[] {
    return this.rows;
  }

  /** Structures that qualified past the cap, so a screen printing `all.length` can say it is a top N. */
  get overflow(): number {
    return this.missed.size;
  }

  /** Drop keys disqualified after the walk, such as structures the dataset turns out to hold. */
  prune(drop: (key: string) => boolean): void {
    for (let i = this.rows.length - 1; i >= 0; i--) {
      if (!drop(this.rows[i].key))
        continue;
      this.byKey.delete(this.rows[i].key);
      this.rows.splice(i, 1);
    }
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

  /** One row per series first, then top up from the rest. */
  take(): SummaryRow[] {
    const out: SummaryRow[] = [];
    const seen = new Set<string>();
    for (const row of this.rows) {
      if (out.length === SUM_ROWS)
        break;
      if (seen.has(row.matrix.id))
        continue;
      seen.add(row.matrix.id);
      out.push(row);
    }
    for (const row of this.rows) {
      if (out.length === SUM_ROWS)
        break;
      if (!out.includes(row))
        out.push(row);
    }
    return out;
  }
}

/** One measured cell of the row being walked, held only while that row's pairs are enumerated. */
export interface RowCell {
  ci: number;
  value: number;
  molIdx: number;
}

/** A measured swap pooled across series. `from`/`to` are the lexically ordered pair, so a swap never
 *  splits into two mirror pools; the display direction is chosen at render. */
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
  /** Compound pairs already counted: one pair recurs in every tier that folded it. */
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
  /** Measured cells behind the column; the tiebreak where there is no reference. */
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
  /** Lineages that tried this substituent: the denominator of "best in 5 of 6". */
  tried: number;
  magnitude: number | null;
  lo: number;
  hi: number;
  magnitudeSeries: number;
}

/** What the R-group card accumulates over the one fit per matrix. */
export interface RGroupAcc {
  wins: Map<string, RGroupWin[]>;
  losses: Map<string, RGroupWin[]>;
  tried: Map<string, Set<string>>;
}

/** What one matrix contributes; its fit arrays are dropped, so the walk never holds them all at once. */
export interface SeriesStat {
  matrix: SarMatrix;
  root: string;
  tier: number;
  /** Distinct measured molecules; not `realCount`, which is captured before pruning. */
  cpd: number;
  realCells: number;
  totalCells: number;
  /** Cells the decomposition cannot express. Zero also where assembly never marked any. */
  impossibleCells: number;
  lo: number;
  hi: number;
  best: {ri: number, ci: number, value: number} | null;
  /** The fit's count-weighted mean: the one fitted quantity comparable between matrices. */
  typical: number;
  /** Most potent prediction that passes the trust gate. */
  bestVirtual: {ri: number, ci: number, value: number} | null;
  /** Unfilled cells with a prediction; with R² it says why a series has nothing to aim at. */
  virtualCells: number;
  trusted: number;
  /** Range of the fitted substituent and core effects, comparable to each other within the matrix. */
  colRange: number | null;
  rowRange: number | null;
  bestRow: {ri: number, effect: number, n: number} | null;
}

export interface StartRow {
  primary: SeriesStat;
  /** Indices into `REASON_WORDS`. */
  reasons: number[];
}

export interface SummaryData {
  compounds: number;
  untested: number;
  measuredCells: number;
  minObserved: number | null;
  maxObserved: number | null;
  /** Predictions that pass the gate but have no structure: an attachment point no picked column fills. */
  trustedNoStructure: number;
  /** Median leave-one-out error across series; null, never zero, when no series has a checked fit. */
  modelError: number | null;
  fitHolds: number;
  unchecked: number;
  lowR2: SarMatrix[];
  axisRole: string | null;
  coreRole: string | null;
  coresAreSeries: boolean;
  /** One fit over every component column; null when the cores came from fragmentation. */
  roleFit: RoleFit | null;
  /** Best measured cell per value of each component column, keyed by column name. */
  roleBest: Map<string, Map<string, MatrixCellRef & {value: number}>>;
  /** Series whose fit stopped short; left out of every pooled comparison, so a card can name the cause. */
  nonConverged: number;
  series: SeriesStat[];
  startHere: StartRow[];
  swaps: SwapPool[];
  /** The same per component column, grouped off the component values, so no component needs a rebuild. */
  swapsByRole: Map<string, SwapPool[]>;
  /** Whether any candidate existed: "none clears the gate" is not "no row has two measured cells". */
  swapCandidates: number;
  rgroups: RGroupRow[];
  rgroupsThin: RGroupRow[];
  rgroupLosers: RGroupRow[];
  test: TopList;
  /** Analogs ranked on gain over their series' leave-one-out error, deduplicated on structure. */
  analogs: TopList;
  /** Same gates except the cell floor: a tiny fit's error is no scale to rank by, so these are listed. */
  analogsThin: TopList;
  /** The best of what the gate turned down, by raw gain; shown only when both lists above are empty. */
  analogsAny: TopList;

  /** Why a prediction is not offered, kept apart because each calls for a different fix. */
  withheldBelowError: number;
  withheldThinSupport: number;
  withheldFitFails: number;
  withheldUnchecked: number;
  withheldNotConverged: number;
  /** Predicted structures the dataset already holds: an assay to run, not a synthesis. */
  alreadyHeld: number;
  /** Qualifying structures past the row cap, so "N worth making" can say it is a top N. */
  analogOverflow: number;
}

/** One side of a measured pair, from a matrix cell or from a row of component columns. */
export interface SwapSide {
  value: string;
  activity: number;
  mol: number;
}

/** Index of the most potent measured cell in one row or column, or -1. */
export function bestMeasured(cells: SarMatrixCell[], dir: number): number {
  let best = -1;
  cells.forEach((cell, i) => {
    if (cell.kind === 'real' && cell.value !== null && (best < 0 || dir * cell.value > dir * cells[best].value!))
      best = i;
  });
  return best;
}
