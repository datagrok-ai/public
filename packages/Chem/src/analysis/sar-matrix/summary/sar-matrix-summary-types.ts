import {RoleFit} from './sar-matrix-role-fit';
import {SarMatrix, SarMatrixCell} from '../sar-matrix-types';
import {MatrixCellRef} from '../sar-matrix-ui-common';

export const SUM_ROWS = 3;

const SUM_POOL = 24;

export const TRUST_R2 = 0.5;

export const MIN_SUPPORT = 3;

export const SWAP_MIN_PAIRS = 3;

export const SWAP_MIN_SERIES = 2;

/** Measured cells kept per row when enumerating swaps; both extremes survive the trim. */
export const SWAP_ROW_CAP = 32;

/** Below three lineages a median is only the mean of two. */
export const RGROUP_MIN_SERIES = 3;

/** The fewest lineages an R-group row is shown from at all: below RGROUP_MIN_SERIES it is listed as thin. */
export const RGROUP_THIN_SERIES = 2;

export const BEST_FIT_MIN_N = 8;

export const LOSER_ROWS = 2;

/** Rows the analog list holds: more than a quarter's synthesis. */
export const ANALOG_LIST_MAX = 200;

/** Sparse series need hundreds of sweeps; under the matrix fill's cap of 50 many would read as unconverged. */
export const SERIES_FIT_SWEEPS = 2000;

export const MIN_LISTED_SUPPORT = 2;

export interface SummaryRow {
  key: string;
  matrix: SarMatrix;
  ri: number;
  ci: number;
  /** Always higher-is-better, so one comparator orders every pool. */
  score: number;
}

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
  readonly hostRowCount: number;
  /** Assayed values the scaling cannot represent: in no matrix, yet not untested either. */
  readonly unscalableCount: number;
  readonly axisRole: string | null;
  /** Null when the cores came from fragmentation. */
  readonly coreRole: string | null;
  readonly coresAreSeries: boolean;
  readonly roleColumns: string[];
  selectRoleValue(role: string, value: string): void;
  selectRows(rows: number[]): void;
  readonly predictVirtual: boolean;
  readonly computing: boolean;
  noMatricesMessage(): string;
  formatActivity(value: number): string;
  cellIdText(cell: SarMatrixCell): string | null;
  observedNeighbours(matrix: SarMatrix, ri: number, ci: number): number;
  revealCell(matrix: SarMatrix, ri: number, ci: number): void;
  revealMatrix(matrix: SarMatrix): void;
  addCellsToMakeList(cells: MatrixCellRef[], emptyMessage: string): void;
}

export function supportOf(row: SummaryRow): number {
  return row.matrix.cells[row.ri][row.ci].support ?? 0;
}

/** Dedupes on insert: a compound recurs in every tier that folded it, so a later dedupe leaves the list short. */
export class TopList {
  private readonly rows: SummaryRow[] = [];
  private readonly byKey = new Map<string, SummaryRow>();
  /** A set, not a counter, so one structure counts once across tiers. */
  private readonly missed = new Set<string>();

  /** `better` picks which occurrence of a key to keep; the score cannot, being the same everywhere. */
  constructor(private readonly better: (a: SummaryRow, b: SummaryRow) => boolean,
    private readonly cap = SUM_POOL) {}

  get all(): readonly SummaryRow[] {
    return this.rows;
  }

  get overflow(): number {
    return this.missed.size;
  }

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

  /** Primitive parameters so a full pool rejects a weak cell before anything is allocated. */
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

export interface RowCell {
  ci: number;
  value: number;
  molIdx: number;
}

/** `from`/`to` are lexically ordered so a swap never splits into two mirror pools. */
export interface SwapPool {
  from: string;
  to: string;
  n: number;
  sum: number;
  min: number;
  max: number;
  nUp: number;
  roots: Set<string>;
  /** One compound pair recurs in every tier that folded it. */
  seen: Set<string>;
  best: {matrix: SarMatrix, ri: number, ciFrom: number, ciTo: number, delta: number} | null;
  sampled: boolean;
}

export interface RGroupWin {
  matrix: SarMatrix;
  ri: number;
  ci: number;
  root: string;
  n: number;
  /** Relative to the reference column, so it is gauge-free and pools across matrices. */
  refDelta: number | null;
  fitHolds: boolean;
}

export interface RGroupRow {
  subst: string;
  win: RGroupWin;
  /** Lineages won after the trust gate, and before it. */
  k: number;
  m: number;
  tried: number;
  magnitude: number | null;
  lo: number;
  hi: number;
  magnitudeSeries: number;
}

export interface RGroupAcc {
  wins: Map<string, RGroupWin[]>;
  losses: Map<string, RGroupWin[]>;
  tried: Map<string, Set<string>>;
}

/** Scalars only, so the walk never holds every matrix's fit arrays at once. */
export interface SeriesStat {
  matrix: SarMatrix;
  tier: number;
  realCells: number;
  totalCells: number;
  /** Zero also where assembly never marked any. */
  impossibleCells: number;
  best: {ri: number, ci: number, value: number} | null;
  /** The one fitted quantity comparable between matrices. */
  typical: number;
  bestVirtual: {ri: number, ci: number, value: number} | null;
  virtualCells: number;
  trusted: number;
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
  /** Null, never zero, when no series has a checked fit. */
  modelError: number | null;
  fitHolds: number;
  unchecked: number;
  lowR2: number;
  roleFit: RoleFit | null;
  roleBest: Map<string, Map<string, MatrixCellRef & {value: number}>>;
  nonConverged: number;
  series: SeriesStat[];
  startHere: StartRow[];
  swaps: SwapPool[];
  swapsByRole: Map<string, SwapPool[]>;
  /** The table column each component reads, by the name the Summary shows: a merged twin has one name for two. */
  roleColumn: Map<string, string>;
  /** "None clears the gate" is not "no row has two measured cells". */
  swapCandidates: number;
  rgroups: RGroupRow[];
  rgroupsThin: RGroupRow[];
  rgroupLosers: RGroupRow[];
  test: TopList;
  analogs: TopList;
  /** Same gates except the cell floor: a tiny fit's error is no scale to rank by. */
  analogsThin: TopList;
  /** Ungated fallback. */
  analogsAny: TopList;

  /** Kept apart because each calls for a different fix. */
  withheldThinSupport: number;
  withheldFitFails: number;
  withheldUnchecked: number;
  analogOverflow: number;
}

