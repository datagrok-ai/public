import * as DG from 'datagrok-api/dg';

/*
 * Every comparator across these modules ends in a content-derived tiebreak: fragments arrive in
 * worker-completion order, so resolving a tie by index or Map order makes the matrices depend on
 * scheduling. Removing one reintroduces that.
 */

export function logSarTime(stage: string, t0: number): void {
  console.log(`SAR Matrix | ${stage}: ${Math.round(performance.now() - t0)} ms`);
}

export interface SeriesMember {
  molIdx: number;
  substSmiles: string;
}

export interface MatchedSeries {
  coreSmiles: string;
  members: SeriesMember[];
}

/** A cluster of structurally related matched series — the rows of one SAR Matrix. */
export interface CoreCluster {
  id: string;
  series: MatchedSeries[];
  /** The remainder every member core reduced to — what this group is keyed on. Empty for a lone
   *  series that no site grouped. Cutting it once more is what relates whole matrices to each other. */
  siteKey: string;
  /** Series are placeholders for compounds that never formed one, so this cluster is worth building
   *  only if its decomposition succeeds. */
  requiresDecomposition?: boolean;
  /** Fragmentation level: 2 for matrices built directly from series, higher for each round that
   *  folds groups of the level below into one. */
  level: number;
  parentId?: string;
  /** The user's own grouping value, used verbatim as the label so the names are the ones they wrote. */
  label?: string;
}

/** `unmeasured` is a compound the dataset already holds whose activity is missing. It is kept apart
 *  from `empty` so the Free-Wilson fill cannot propose synthesizing something that already exists. */
export type SarMatrixCellKind = 'real' | 'virtual' | 'empty' | 'unmeasured';

export interface SarMatrixCell {
  kind: SarMatrixCellKind;
  /** Observed (real) or predicted (virtual) scaled activity; null when empty or unmeasured. */
  value: number | null;
  molIdx: number | null;
  smiles: string | null;
  /** For virtual cells: smaller of the row's and column's real-cell counts backing the Free-Wilson
   *  prediction. Higher means more trustworthy. */
  support?: number;
  /** For virtual cells: sum of the row's and column's real-cell counts. */
  references?: number;
  /** For real cells: the additive (Free-Wilson) fitted value; a large gap to `value` flags a
   *  non-additive cell. Undefined when predictions weren't computed. */
  fit?: number;
}

export interface SarMatrixRow {
  coreSmiles: string;
  /** The core carrying this row's folded substituents, with the column position still open. Falls
   *  back to the core when nothing is folded or the fragments cannot be linked. */
  keySmiles: string;
  label: string;
  /** Substituents fixed by this row, for every position that is not the column axis. Empty for
   *  single-position matrices. */
  foldedValues: {[position: string]: string};
}

export interface SarMatrixColumn {
  /** e.g. 'R1', 'R2'. Single-position fallback matrices use the constant 'R1'. */
  position: string;
  substSmiles: string;
}

/** A single SAR Matrix: related cores (rows) × substituents (columns). */
export interface SarMatrix {
  id: string;
  /** Stable human label ("Series A", "Series B", ...), assigned before ranking. */
  label: string;
  rows: SarMatrixRow[];
  columns: SarMatrixColumn[];
  /** cells[rowIndex][columnIndex]. */
  cells: SarMatrixCell[][];
  minActivity: number;
  maxActivity: number;
  realCount: number;
  virtualCount: number;
  scores: {[scheme: string]: number};
  /** Active R-group positions, richest first. Length 1 means a single-position (fallback) view. */
  positions: string[];
  /** Reference (most-frequent-observed) substituent per position — pinned while one position varies. */
  refValues: {[position: string]: string};
  siteKey: string;
  /** Fragmentation level (2 = built directly from series). Higher levels hold the same compounds
   *  over cores that agree one cut deeper — fewer rows covering more chemistry. */
  level: number;
  parentId?: string;
  /** LOO cross-validated quality of the Free-Wilson fit; null when too few observations to
   *  cross-validate. `n` cross-validatable observed cells out of `total`. */
  confidence?: {r2: number, rmse: number, n: number, total: number} | null;
}

/** Close a grid, tolerating one that cannot: a view-less grid may not support close, and dropping
 *  the reference is enough. */
export function closeGridQuietly(grid: DG.Grid | null | undefined): void {
  try {
    grid?.close?.();
  } catch (e) {
    // Nothing to do; the caller drops its reference either way.
  }
}

/** A score fit to sit in a numeric column: anything non-finite becomes empty, since a sentinel like
 *  -Infinity would stretch the column's histogram across the whole axis. */
export function finiteOrNaN(value: number | undefined): number {
  return value === undefined || !Number.isFinite(value) ? NaN : value;
}

