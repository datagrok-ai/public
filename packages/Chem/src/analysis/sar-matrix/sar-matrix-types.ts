/**
 * Shared type definitions for the SAR Matrix analysis.
 *
 * A SAR Matrix groups structurally related cores (rows) against their shared
 * substituents (columns), colors each real compound by potency, and fills the
 * unmade core-substituent combinations with predicted virtual analogs.
 *
 * v1 is single-position: each matched series varies at one MMP cut site, so the
 * varying substituent is the column directly (no R-group decomposition). A
 * multi-position (R1/R2) view is a planned follow-up.
 */

/** One member of a matched series: a molecule and its varying substituent. */
export interface SeriesMember {
  /** Index of the molecule in the parent dataframe. */
  molIdx: number;
  /** SMILES of the varying substituent. */
  substSmiles: string;
}

/** A matched series: molecules sharing the same MMP core, differing at one site. */
export interface MatchedSeries {
  /** SMILES of the shared core. */
  coreSmiles: string;
  members: SeriesMember[];
}

/** A cluster of structurally related matched series — the rows of one SAR Matrix. */
export interface CoreCluster {
  id: string;
  /** One series per related core; each becomes a matrix row. */
  series: MatchedSeries[];
  /** The remainder every member core reduced to — what this group is keyed on. Empty for a lone
   *  series that no site grouped. Cutting it once more is what relates whole matrices to each other. */
  siteKey: string;
  /** Fragmentation level that produced this group: 2 for the matrices built directly from series,
   *  higher for each round that folds groups of the level below into one. */
  level: number;
  /** The coarser group this one was folded into, when a further level grouped it. */
  parentId?: string;
}

export type SarMatrixCellKind = 'real' | 'virtual' | 'empty';

/** One cell of a SAR Matrix. */
export interface SarMatrixCell {
  kind: SarMatrixCellKind;
  /** Observed (real) or predicted (virtual) scaled activity; null when empty. */
  value: number | null;
  /** Parent-df molecule index for real cells; null otherwise. */
  molIdx: number | null;
  /** Compound SMILES for real cells; null otherwise. */
  smiles: string | null;
  /** For virtual cells: how many observations back the Free-Wilson prediction — the smaller of the
   *  row's and column's real-cell counts. Higher means more trustworthy. Undefined for real/empty. */
  support?: number;
  /** For virtual cells: how many measured compounds went into the prediction at all — the row's plus
   *  the column's real-cell counts, which is exactly what the context panel lists under "Measured with
   *  this core" and "Measured with this substituent". Undefined for real/empty. */
  references?: number;
  /** For real cells: the additive (Free-Wilson) model's fitted value at this cell, so the observed
   *  potency can be compared with what the model expects. A large gap flags a non-additive cell.
   *  Undefined for virtual/empty, and for real cells when predictions weren't computed. */
  fit?: number;
}

/** A matrix row: one core. */
export interface SarMatrixRow {
  coreSmiles: string;
  /** What the row stands for: the core carrying this row's folded substituents, with the column
   *  position still open. Rows of one matrix routinely share a core byte for byte — every member was
   *  decomposed against the same anchor — so the core alone draws every row identically. Falls back to
   *  the core when nothing is folded or the fragments cannot be linked. */
  keySmiles: string;
  label: string;
  /** Substituents fixed by this row, for every position that is not the column axis. Two rows sharing
   *  a core differ only here, so a compound varying at more than one position still has a cell of its
   *  own. Empty for single-position matrices. */
  foldedValues: {[position: string]: string};
}

/** A matrix column: one varying substituent at one R-group position. */
export interface SarMatrixColumn {
  /** R-group position this column belongs to, e.g. 'R1', 'R2'. Single-position fallback
   *  matrices use the constant 'R1' for every column. */
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
  /** Ranking scores keyed by scheme name. */
  scores: {[scheme: string]: number};
  /** Active R-group positions, richest first (e.g. ['R1', 'R2']). Length 1 means this matrix
   *  is a single-position (fallback) view. */
  positions: string[];
  /** Reference (most-frequent-observed) substituent SMILES per position — the value every
   *  other active position is pinned to while one position's column group varies. */
  refValues: {[position: string]: string};
  /** The remainder this matrix's cores were grouped on (`CoreCluster.siteKey`). */
  siteKey: string;
  /** Fragmentation level this matrix was built at (2 = built directly from series). A higher level
   *  holds the same compounds as the matrices below it, over cores that agree one cut deeper — fewer
   *  rows covering more chemistry. */
  level: number;
  /** Id of the matrix one level up, which covers this one's compounds among others. */
  parentId?: string;
  /** Leave-one-out cross-validated quality of the Free-Wilson fit over this matrix's observed cells;
   *  null when there are too few observations to cross-validate. Tells the user how far to trust the
   *  virtual predictions. `n` is how many observed cells were cross-validatable out of `total`. */
  confidence?: {r2: number, rmse: number, n: number, total: number} | null;
}
