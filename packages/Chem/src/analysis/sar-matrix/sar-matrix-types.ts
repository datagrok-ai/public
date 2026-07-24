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
  /** Fragment id of the varying substituent. */
  substFragId: number;
  /** SMILES of the varying substituent. */
  substSmiles: string;
}

/** A matched series: molecules sharing the same MMP core, differing at one site. */
export interface MatchedSeries {
  /** Fragment id of the shared core. */
  coreFragId: number;
  /** SMILES of the shared core. */
  coreSmiles: string;
  members: SeriesMember[];
}

/** A cluster of structurally related matched series — the rows of one SAR Matrix. */
export interface CoreCluster {
  id: string;
  /** One series per related core; each becomes a matrix row. */
  series: MatchedSeries[];
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
  /** For real cells: the additive (Free-Wilson) model's fitted value at this cell, so the observed
   *  potency can be compared with what the model expects. A large gap flags a non-additive cell.
   *  Undefined for virtual/empty, and for real cells when predictions weren't computed. */
  fit?: number;
}

/** A matrix row: one core. */
export interface SarMatrixRow {
  coreFragId: number;
  coreSmiles: string;
  label: string;
}

/** A matrix column: one varying substituent at one R-group position. */
export interface SarMatrixColumn {
  /** R-group position this column belongs to, e.g. 'R1', 'R2'. Single-position fallback
   *  matrices use the constant 'R1' for every column. */
  position: string;
  substSmiles: string;
  /** Number of real compounds observed in this column. */
  count: number;
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
  /** Leave-one-out cross-validated quality of the Free-Wilson fit over this matrix's observed cells;
   *  null when there are too few observations to cross-validate. Tells the user how far to trust the
   *  virtual predictions. */
  confidence?: {r2: number, rmse: number, n: number} | null;
}
