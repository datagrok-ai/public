import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from '../utils/proteomics-types';
import {addPrimaryColumnIfNeeded, detectDelimiter} from './shared-utils';
import {resolveGeneLabels} from '../utils/gene-label-resolver';
import {setGroups} from '../analysis/experiment-setup';

/** Column name variants for protein-level identifiers in Spectronaut Candidates output.
 * Modern Spectronaut prefixes with `PG.`; older / re-exported reports may strip it. */
const PROTEIN_ID_COLUMNS = ['PG.ProteinGroups', 'ProteinGroups', 'Protein Groups'];

/** Column name variants for gene symbols. */
const GENE_NAME_COLUMNS = ['PG.Genes', 'Genes', 'Gene Names'];

/** Column name variants for log2 fold change. Spectronaut's canonical name is
 * `AVG Log2 Ratio`; some export profiles abbreviate or rename it. */
const LOG2FC_COLUMNS = ['AVG Log2 Ratio', 'Log2 Ratio', 'AVG Log2 Fold Change'];

/** Column name variants for the FDR-adjusted significance metric. Spectronaut emits
 * `Qvalue` from its built-in Storey FDR; some export profiles use `Q-Value` or
 * `AdjustedPValue`. */
const ADJ_P_COLUMNS = ['Qvalue', 'Q-Value', 'Q.Value', 'AdjustedPValue', 'Adjusted P-Value'];

/** Column name variants for the raw p-value (optional in Candidates reports). */
const P_VALUE_COLUMNS = ['Pvalue', 'PValue', 'P-Value', 'P.Value'];

/** The declared contrast string, kept verbatim. Form: "group1 / group2". */
const COMPARISON_COLUMNS = ['Comparison (group1/group2)', 'Comparison'];

/** Per-group mean abundance, swapped when a row's sign is flipped (A6: only
 * when BOTH are present). */
const AVG_GROUP_QTY_NUM_COLUMNS = ['AVG Group Quantity Numerator'];
const AVG_GROUP_QTY_DEN_COLUMNS = ['AVG Group Quantity Denominator'];

/** Declared numerator/denominator condition labels, relabeled on flip. */
const CONDITION_NUM_COLUMNS = ['Condition Numerator'];
const CONDITION_DEN_COLUMNS = ['Condition Denominator'];

/** Decoy / contaminant prefixes to filter at the row level. Mirrors the FragPipe
 * parser — Spectronaut libraries may include MaxQuant-style CON__/REV__ entries
 * when the FASTA was re-used across pipelines. */
const CONTAMINANT_PREFIXES = ['contam_', 'rev_', 'CON__', 'REV__'];

/** Default significance thresholds when the input does not pre-compute a boolean
 * `significant` column. Matches the defaults in `runDifferentialExpression`. */
const DEFAULT_FC_THRESHOLD = 1.0;
const DEFAULT_P_THRESHOLD = 0.05;

/** Finds the first existing column from a list of name candidates. */
function findCol(df: DG.DataFrame, names: readonly string[]): DG.Column | null {
  for (const name of names) {
    const col = df.col(name);
    if (col)
      return col;
  }
  return null;
}

/** Marks rows for exclusion when the protein-group identifier starts with any
 * contaminant/decoy prefix. */
function filterContaminantRows(df: DG.DataFrame, idColName: string): void {
  const col = df.col(idColName);
  if (!col) return;
  for (let i = 0; i < df.rowCount; i++) {
    const val = col.get(i);
    if (typeof val !== 'string') continue;
    for (const prefix of CONTAMINANT_PREFIXES) {
      if (val.startsWith(prefix)) {
        df.filter.set(i, false, false);
        break;
      }
    }
  }
}

/** Renames a column to the canonical downstream name. Skips the rename if the
 * canonical name is already taken by a different column — defensive against
 * malformed inputs that already contain both `Qvalue` and `adj.p-value`. */
function renameToCanonical(df: DG.DataFrame, src: DG.Column, canonical: string): void {
  if (src.name === canonical) return;
  const existing = df.col(canonical);
  if (existing && existing !== src) return;
  src.name = canonical;
}

/**
 * Normalizes each row's sign so positive log2FC = enriched in the canonical
 * group1. Ported from CK-omics `create_subset_data` (CKomics_tool2.py:1625-1662)
 * and locked by 13-WAVE0-FINDINGS.md A2: the canonical orientation is the first
 * parseable declared "g1 / g2"; ONLY rows whose declared comparison is the exact
 * reverse ("g2 / g1") are flipped — never unconditional. A flip negates log2FC,
 * swaps AVG Group Quantity Numerator/Denominator (only if BOTH present — A6),
 * and relabels Comparison + Condition Numerator/Denominator to canonical.
 * Pure (no shell/UI side effects). Single-orientation / unparseable input is
 * left byte-identical (Pitfall 4 — never invert on ambiguity).
 */
function normalizeCandidatesSign(df: DG.DataFrame): void {
  const cmpCol = findCol(df, COMPARISON_COLUMNS);
  if (!cmpCol) return;
  const n = df.rowCount;

  // Canonical orientation = first parseable "g1 / g2".
  let canonical: string | null = null;
  let canonG1 = '';
  let canonG2 = '';
  for (let i = 0; i < n; i++) {
    const v = cmpCol.get(i);
    if (typeof v === 'string' && v.includes(' / ')) {
      const parts = v.split(' / ');
      if (parts.length === 2 && parts[0] && parts[1]) {
        canonical = v;
        canonG1 = parts[0];
        canonG2 = parts[1];
        break;
      }
    }
  }
  if (canonical === null) return; // no parseable comparison → all canonical
  const reversed = `${canonG2} / ${canonG1}`;

  // Flip ONLY rows whose declared comparison is the exact reverse of canonical.
  const flip = new Uint8Array(n);
  let anyFlip = false;
  for (let i = 0; i < n; i++) {
    if (cmpCol.get(i) === reversed) {
      flip[i] = 1;
      anyFlip = true;
    }
  }
  if (!anyFlip) return; // single-orientation / no reversed rows → unchanged

  // log2FC sign flip. Bulk getRawData → number[] → init (memory
  // feedback_dg_column_bulk_init). FLOAT_NULL is preserved by writing the
  // sentinel number back, never returning JS null (memory
  // feedback_dg_column_init_null_sentinel). number[] (not Float32Array) keeps
  // full double precision for large AVG Group Quantity values.
  const fcCol = df.col('log2FC')!;
  const fcRaw = fcCol.getRawData() as Float32Array | Float64Array;
  const newFc: number[] = new Array(n);
  for (let i = 0; i < n; i++) {
    const v = fcRaw[i];
    newFc[i] = (v === DG.FLOAT_NULL) ? DG.FLOAT_NULL : (flip[i] ? -v : v);
  }
  fcCol.init((i) => newFc[i]);

  // AVG Group Quantity swap — only when BOTH columns exist (A6 guard).
  const numCol = findCol(df, AVG_GROUP_QTY_NUM_COLUMNS);
  const denCol = findCol(df, AVG_GROUP_QTY_DEN_COLUMNS);
  if (numCol && denCol) {
    const numRaw = numCol.getRawData() as Float32Array | Float64Array;
    const denRaw = denCol.getRawData() as Float32Array | Float64Array;
    const newNum: number[] = new Array(n);
    const newDen: number[] = new Array(n);
    for (let i = 0; i < n; i++) {
      const nv = numRaw[i];
      const dv = denRaw[i];
      newNum[i] = flip[i] ? dv : nv;
      newDen[i] = flip[i] ? nv : dv;
    }
    numCol.init((i) => newNum[i]);
    denCol.init((i) => newDen[i]);
  }

  // Relabel Comparison + Condition for flipped rows (string init, not per-row
  // set — memory feedback_dg_column_bulk_init).
  const cmpVals: string[] = new Array(n);
  for (let i = 0; i < n; i++)
    cmpVals[i] = flip[i] ? canonical! : (cmpCol.get(i) as string);
  cmpCol.init((i) => cmpVals[i]);

  const condNumCol = findCol(df, CONDITION_NUM_COLUMNS);
  const condDenCol = findCol(df, CONDITION_DEN_COLUMNS);
  if (condNumCol && condDenCol) {
    const cnVals: string[] = new Array(n);
    const cdVals: string[] = new Array(n);
    for (let i = 0; i < n; i++) {
      cnVals[i] = flip[i] ? canonG1 : (condNumCol.get(i) as string);
      cdVals[i] = flip[i] ? canonG2 : (condDenCol.get(i) as string);
    }
    condNumCol.init((i) => cnVals[i]);
    condDenCol.init((i) => cdVals[i]);
  }
}

/** Parses a Spectronaut Candidates report into a DataFrame shaped like the output
 * of `runDifferentialExpression` — pre-computed `log2FC`, `p-value` (optional),
 * `adj.p-value`, and `significant` columns. Sets `proteomics.de_complete` so the
 * downstream volcano/heatmap/enrichment viewers light up without running DE.
 *
 * Required input columns (any spelling):
 *   - Protein group: `PG.ProteinGroups` / `ProteinGroups`
 *   - log2FC:        `AVG Log2 Ratio` / `Log2 Ratio` / `AVG Log2 Fold Change`
 *   - Adj. p-value:  `Qvalue` / `Q-Value` / `AdjustedPValue` (Spectronaut's
 *                    Storey-FDR output is the canonical Qvalue)
 *
 * Optional:
 *   - Gene symbol:   `PG.Genes` / `Genes` / `Gene Names`
 *   - Raw p-value:   `Pvalue` / `PValue` / `P-Value`
 *   - `Comparison` / `Comparison (group1/group2)` — kept as-is for reference.
 *
 * Notes:
 *   - Multi-comparison files (one row per protein × per comparison) are accepted
 *     verbatim; users filter the resulting DataFrame by the `Comparison` column.
 *   - Significance is computed from |log2FC| ≥ 1.0 and adj.p-value ≤ 0.05 by
 *     default, matching `runDifferentialExpression`'s defaults. */
export async function parseSpectronautCandidatesText(text: string): Promise<DG.DataFrame> {
  const delimiter = detectDelimiter(text);
  const raw = DG.DataFrame.fromCsv(text, {delimiter});

  const idCol = findCol(raw, PROTEIN_ID_COLUMNS);
  if (!idCol)
    throw new Error(`Spectronaut Candidates: missing protein-group column ` +
      `(expected one of ${PROTEIN_ID_COLUMNS.join(', ')})`);

  if (!findCol(raw, LOG2FC_COLUMNS))
    throw new Error(`Spectronaut Candidates: missing log2 ratio column ` +
      `(expected one of ${LOG2FC_COLUMNS.join(', ')})`);

  if (!findCol(raw, ADJ_P_COLUMNS))
    throw new Error(`Spectronaut Candidates: missing q-value column ` +
      `(expected one of ${ADJ_P_COLUMNS.join(', ')})`);

  raw.filter.init((_i) => true);
  filterContaminantRows(raw, idCol.name);
  raw.filter.fireChanged();

  const df = raw.clone(raw.filter);

  // Re-resolve columns on the cloned frame (clone preserves names but not refs)
  const idCol2 = findCol(df, PROTEIN_ID_COLUMNS)!;
  const log2fcCol = findCol(df, LOG2FC_COLUMNS)!;
  const adjPCol = findCol(df, ADJ_P_COLUMNS)!;
  const geneCol = findCol(df, GENE_NAME_COLUMNS);
  const pValCol = findCol(df, P_VALUE_COLUMNS);

  renameToCanonical(df, log2fcCol, 'log2FC');
  renameToCanonical(df, adjPCol, 'adj.p-value');
  if (pValCol)
    renameToCanonical(df, pValCol, 'p-value');

  idCol2.semType = SEMTYPE.PROTEIN_ID;
  if (geneCol) geneCol.semType = SEMTYPE.GENE_SYMBOL;
  df.col('log2FC')!.semType = SEMTYPE.LOG2FC;
  df.col('adj.p-value')!.semType = SEMTYPE.P_VALUE;
  if (df.col('p-value'))
    df.col('p-value')!.semType = SEMTYPE.P_VALUE;

  addPrimaryColumnIfNeeded(df, idCol2.name, 'Primary Protein ID', SEMTYPE.PROTEIN_ID);
  if (geneCol)
    addPrimaryColumnIfNeeded(df, geneCol.name, 'Primary Gene Name', SEMTYPE.GENE_SYMBOL);

  // R3/D-08: per-row sign normalization BEFORE significance (|log2FC| is
  // unchanged by a sign flip, so significance is orientation-invariant).
  normalizeCandidatesSign(df);

  const sigCol = df.columns.addNewBool('significant');
  const fcRaw = df.col('log2FC')!.getRawData() as Float32Array | Float64Array;
  const adjPRaw = df.col('adj.p-value')!.getRawData() as Float32Array | Float64Array;
  sigCol.init((i) => {
    const fc = fcRaw[i];
    const adjP = adjPRaw[i];
    if (fc === DG.FLOAT_NULL || adjP === DG.FLOAT_NULL) return false;
    return Math.abs(fc) >= DEFAULT_FC_THRESHOLD && adjP <= DEFAULT_P_THRESHOLD;
  });

  df.setTag('proteomics.source', 'spectronaut-candidates');
  df.setTag('proteomics.de_complete', 'true');
  df.setTag('proteomics.de_method', 'spectronaut');

  // Name the contrast groups (e.g. DMT / WT) from the report's Condition
  // Numerator / Denominator (or the "num / den" Comparison string) so the
  // volcano legend + title show real names instead of generic group1/group2.
  // Candidates carries no per-sample intensities, so the column lists stay
  // empty — downstream Annotate/Normalize/Impute/DE are skipped via
  // de_complete, and the volcano reads only the names. group1 = Numerator
  // (positive log2FC, magenta), group2 = Denominator (cyan).
  let g1Name = firstNonEmpty(findCol(df, CONDITION_NUM_COLUMNS));
  let g2Name = firstNonEmpty(findCol(df, CONDITION_DEN_COLUMNS));
  if (!g1Name || !g2Name) {
    const cmp = firstNonEmpty(findCol(df, COMPARISON_COLUMNS));
    if (cmp && cmp.includes('/')) {
      const [a, b] = cmp.split('/').map((s) => s.trim());
      g1Name = g1Name ?? (a || null);
      g2Name = g2Name ?? (b || null);
    }
  }
  if (g1Name && g2Name) {
    setGroups(df, {
      group1: {name: g1Name, columns: []},
      group2: {name: g2Name, columns: []},
    });
  }

  await resolveGeneLabels(df);

  return df;
}

/** First non-empty trimmed string value in a column, or null. */
function firstNonEmpty(col: DG.Column | null): string | null {
  if (!col) return null;
  for (let i = 0; i < col.length; i++) {
    const v = col.get(i);
    if (v != null && String(v).trim().length > 0) return String(v).trim();
  }
  return null;
}
