import * as DG from 'datagrok-api/dg';
import {SEMTYPE} from './proteomics-types';

/** Finds a column by semantic type, falling back to column name heuristics.
 * Lowercases each hint so callers passing mixed-case strings still match. */
export function findColumn(df: DG.DataFrame, semType: string, nameHints: string[]): DG.Column | null {
  const bySemType = df.columns.toList().find((c) => c.semType === semType);
  if (bySemType)
    return bySemType;

  for (const hint of nameHints) {
    const hintLower = hint.toLowerCase();
    const byName = df.columns.toList().find((c) => c.name.toLowerCase().includes(hintLower));
    if (byName)
      return byName;
  }
  return null;
}

/** Finds common proteomics columns in a dataframe. */
export function findProteomicsColumns(df: DG.DataFrame) {
  return {
    proteinId: findColumn(df, SEMTYPE.PROTEIN_ID, ['protein id', 'majority protein id', 'accession', 'uniprot']),
    geneName: findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol']),
    log2fc: findColumn(df, SEMTYPE.LOG2FC, ['log2fc', 'log2 fold', 'logfc']),
    pValue: findColumn(df, SEMTYPE.P_VALUE, ['p-value', 'pvalue', 'adj.p', 'fdr', 'q-value']),
  };
}
