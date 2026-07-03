export const SEMTYPE = {
  PROTEIN_ID: 'Proteomics-ProteinId',
  GENE_SYMBOL: 'Proteomics-GeneSymbol',
  LOG2FC: 'Proteomics-Log2FC',
  P_VALUE: 'Proteomics-PValue',
  INTENSITY: 'Proteomics-Intensity',
  SUBCELLULAR_LOCATION: 'Proteomics-SubcellularLocation',
  DISPLAY_NAME: 'Proteomics-DisplayName',
  SOURCE_ID: 'Proteomics-SourceId',
  NUMERATOR_MEAN: 'Proteomics-NumeratorMean',
  DENOMINATOR_MEAN: 'Proteomics-DenominatorMean',
} as const;

/** Default significance cutoffs — the SINGLE source for the |log2FC| ≥ 1 and
 * adj.p ≤ 0.05 thresholds used as dialog defaults, function-parameter defaults,
 * and `?? ` fallbacks across DE, enrichment, the volcano, and the demos. Never
 * re-inline the literals; import these so the pipeline stays internally
 * consistent when a default is retuned. */
export const DEFAULT_FC_THRESHOLD = 1.0;
export const DEFAULT_P_THRESHOLD = 0.05;

/** Locked direction palette (D-04): numerator/group1 = magenta, denominator/
 * group2 = cyan, not-significant = gray, as ARGB ints. Single source shared by
 * the volcano direction column and the UniProt panel's per-group bar chart —
 * lives here (a leaf module) so both can import without a circular dependency. */
export const DIRECTION_COLORS_BASE = {
  enrichedG1: 0xFFFF00FF, // magenta (ARGB)
  enrichedG2: 0xFF00FFFF, // cyan (ARGB)
  notSig: 0xFFAAAAAA, // gray (ARGB)
} as const;
