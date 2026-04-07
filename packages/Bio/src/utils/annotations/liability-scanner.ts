/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {ISeqHandler} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import {
  SeqAnnotation, SeqAnnotationHit, RowAnnotationData,
  AnnotationVisualType, AnnotationCategory, LiabilitySeverity,
  ANNOTATION_COLORS,
} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';
import {
  getOrCreateAnnotationColumn, setColumnAnnotations, setRowAnnotations,
  getColumnAnnotations, getRowAnnotations, mergeRowHits,
} from './annotation-manager';

/** A single liability scanning rule. */
export interface LiabilityRule {
  id: string;
  name: string;
  pattern: RegExp;
  length: number;
  severity: LiabilitySeverity;
  /** Sub-category for grouping (e.g. "deamidation", "oxidation") */
  ruleCategory: string;
  color: string;
  enabled: boolean;
}

/** Built-in liability rules for antibody engineering. */
export const BUILTIN_LIABILITY_RULES: LiabilityRule[] = [
  {id: 'deamid-ng', name: 'Deamidation (NG)', pattern: /NG/g, length: 2, severity: LiabilitySeverity.High, ruleCategory: 'deamidation', color: ANNOTATION_COLORS.liability.deamidation, enabled: true},
  {id: 'deamid-ns', name: 'Deamidation (NS)', pattern: /NS/g, length: 2, severity: LiabilitySeverity.Medium, ruleCategory: 'deamidation', color: ANNOTATION_COLORS.liability.deamidation, enabled: true},
  {id: 'deamid-na', name: 'Deamidation (NA)', pattern: /NA/g, length: 2, severity: LiabilitySeverity.Low, ruleCategory: 'deamidation', color: ANNOTATION_COLORS.liability.deamidation, enabled: true},
  {id: 'deamid-nd', name: 'Deamidation (ND)', pattern: /ND/g, length: 2, severity: LiabilitySeverity.Low, ruleCategory: 'deamidation', color: ANNOTATION_COLORS.liability.deamidation, enabled: true},
  {id: 'deamid-nt', name: 'Deamidation (NT)', pattern: /NT/g, length: 2, severity: LiabilitySeverity.Low, ruleCategory: 'deamidation', color: ANNOTATION_COLORS.liability.deamidation, enabled: true},
  {id: 'isom-dg', name: 'Isomerization (DG)', pattern: /DG/g, length: 2, severity: LiabilitySeverity.High, ruleCategory: 'isomerization', color: ANNOTATION_COLORS.liability.isomerization, enabled: true},
  {id: 'isom-ds', name: 'Isomerization (DS)', pattern: /DS/g, length: 2, severity: LiabilitySeverity.Medium, ruleCategory: 'isomerization', color: ANNOTATION_COLORS.liability.isomerization, enabled: true},
  {id: 'oxid-m', name: 'Oxidation (Met)', pattern: /M/g, length: 1, severity: LiabilitySeverity.Medium, ruleCategory: 'oxidation', color: ANNOTATION_COLORS.liability.oxidation, enabled: true},
  {id: 'oxid-w', name: 'Oxidation (Trp)', pattern: /W/g, length: 1, severity: LiabilitySeverity.Low, ruleCategory: 'oxidation', color: ANNOTATION_COLORS.liability.oxidation, enabled: true},
  {id: 'glyco-nxst', name: 'N-glycosylation', pattern: /N[^P][ST]/g, length: 3, severity: LiabilitySeverity.High, ruleCategory: 'glycosylation', color: ANNOTATION_COLORS.liability.glycosylation, enabled: true},
  {id: 'free-cys', name: 'Free Cysteine', pattern: /C/g, length: 1, severity: LiabilitySeverity.Info, ruleCategory: 'freeCysteine', color: ANNOTATION_COLORS.liability.freeCysteine, enabled: false},
];

/** Extracts a canonical single-letter string from a sequence handler for a given row. */
function getCanonicalString(sh: ISeqHandler, rowIdx: number): string {
  const splitted = sh.getSplitted(rowIdx);
  const chars: string[] = new Array(splitted.length);
  for (let i = 0; i < splitted.length; i++)
    chars[i] = splitted.getOriginal(i);
  return chars.join('');
}

export interface ScanLiabilitiesResult {
  annotations: SeqAnnotation[];
  rowData: RowAnnotationData[];
  totalHits: number;
}

/** Scans all rows of a macromolecule column for liability motifs.
 *  Returns column-level SeqAnnotation entries + per-row SeqAnnotationHit arrays. */
export function scanLiabilities(
  col: DG.Column<string>,
  sh: ISeqHandler,
  rules: LiabilityRule[],
): ScanLiabilitiesResult {
  const enabledRules = rules.filter((r) => r.enabled);
  const posList = sh.posList;

  // Track which rules had hits
  const ruleHitCounts = new Map<string, number>();

  const rowData: RowAnnotationData[] = new Array(col.length);
  let totalHits = 0;

  for (let rowIdx = 0; rowIdx < col.length; rowIdx++) {
    const seq = getCanonicalString(sh, rowIdx);
    const hits: SeqAnnotationHit[] = [];

    for (const rule of enabledRules) {
      // Reset regex lastIndex for global patterns
      rule.pattern.lastIndex = 0;
      let match: RegExpExecArray | null;
      while ((match = rule.pattern.exec(seq)) !== null) {
        hits.push({
          annotationId: rule.id,
          positionIndex: match.index,
          positionName: match.index < posList.length ? posList[match.index] : undefined,
          matchedMonomers: match[0],
        });
        ruleHitCounts.set(rule.id, (ruleHitCounts.get(rule.id) ?? 0) + 1);
        totalHits++;
      }
    }
    rowData[rowIdx] = hits;
  }

  // Build column-level annotations only for rules that had hits
  const annotations: SeqAnnotation[] = enabledRules
    .filter((r) => ruleHitCounts.has(r.id))
    .map((r) => ({
      id: r.id,
      name: r.name,
      description: `${r.ruleCategory} liability pattern (${ruleHitCounts.get(r.id)} hits)`,
      start: null,
      end: null,
      visualType: r.length === 1 ? AnnotationVisualType.Point : AnnotationVisualType.Motif,
      category: AnnotationCategory.Liability,
      color: r.color,
      severity: r.severity,
      motifPattern: r.pattern.source,
      autoGenerated: true,
    }));

  return {annotations, rowData, totalHits};
}

/** Applies liability scan results to the DataFrame (writes tags + companion column). */
export function applyLiabilityScanResults(
  df: DG.DataFrame,
  seqCol: DG.Column<string>,
  result: ScanLiabilitiesResult,
): void {
  // Merge with existing annotations, removing old liability entries
  const existing = getColumnAnnotations(seqCol)
    .filter((a) => a.category !== AnnotationCategory.Liability);
  setColumnAnnotations(seqCol, [...existing, ...result.annotations]);

  // Write per-row data to hidden companion column, preserving region hits from numbering
  const annotCol = getOrCreateAnnotationColumn(df, seqCol);
  for (let i = 0; i < result.rowData.length; i++) {
    const existingHits = getRowAnnotations(annotCol, i) ?? [];
    setRowAnnotations(annotCol, i, mergeRowHits(existingHits, result.rowData[i], false, true));
  }
}

/** Creates a liability summary count column (total hits per row). */
export function createLiabilitySummaryColumn(
  df: DG.DataFrame,
  seqCol: DG.Column<string>,
  result: ScanLiabilitiesResult,
): DG.Column<number> {
  const colName = `${seqCol.name}_liability_count`;
  const counts = result.rowData.map((hits) => hits.length);
  const col = df.columns.addNewInt(colName);
  for (let i = 0; i < counts.length; i++)
    col.set(i, counts[i]);
  return col;
}
