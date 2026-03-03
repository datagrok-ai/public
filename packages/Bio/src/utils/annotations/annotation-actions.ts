import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {AnnotationCategory, SeqAnnotationHit} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';
import {ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getAnnotationColumnName, getColumnAnnotations, cacheAllRowAnnotations} from './annotation-manager';

/** Filters the DataFrame to show only rows that have at least one liability hit. */
export function filterByLiabilityHits(df: DG.DataFrame, seqCol: DG.Column<string>): void {
  const annotColName = getAnnotationColumnName(seqCol.name);
  let annotCol: DG.Column<string> | null = null;
  try {
    annotCol = df.columns.byName(annotColName) as DG.Column<string>;
  } catch { /* not found */ }

  if (!annotCol) {
    grok.shell.warning('No annotation data found. Run liability scanning first.');
    return;
  }

  const rowData = cacheAllRowAnnotations(annotCol);
  const bs = DG.BitSet.create(df.rowCount);
  for (let i = 0; i < df.rowCount; i++) {
    if (rowData[i] && rowData[i]!.length > 0)
      bs.set(i, true);
  }

  df.filter.copyFrom(bs);
  grok.shell.info(`Filtered to ${bs.trueCount} rows with liability hits`);
}

/** Selects all rows that contain a specific annotation hit. */
export function selectRowsWithAnnotation(df: DG.DataFrame, seqCol: DG.Column<string>, annotationId: string): void {
  const annotColName = getAnnotationColumnName(seqCol.name);
  let annotCol: DG.Column<string> | null = null;
  try {
    annotCol = df.columns.byName(annotColName) as DG.Column<string>;
  } catch { /* not found */ }

  if (!annotCol) {
    grok.shell.warning('No annotation data found.');
    return;
  }

  const rowData = cacheAllRowAnnotations(annotCol);
  const bs = DG.BitSet.create(df.rowCount);
  for (let i = 0; i < df.rowCount; i++) {
    if (rowData[i]?.some((h) => h.annotationId === annotationId))
      bs.set(i, true);
  }

  df.selection.copyFrom(bs);
  grok.shell.info(`Selected ${bs.trueCount} rows with ${annotationId} hits`);
}

/** Extracts a named region annotation as a new column.
 *  Uses per-row region spans from the companion column when available (unaligned data),
 *  falls back to column-level position names (aligned/MSA data). */
export function extractAnnotatedRegion(
  df: DG.DataFrame,
  seqCol: DG.Column<string>,
  annotationName: string,
  seqHelper: ISeqHelper,
): DG.Column<string> | null {
  const annotations = getColumnAnnotations(seqCol);
  const annot = annotations.find((a) =>
    a.name === annotationName && a.category === AnnotationCategory.Structure);

  if (!annot) {
    grok.shell.warning(`Region annotation "${annotationName}" not found.`);
    return null;
  }

  const sh = seqHelper.getSeqHandler(seqCol);
  const colName = `${seqCol.name}(${annotationName})`;

  // Try per-row extraction using companion column region spans
  const annotColName = getAnnotationColumnName(seqCol.name);
  let annotCol: DG.Column<string> | null = null;
  try { annotCol = df.columns.byName(annotColName) as DG.Column<string>; } catch { /* not found */ }

  if (annotCol) {
    const allRowData = cacheAllRowAnnotations(annotCol);
    const hasPerRowRegions = allRowData.some((rd) =>
      rd?.some((h: SeqAnnotationHit) => h.annotationId === annot.id && h.endPositionIndex != null));

    if (hasPerRowRegions) {
      const regCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, colName, df.rowCount);
      for (let i = 0; i < df.rowCount; i++) {
        const rowHits = allRowData[i];
        const regionHit = rowHits?.find((h: SeqAnnotationHit) =>
          h.annotationId === annot.id && h.endPositionIndex != null);
        if (regionHit) {
          const splitted = sh.getSplitted(i);
          const parts: string[] = [];
          for (let p = regionHit.positionIndex; p <= regionHit.endPositionIndex!; p++) {
            if (p < splitted.length)
              parts.push(splitted.getOriginal(p));
          }
          regCol.set(i, parts.join(sh.separator || ''));
        } else
          regCol.set(i, '');
      }
      df.columns.add(regCol);
      grok.data.detectSemanticTypes(df);
      grok.shell.info(`Extracted region ${annotationName} as column "${colName}"`);
      return regCol;
    }
  }

  // Fall back to column-level position names (aligned/MSA data)
  if (annot.start == null || annot.end == null) {
    grok.shell.warning(`Region annotation "${annotationName}" has no position range.`);
    return null;
  }

  const startIdx = sh.posList.indexOf(annot.start);
  const endIdx = sh.posList.indexOf(annot.end);

  if (startIdx < 0 || endIdx < 0) {
    grok.shell.warning(`Position names "${annot.start}" or "${annot.end}" not found in position list.`);
    return null;
  }

  const regCol = sh.getRegion(startIdx, endIdx, colName);
  df.columns.add(regCol);
  grok.data.detectSemanticTypes(df);
  grok.shell.info(`Extracted region ${annotationName} as column "${colName}"`);
  return regCol;
}
