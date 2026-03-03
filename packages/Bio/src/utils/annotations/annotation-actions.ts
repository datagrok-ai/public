import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {AnnotationCategory} from '@datagrok-libraries/bio/src/utils/macromolecule/annotations';
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

/** Extracts a named region annotation as a new column using ISeqHandler.getRegion(). */
export function extractAnnotatedRegion(
  df: DG.DataFrame,
  seqCol: DG.Column<string>,
  annotationName: string,
  seqHelper: any,
): DG.Column<string> | null {
  const annotations = getColumnAnnotations(seqCol);
  const annot = annotations.find((a) =>
    a.name === annotationName && a.category === AnnotationCategory.Structure);

  if (!annot || annot.start == null || annot.end == null) {
    grok.shell.warning(`Region annotation "${annotationName}" not found or has no position range.`);
    return null;
  }

  const sh = seqHelper.getSeqHandler(seqCol);
  const startIdx = sh.posList.indexOf(annot.start);
  const endIdx = sh.posList.indexOf(annot.end);

  if (startIdx < 0 || endIdx < 0) {
    grok.shell.warning(`Position names "${annot.start}" or "${annot.end}" not found in position list.`);
    return null;
  }

  const colName = `${seqCol.name}(${annotationName})`;
  const regCol = sh.getRegion(startIdx, endIdx, colName);
  df.columns.add(regCol);
  grok.data.detectSemanticTypes(df);
  grok.shell.info(`Extracted region ${annotationName} as column "${colName}"`);
  return regCol;
}
