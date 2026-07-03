import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PUBLISHED_TAGS, isPublished} from './publish-state';
import {applyVolcanoFormulaLines} from './publish-project';
import {wireEnrichmentToVolcano} from '../viewers/enrichment-viewers';

const ENRICHMENT_WIRED_TAG = 'proteomics.enrichment_wired';

/**
 * Post-open recovery hook for reopened published Projects. Idempotent — safe
 * to call multiple times on the same DataFrame.
 *
 * Self-heals two failure modes:
 *
 *  - **B-2 / W-7 — volcano formula lines stripped by serializer.** Reads FC
 *    and p threshold values from the `proteomics.published_fc_threshold` /
 *    `proteomics.published_p_threshold` tags (Plan 02's belt-and-braces #2
 *    write via `setPublishedTags`) and re-applies via the exported
 *    {@link applyVolcanoFormulaLines}. No-op when formula lines are already
 *    present.
 *
 *  - **W-5 — enrichment cross-DF highlight subscription not re-established
 *    on reopen.** When the published Project carries both the protein DF
 *    AND a published enrichment DF (both tagged `proteomics.published=true`),
 *    wires the existing `wireEnrichmentToVolcano` subscription pattern from
 *    `src/viewers/enrichment-viewers.ts`. Sentinel-tag duplicate-guard
 *    (`proteomics.enrichment_wired`) prevents double subscriptions when the
 *    user reopens the project repeatedly in one session.
 */
export async function recoverPublishedProject(df: DG.DataFrame): Promise<void> {
  if (!isPublished(df)) return;

  // PART 1 — Re-apply formula lines on the volcano (B-2 / W-7)
  const tv = findTableViewFor(df);
  if (tv) {
    const volcano = Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.SCATTER_PLOT);
    if (volcano) {
      const fcRaw = df.getTag(PUBLISHED_TAGS.PUBLISHED_FC_THRESHOLD);
      const pRaw = df.getTag(PUBLISHED_TAGS.PUBLISHED_P_THRESHOLD);
      if (fcRaw && pRaw) {
        const fc = parseFloat(fcRaw);
        const p = parseFloat(pRaw);
        if (Number.isFinite(fc) && Number.isFinite(p)) {
          if (!hasFormulaLines(df, volcano)) {
            try { applyVolcanoFormulaLines(volcano, fc, p); } catch { /* best-effort */ }
          }
        }
      }
    }
  }

  // PART 2 — Re-establish enrichment cross-DF subscription (W-5)
  const enrichDf = findPublishedEnrichmentSibling(df);
  if (enrichDf && enrichDf.getTag(ENRICHMENT_WIRED_TAG) !== 'true') {
    try {
      wireEnrichmentToVolcano(enrichDf, df);
      enrichDf.setTag(ENRICHMENT_WIRED_TAG, 'true');
    } catch { /* best-effort */ }
  }
}

function findTableViewFor(df: DG.DataFrame): DG.TableView | null {
  try {
    const views = (grok.shell as any).tableViews as DG.TableView[] | undefined;
    if (Array.isArray(views)) {
      const match = views.find((v) => v?.dataFrame === df);
      if (match) return match;
    }
  } catch { /* fall through */ }
  try {
    const tv = grok.shell.tv;
    if (tv && tv.dataFrame === df) return tv;
  } catch { /* fall through */ }
  return null;
}

function hasFormulaLines(df: DG.DataFrame, viewer: DG.Viewer): boolean {
  try {
    const items = ((df as any).meta?.formulaLines?.items ?? []) as Array<{formula?: string}>;
    if (Array.isArray(items) && items.length > 0) return true;
  } catch { /* fall through */ }
  try {
    const opts: any = (viewer as any).getOptions?.();
    const lookLines = opts?.look?.formulaLines;
    if (typeof lookLines === 'string' && lookLines !== '' && lookLines !== '[]') return true;
    if (Array.isArray(lookLines) && lookLines.length > 0) return true;
  } catch { /* fall through */ }
  return false;
}

function findPublishedEnrichmentSibling(proteinDf: DG.DataFrame): DG.DataFrame | null {
  try {
    const tables = (grok.shell as any).tables as DG.DataFrame[] | undefined;
    if (!Array.isArray(tables)) return null;
    return tables.find((t) =>
      t !== proteinDf &&
      t?.getTag?.('proteomics.enrichment') === 'true' &&
      t?.getTag?.(PUBLISHED_TAGS.PUBLISHED) === 'true',
    ) ?? null;
  } catch {
    return null;
  }
}
