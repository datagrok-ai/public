import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  META_COLUMNS, PUBLISHED_TAGS, PublishedMetadata,
  isPublished, getPublishedMetadata, buildMailtoUrl,
} from '../publishing/publish-state';
import {findHostDataFrameForProtein} from './uniprot-panel';
import {getGroups} from '../analysis/experiment-setup';

/**
 * Reviewer-side audit context panel. Opens when a reviewer clicks a protein
 * row in a published Project. Surfaces the metadata the analyst stamped at
 * publish time (DE method, thresholds, group names, target, share date, sharer
 * name) plus PUB-13 mailto button. First-line `isPublished(df)` guard ensures
 * it does NOT render on the analyst's live working DataFrame.
 *
 * Reads metadata column FIRST, tag SECOND per Pitfall 3 — the column survives
 * the serializer; the tag is a defensive secondary. All user-facing strings
 * pass the Pitfall 14 biologist-jargon audit. Imports {@link buildMailtoUrl}
 * from `../publishing/publish-state` rather than `../publishing/share-dialog`
 * (B-1 fix — avoids the wave-2 → wave-4 dependency cycle).
 */

/** Human-friendly expansion of the DE method ID for the audit "Method" row.
 *  CONTEXT.md domain pin: narrate what an expert would assume known. */
function methodLabel(method: string): string {
  switch (method) {
    case 'limma': return 'limma — moderated t-test';
    case 'deqms': return 'DEqMS — peptide-count-aware moderated t-test';
    case 't-test': return 'Welch t-test';
    case 'spectronaut': return 'Spectronaut Candidates (pre-computed DE)';
    default: return method || '(unknown method)';
  }
}

function dateSlice(d: Date | string | number | null | undefined): string {
  if (d == null) return 'unknown date';
  try {
    if (d instanceof Date) {
      if (Number.isNaN(d.getTime())) return 'unknown date';
      return d.toISOString().slice(0, 10);
    }
    const parsed = new Date(d);
    if (Number.isNaN(parsed.getTime())) return 'unknown date';
    return parsed.toISOString().slice(0, 10);
  } catch {
    return 'unknown date';
  }
}

function fieldRow(label: string, value: string): HTMLDivElement {
  const labelEl = ui.divText(label);
  labelEl.style.fontWeight = '600';
  labelEl.style.fontSize = '12px';
  labelEl.style.color = '#666';
  const valueEl = ui.divText(value);
  valueEl.style.marginBottom = '8px';
  const row = ui.div([labelEl, valueEl]);
  row.style.marginBottom = '6px';
  return row;
}

function readSharerEmail(df: DG.DataFrame, meta: PublishedMetadata): string | null {
  if (meta.publishedByEmail != null && meta.publishedByEmail !== '')
    return meta.publishedByEmail;
  try {
    const col = df.col(META_COLUMNS.PUBLISHED_BY_EMAIL);
    if (col != null) {
      const v = col.get(0);
      if (v != null && v !== '') return String(v);
    }
  } catch { /* fall through */ }
  try {
    const t = df.getTag(PUBLISHED_TAGS.PUBLISHED_BY_EMAIL);
    if (t) return t;
  } catch { /* fall through */ }
  return null;
}

export function publishedAnalysisPanel(proteinId: string): DG.Widget {
  const df = findHostDataFrameForProtein(proteinId);
  if (df == null) return new DG.Widget(ui.div());
  if (!isPublished(df)) return new DG.Widget(ui.div());

  const meta = getPublishedMetadata(df);
  if (meta == null)
    return new DG.Widget(ui.divText('Shared analysis metadata unavailable. Ask the sharer to re-share.'));

  const body = ui.div();
  body.style.padding = '8px 4px';

  if (meta.supersededBy != null && meta.supersededBy !== '') {
    const newerLink = ui.link('Newer version available — open it', async () => {
      try {
        const newer = await grok.dapi.projects.find(meta.supersededBy!);
        await newer.open();
      } catch (e) {
        grok.shell.warning(`Could not open the newer version: ${(e as Error)?.message ?? e}`);
      }
    });
    const banner = ui.div([newerLink]);
    banner.style.padding = '8px';
    banner.style.marginBottom = '10px';
    banner.style.backgroundColor = '#fff8dc';
    banner.style.borderRadius = '4px';
    body.appendChild(banner);
  }

  body.appendChild(ui.h2('Shared analysis details'));

  const groups = getGroups(df);
  const comparison = groups != null
    ? `${groups.group2.name} vs ${groups.group1.name}`
    : 'unavailable';

  const fcFold = Number.isFinite(meta.fcThreshold) ? Math.pow(2, meta.fcThreshold) : NaN;
  const fcLabel = Number.isFinite(meta.fcThreshold)
    ? `log2(${meta.fcThreshold}) (${Number.isFinite(fcFold) ? fcFold.toFixed(2) : 'NaN'}-fold)`
    : '(unknown)';
  const pLabel = Number.isFinite(meta.pThreshold)
    ? `${meta.pThreshold} (FDR-adjusted)`
    : '(unknown)';

  body.appendChild(fieldRow('Target', meta.target || '(unknown)'));
  body.appendChild(fieldRow('Shared', dateSlice(meta.publishedAt)));
  body.appendChild(fieldRow('Shared by', meta.publishedBy || '(unknown)'));
  body.appendChild(fieldRow('Method', methodLabel(meta.deMethod)));
  body.appendChild(fieldRow('Comparison', comparison));
  body.appendChild(fieldRow('Fold-change cutoff', fcLabel));
  body.appendChild(fieldRow('p-value cutoff', pLabel));

  const sharerEmail = readSharerEmail(df, meta);
  if (sharerEmail != null) {
    const projectName = df.name || 'shared analysis';
    const mailtoUrl = buildMailtoUrl({
      sharerEmail,
      sharerName: meta.publishedBy || 'colleague',
      projectName,
      publishedDateStr: dateSlice(meta.publishedAt),
    });
    const mailLink = ui.link('Request re-run with different parameters', mailtoUrl);
    (mailLink as HTMLAnchorElement).href = mailtoUrl;
    const mailRow = ui.div([mailLink]);
    mailRow.style.marginTop = '12px';
    body.appendChild(mailRow);
  }

  return new DG.Widget(body);
}
