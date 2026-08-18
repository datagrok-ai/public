import * as grok from 'datagrok-api/grok';
import {AlignmentRow, HistoryRow, T_ALIGNMENT, T_HISTORY} from '../domain/constants';

export interface DriftFinding {
  kind: 'multiple-in-review' | 'no-live-row' | 'interrupted-flip' | 'missing-frozen-artifact' |
    'table-grant-missing';
  publicationId: string;
  detail: string;
}

/** Scans the catalog invariants the server cannot enforce by itself:
 * - at most one in-review (pending/rejected) live row per publication — the business
 *   key enforces per-status singletons, but pending + rejected may coexist;
 * - a publication with history but no live row (interrupted republish);
 * - a publication whose only live row is in-review while an approved version exists
 *   in history newer than any approval (interrupted approve flip);
 * - frozen artifact funccalls that no longer resolve (`deep` only). */
export async function driftCheck(options?: {deep?: boolean}): Promise<DriftFinding[]> {
  const findings: DriftFinding[] = [];
  const live = await grok.dapi.domains.table<AlignmentRow>(T_ALIGNMENT).query({limit: 10000});
  const hist = await grok.dapi.domains.table<HistoryRow>(T_HISTORY)
    .query({columns: ['publication_id', 'revision', 'final_status'] as any, limit: 10000});

  const liveByPub = new Map<string, AlignmentRow[]>();
  for (const row of live)
    liveByPub.set(row.publication_id, [...liveByPub.get(row.publication_id) ?? [], row]);

  for (const [pub, rows] of liveByPub) {
    const inReview = rows.filter((r) => r.status !== 'approved');
    if (inReview.length > 1) {
      findings.push({kind: 'multiple-in-review', publicationId: pub,
        detail: `statuses: ${rows.map((r) => `v${r.revision}:${r.status}`).join(', ')}`});
    }
  }

  const histPubs = new Set(hist.map((h) => h.publication_id));
  for (const pub of histPubs) {
    if (!liveByPub.has(pub)) {
      findings.push({kind: 'no-live-row', publicationId: pub,
        detail: 'every version is archived; the last archived revision should be ' +
          're-inserted or the publication retired deliberately'});
      continue;
    }
    const rows = liveByPub.get(pub)!;
    if (rows.every((r) => r.status !== 'approved') &&
        hist.some((h) => h.publication_id === pub && h.final_status === 'approved')) {
      findings.push({kind: 'interrupted-flip', publicationId: pub,
        detail: 'no approved live row, but an approved version is archived — ' +
          'an approve flip may have been interrupted'});
    }
  }

  if (options?.deep === true) {
    for (const row of live) {
      try {
        await grok.dapi.functions.calls.allPackageVersions().find(row.artifact_id);
      } catch (_) {
        findings.push({kind: 'missing-frozen-artifact', publicationId: row.publication_id,
          detail: `frozen meta call ${row.artifact_id} of v${row.revision} does not resolve`});
      }
    }
  }
  return findings;
}
