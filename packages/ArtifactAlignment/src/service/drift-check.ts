import * as grok from 'datagrok-api/grok';
import {AlignmentRow, T_ALIGNMENT} from '../domain/constants';

export interface DriftFinding {
  kind: 'missing-frozen-artifact';
  publicationId: string;
  detail: string;
}

/** Verifies the one invariant the platform cannot: every live version row's frozen
 * artifact FuncCall still resolves. Row-level invariants are guaranteed by the
 * transactional publish/approve flows; this referential scan remains only because
 * FuncCalls are not entities (no FK, no ACL) — deletable until that core work lands. */
export async function driftCheck(): Promise<DriftFinding[]> {
  const findings: DriftFinding[] = [];
  const live = await grok.dapi.domains.table<AlignmentRow>(T_ALIGNMENT).query({limit: 10000});
  for (const row of live) {
    try {
      await grok.dapi.functions.calls.allPackageVersions().find(row.artifact_id);
    } catch (_) {
      findings.push({kind: 'missing-frozen-artifact', publicationId: row.publication_id,
        detail: `frozen meta call ${row.artifact_id} of v${row.revision} does not resolve`});
    }
  }
  return findings;
}
