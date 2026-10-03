/** A program's data package manifest against the platform (F-43): reads the records the manifest
 * lists (criteria and their state, documents, evidence tables, saved runs and the criteria they
 * were selected for, datasets, snapshots) and hands them to ./model/data-package. Read only. */
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {portalDb} from './generated/db';
import {programContext} from './portal';
import {documentsOf} from './documents';
import {savedRunsOf} from './model-runs';
import {snapshotsOf} from './snapshots';
import {datasetsOf} from './program-data';
import {generatedOf} from './evidence-tables';
import {userNames} from './approvals';
import {EvidenceCheck} from './model/evidence-tables';
import {buildManifest, MANIFEST_COLUMNS, manifestCell, ManifestLine} from './model/data-package';
import {filterValue, ROW_LIMIT} from './rows';

const iso = (d: any) => d ? dayjs(d).toISOString() : null;

export interface DataPackage {
  program: {id: string, code: string, name: string};
  portfolio: string;
  lines: ManifestLine[];
}

/** The manifest of one program, as its records stand now. */
export async function dataPackageOf(programId: string): Promise<DataPackage> {
  const ctx = await programContext(programId);
  const [documents, runs, snapshots, artifacts, statuses] = await Promise.all([
    documentsOf(programId),
    savedRunsOf(programId).catch(() => []),
    snapshotsOf(programId).catch(() => []),
    portalDb.criterionArtifacts.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT}).catch(() => []),
    portalDb.criterionStatuses.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT}).catch(() => []),
  ]);
  const criterionOfStatus = new Map(statuses.map((s) => [s.id, s.gate_criterion_id]));
  const selectedFor = new Map<string, string[]>();
  for (const a of artifacts.filter((x) => x.role === 'selected')) {
    const c = criterionOfStatus.get(a.criterion_status_id);
    if (c)
      selectedFor.set(a.model_run_id, [...selectedFor.get(a.model_run_id) ?? [], c]);
  }
  const names = await userNames([...documents.flatMap((d) => [d.linked_by, d.approved_by]),
    ...runs.map((r) => r.run.run_author ?? r.run.author_id), ...snapshots.map((s) => s.taken_by)]).catch(() => new Map());
  const lines = buildManifest({
    stages: ctx.stages.map((s) => ({id: s.id, code: s.code, ordinal: s.ordinal})),
    stageId: ctx.program.stage_id ?? null,
    criteria: ctx.criteria.map((c) => ({id: c.id, stageId: c.stageId, code: c.code, name: c.name, kind: c.kind, required: c.required})),
    states: ctx.states,
    evaluated: new Set((ctx.readiness?.criteria ?? []).filter((c) => c.basis === 'evaluated').map((c) => c.id)),
    documents: documents.map((d) => {
      const g = generatedOf(d);
      return {id: d.id, title: d.title, kind: d.kind, version: d.doc_version ?? null, status: d.status, path: d.path,
        criterionId: d.gate_criterion_id ?? null, linkedBy: d.linked_by ?? null, linkedOn: iso(d.linked_on ?? d.created_on),
        approvedBy: d.approved_by ?? null, approvedOn: iso(d.approved_on), writtenBy: d.written_by ?? null,
        generated: g ? {sources: g.sources.map((s) => s.kind.replace(/_/g, ' ')),
          state: (d.evidence_check as EvidenceCheck | undefined)?.state ?? null} : null};
    }),
    runs: runs.map((r) => ({id: r.run.id, runId: r.run.run_id ?? null, title: r.run.title ?? r.run.model_name,
      model: r.model?.name ?? r.run.model_name, modelVersion: r.run.model_version ?? null, stageId: r.run.stage_id ?? null,
      unitOperation: r.unitOperation, runType: r.runType, purpose: r.run.purpose ?? null,
      by: r.run.run_author ?? r.run.author_id, on: iso(r.run.run_started ?? r.run.created_on), isReference: r.run.is_reference ?? null,
      selectedFor: selectedFor.get(r.run.id) ?? []})),
    datasets: datasetsOf(ctx.portfolio).map((d) => ({code: d.code, name: d.name, connection: d.connection})),
    snapshots: snapshots.map((s) => ({id: s.id, label: s.label, by: s.taken_by ?? null, on: iso(s.taken_on)})),
    names,
  });
  return {program: {id: ctx.program.id, code: ctx.program.code, name: ctx.program.name}, portfolio: ctx.portfolio.code, lines};
}

/** The manifest as a table: one string column per manifest column. */
export function dataPackageFrame(p: DataPackage): DG.DataFrame {
  const df = DG.DataFrame.fromColumns(MANIFEST_COLUMNS.map((c) =>
    DG.Column.fromStrings(c.name, p.lines.map((l) => manifestCell(l, c.key)))));
  df.name = `${p.program.code} data package`;
  return df;
}
