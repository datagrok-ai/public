import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {
  AlignmentRow, HistoryRow, PublicationStatus,
  SCHEMA, T_ALIGNMENT, T_COMPOUND, T_HISTORY, T_TAG,
} from '../domain/constants';
import {currentUserIn, getProgramGroups} from '../domain/security';
import {cloneRun} from './clone-service';

const alignment = () => grok.dapi.domains.table<AlignmentRow>(T_ALIGNMENT);
const history = () => grok.dapi.domains.table<HistoryRow>(T_HISTORY);

const keyFilter = (programId: string, studyId: string | null | undefined, name: string): DG.DomainConditionTree =>
  DG.and(
    DG.cond('program_id', '=', programId),
    DG.cond('study_id', studyId == null ? 'is' : '=', studyId ?? null),
    DG.cond('name', '=', name),
  );

async function ensureRegistryRows(table: string, names: string[], column: string): Promise<string[]> {
  const client = grok.dapi.domains.table(table);
  const ids: string[] = [];
  for (const name of names) {
    const existing = await client.getByKey({[column]: name});
    ids.push(existing?.id ?? (await client.insert({[column]: name}))[0].id);
  }
  return ids;
}

const joinNames = (links?: {name: string}[]) =>
  (links ?? []).map((l) => l.name).join(',');

/** Transaction ops copying a live version row into alignment_history and removing
 * it — always combined with the write they make room for, so the catalog never
 * exposes an intermediate state. */
function archiveOps(row: AlignmentRow, supersededOn: dayjs.Dayjs): DG.DomainTransactionOp[] {
  return [
    {op: 'insert', table: 'alignment_history', values: {
      publication_id: row.publication_id,
      revision: row.revision,
      name: row.name,
      artifact_id: row.artifact_id,
      program_id: row.program_id,
      study_id: row.study_id,
      workstream: row.workstream,
      artifact_type: row.artifact_type,
      artifact_author: row.artifact_author,
      artifact_created_on: row.artifact_created_on,
      description: row.description,
      final_status: row.status,
      approved_by: row.approved_by,
      approved_on: row.approved_on,
      reject_reason: row.reject_reason,
      path: row.path,
      tags_snapshot: joinNames(row.tags),
      compounds_snapshot: joinNames(row.compounds),
      superseded_on: supersededOn,
      original_row_id: row.id,
    }},
    {op: 'delete', table: 'alignment', id: row.id!},
  ];
}

async function liveRows(publicationId: string): Promise<AlignmentRow[]> {
  return alignment().query({
    filter: DG.cond('publication_id', '=', publicationId),
    expand: ['tags', 'compounds'],
  });
}

export interface PublishRequest {
  /** A saved run's meta FuncCall id; required for workflow runs. */
  sourceMetaCallId?: string;
  /** A live in-memory function run — publishes without a prior history save. */
  sourceCall?: DG.FuncCall;
  programId: string;
  studyId?: string | null;
  name: string;
  workstream?: string;
  description?: string;
  /** Compound registration codes; created in the registry when missing. */
  compounds?: string[];
  /** Free tag names; created when missing. Defaults to the previous version's tags. */
  tags?: string[];
  path?: string;
  /** Leaves the new version pending even though the gate is off — the explicit-review
   * path (and the way tests exercise the pending state). */
  skipAutoApprove?: boolean;
}

export interface PublishResult {
  publicationId: string;
  revision: number;
  rowId: string;
  /** The frozen meta FuncCall id. */
  artifactId: string;
  status: PublicationStatus;
}

/** Publishes a saved Compute2 run — a workflow run, or a single function run (an RFV
 * model run or an individually saved step) — into a program: freezes it into a
 * stamped catalog copy, then inserts the next version row of the publication addressed
 * by the (program, study, name) republish key. The previous approved version stays
 * live until this version is approved; a stale in-review row (rejected or superseded
 * draft) is archived first. The review gate is OFF for workflow runs, so the service
 * approves in the same breath — unless the caller lacks approval rights, in which
 * case the version stays pending. */
export async function publishWorkflowRun(req: PublishRequest): Promise<PublishResult> {
  const source = req.sourceCall ?? req.sourceMetaCallId;
  if (source == null)
    throw new Error('Either sourceCall or sourceMetaCallId is required');
  const groups = await getProgramGroups(req.programId);
  await requireProgramWriter(groups, 'publish into');
  const key = keyFilter(req.programId, req.studyId, req.name);

  const liveByKey = await alignment().query({filter: key, expand: ['tags', 'compounds']});
  const publicationId = liveByKey[0]?.publication_id ?? crypto.randomUUID();
  const revision = Math.max(0, ...liveByKey.map((r) => r.revision)) + 1;

  const now = dayjs();
  const staleInReview = liveByKey.find((r) => r.status !== 'approved');
  const previous = liveByKey.find((r) => r.status === 'approved') ?? staleInReview;

  const clone = await cloneRun(source, {
    audience: groups.viewers ?? undefined,
    publicationId,
    title: `${req.name} v${revision}`,
    description: req.description,
  });

  const tagNames = req.tags ?? previous?.tags?.map((t) => t.name) ?? [];
  const compoundCodes = req.compounds ?? previous?.compounds?.map((c) => c.name) ?? [];
  // archive the superseded draft and insert the next version atomically
  const results = await grok.dapi.domains.transaction(SCHEMA, [
    ...(staleInReview != null ? archiveOps(staleInReview, now) : []),
    {op: 'insert', table: 'alignment', values: {
      publication_id: publicationId,
      revision,
      name: req.name,
      artifact_id: clone.metaCall.id,
      program_id: req.programId,
      study_id: req.studyId ?? undefined,
      workstream: req.workstream ?? previous?.workstream ?? undefined,
      artifact_type: clone.artifactType,
      artifact_author: clone.sourceAuthorId,
      artifact_created_on: clone.sourceStartedOn,
      description: req.description ?? previous?.description ?? undefined,
      path: req.path ?? previous?.path ?? undefined,
      tags: await ensureRegistryRows(T_TAG, tagNames, 'name'),
      compounds: await ensureRegistryRows(T_COMPOUND, compoundCodes, 'registration_code'),
    }},
  ]);
  const ins = results[results.length - 1] as DG.DomainInsertResult;
  if (!ins.created)
    throw new Error(`Version row was not created (${ins.status})`);

  // Gate off for workflow runs: auto-approve with the same audit trail as a human
  // approval. Under the emulated service identity this runs as the publishing user,
  // so it succeeds only when they hold approval rights; otherwise the version stays
  // pending for a real approver.
  let status: PublicationStatus = 'pending';
  if (req.skipAutoApprove !== true) {
    try {
      await approvePublication(ins.id, {auto: true});
      status = 'approved';
    } catch (_) {/* stays pending */}
  }

  return {publicationId, revision, rowId: ins.id, artifactId: clone.metaCall.id, status};
}

/** Flips a pending version to approved: archives the previously approved version and
 * stamps the approval columns in one transaction, so the publication always shows
 * exactly one approved row. Human approvals enforce approvers-group membership and
 * author != approver service-side; `auto` marks the gate-off same-breath approval. */
export async function approvePublication(rowId: string, options?: {auto?: boolean}): Promise<void> {
  await DG.retryOnVersionConflict(async () => {
    const row = await alignment().get(rowId);
    if (row == null)
      throw new Error(`Version row ${rowId} not found`);
    if (row.status !== 'pending')
      throw new Error(`Only a pending version can be approved (status: ${row.status})`);
    const user = await grok.dapi.users.current();
    if (options?.auto !== true) {
      const groups = await getProgramGroups(row.program_id);
      if (!await currentUserIn(groups.approvers))
        throw new Error('Only members of the program approvers group can approve');
      if (row.artifact_author != null && row.artifact_author === user.id)
        throw new Error('Authors cannot approve their own publication');
    }
    const now = dayjs();
    const currentApproved = (await liveRows(row.publication_id)).find((r) => r.status === 'approved');
    await grok.dapi.domains.transaction(SCHEMA, [
      ...(currentApproved != null && currentApproved.id !== row.id ? archiveOps(currentApproved, now) : []),
      {op: 'update', table: 'alignment', id: rowId,
        values: {status: 'approved', approved_by: user.id, approved_on: now},
        expectedVersion: row.version},
    ]);
  });
}

/** Rejects a pending version with a reason. The rejected row stays live (occupying
 * the single in-review slot) until the author resubmits, which archives it; the
 * audience keeps seeing the last approved version throughout. */
export async function rejectPublication(rowId: string, reason: string,
  options?: {auto?: boolean}): Promise<void> {
  await DG.retryOnVersionConflict(async () => {
    const row = await alignment().get(rowId);
    if (row == null)
      throw new Error(`Version row ${rowId} not found`);
    if (row.status !== 'pending')
      throw new Error(`Only a pending version can be rejected (status: ${row.status})`);
    const user = await grok.dapi.users.current();
    if (options?.auto !== true) {
      const groups = await getProgramGroups(row.program_id);
      if (!await currentUserIn(groups.approvers))
        throw new Error('Only members of the program approvers group can reject');
      if (row.artifact_author != null && row.artifact_author === user.id)
        throw new Error('Authors cannot review their own publication');
    }
    await alignment().update(rowId, {
      status: 'rejected', reject_reason: reason, approved_by: user.id, approved_on: dayjs(),
    } as any, {version: row.version});
  });
}

/** Per-program write authority: transitive membership in the program's
 * contributors or approvers group — the same groups the program row's Edit is
 * shared to. Resolved through the acting user's session (holds under
 * impersonation, unlike DomainRow.permissions(), whose client-side admin
 * fast-path always answers true in an admin browser session). The table-level
 * umbrella grants are wider (a master-row-scoped write predicate is a core
 * ask), so the services enforce per-program isolation themselves. */
async function requireProgramWriter(
  groups: {contributors: DG.Group | null, approvers: DG.Group | null}, action: string): Promise<void> {
  if (!await currentUserIn(groups.contributors) && !await currentUserIn(groups.approvers))
    throw new Error(`Only members of the program contributors or approvers groups can ${action} it`);
}

/** Post-publish curation: retags/moves a live version row. Requires Edit on the
 * curation columns (contributors umbrella) plus program-group membership. */
export async function updateCuration(rowId: string,
  curation: {tags?: string[], path?: string}): Promise<void> {
  const row = await alignment().get(rowId);
  if (row == null)
    throw new Error(`Version row ${rowId} not found`);
  await requireProgramWriter(await getProgramGroups(row.program_id), 'curate rows of');
  const values: any = {};
  if (curation.tags != null)
    values.tags = await ensureRegistryRows(T_TAG, curation.tags, 'name');
  if (curation.path != null)
    values.path = curation.path;
  await alignment().update(rowId, values);
}

export async function getLiveVersions(publicationId: string): Promise<AlignmentRow[]> {
  return liveRows(publicationId);
}

export async function getPublicationHistory(publicationId: string): Promise<HistoryRow[]> {
  return history().query({filter: DG.cond('publication_id', '=', publicationId), sort: '!revision'});
}

/** The discovery slice: approved rows only, security-trimmed by the caller's
 * program-row access; exactly one row per publication by construction. */
export async function discover(programId?: string): Promise<AlignmentRow[]> {
  const conds: DG.DomainCondition[] = [DG.cond('status', '=', 'approved')];
  if (programId != null)
    conds.push(DG.cond('program_id', '=', programId));
  return alignment().query({filter: DG.and(...conds), expand: ['tags', 'compounds'], sort: '!updated_on'});
}
