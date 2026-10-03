/** The optional approval workflow against the platform (F-30, FR-57): proposing a gate decision,
 * submitting a document or a model card, reviewing a request, and the lists of what waits. An
 * approved decision goes through recordDecision (evidence, workspace growth), an approved document
 * through linkDocumentTo (criterion met, older version superseded). The rules live in
 * ./model/approvals.
 *
 * `approval_request` is table-secured, like program_role. The server lets only the portfolios'
 * contributors groups (approvers are nested in them) write it (./write-grants); which contributor
 * may review which request is the app's check. Master security through the program would hide a
 * request from every approver but the program's author until GROK-21041 is fixed. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {ApprovalRequestRow, portalDb, PortfolioRow, ProgramRow} from './generated/db';
import {isApprover, isContributor, programContext, recordDecision, RecordedDecision, verdictFor} from './portal';
import {programAccess, rolesOf} from './program-roles';
import {notifyUser} from './notify';
// evidence.ts and snapshots.ts reach portal.ts, which does not import this file
import {suggestedEvidence} from './evidence';
import {summaryOf, takeSnapshot} from './snapshots';
import {snapshotEvidence} from './model/snapshots';
import {linkDocumentTo} from './documents';
import {describeDocument} from './model/documents';
import {EvidenceItem} from './model/evidence';
import {Decision} from './model/gates';
import {
  approvalOf, approvalsNeeded, approversToTell, checkProposal, checkReview, checkWithdraw, decisionTitle, describeProgress,
  RequestKind, RequestRef, requestNotice, Review, reviewsOf, ReviewVerdict, statusAfter, waitingFor,
} from './model/approvals';
import {newestFirstBy, oldestFirstBy} from './model/time';
import {filterValue, ROW_LIMIT} from './rows';

export const requestRef = (r: ApprovalRequestRow): RequestRef => ({id: r.id, kind: r.kind as RequestKind,
  status: r.status as RequestRef['status'], title: r.title, proposedBy: r.proposed_by ?? null,
  proposedOn: r.proposed_on ? dayjs(r.proposed_on as any).toISOString() : null, source: r.source === 'data' ? 'data' : 'person',
  needed: r.approvals_needed ?? 1, reviews: reviewsOf(r.reviews), programId: r.program_id ?? null, portfolioId: r.portfolio_id});

const newestFirst = (rows: ApprovalRequestRow[]) =>
  rows.sort(newestFirstBy((x) => x.created_on));

/** Every request of a program, newest first. */
export async function requestsOf(programId: string): Promise<ApprovalRequestRow[]> {
  return newestFirst(await portalDb.approvalRequests.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT}));
}

/** The requests about one document or model card, newest first. */
async function requestsFor(kind: RequestKind, target: string): Promise<ApprovalRequestRow[]> {
  return newestFirst(await portalDb.approvalRequests.query({filter: `kind = ${filterValue(kind)} and target = ${filterValue(target)}`, limit: ROW_LIMIT}));
}

export interface Proposal {
  request: ApprovalRequestRow;
  line: string;
  /** How telling the approvers went: `2 of 2 approver(s) notified`. */
  notified?: string;
}

/** Groups are walked this deep for their users: a group of groups of users. */
const GROUP_DEPTH = 4;

/** The users in [groupIds], nested groups included, by user id. */
async function usersOfGroups(groupIds: string[]): Promise<string[]> {
  const users: string[] = [];
  const seen = new Set<string>();
  let frontier = groupIds;
  for (let depth = 0; depth < GROUP_DEPTH && frontier.length; depth++) {
    const next: string[] = [];
    for (const id of frontier) {
      if (seen.has(id))
        continue;
      seen.add(id);
      const g = await grok.dapi.groups.find(id).catch(() => null);
      if (!g)
        continue;
      if (g.personal) {
        const u = await grok.dapi.groups.getUser(g).catch(() => null);
        if (u)
          users.push(u.id);
      } else
        next.push(...(g.members ?? []).map((m) => m.id));
    }
    frontier = next;
  }
  return users;
}

/** The users who review a request about [programId]: the portfolio's approvers group and the
 * program's approver roles; for a model card (no program), the portfolio's approvers. */
export async function approverUsers(portfolio: PortfolioRow, programId: string | null): Promise<string[]> {
  const roles = programId ? (await rolesOf(programId).catch(() => [])).filter((r) => r.role === 'approver') : [];
  return usersOfGroups([portfolio.approvers, ...roles.map((r) => r.group_id)].filter((x): x is string => !!x));
}

/** Tells the approvers of [request], in the platform's notifications, that it waits for them. The
 * proposer is not told. Never throws: says how it went. The server lets a user notify others only
 * when they may edit users; otherwise each approver finds it under Waiting for me. */
export async function tellApprovers(request: ApprovalRequestRow, portfolio: PortfolioRow): Promise<string> {
  try {
    const me = DG.User.current();
    const ids = approversToTell(await approverUsers(portfolio, request.program_id ?? null), me.id);
    if (!ids.length)
      return 'no other approver to tell';
    const text = requestNotice(request, me.friendlyName);
    const results = await Promise.all(ids.map((id) => notifyUser(id, text)));
    const failed = results.find((r) => !r.ok);
    return `${results.filter((r) => r.ok).length} of ${ids.length} approver(s) notified` + (failed ? `; ${failed.detail}` : '');
  } catch (e: any) {
    return `approvers not notified (${e?.message ?? e})`;
  }
}

/** Proposes a gate decision in a portfolio that asks for approval of it. Checked like the decision
 * itself (the verdict as the proposer, an approver), and refused while another decision on the
 * program waits. The evidence is frozen now: the items given, else the gate panel's suggestion,
 * and a snapshot taken now when asked. Nothing changes on the program until it is approved. */
export async function proposeDecision(programId: string, decision: Decision, rationale: string,
  by: {source?: 'person' | 'data', items?: EvidenceItem[], snapshot?: boolean, evidence?: string, toStageId?: string | null} = {}):
  Promise<Proposal> {
  const ctx = await programContext(programId);
  const config = approvalOf(ctx.portfolio.approval);
  const verdict = verdictFor(ctx, decision, by.toStageId ?? null, rationale);
  const pending = (await requestsOf(programId)).filter((r) => r.kind === 'decision').map(requestRef);
  const refused = [...checkProposal({kind: 'decision', config, decision, mayPropose: ctx.approver, pending}),
    ...(verdict.allowed ? [] : verdict.reasons)];
  if (refused.length)
    throw new Error(refused.join('; '));
  let items = by.items ?? await suggestedEvidence(programId, ctx.readiness,
    ctx.stages.find((s) => s.id === ctx.readiness?.stageId)?.name ?? '').catch(() => []);
  if (by.snapshot) {
    const snap = await takeSnapshot(ctx.program, {reason: `at the ${decision} proposal`});
    items = [...items, snapshotEvidence(snap.snapshot.id, summaryOf(snap.snapshot))];
  }
  const to = ctx.stages.find((s) => s.id === verdict.toStageId);
  const [r] = await portalDb.approvalRequests.insert({portfolio_id: ctx.portfolio.id, program_id: programId, kind: 'decision',
    decision: decision as any, to_stage_id: verdict.toStageId ?? undefined, title: decisionTitle(ctx.program.code, decision, to?.name),
    rationale, source: by.source ?? 'person', evidence: {items, ...(by.evidence ? {text: by.evidence} : {})},
    proposed_by: DG.User.current().id, proposed_on: dayjs(), approvals_needed: approvalsNeeded(config), status: 'pending'});
  const request = await portalDb.approvalRequests.get(r.id);
  return {request, line: `${request.title}: proposed; waits for ${describeProgress(requestRef(request)).replace(/^0 of /, '')}`,
    notified: await tellApprovers(request, ctx.portfolio)};
}

/** What a gate decision does: recorded at once, or proposed when the portfolio asks for approval of
 * it. The one entry point for the gate panel and for advance by data. */
export async function decideOrPropose(programId: string, decision: Decision, rationale: string,
  by: {source?: 'approver' | 'data', items?: EvidenceItem[], snapshot?: boolean, evidence?: string, toStageId?: string | null} = {}):
  Promise<{recorded?: RecordedDecision, proposed?: Proposal}> {
  const program = await portalDb.programs.get(programId);
  const portfolio = await portalDb.portfolios.get(program.portfolio_id);
  const config = approvalOf(portfolio.approval);
  if ((config?.decisions ?? []).includes(decision as any)) {
    return {proposed: await proposeDecision(programId, decision, rationale, {source: by.source === 'data' ? 'data' : 'person',
      items: by.items, snapshot: by.snapshot, evidence: by.evidence, toStageId: by.toStageId})};
  }
  return {recorded: await recordDecision(programId, decision, by.toStageId ?? null, rationale,
    {source: by.source ?? 'approver', items: by.items, snapshot: by.snapshot, evidence: by.evidence})};
}

/** Submits a document for approval (draft or rejected → submitted). Contributors and approvers. */
export async function submitDocument(program: ProgramRow, documentId: string, rationale?: string): Promise<Proposal> {
  const ctx = await programContext(program.id);
  const access = await programAccess(ctx.portfolio, program.id);
  const doc = await portalDb.documents.get(documentId);
  if (doc.program_id !== program.id)
    throw new Error(`${program.code} has no document ${documentId}`);
  const refused = checkProposal({kind: 'document', config: approvalOf(ctx.portfolio.approval),
    mayPropose: access.approver || access.contributor, pending: (await requestsFor('document', doc.id)).map(requestRef)});
  if (!['draft', 'rejected'].includes(doc.status))
    refused.push(`the document is ${doc.status}; only a draft or a rejected document is submitted`);
  if (refused.length)
    throw new Error(refused.join('; '));
  const [r] = await portalDb.approvalRequests.insert({portfolio_id: ctx.portfolio.id, program_id: program.id, kind: 'document',
    target: doc.id, title: `${program.code}: approve ${describeDocument({title: doc.title, version: doc.doc_version, status: 'draft'})
      .replace(/ · draft$/, '')}`, rationale: rationale ?? undefined, source: 'person',
    proposed_by: DG.User.current().id, proposed_on: dayjs(), approvals_needed: approvalsNeeded(approvalOf(ctx.portfolio.approval)),
    status: 'pending'});
  // the request first, then the document: a failed insert leaves no document stuck as submitted
  await markSubmitted(r.id, () => portalDb.documents.update(doc.id, {status: 'submitted'}));
  const request = await portalDb.approvalRequests.get(r.id);
  return {request, line: `${request.title}: submitted for approval`, notified: await tellApprovers(request, ctx.portfolio)};
}

/** Marks the target of request [id] as waiting; if that fails, the request is removed again. */
async function markSubmitted(id: string, mark: () => Promise<unknown>): Promise<void> {
  try {
    await mark();
  } catch (e: any) {
    await portalDb.approvalRequests.delete(id).catch(() => undefined);
    throw e;
  }
}

/** Submits a model card for approval (draft or rejected → in review), in a portfolio whose
 * catalogue holds the model. Contributors and approvers of that portfolio. */
export async function submitModelCard(portfolio: PortfolioRow, modelCode: string, version?: string | null,
  rationale?: string): Promise<Proposal> {
  const model = await portalDb.models.getByKey({code: modelCode});
  if (!model)
    throw new Error(`${modelCode} is not in the model catalogue`);
  if (!(model.portfolios ?? '').split('|').includes(portfolio.code))
    throw new Error(`${modelCode} is not in ${portfolio.code}'s catalogue`);
  const v = version || model.current_version;
  const card = v ? await portalDb.modelCards.getByKey({model_id: model.id, model_version: v}) : null;
  if (!card)
    throw new Error(`${modelCode} has no card${v ? ` for version ${v}` : ''}`);
  const may = await isApprover(portfolio) || await isContributor(portfolio);
  const refused = checkProposal({kind: 'model_card', config: approvalOf(portfolio.approval), mayPropose: may,
    pending: (await requestsFor('model_card', card.id)).map(requestRef)});
  if (!['draft', 'rejected'].includes(card.approval))
    refused.push(`the card is ${card.approval.replace('_', ' ')}; only a draft or a rejected card is submitted`);
  if (refused.length)
    throw new Error(refused.join('; '));
  const [r] = await portalDb.approvalRequests.insert({portfolio_id: portfolio.id, kind: 'model_card', target: card.id,
    title: `${model.name} ${card.model_version}: approve model card`, rationale: rationale ?? undefined, source: 'person',
    proposed_by: DG.User.current().id, proposed_on: dayjs(), approvals_needed: approvalsNeeded(approvalOf(portfolio.approval)),
    status: 'pending'});
  await markSubmitted(r.id, () => portalDb.modelCards.update(card.id, {approval: 'in_review'}));
  const request = await portalDb.approvalRequests.get(r.id);
  return {request, line: `${request.title}: submitted for approval`, notified: await tellApprovers(request, portfolio)};
}

/** Whether the current user reviews [r]: an approver of its program (portfolio or program role),
 * or of its portfolio for a model card. */
async function mayReview(r: ApprovalRequestRow, portfolios: Map<string, PortfolioRow> = new Map()): Promise<boolean> {
  const portfolio = portfolios.get(r.portfolio_id) ?? await portalDb.portfolios.get(r.portfolio_id);
  portfolios.set(portfolio.id, portfolio);
  if (r.program_id)
    return (await programAccess(portfolio, r.program_id)).approver;
  return isApprover(portfolio);
}

export interface ReviewResult {
  request: ApprovalRequestRow;
  line: string;
  /** What approving did: the decision recorded, the document approved. */
  effect?: RecordedDecision;
  notified?: string;
}

/** Approves or rejects a request with a rationale. On the last approval it takes effect: the
 * decision is recorded by the current user (source as proposed), the document approved, the card
 * approved. If that fails, the review is not kept and the request stays pending. A rejection is
 * logged on the request; the document or card is marked rejected. The proposer is told. */
export async function reviewRequest(requestId: string, verdict: ReviewVerdict, rationale: string): Promise<ReviewResult> {
  const row = await portalDb.approvalRequests.get(requestId);
  const ref = requestRef(row);
  const me = DG.User.current().id;
  const refused = checkReview(ref, {reviewer: me, isApprover: await mayReview(row), verdict, rationale});
  if (refused.length)
    throw new Error(refused.join('; '));
  const reviews: Review[] = [...ref.reviews, {by: me, on: dayjs().toISOString(), verdict, rationale: rationale.trim()}];
  const status = statusAfter(reviews, ref.needed);
  // claim the request first, with its new status, against the version read: a second reviewer (or a
  // second click) at the same time fails here, or then finds it closed, instead of losing a review
  // or applying the decision twice. If the effect fails, the claim is undone.
  try {
    await portalDb.approvalRequests.update(row.id, {reviews: {items: reviews}, status} as any, {version: row.version});
  } catch (e: any) {
    throw new Error(`${row.title}: reviewed by someone else meanwhile; reload and try again (${e?.message ?? e})`);
  }
  let result = '';
  let effect: RecordedDecision | undefined;
  let transition: string | undefined;
  if (status === 'approved') {
    try {
      if (row.kind === 'decision') {
        const ev = (row.evidence ?? {}) as {items?: EvidenceItem[], text?: string};
        effect = await recordDecision(row.program_id!, row.decision as Decision, row.to_stage_id ?? null, row.rationale ?? '',
          {source: row.source === 'data' ? 'data' : 'approver', items: ev.items ?? [], evidence: ev.text, approvedRequest: row.id});
        transition = effect.transitionId;
        result = `${row.decision} recorded` + (effect.evidence?.written ? ` with ${effect.evidence.written} piece(s) of evidence` : '') +
          (effect.workspace ? `; ${effect.workspace.split('\n')[0]}` : '');
      } else if (row.kind === 'document') {
        const program = await portalDb.programs.get(row.program_id!);
        const doc = await portalDb.documents.get(row.target!);
        const linked = await linkDocumentTo(program, {path: doc.path, status: 'approved'}, {approvedRequest: row.id});
        result = `document approved${linked.criterion ? `; ${linked.criterion.code} ${linked.criterion.state}` : ''}` +
          (linked.superseded ? `; ${linked.superseded} older version(s) superseded` : '');
      } else {
        await portalDb.modelCards.update(row.target!, {approval: 'approved', approved_on: dayjs()});
        result = 'model card approved';
      }
    } catch (e: any) {
      await portalDb.approvalRequests.update(row.id, {reviews: {items: ref.reviews}, status: 'pending'} as any)
        .catch(() => undefined);
      throw new Error(`Approved, but it could not take effect, so the request stays pending: ${e?.message ?? e}`);
    }
  } else if (status === 'rejected') {
    if (row.kind === 'document')
      await portalDb.documents.update(row.target!, {status: 'rejected'}).catch(() => undefined);
    else if (row.kind === 'model_card')
      await portalDb.modelCards.update(row.target!, {approval: 'rejected'}).catch(() => undefined);
    result = 'rejected';
  }
  await portalDb.approvalRequests.update(row.id, {reviews: {items: reviews}, status,
    ...(status !== 'pending' ? {closed_on: dayjs(), result} : {}),
    ...(transition ? {stage_transition_id: transition} : {})} as any);
  const request = await portalDb.approvalRequests.get(row.id);
  const line = `${request.title}: ${status === 'pending' ? `approved by you; ${describeProgress(requestRef(request))}` : result}`;
  let notified: string | undefined;
  if (status !== 'pending' && row.source !== 'data' && row.proposed_by && row.proposed_by !== me) {
    notified = (await notifyUser(row.proposed_by, `${request.title}: ${status} by ${DG.User.current().friendlyName}` +
      ` — ${rationale.trim()}`)).detail;
  }
  return {request, line, effect, notified};
}

/** Withdraws a pending request; only its proposer may. A document or card goes back to draft. */
export async function withdrawRequest(requestId: string): Promise<ApprovalRequestRow> {
  const row = await portalDb.approvalRequests.get(requestId);
  const refused = checkWithdraw(requestRef(row), DG.User.current().id);
  if (refused.length)
    throw new Error(refused.join('; '));
  if (row.kind === 'document')
    await portalDb.documents.update(row.target!, {status: 'draft'}).catch(() => undefined);
  else if (row.kind === 'model_card')
    await portalDb.modelCards.update(row.target!, {approval: 'draft'}).catch(() => undefined);
  await portalDb.approvalRequests.update(row.id, {status: 'withdrawn', closed_on: dayjs(), result: 'withdrawn by the proposer'} as any);
  return portalDb.approvalRequests.get(row.id);
}

/** The pending requests the current user may review, oldest first. */
export async function waitingForMe(): Promise<ApprovalRequestRow[]> {
  const pending = await portalDb.approvalRequests.query({filter: 'status = "pending"', limit: ROW_LIMIT});
  const portfolios = new Map<string, PortfolioRow>();
  const allowed = new Set<string>();
  for (const r of pending) {
    if (await mayReview(r, portfolios).catch(() => false))
      allowed.add(r.id);
  }
  const me = DG.User.current().id;
  const ids = new Set(waitingFor(pending.map(requestRef), me, (r) => allowed.has(r.id)).map((r) => r.id));
  return pending.filter((r) => ids.has(r.id))
    .sort(oldestFirstBy((x) => x.proposed_on ?? x.created_on));
}

/** The names of the users behind [ids], by id. */
export async function userNames(ids: (string | null | undefined)[]): Promise<Map<string, string>> {
  const out = new Map<string, string>();
  for (const id of new Set(ids.filter((x): x is string => !!x)))
    out.set(id, (await grok.dapi.users.find(id).catch(() => null))?.friendlyName ?? id);
  return out;
}

/** `advance to Lead Generation · proposed by Ann 2026-10-02 · 0 of 1 approvals` */
export function describeRequest(r: ApprovalRequestRow, names: Map<string, string>): string {
  const by = r.source === 'data' ? 'by data' : `by ${names.get(r.proposed_by ?? '') ?? 'someone'}`;
  const on = r.proposed_on ? dayjs(r.proposed_on as any).format('YYYY-MM-DD') : '';
  return `${r.title} · proposed ${by} ${on} · ${describeProgress(requestRef(r))}`;
}

/** One row per request: what, status, progress, who proposed it and when, the rationale, the reviews. */
export function requestsFrame(rows: ApprovalRequestRow[], names: Map<string, string>, name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('id', rows.map((r) => r.id)),
    DG.Column.fromStrings('kind', rows.map((r) => r.kind)),
    DG.Column.fromStrings('title', rows.map((r) => r.title)),
    DG.Column.fromStrings('status', rows.map((r) => r.status)),
    DG.Column.fromStrings('progress', rows.map((r) => describeProgress(requestRef(r)))),
    DG.Column.fromStrings('proposed by', rows.map((r) => r.source === 'data' ? 'data' : names.get(r.proposed_by ?? '') ?? '')),
    DG.Column.fromStrings('proposed on', rows.map((r) => r.proposed_on ? dayjs(r.proposed_on as any).toISOString() : '')),
    DG.Column.fromStrings('rationale', rows.map((r) => r.rationale ?? '')),
    DG.Column.fromStrings('reviews', rows.map((r) => reviewsOf(r.reviews)
      .map((x) => `${names.get(x.by) ?? x.by} ${x.verdict}d: ${x.rationale}`).join('; '))),
    DG.Column.fromStrings('result', rows.map((r) => r.result ?? '')),
  ]);
  df.name = name;
  return df;
}
