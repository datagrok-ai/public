/** The optional approval workflow (F-30, FR-57): the ArtifactAlignment pattern rebuilt on the
 * portal's own tables. A portfolio that opts in names the gate decisions and the artifacts that
 * need a second person. Such a decision is proposed instead of taking effect; a document or model
 * card is submitted instead of approved. It takes effect once `approvers` approvers other than the
 * proposer approve it, each with a rationale; one rejection ends it. Off by default: without
 * `approval` every decision and approval takes effect at once, as before.
 *
 * Pure: the caller reads the request, the portfolio's config and who the user is. */
import type {Decision} from './gates';

export type ApprovalArtifact = 'document' | 'model_card';
export type ApprovalDecision = Exclude<Decision, 'sync'>;
export type RequestKind = 'decision' | ApprovalArtifact;
export type RequestStatus = 'pending' | 'approved' | 'rejected' | 'withdrawn';
export type ReviewVerdict = 'approve' | 'reject';

export const APPROVAL_DECISIONS: ApprovalDecision[] = ['advance', 'recycle', 'hold', 'terminate'];
export const APPROVAL_ARTIFACTS: ApprovalArtifact[] = ['document', 'model_card'];
/** More than this many approvers per request is a committee, not a workflow. */
export const MAX_APPROVERS = 5;

export interface ApprovalConfig {
  /** Gate decisions that are proposed and need approval: advance, recycle, hold, terminate. */
  decisions?: ApprovalDecision[];
  /** Artifacts submitted for approval instead of approved at once. */
  artifacts?: ApprovalArtifact[];
  /** How many approvers other than the proposer must approve; 1 when left out. */
  approvers?: number;
}

/** Every problem with a portfolio's `approval`, as sentences. */
export function validateApproval(a: ApprovalConfig | undefined | null, stageAuthority: string): string[] {
  if (a === undefined || a === null)
    return [];
  if (typeof a !== 'object' || Array.isArray(a))
    return ['approval must be an object: {decisions, artifacts, approvers}'];
  const errors: string[] = [];
  if (a.decisions !== undefined && !Array.isArray(a.decisions))
    errors.push('approval.decisions must be a list');
  for (const d of Array.isArray(a.decisions) ? a.decisions : []) {
    if (!APPROVAL_DECISIONS.includes(d))
      errors.push(`approval.decisions: "${d}" is not one of ${APPROVAL_DECISIONS.join(', ')}`);
  }
  if (Array.isArray(a.decisions) && a.decisions.length && stageAuthority !== 'portal')
    errors.push('approval.decisions needs stageAuthority "portal": where the registry owns the stage, nobody decides here');
  if (a.artifacts !== undefined && !Array.isArray(a.artifacts))
    errors.push('approval.artifacts must be a list');
  for (const x of Array.isArray(a.artifacts) ? a.artifacts : []) {
    if (!APPROVAL_ARTIFACTS.includes(x))
      errors.push(`approval.artifacts: "${x}" is not one of ${APPROVAL_ARTIFACTS.join(', ')}`);
  }
  if (a.approvers !== undefined && (!Number.isInteger(a.approvers) || a.approvers < 1 || a.approvers > MAX_APPROVERS))
    errors.push(`approval.approvers must be a whole number from 1 to ${MAX_APPROVERS}`);
  if (!(a.decisions?.length) && !(a.artifacts?.length))
    errors.push('approval names no decisions and no artifacts; leave it out to turn the workflow off');
  return errors;
}

/** The config as stored on the portfolio row, or null when the workflow is off. */
export function approvalOf(stored: unknown): ApprovalConfig | null {
  const a = stored as ApprovalConfig | null | undefined;
  return a && typeof a === 'object' && ((a.decisions?.length ?? 0) > 0 || (a.artifacts?.length ?? 0) > 0) ? a : null;
}

/** Whether a decision, or an artifact's approval, of this kind needs a second person. */
export function needsApproval(config: ApprovalConfig | null, kind: RequestKind, decision?: string | null): boolean {
  if (!config)
    return false;
  if (kind === 'decision')
    return !!decision && (config.decisions ?? []).includes(decision as ApprovalDecision);
  return (config.artifacts ?? []).includes(kind);
}

export const approvalsNeeded = (config: ApprovalConfig | null): number => config?.approvers ?? 1;

export interface Review {
  /** The reviewer's user id. */
  by: string;
  on: string;
  verdict: ReviewVerdict;
  rationale: string;
}

/** A request as the rules see it. */
export interface RequestRef {
  id: string;
  kind: RequestKind;
  status: RequestStatus;
  title: string;
  proposedBy: string | null;
  proposedOn?: string | null;
  /** `data`: proposed by the program's data, so no person is excluded from approving it. */
  source: 'person' | 'data';
  needed: number;
  reviews: Review[];
  programId?: string | null;
  portfolioId: string;
}

/** The reviews stored on a request: {items: [...]}, since a json column holds an object. */
export const reviewsOf = (stored: unknown): Review[] =>
  Array.isArray((stored as any)?.items) ? (stored as any).items : [];

/** Why a proposal or submission is refused, empty when it may be made. */
export function checkProposal(c: {kind: RequestKind, config: ApprovalConfig | null, decision?: string | null,
  mayPropose: boolean, pending: RequestRef[]}): string[] {
  const refused: string[] = [];
  if (!needsApproval(c.config, c.kind, c.decision))
    refused.push(c.kind === 'decision' ? `the portfolio does not ask for approval of ${c.decision ?? 'this decision'}: record it` :
      `the portfolio does not ask for approval of a ${c.kind.replace('_', ' ')}: approve it`);
  if (!c.mayPropose)
    refused.push(c.kind === 'decision' ? 'Only approvers propose gate decisions' :
      `Only contributors and approvers submit a ${c.kind.replace('_', ' ')}`);
  const open = c.pending.find((r) => r.status === 'pending');
  if (open)
    refused.push(`"${open.title}" is waiting for approval: approve, reject or withdraw it first`);
  return refused;
}

/** Why [reviewer] may not review [r] now, empty when they may. The proposer never approves their
 * own proposal; a proposal by data may be approved by any approver. One review per person. */
export function checkReview(r: RequestRef, c: {reviewer: string, isApprover: boolean, verdict: string,
  rationale: string | null | undefined}): string[] {
  const refused: string[] = [];
  if (r.status !== 'pending')
    refused.push(`"${r.title}" is ${r.status}, not waiting for approval`);
  if (!['approve', 'reject'].includes(c.verdict))
    refused.push(`verdict "${c.verdict}" is not approve or reject`);
  if (!c.isApprover)
    refused.push('Only approvers of the portfolio or the program review a request');
  if (r.source !== 'data' && r.proposedBy && r.proposedBy === c.reviewer)
    refused.push('You proposed it: another approver must review it');
  if (r.reviews.some((x) => x.by === c.reviewer))
    refused.push('You have reviewed it already');
  if (!c.rationale?.trim())
    refused.push('A review needs a rationale');
  return refused;
}

/** What a request becomes after [reviews]: rejected on any rejection, approved once [needed]
 * different people approved, else still pending. */
export function statusAfter(reviews: Review[], needed: number): RequestStatus {
  if (reviews.some((x) => x.verdict === 'reject'))
    return 'rejected';
  const approvers = new Set(reviews.filter((x) => x.verdict === 'approve').map((x) => x.by));
  return approvers.size >= Math.max(1, needed) ? 'approved' : 'pending';
}

/** Why [user] may not withdraw [r], empty when they may: the proposer, while it is pending. */
export function checkWithdraw(r: RequestRef, user: string): string[] {
  const refused: string[] = [];
  if (r.status !== 'pending')
    refused.push(`"${r.title}" is ${r.status}, not waiting for approval`);
  if (r.source === 'data' || r.proposedBy !== user)
    refused.push('Only the person who proposed it withdraws it; an approver rejects it');
  return refused;
}

/** The requests waiting for [user]: pending, not theirs (unless proposed by data), not reviewed by
 * them yet, and where [mayReview] says they are an approver. Oldest first. */
export function waitingFor<T extends RequestRef>(requests: T[], user: string, mayReview: (r: T) => boolean): T[] {
  return requests.filter((r) => r.status === 'pending' && (r.source === 'data' || r.proposedBy !== user) &&
    !r.reviews.some((x) => x.by === user) && mayReview(r))
    .sort((a, b) => String(a.proposedOn ?? '').localeCompare(String(b.proposedOn ?? '')));
}

/** `1 of 2 approvals` */
export function describeProgress(r: Pick<RequestRef, 'reviews' | 'needed' | 'status'>): string {
  const approvals = new Set(r.reviews.filter((x) => x.verdict === 'approve').map((x) => x.by)).size;
  return r.status === 'pending' ? `${approvals} of ${r.needed} approval${r.needed === 1 ? '' : 's'}` : r.status;
}

/** The one-line title of a proposed gate decision: `EXB-101: advance to Lead Generation`. */
export function decisionTitle(program: string, decision: string, toStage?: string | null): string {
  return `${program}: ${decision}${toStage && (decision === 'advance' || decision === 'recycle') ? ` to ${toStage}` : ''}`;
}

/** The setting in words, for the program page and the sync history. */
export function describeApproval(config: ApprovalConfig | null): string {
  if (!config)
    return 'Decisions and approvals take effect at once';
  const what = [...(config.decisions ?? []), ...(config.artifacts ?? []).map((a) => a.replace('_', ' '))];
  const n = approvalsNeeded(config);
  return `Needs ${n} approver${n === 1 ? '' : 's'} besides the proposer: ${what.join(', ')}`;
}

/** Whether an approval may take effect now on [r]: it waits, or its last approval has claimed it
 * (status approved, not yet closed) and the effect is being applied. */
export const takesEffect = (r: {status: string, closed_on?: unknown}): boolean =>
  r.status === 'pending' || (r.status === 'approved' && !r.closed_on);

/** At most this many approvers are told about one new request. */
export const MAX_APPROVERS_TOLD = 25;

/** The notice an approver gets in the platform's notifications when a request starts waiting:
 * what, who proposed it, and why. */
export function requestNotice(r: {title: string, source?: string | null, rationale?: string | null},
  proposer: string | null): string {
  const by = r.source === 'data' ? 'by data' : `by ${proposer || 'someone'}`;
  return `${r.title}: waits for your approval, proposed ${by}${r.rationale?.trim() ? ` — ${r.rationale.trim()}` : ''}`;
}

/** Who is told about a new request: each approver once, never the proposer, at most [max]. */
export function approversToTell(approvers: (string | null | undefined)[], proposer: string | null,
  max = MAX_APPROVERS_TOLD): string[] {
  return [...new Set(approvers.filter((x): x is string => !!x && x !== proposer))].slice(0, max);
}
