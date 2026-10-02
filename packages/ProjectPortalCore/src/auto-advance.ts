/** Stage advance driven by data against the platform: each program of a portfolio that opts in is
 * checked against its next gate, and one that the data alone meets is advanced like a manual
 * decision (gate log, criterion evidence, workspace folders), with source `data`. The rules live in
 * ./model/auto-advance. */
import {portalDb, PortfolioRow} from './generated/db';
import {isApprover, portfolioByCode, programContext, recordDecision} from './portal';
import {AutoAdvanceConfig, AutoAdvanceOutcome, planAutoAdvance, summarizeAutoAdvance} from './model/auto-advance';
import {ProgramStatus} from './model/config';
import {approvalOf, needsApproval} from './model/approvals';
import {proposeDecision, requestsOf} from './approvals';

export interface AutoAdvanceResult {
  portfolio: string;
  config: AutoAdvanceConfig | null;
  dryRun: boolean;
  outcomes: (AutoAdvanceOutcome & {
    /** Set when the move was planned but could not be recorded. */
    error?: string;
    /** The workspace report of a move that grew the program's folders. */
    workspace?: string;
  })[];
}

export const autoAdvanceOf = (p: PortfolioRow): AutoAdvanceConfig | null => {
  const a = p.auto_advance as AutoAdvanceConfig | undefined | null;
  return a?.through ? a : null;
};

/** Checks every program of [code] against its next gate and, unless [dryRun], advances those the
 * data meets. Off (nothing read, nothing written) when the portfolio declares no autoAdvance.
 * Recording needs the caller to be one of the portfolio's approvers, as a manual advance does. */
export async function autoAdvance(code: string, dryRun: boolean): Promise<AutoAdvanceResult> {
  const portfolio = await portfolioByCode(code);
  const config = autoAdvanceOf(portfolio);
  const result: AutoAdvanceResult = {portfolio: code, config, dryRun, outcomes: []};
  if (!config)
    return result;
  if (portfolio.stage_authority !== 'portal')
    throw new Error(`${code}: autoAdvance needs stage authority "portal"; the registry owns these stages`);
  if (!dryRun && !await isApprover(portfolio))
    throw new Error(`${code}: only the portfolio's approvers may let data advance its programs; preview with dryRun`);
  const programs = (await portalDb.programs.query({filter: `portfolio_id = "${portfolio.id}"`}))
    .sort((a, b) => a.code.localeCompare(b.code));
  for (const p of programs) {
    const status = (p.status ?? 'active') as ProgramStatus;
    // the cheap refusals first: an ended or unstaged program never reads its data
    const quick = status !== 'active' || !p.stage_id ?
      planAutoAdvance({program: {code: p.code, status, stageId: p.stage_id ?? null}, stages: [], config, readiness: null}) : null;
    if (quick) {
      result.outcomes.push(quick);
      continue;
    }
    const ctx = await programContext(p.id);
    const outcome = planAutoAdvance({program: {code: p.code, status, stageId: p.stage_id ?? null},
      stages: ctx.stages.map((s) => ({id: s.id, code: s.code, name: s.name, ordinal: s.ordinal})), config,
      readiness: ctx.readiness});
    // where an advance needs approval, data proposes it and an approver approves (F-30)
    if (outcome.to && needsApproval(approvalOf(portfolio.approval), 'decision', 'advance')) {
      const waiting = (await requestsOf(p.id).catch(() => [])).some((r) => r.kind === 'decision' && r.status === 'pending');
      if (dryRun || waiting) {
        result.outcomes.push({...outcome, proposal: waiting ? 'pending' : 'would'});
        continue;
      }
      try {
        await proposeDecision(p.id, 'advance',
          `Advanced by data: every required criterion of ${outcome.to.name} met from the program's data`,
          {source: 'data', evidence: outcome.evidence, toStageId: outcome.to.id});
        result.outcomes.push({...outcome, proposal: 'new'});
      } catch (e: any) {
        result.outcomes.push({...outcome, to: null, reason: `not proposed: ${e?.message ?? e}`, error: e?.message ?? String(e)});
      }
      continue;
    }
    if (outcome.to && !dryRun) {
      try {
        const r = await recordDecision(p.id, 'advance', outcome.to.id,
          `Advanced by data: every required criterion of ${outcome.to.name} met from the program's data`,
          {source: 'data', evidence: outcome.evidence});
        result.outcomes.push({...outcome, workspace: r.workspace});
      } catch (e: any) {
        result.outcomes.push({...outcome, to: null, reason: `not recorded: ${e?.message ?? e}`, error: e?.message ?? String(e)});
      }
    } else
      result.outcomes.push(outcome);
  }
  return result;
}

export const summarizeAutoAdvanceResult = (r: AutoAdvanceResult): string =>
  summarizeAutoAdvance(r.portfolio, r.config, r.outcomes, r.dryRun);
