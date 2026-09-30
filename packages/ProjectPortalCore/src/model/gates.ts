/** Whether a stage decision may be recorded, and why not. Pure; the app asks before it
 * enables a button, and the recording function asks again before it writes. */
import {ProgramStatus, StageAuthority} from './config';

export type Decision = 'advance' | 'recycle' | 'hold' | 'terminate';
export type CriterionState = 'not_started' | 'in_progress' | 'met' | 'not_met' | 'waived';

export interface GateStage {
  id: string;
  code: string;
  name: string;
  ordinal: number;
  isTerminal?: boolean;
}

export interface GateCriterion {
  id: string;
  stageId: string;
  code: string;
  name: string;
  required: boolean;
}

export interface GateInput {
  decision: Decision;
  stageAuthority: StageAuthority;
  programStatus: ProgramStatus;
  /** The program's own model: only these stages are reachable. */
  stages: GateStage[];
  currentStageId: string | null;
  targetStageId?: string | null;
  /** Criteria of the target stage. */
  criteria: GateCriterion[];
  criterionStates: Record<string, CriterionState>;
  rationale?: string;
  /** Whether the user is in the portfolio's approvers group. */
  isApprover: boolean;
}

export interface GateVerdict {
  allowed: boolean;
  reasons: string[];
  /** The stage the decision would move to; null for hold / terminate. */
  toStageId: string | null;
  /** Required criteria still open or not met, by code. */
  outstanding: string[];
}

/** The stage after [currentStageId]; for an unstaged program, the first stage. */
export function nextStage(stages: GateStage[], currentStageId: string | null): GateStage | null {
  const ordered = [...stages].sort((a, b) => a.ordinal - b.ordinal);
  if (currentStageId === null)
    return ordered[0] ?? null;
  const i = ordered.findIndex((s) => s.id === currentStageId);
  return i < 0 ? null : ordered[i + 1] ?? null;
}

export function checkDecision(g: GateInput): GateVerdict {
  const reasons: string[] = [];
  const outstanding: string[] = [];
  const byId = new Map(g.stages.map((s) => [s.id, s]));
  const current = g.currentStageId === null ? null : byId.get(g.currentStageId) ?? null;
  let toStageId: string | null = null;

  if (g.stageAuthority === 'feed')
    reasons.push('The registry owns this portfolio\'s stages; record the decision there and sync');
  if (!g.isApprover)
    reasons.push('Only the portfolio\'s approvers record gate decisions');
  if (['terminated', 'completed', 'transitioned'].includes(g.programStatus))
    reasons.push(`The program is ${g.programStatus}`);
  if (g.currentStageId !== null && current === null)
    reasons.push('The program\'s stage is not part of its stage model');

  switch (g.decision) {
  case 'advance': {
    const next = nextStage(g.stages, g.currentStageId);
    const target = g.targetStageId ? byId.get(g.targetStageId) ?? null : next;
    if (g.targetStageId && target === null)
      reasons.push('The target stage belongs to a different stage model');
    else if (target === null)
      reasons.push(current?.isTerminal ? 'The program is at its terminal stage; hand it off instead' :
        'There is no next stage');
    else if (next !== null && target.id !== next.id)
      reasons.push(`Advancing moves one stage at a time: the next stage is ${next.name}`);
    else {
      toStageId = target.id;
      for (const c of g.criteria.filter((c) => c.stageId === target.id && c.required)) {
        const state = g.criterionStates[c.id] ?? 'not_started';
        if (state !== 'met' && state !== 'waived')
          outstanding.push(c.code);
      }
      if (outstanding.length > 0)
        reasons.push(`${outstanding.length} required criteria of ${target.name} are not met or waived`);
    }
    break;
  }
  case 'recycle': {
    const target = g.targetStageId ? byId.get(g.targetStageId) ?? null : null;
    if (target === null)
      reasons.push('Recycling needs an earlier stage of the same model');
    else if (current === null || target.ordinal >= current.ordinal)
      reasons.push('Recycling goes back to an earlier stage');
    else
      toStageId = target.id;
    if (!g.rationale?.trim())
      reasons.push('Recycling needs a rationale');
    break;
  }
  case 'hold':
  case 'terminate':
    if (!g.rationale?.trim())
      reasons.push(`${g.decision === 'hold' ? 'Holding' : 'Terminating'} needs a rationale`);
    if (g.decision === 'hold' && g.programStatus === 'on_hold')
      reasons.push('The program is already on hold');
    break;
  }
  return {allowed: reasons.length === 0, reasons, toStageId, outstanding};
}

/** The status a decision leaves the program in. */
export function statusAfter(decision: Decision, current: ProgramStatus): ProgramStatus {
  switch (decision) {
  case 'hold': return 'on_hold';
  case 'terminate': return 'terminated';
  default: return current === 'on_hold' ? 'active' : current;
  }
}
