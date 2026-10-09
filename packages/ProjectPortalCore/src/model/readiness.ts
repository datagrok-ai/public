/** Readiness: computed, never typed (FR-13). One function over every criterion of a stage:
 * documents and decisions read their recorded status; data conditions with a structured rule are
 * evaluated over the program's data, and a recorded waiver still wins. Pure, so the gate panel,
 * the readiness report and the decision check all see the same answer. */
import {CriterionState} from './gates';
import {ConditionResult, DataCondition, DataRow, evaluateCondition} from './conditions';

export interface ReadinessCriterion {
  id: string;
  code: string;
  name: string;
  kind: string;
  required: boolean;
  condition?: DataCondition | null;
  /** Free-text rule, shown when there is no structured condition. */
  rule?: string | null;
}

export type ReadinessState = 'met' | 'waived' | 'in_progress' | 'not_met' | 'missing';

export interface CriterionReadiness {
  id: string;
  code: string;
  name: string;
  kind: string;
  required: boolean;
  state: ReadinessState;
  /** `evaluated`: computed from data now; `recorded`: read from criterion_status. */
  basis: 'evaluated' | 'recorded';
  result?: ConditionResult;
  /** What is still needed, in one sentence; empty when met or waived. */
  gap: string;
}

export interface Readiness {
  stageId: string | null;
  criteria: CriterionReadiness[];
  /** Required criteria met or waived, over required criteria. */
  requiredMet: number;
  requiredTotal: number;
  /** 0..1; null when the gate has no required criteria ([empty]). */
  score: number | null;
  /** The gate has no required criteria: nothing to measure, so it is not ready and never 100%.
   * Shown as "no criteria"; auto-advance never moves a program through it. */
  empty: boolean;
  /** The open required criteria, one sentence each, in declaration order. */
  gaps: string[];
  /** Criterion states for the gate rules: evaluated data conditions folded in. */
  states: Record<string, CriterionState>;
  /** Whether the program dataset could be read; when not, evaluated criteria fall back to recorded. */
  dataError?: string;
}

const MAX_GAP_ROWS = 5;

function conditionGap(r: ConditionResult, c: DataCondition): string {
  if (r.met)
    return '';
  if (r.inScope === 0)
    return c.scope === 'lead_set' ? 'no lead-set data yet' : 'no data for this program yet';
  if (r.gaps.length > 0) {
    const shown = r.gaps.slice(0, MAX_GAP_ROWS).map((g) => `${g.id} lacks ${g.missing.join(', ')}`).join('; ');
    const more = r.gaps.length > MAX_GAP_ROWS ? `; ${r.gaps.length - MAX_GAP_ROWS} more` : '';
    return `${r.summary}: ${shown}${more}`;
  }
  return r.summary;
}

export function computeReadiness(stageId: string | null, criteria: ReadinessCriterion[],
  recorded: Record<string, CriterionState>, rows: DataRow[] | null, leadRoles?: string[], dataError?: string): Readiness {
  const out: CriterionReadiness[] = [];
  const states: Record<string, CriterionState> = {...recorded};
  for (const c of criteria) {
    const stored = recorded[c.id] ?? 'not_started';
    const base = {id: c.id, code: c.code, name: c.name, kind: c.kind, required: c.required};
    if (stored === 'waived') {
      out.push({...base, state: 'waived', basis: 'recorded', gap: ''});
      continue;
    }
    if (c.kind === 'data_condition' && c.condition && rows !== null) {
      const result = evaluateCondition(c.condition, rows, leadRoles);
      states[c.id] = result.met ? 'met' : 'not_met';
      out.push({...base, state: result.met ? 'met' : 'not_met', basis: 'evaluated', result,
        gap: conditionGap(result, c.condition)});
      continue;
    }
    const state: ReadinessState = stored === 'not_started' ? 'missing' : stored;
    const gap = state === 'met' ? '' :
      c.kind === 'document' ? (state === 'in_progress' ? 'document in progress' : 'document not provided') :
      c.kind === 'decision' ? 'decision not recorded' :
      c.condition ? `program data unavailable${dataError ? ` (${dataError})` : ''}; not recorded either` :
      `not recorded${c.rule ? ` (rule: ${c.rule})` : ''}`;
    out.push({...base, state, basis: 'recorded', gap});
  }
  const required = out.filter((c) => c.required);
  const requiredMet = required.filter((c) => c.state === 'met' || c.state === 'waived').length;
  return {
    stageId, criteria: out, requiredMet, requiredTotal: required.length,
    score: required.length === 0 ? null : requiredMet / required.length, empty: required.length === 0,
    gaps: required.filter((c) => c.gap).map((c) => `${c.name}: ${c.gap}`),
    states, dataError,
  };
}

/** Ready: the gate has required criteria and every one is met or waived. An empty gate is not. */
export const isReady = (r: Pick<Readiness, 'score'> | null | undefined): boolean => !!r && r.score === 1;

/** How an empty gate is shown, wherever a percentage would be. */
export const NO_CRITERIA = 'no criteria';

/** `75%`, or `no criteria` for a gate without required criteria. */
export const readinessLabel = (r: Pick<Readiness, 'score'>): string =>
  r.score === null ? NO_CRITERIA : `${Math.round(r.score * 100)}%`;

export function summarizeReadiness(programCode: string, stageName: string, r: Readiness): string {
  const head = `${programCode} → ${stageName}: ` + (r.empty ?
    `no ${r.criteria.length ? 'required ' : ''}criteria defined; not ready` :
    `${r.requiredMet} of ${r.requiredTotal} required criteria met or waived (${readinessLabel(r)})`);
  const lines = [head];
  for (const c of r.criteria) {
    const detail = c.result ? ` — ${c.result.summary}` : '';
    lines.push(`  ${c.state.padEnd(11)} ${c.code}${c.required ? '' : ' (optional)'} [${c.basis}]${detail}`);
  }
  if (r.gaps.length)
    lines.push('  Gaps:', ...r.gaps.map((g) => `    - ${g}`));
  if (r.dataError)
    lines.push(`  Program data could not be read: ${r.dataError}`);
  return lines.join('\n');
}
