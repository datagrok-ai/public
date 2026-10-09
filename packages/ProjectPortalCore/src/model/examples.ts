/** The example history (files/examples/history.json): what a registry feed cannot carry for the EX-*
 * portfolios. Gate decisions with rationale and evidence, criteria recorded met, observations and a
 * follow, all dated in the past. The feed places nothing in a portfolio whose portal owns the stage,
 * so the decisions are also what puts the example programs on their tracks. Pure: the checks and the
 * plan of what is still missing; ../examples writes it. */
import {CriterionState} from './gates';
import {DecisionOutcome, DECISION_OUTCOMES} from './decision-records';
import {ObservationSubject} from './observations';
import {PortfolioConfig} from './config';

export type ExampleDecisionKind = 'advance' | 'hold' | 'terminate' | 'record';

export interface ExampleDecision {
  portfolio: string;
  program: string;
  decision: ExampleDecisionKind;
  /** advance: the stage entered; record: the gate's stage. */
  to?: string;
  /** record only: the stage the program was at when the registry passed the gate. */
  from?: string;
  /** record only: go, hold or no_go. */
  outcome?: DecisionOutcome;
  /** YYYY-MM-DD. */
  on: string;
  rationale: string;
  links?: {name: string, url: string}[];
}

export interface ExampleCriterion {
  portfolio: string;
  program: string;
  stage: string;
  criterion: string;
  status: CriterionState;
  on: string;
  evidence: string;
}

export interface ExampleObservation {
  portfolio: string;
  program: string;
  subject: ObservationSubject;
  code?: string;
  on: string;
  title: string;
  finding: string;
}

export interface ExampleHistory {
  criteria: ExampleCriterion[];
  decisions: ExampleDecision[];
  observations: ExampleObservation[];
  follows: {portfolio: string, program: string}[];
}

export const isExamplePortfolio = (code: string): boolean => code.startsWith('EX-');

const DAY = /^\d{4}-\d{2}-\d{2}$/;

/** Every problem with [h] against the portfolio files, as sentences. Each program's advances must
 * walk its model one stage at a time from the first stage; records only where the registry owns the
 * stage; every stage and criterion must exist. */
export function validateExampleHistory(h: ExampleHistory, portfolios: PortfolioConfig[]): string[] {
  const errors: string[] = [];
  const byCode = new Map(portfolios.map((p) => [p.code, p]));
  const where = (x: {portfolio: string, program: string}) => `${x.portfolio} ${x.program}`;
  const portfolioOf = (x: {portfolio: string, program: string}): PortfolioConfig | null => {
    if (!isExamplePortfolio(x.portfolio)) {
      errors.push(`${where(x)}: examples go to EX-* portfolios only`);
      return null;
    }
    const p = byCode.get(x.portfolio) ?? null;
    if (!p)
      errors.push(`${where(x)}: no such portfolio`);
    return p;
  };
  const stage = (p: PortfolioConfig, code: string | undefined) => p.stageModel.stages.find((s) => s.code === code) ?? null;
  for (const c of h.criteria) {
    const p = portfolioOf(c);
    const s = p ? stage(p, c.stage) : null;
    if (p && !s)
      errors.push(`${where(c)}: no stage ${c.stage}`);
    else if (s && !(s.criteria ?? []).some((x) => x.code === c.criterion))
      errors.push(`${where(c)}: stage ${c.stage} has no criterion ${c.criterion}`);
    if (!DAY.test(c.on))
      errors.push(`${where(c)}: "${c.on}" is not a YYYY-MM-DD date`);
  }
  for (const x of [...h.observations, ...h.follows])
    portfolioOf(x);
  const at = new Map<string, number>();
  const sorted = [...h.decisions].sort((a, b) => a.on.localeCompare(b.on));
  for (const d of sorted) {
    const p = portfolioOf(d);
    if (!DAY.test(d.on))
      errors.push(`${where(d)}: "${d.on}" is not a YYYY-MM-DD date`);
    if (!d.rationale?.trim())
      errors.push(`${where(d)}: ${d.decision} needs a rationale`);
    if (!p)
      continue;
    const ordered = [...p.stageModel.stages].sort((a, b) => a.ordinal - b.ordinal);
    if (d.decision === 'record') {
      if (p.stageAuthority !== 'feed')
        errors.push(`${where(d)}: a record is for a portfolio whose registry owns the stage`);
      if (!DECISION_OUTCOMES.includes(d.outcome as DecisionOutcome))
        errors.push(`${where(d)}: a record needs an outcome (${DECISION_OUTCOMES.join(', ')})`);
      if (!stage(p, d.to) || !stage(p, d.from))
        errors.push(`${where(d)}: a record needs from and to stages of ${p.code}`);
      continue;
    }
    if (p.stageAuthority !== 'portal') {
      errors.push(`${where(d)}: the registry owns ${p.code}'s stages; use a record`);
      continue;
    }
    const key = where(d);
    const now = at.get(key) ?? -1;
    if (d.decision === 'advance') {
      const next = ordered[now + 1];
      if (!next || d.to !== next.code)
        errors.push(`${key}: advances one stage at a time; next is ${next?.code ?? 'none'}, not ${d.to}`);
      else
        at.set(key, now + 1);
    } else if (now < 0)
      errors.push(`${key}: ${d.decision} before the program is on its track`);
  }
  return errors;
}

/** What the platform holds now, for the programs the history names. */
export interface SeedState {
  /** By `portfolio program`; absent when the program is not there (not synced yet). */
  programs: Map<string, {id: string, stage: string | null, status: string}>;
  /** The gate log of those programs. */
  transitions: {programId: string, decision: string, rationale: string | null}[];
  /** Recorded criterion statuses: `programId stage criterion`. */
  recorded: Set<string>;
  /** Observation titles by program id. */
  observations: {programId: string, title: string}[];
  /** Program ids the current user follows. */
  follows: Set<string>;
}

export interface SeedPlan {
  criteria: ExampleCriterion[];
  decisions: ExampleDecision[];
  /** The stage and status each program is left in, by `portfolio program`, where decisions move it. */
  moves: {key: string, stage: string | null, status: string | null}[];
  observations: ExampleObservation[];
  follows: {portfolio: string, program: string}[];
  /** What is left alone, and why, one sentence each. */
  skipped: string[];
}

export const exampleKey = (x: {portfolio: string, program: string}): string => `${x.portfolio} ${x.program}`;

/** What is still missing of [h] given [s]: idempotent, so a second run plans nothing. A decision is
 * there when the program's gate log holds one of its kind with its rationale. A program whose stage
 * is not where the written part of its history leaves it has been moved by someone since: its
 * remaining decisions are skipped, never forced. A criterion already recorded is never overwritten. */
export function planExampleSeed(h: ExampleHistory, s: SeedState): SeedPlan {
  const plan: SeedPlan = {criteria: [], decisions: [], moves: [], observations: [], follows: [], skipped: []};
  const missing = (x: {portfolio: string, program: string}) => {
    if (!isExamplePortfolio(x.portfolio)) {
      plan.skipped.push(`${exampleKey(x)}: not an EX-* portfolio`);
      return null;
    }
    const p = s.programs.get(exampleKey(x));
    if (!p)
      plan.skipped.push(`${exampleKey(x)}: no such program (sync the portfolio first)`);
    return p ?? null;
  };
  const byProgram = new Map<string, ExampleDecision[]>();
  for (const d of [...h.decisions].sort((a, b) => a.on.localeCompare(b.on)))
    byProgram.set(exampleKey(d), [...byProgram.get(exampleKey(d)) ?? [], d]);
  for (const [key, list] of byProgram) {
    const p = missing(list[0]);
    if (!p)
      continue;
    const done = (d: ExampleDecision) => s.transitions.some((t) => t.programId === p.id && t.decision === d.decision &&
      (t.rationale ?? '') === d.rationale);
    const todo = list.filter((d) => !done(d));
    if (!todo.length)
      continue;
    // where the written part leaves the program, and where the whole history leaves it
    const stageAfter = (ds: ExampleDecision[]) => ds.filter((d) => d.decision === 'advance').map((d) => d.to!).pop() ?? null;
    const written = list.filter(done);
    const portalOwned = list.some((d) => d.decision !== 'record');
    if (portalOwned && p.stage !== stageAfter(written)) {
      plan.skipped.push(`${key}: at ${p.stage ?? 'no stage'}, not ${stageAfter(written) ?? 'unstaged'} as its example history ` +
        `leaves it; ${todo.length} decision(s) left alone`);
      continue;
    }
    plan.decisions.push(...todo);
    if (portalOwned) {
      const last = todo[todo.length - 1];
      plan.moves.push({key, stage: stageAfter(list),
        status: last.decision === 'hold' ? 'on_hold' : last.decision === 'terminate' ? 'terminated' : null});
    }
  }
  for (const c of h.criteria) {
    const p = missing(c);
    if (p && !s.recorded.has(`${p.id} ${c.stage} ${c.criterion}`))
      plan.criteria.push(c);
  }
  for (const o of h.observations) {
    const p = missing(o);
    if (p && !s.observations.some((x) => x.programId === p.id && x.title === o.title))
      plan.observations.push(o);
  }
  for (const f of h.follows) {
    const p = missing(f);
    if (p && !s.follows.has(p.id))
      plan.follows.push(f);
  }
  plan.decisions.sort((a, b) => a.on.localeCompare(b.on));
  return plan;
}

export function describeSeedPlan(plan: SeedPlan, dryRun: boolean): string {
  const verb = dryRun ? 'would write' : 'wrote';
  const n = plan.decisions.length + plan.criteria.length + plan.observations.length + plan.follows.length;
  const lines = [n === 0 ? 'Example history: nothing to write; it is all there' :
    `Example history: ${verb} ${plan.decisions.length} decision(s), ${plan.criteria.length} criterion status(es), ` +
    `${plan.observations.length} observation(s), ${plan.follows.length} follow(s)`];
  for (const d of plan.decisions) {
    const what = d.decision === 'record' ? `record ${d.outcome} for ${d.to}` : d.decision === 'advance' ? `advance to ${d.to}` : d.decision;
    lines.push(`  ${d.on} ${exampleKey(d)}: ${what}`);
  }
  for (const c of plan.criteria)
    lines.push(`  ${c.on} ${exampleKey(c)}: ${c.stage} ${c.criterion} ${c.status}`);
  for (const o of plan.observations)
    lines.push(`  ${o.on} ${exampleKey(o)}: observation "${o.title}"`);
  for (const f of plan.follows)
    lines.push(`  ${exampleKey(f)}: followed`);
  lines.push(...plan.skipped.map((x) => `  skipped ${x}`));
  return lines.join('\n');
}
