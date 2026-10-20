/** Stage advance driven by data (opt-in, off by default). A portfolio whose stages the portal owns
 * may let a program move one stage when the data alone meets the next gate, up to a limit stage, so
 * the top gates still need an approver. Pure: the sync and `autoAdvancePortfolio` ask the same rules,
 * over the readiness the gate panel shows. */
import type {ProgramStatus, StageAuthority} from './config';
import type {Readiness} from './readiness';

export interface AutoAdvanceConfig {
  /** The last stage a program may enter by data. Entering a later stage takes an approver. */
  through: string;
  /** Whether a waived criterion counts as met. Off by default: a waiver is a person's call. */
  countWaivers?: boolean;
}

/** Every problem with a portfolio's autoAdvance, as sentences. */
export function validateAutoAdvance(a: AutoAdvanceConfig | undefined, stageAuthority: StageAuthority,
  stageCodes: string[]): string[] {
  if (a === undefined || a === null)
    return [];
  const errors: string[] = [];
  if (stageAuthority !== 'portal')
    errors.push('autoAdvance needs stageAuthority "portal": where the registry owns the stage, the sync already moves it');
  if (!a.through)
    errors.push('autoAdvance: through is required (the last stage a program may enter by data)');
  else if (!stageCodes.includes(a.through))
    errors.push(`autoAdvance: through "${a.through}" is not a stage of the portfolio's model`);
  if (a.countWaivers !== undefined && typeof a.countWaivers !== 'boolean')
    errors.push('autoAdvance: countWaivers is true or false');
  return errors;
}

export interface AutoAdvanceStage {
  id: string;
  code: string;
  name: string;
  ordinal: number;
}

export interface AutoAdvanceInput {
  program: {code: string, status: ProgramStatus | null, stageId: string | null};
  /** The program's model, any order. */
  stages: AutoAdvanceStage[];
  config: AutoAdvanceConfig;
  /** Readiness for the program's next stage, as programContext computes it; null at the last stage. */
  readiness: Readiness | null;
}

export interface AutoAdvanceOutcome {
  program: string;
  from: AutoAdvanceStage | null;
  /** The stage the program moves to; null when it stays. */
  to: AutoAdvanceStage | null;
  /** Why it stays, in one sentence; empty when it moves. */
  reason: string;
  /** What the data said, a line per required criterion; set when it moves. */
  evidence: string;
  /** In a portfolio that asks for approval of an advance (F-30) the move is proposed, not made:
   * `new` proposed in this pass, `pending` already waiting, `would` in a dry run. */
  proposal?: 'new' | 'pending' | 'would';
}

/** Whether one program moves by data, and to where. A program moves one stage at most, and only
 * when it is active, already in its model, the next stage is not past `through`, the stage has
 * required criteria, and every one of them was evaluated met from the data (or waived, when the
 * portfolio counts waivers). A criterion recorded met by hand is a person's call, not data. */
export function planAutoAdvance(input: AutoAdvanceInput): AutoAdvanceOutcome {
  const ordered = [...input.stages].sort((a, b) => a.ordinal - b.ordinal);
  const from = ordered.find((s) => s.id === input.program.stageId) ?? null;
  const stay = (reason: string): AutoAdvanceOutcome => ({program: input.program.code, from, to: null, reason, evidence: ''});
  const status = input.program.status ?? 'active';
  if (status !== 'active')
    return stay(`${status === 'on_hold' ? 'on hold' : status.replace(/_/g, ' ')}: only active programs advance by data`);
  if (input.program.stageId === null)
    return stay('unstaged: a program enters its model by a decision');
  if (from === null)
    return stay('its stage is not part of its model');
  const next = ordered[ordered.indexOf(from) + 1] ?? null;
  if (next === null)
    return stay('at the last stage of its model');
  const through = ordered.find((s) => s.code === input.config.through);
  if (!through)
    return stay(`autoAdvance.through "${input.config.through}" is not a stage of its model`);
  if (next.ordinal > through.ordinal)
    return stay(`${next.name} is past ${through.name}, the last stage data may enter: an approver decides`);
  const r = input.readiness;
  if (!r || r.stageId !== next.id)
    return stay(`readiness for ${next.name} was not computed`);
  if (r.dataError)
    return stay(`program data could not be read (${r.dataError})`);
  const required = r.criteria.filter((c) => c.required);
  if (required.length === 0)
    return stay(`${next.name} has no required criteria for data to decide`);
  const open: string[] = [];
  for (const c of required) {
    if (c.state === 'waived') {
      if (!input.config.countWaivers)
        open.push(`${c.code} is met only by waiver`);
    } else if (c.state !== 'met')
      open.push(`${c.code} is ${c.state.replace(/_/g, ' ')}`);
    else if (c.basis !== 'evaluated')
      open.push(`${c.code} is recorded met, not evaluated from data`);
  }
  if (open.length)
    return stay(`${next.name} is not met by data: ${open.join('; ')}`);
  const evidence = [`${r.requiredMet} of ${r.requiredTotal} required criteria of ${next.name} met from data:`,
    ...required.map((c) => `${c.code}: ${c.state === 'waived' ? 'waived (counted)' : c.result?.summary ?? 'met'}`)]
    .join('\n');
  return {program: input.program.code, from, to: next, reason: '', evidence};
}

/** One line per program: `EXS-GPC3: HF → H2L by data` or `EXS-KIN1 stays: unstaged: …`. */
export function describeAutoAdvance(o: AutoAdvanceOutcome): string {
  const proposal = o.proposal ? {new: ', proposed for approval', pending: ', already proposed and waiting for approval',
    would: ', to be proposed for approval'}[o.proposal] : '';
  return o.to ? `${o.program}: ${o.from?.code ?? '—'} → ${o.to.code} by data${proposal}` : `${o.program} stays: ${o.reason}`;
}

/** The report of one pass over a portfolio: the moves first, then the programs that stay. */
export function summarizeAutoAdvance(portfolio: string, config: AutoAdvanceConfig | null, outcomes: AutoAdvanceOutcome[],
  dryRun: boolean): string {
  if (!config)
    return `${portfolio}: stages advance only by decision (no autoAdvance)`;
  const moved = outcomes.filter((o) => o.to);
  const proposed = moved.filter((o) => o.proposal).length;
  const head = `${portfolio}: ${moved.length} of ${outcomes.length} program(s) ${dryRun ? 'would advance' : 'advanced'} by data` +
    ` (through ${config.through}${config.countWaivers ? ', waivers count' : ''})` +
    (proposed ? `; ${proposed} of them ${dryRun ? 'would be' : 'are'} proposals waiting for approval` : '');
  return [head, ...moved.map(describeAutoAdvance), ...outcomes.filter((o) => !o.to).map(describeAutoAdvance)]
    .map((l, i) => i === 0 ? l : `  ${l}`).join('\n');
}

/** The setting in words, for the sync history. */
export function describeAutoAdvanceConfig(config: AutoAdvanceConfig | null, throughName?: string): string {
  if (!config)
    return 'Stages advance only by decision';
  return `Stages advance by data through ${throughName ?? config.through}` +
    `; ${config.countWaivers ? 'waivers count as met' : 'waivers do not count'}; later gates need an approver`;
}
