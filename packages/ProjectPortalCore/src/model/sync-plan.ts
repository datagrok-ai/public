/** Turning one portfolio's feed rows into a plan against what the portal already holds.
 * Pure: the caller reads the rows and the programs, and applies the plan. The plan is what
 * a dry run prints and what a test asserts on. */
import {isProgramField, normalizeStageSpelling, PROGRAM_STATUSES, ProgramStatus, StageAuthority} from './config';

export interface StageRef {
  id: string;
  code: string;
  name: string;
  ordinal: number;
  aliases?: string[];
}

export interface ProgramRef {
  id: string;
  code: string;
  portfolioCode: string;
  stageId?: string | null;
  status?: ProgramStatus | null;
  source?: 'feed' | 'manual' | null;
  externalRef?: string | null;
  predecessorId?: string | null;
  name?: string;
  modality?: string | null;
  target?: string | null;
  therapeuticArea?: string | null;
  sourceStage?: string | null;
}

export type FeedRow = Record<string, unknown>;

/** The program columns a feed may set. */
export interface ProgramValues {
  code: string;
  name: string;
  status: ProgramStatus;
  modality: string | null;
  target: string | null;
  therapeutic_area: string | null;
  source_stage: string | null;
  stage_id: string | null;
  external_ref: string;
  predecessor_id?: string | null;
  details: Record<string, unknown> | null;
}

export interface StageChange {
  code: string;
  programId: string | null;
  fromStageId: string | null;
  toStageId: string | null;
}

export interface SyncPlan {
  portfolio: string;
  system: string;
  inserts: ProgramValues[];
  updates: {id: string, code: string, changes: Partial<ProgramValues>}[];
  unchanged: string[];
  /** Stage moves the feed made; recorded as `sync` transitions (feed authority only). */
  stageChanges: StageChange[];
  /** Feed spellings no stage of the portfolio's model answers to. */
  unresolvedStages: {code: string, spelling: string}[];
  /** Programs this feed used to deliver and no longer does. */
  missing: {id: string, code: string, status: ProgramStatus}[];
  /** New programs that continue a program of another portfolio (same system, same key). */
  lineage: {code: string, predecessorId: string, fromPortfolio: string}[];
  /** Keys another portfolio's feed still delivers as active: the two feeds overlap. */
  overlaps: {code: string, otherPortfolio: string}[];
  /** Rows that could not be used at all. */
  rejected: {row: number, reason: string}[];
  /** Values the portfolio's mapping does not cover; the row is still used. */
  warnings: string[];
}

export interface PlanInput {
  portfolioCode: string;
  system: string;
  stageAuthority: StageAuthority;
  stages: StageRef[];
  statusMap?: Record<string, ProgramStatus>;
  /** Programs of THIS portfolio. */
  existing: ProgramRef[];
  /** Programs of OTHER portfolios whose external_ref is in the same system. */
  elsewhere: ProgramRef[];
  rows: FeedRow[];
}

export const externalRef = (system: string, code: string): string => `${system}:${code}`;

const str = (v: unknown): string | null => {
  if (v === undefined || v === null)
    return null;
  const s = String(v).trim();
  return s === '' || s.toUpperCase() === 'NULL' ? null : s;
};

/** The stage a feed spelling means, by code first, then name, then alias; null when none. */
export function resolveStage(stages: StageRef[], spelling: string | null): StageRef | null {
  if (spelling === null)
    return null;
  const key = normalizeStageSpelling(spelling);
  if (key === '')
    return null;
  return stages.find((s) => normalizeStageSpelling(s.code) === key) ??
    stages.find((s) => normalizeStageSpelling(s.name) === key) ??
    stages.find((s) => (s.aliases ?? []).some((a) => normalizeStageSpelling(a) === key)) ?? null;
}

export function planSync(input: PlanInput): SyncPlan {
  const plan: SyncPlan = {portfolio: input.portfolioCode, system: input.system, inserts: [], updates: [],
    unchanged: [], stageChanges: [], unresolvedStages: [], missing: [], lineage: [], overlaps: [], rejected: [], warnings: []};
  const byCode = new Map(input.existing.map((p) => [p.code, p]));
  const seen = new Set<string>();
  const feedsStage = input.stageAuthority === 'feed';

  input.rows.forEach((row, i) => {
    const code = str(row.code);
    if (code === null) {
      plan.rejected.push({row: i, reason: 'no code'});
      return;
    }
    if (seen.has(code)) {
      plan.rejected.push({row: i, reason: `code ${code} is delivered twice; the first row wins`});
      return;
    }
    seen.add(code);

    const rawStatus = str(row.status);
    let status: ProgramStatus = 'active';
    if (rawStatus !== null) {
      const mapped = input.statusMap?.[rawStatus] ?? rawStatus.toLowerCase();
      if ((PROGRAM_STATUSES as string[]).includes(mapped))
        status = mapped as ProgramStatus;
      else
        plan.warnings.push(`${code}: status "${rawStatus}" is not mapped; taken as active`);
    }
    // the feed's own stage places the program only where the feed is the stage authority;
    // otherwise it is kept verbatim as source_stage and the portal's gate log decides
    const spelling = feedsStage ? str(row.stage) ?? str(row.source_stage) : null;
    const stage = resolveStage(input.stages, spelling);
    if (spelling !== null && stage === null)
      plan.unresolvedStages.push({code, spelling});

    const details: Record<string, unknown> = {};
    for (const [k, v] of Object.entries(row)) {
      if (!isProgramField(k) && v !== null && v !== undefined && v !== '')
        details[k] = v;
    }
    const values: ProgramValues = {
      code,
      name: str(row.name) ?? code,
      status,
      modality: str(row.modality),
      target: str(row.target),
      therapeutic_area: str(row.therapeutic_area),
      source_stage: str(row.source_stage) ?? str(row.stage),
      stage_id: stage?.id ?? null,
      external_ref: externalRef(input.system, code),
      details: Object.keys(details).length > 0 ? details : null,
    };

    const current = byCode.get(code);
    if (current === undefined) {
      const other = input.elsewhere.find((p) => p.externalRef === values.external_ref);
      if (other) {
        if (other.status === 'active' || other.status === 'on_hold')
          plan.overlaps.push({code, otherPortfolio: other.portfolioCode});
        values.predecessor_id = other.id;
        plan.lineage.push({code, predecessorId: other.id, fromPortfolio: other.portfolioCode});
      }
      if (!feedsStage)
        values.stage_id = null;   // the portal places it; see the unstaged column
      plan.inserts.push(values);
      if (feedsStage && values.stage_id !== null)
        plan.stageChanges.push({code, programId: null, fromStageId: null, toStageId: values.stage_id});
      return;
    }

    const changes: Partial<ProgramValues> = {};
    const compare: [keyof ProgramValues, unknown][] = [
      ['name', current.name], ['status', current.status], ['modality', current.modality],
      ['target', current.target], ['therapeutic_area', current.therapeuticArea],
      ['source_stage', current.sourceStage], ['external_ref', current.externalRef]];
    for (const [field, was] of compare) {
      if ((values[field] ?? null) !== (was ?? null))
        (changes as any)[field] = values[field];
    }
    if (feedsStage && (values.stage_id ?? null) !== (current.stageId ?? null)) {
      // an unresolvable spelling never clears a stage the program already has
      if (values.stage_id !== null || spelling === null) {
        changes.stage_id = values.stage_id;
        plan.stageChanges.push({code, programId: current.id, fromStageId: current.stageId ?? null,
          toStageId: values.stage_id});
      }
    }
    if (Object.keys(changes).length === 0)
      plan.unchanged.push(code);
    else
      plan.updates.push({id: current.id, code, changes});
  });

  for (const p of input.existing) {
    if (p.source !== 'feed' || seen.has(p.code) || p.status === 'transitioned' || p.status === 'missing_from_feed')
      continue;
    // gone from this feed: if another portfolio now carries it, it moved on; otherwise flag it
    const successor = input.elsewhere.find((o) => o.predecessorId === p.id ||
      (p.externalRef && o.externalRef === p.externalRef));
    plan.missing.push({id: p.id, code: p.code, status: successor ? 'transitioned' : 'missing_from_feed'});
  }
  return plan;
}

export function summarizePlan(plan: SyncPlan): string {
  const parts = [
    `${plan.portfolio}: ${plan.inserts.length} new, ${plan.updates.length} changed, ${plan.unchanged.length} unchanged`,
  ];
  if (plan.stageChanges.length)
    parts.push(`${plan.stageChanges.length} stage change(s)`);
  if (plan.lineage.length)
    parts.push(`${plan.lineage.length} continued from another portfolio`);
  if (plan.missing.length)
    parts.push(`${plan.missing.length} no longer in the feed`);
  if (plan.unresolvedStages.length)
    parts.push(`${plan.unresolvedStages.length} unresolved stage spelling(s): ` +
      plan.unresolvedStages.map((u) => `${u.code}="${u.spelling}"`).join(', '));
  if (plan.overlaps.length)
    parts.push(`OVERLAP with ${[...new Set(plan.overlaps.map((o) => o.otherPortfolio))].join(', ')}: ` +
      plan.overlaps.map((o) => o.code).join(', '));
  if (plan.rejected.length)
    parts.push(`${plan.rejected.length} row(s) rejected`);
  if (plan.warnings.length)
    parts.push(`${plan.warnings.length} warning(s)`);
  return parts.join('; ');
}
