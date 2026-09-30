/** A portfolio as its configuration file declares it (files/portfolios/*.json). Pure data and
 * pure checks: nothing here touches the platform, so the rules hold wherever the file is read. */

export type Setting = 'discovery' | 'development' | 'rollout';
export type StageAuthority = 'feed' | 'portal';
export type CriterionKind = 'document' | 'data_condition' | 'decision';
export type ProgramStatus = 'active' | 'on_hold' | 'terminated' | 'completed' | 'transitioned' | 'missing_from_feed';

export const PROGRAM_STATUSES: ProgramStatus[] =
  ['active', 'on_hold', 'terminated', 'completed', 'transitioned', 'missing_from_feed'];

export interface CriterionConfig {
  code: string;
  name: string;
  kind: CriterionKind;
  rule?: string;
  required?: boolean;
}

export interface StageConfig {
  code: string;
  name: string;
  ordinal: number;
  phase?: string;
  isGate?: boolean;
  isTerminal?: boolean;
  governance?: string;
  aliases?: string[];
  description?: string;
  criteria?: CriterionConfig[];
}

export interface StageModelConfig {
  code: string;
  name: string;
  version?: string;
  setting: Setting;
  description?: string;
  stages: StageConfig[];
}

export interface FeedConfig {
  code: string;
  /** The source system. Two feeds naming the same system read the same key space. */
  system: string;
  connection: string;
  /** Returns one row per program: `code` and `name` are required; `stage`, `source_stage`, `status`,
   * `modality`, `target`, `therapeutic_area` are mapped when present; anything else lands in `details`. */
  sql: string;
  mapping?: {status?: Record<string, ProgramStatus>};
}

export interface PortfolioConfig {
  code: string;
  name: string;
  setting: Setting;
  modalityClass?: 'small_molecule' | 'biologic' | 'any';
  stageAuthority: StageAuthority;
  programNoun?: string;
  /** The root space the workspace builder creates this portfolio's tree under. */
  workspaceRoot?: string;
  /** Folders created under every program's folder: workstreams, never compounds (FR-26). */
  workstreams?: string[];
  description?: string;
  groups?: {viewers?: string, contributors?: string, approvers?: string};
  /** Where a program goes when it leaves this portfolio's terminal stage. */
  handoff?: {to: string, fromStage: string};
  stageModel: StageModelConfig;
  feeds: FeedConfig[];
}

export interface SourceConfig {
  system: string;
  connection: string;
  description?: string;
  /** Every key in the system: what the routing audit checks the portfolios' feeds against. */
  sql: string;
}

const CODE = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;
const PROGRAM_FIELDS = new Set(['code', 'name', 'stage', 'source_stage', 'status', 'modality', 'target',
  'therapeutic_area']);

/** Every problem with one portfolio file, as sentences; empty when it can be loaded. */
export function validatePortfolio(c: PortfolioConfig): string[] {
  const errors: string[] = [];
  const need = (value: unknown, what: string) => {
    if (value === undefined || value === null || value === '')
      errors.push(`${what} is required`);
  };
  need(c.code, 'code');
  need(c.name, 'name');
  if (c.code && !CODE.test(c.code))
    errors.push(`code "${c.code}" may hold only letters, digits, '_', '.', '-'`);
  if (!['discovery', 'development', 'rollout'].includes(c.setting))
    errors.push(`setting "${c.setting}" is not one of discovery, development, rollout`);
  if (!['feed', 'portal'].includes(c.stageAuthority))
    errors.push(`stageAuthority "${c.stageAuthority}" is not one of feed, portal`);
  if (!c.stageModel)
    errors.push('stageModel is required');
  else
    errors.push(...validateStageModel(c.stageModel));

  const feedCodes = new Set<string>();
  for (const f of c.feeds ?? []) {
    need(f.code, 'feed code');
    need(f.system, `feed ${f.code}: system`);
    need(f.connection, `feed ${f.code}: connection`);
    need(f.sql, `feed ${f.code}: sql`);
    if (feedCodes.has(f.code))
      errors.push(`feed code "${f.code}" is declared twice`);
    feedCodes.add(f.code);
    for (const [raw, status] of Object.entries(f.mapping?.status ?? {})) {
      if (!(PROGRAM_STATUSES as string[]).includes(status))
        errors.push(`feed ${f.code}: status "${raw}" maps to "${status}", which is not a program status`);
    }
  }
  if (c.workstreams !== undefined) {
    if (!c.workspaceRoot)
      errors.push('workstreams need a workspaceRoot to hang under');
    const seen = new Set<string>();
    for (const w of c.workstreams) {
      const key = normalizeSpaceName(w ?? '').toLowerCase();
      if (!key)
        errors.push(`workstream "${w}" has no letters or digits, so it cannot name a space`);
      else if (seen.has(key))
        errors.push(`workstream "${w}" is declared twice (space names ignore case and punctuation)`);
      else if (w.includes('|'))
        errors.push(`workstream "${w}" may not contain '|'`);
      seen.add(key);
    }
  }
  if (c.stageAuthority === 'feed' && (c.feeds ?? []).length === 0)
    errors.push('stageAuthority "feed" needs at least one feed');
  if (c.handoff) {
    const from = c.stageModel?.stages?.find((s) => s.code === c.handoff!.fromStage);
    if (!from)
      errors.push(`handoff.fromStage "${c.handoff.fromStage}" is not a stage of ${c.stageModel?.code}`);
    else if (!from.isTerminal)
      errors.push(`handoff.fromStage "${from.code}" must be a terminal stage`);
    if (c.handoff.to === c.code)
      errors.push('a portfolio cannot hand off to itself');
  }
  return errors;
}

export function validateStageModel(m: StageModelConfig): string[] {
  const errors: string[] = [];
  const where = `stage model ${m.code ?? '?'}`;
  if (!m.code || !CODE.test(m.code))
    errors.push(`${where}: code is required and may hold only letters, digits, '_', '.', '-'`);
  if (!m.stages?.length) {
    errors.push(`${where}: has no stages`);
    return errors;
  }
  const codes = new Map<string, string>();
  const ordinals = new Set<number>();
  const spellings = new Map<string, string>();
  for (const s of m.stages) {
    if (!s.code)
      errors.push(`${where}: a stage has no code`);
    else if (codes.has(s.code))
      errors.push(`${where}: stage code "${s.code}" is declared twice`);
    codes.set(s.code, s.code);
    if (!Number.isInteger(s.ordinal) || s.ordinal < 1)
      errors.push(`${where}: stage ${s.code} has ordinal ${s.ordinal}; ordinals are positive integers`);
    else if (ordinals.has(s.ordinal))
      errors.push(`${where}: ordinal ${s.ordinal} is used by more than one stage`);
    ordinals.add(s.ordinal);
    // a feed spelling must resolve to exactly one stage
    for (const spelling of [s.code, s.name, ...(s.aliases ?? [])]) {
      const key = normalizeStageSpelling(spelling ?? '');
      const owner = spellings.get(key);
      if (key && owner && owner !== s.code)
        errors.push(`${where}: "${spelling}" would resolve to both ${owner} and ${s.code}`);
      spellings.set(key, s.code);
    }
    const criteria = new Set<string>();
    for (const cr of s.criteria ?? []) {
      if (criteria.has(cr.code))
        errors.push(`${where}: stage ${s.code} declares criterion "${cr.code}" twice`);
      criteria.add(cr.code);
      if (!['document', 'data_condition', 'decision'].includes(cr.kind))
        errors.push(`${where}: criterion ${s.code}/${cr.code} has kind "${cr.kind}"`);
      if (cr.kind === 'data_condition' && !cr.rule)
        errors.push(`${where}: data condition ${s.code}/${cr.code} has no rule`);
    }
  }
  return errors;
}

/** Problems across the whole set of portfolios one instance runs. */
export function validatePortfolioSet(portfolios: PortfolioConfig[]): string[] {
  const errors: string[] = [];
  const codes = new Set<string>();
  const models = new Map<string, string>();
  for (const p of portfolios) {
    errors.push(...validatePortfolio(p).map((e) => `${p.code}: ${e}`));
    if (codes.has(p.code))
      errors.push(`portfolio code "${p.code}" is used twice`);
    codes.add(p.code);
    // two portfolios may share a stage model, but then they must declare the same one
    const m = p.stageModel;
    if (m) {
      const fingerprint = JSON.stringify(m);
      const seen = models.get(m.code);
      if (seen !== undefined && seen !== fingerprint)
        errors.push(`stage model ${m.code} is declared differently by two portfolios; bump its code to version it`);
      models.set(m.code, fingerprint);
    }
  }
  for (const p of portfolios) {
    if (p.handoff && !codes.has(p.handoff.to))
      errors.push(`${p.code}: hands off to "${p.handoff.to}", which is not configured`);
  }
  return errors;
}

/** How a feed spelling and a stage are compared: case, spaces and punctuation do not matter. */
export function normalizeStageSpelling(s: string): string {
  return s.toLowerCase().replace(/[^a-z0-9]/g, '');
}

export function isProgramField(column: string): boolean {
  return PROGRAM_FIELDS.has(column);
}

/** The form the server stores a space name under: every character but letters and digits dropped,
 * case kept. Two names with the same form cannot be siblings. */
export function normalizeSpaceName(name: string): string {
  return String(name ?? '').replace(/[^A-Za-z0-9]/g, '');
}
