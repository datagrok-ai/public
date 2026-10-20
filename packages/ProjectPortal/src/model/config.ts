/** A portfolio as its configuration file declares it (files/portfolios/*.json). Pure data and
 * pure checks: nothing here touches the platform, so the rules hold wherever the file is read. */
import {DataCondition, ProgramDataConfig, validateCondition} from './conditions';
import {validateWorkUnits, WorkUnitFeedConfig, WorkUnitKind} from './work-units';
import {validateStageFolders} from './workspace-plan';
import {CompoundFeedConfig, validateCompoundFeed} from './compounds';
import {mergeCatalogues, ModelCatalogueConfig} from './catalogue';
import {CriterionArtifactKind, RunFilter, validateRunCriterion} from './run-selection';
import {MAIN_DATASET, ProgramDatasetConfig, validateProgramDatasets} from './program-data';
import {DmtaConfig, validateDmta} from './dmta';
import {AskSmeConfig, validateAskSme} from './discussion';
import {AttentionRule, validateAttentionRules} from './attention';
import {validateSchedule} from './sync-runs';
import {AutoAdvanceConfig, validateAutoAdvance} from './auto-advance';
import {PersonaConfig, validatePersonas} from './personas';
import {DocumentKindConfig, validateDocumentKinds} from './documents';
import {ApprovalConfig, validateApproval} from './approvals';
import {FlowFoldersConfig, validateFlowFolders, validateFlowSet} from './flow-folders';

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
  /** The rule as a person reads it. */
  rule?: string;
  /** The rule as the portal evaluates it over the program dataset; data conditions only. Without
   * it, a data condition is recorded met or waived by hand. */
  condition?: DataCondition;
  required?: boolean;
  /** Applies only to programs of this modality; empty means every modality. */
  modality?: string;
  /** What the criterion rests on: `model_run` means a run saved to the program is chosen for it
   * from the gate panel. */
  artifact?: CriterionArtifactKind;
  /** Which of the program's runs it offers, by unit operation, run type or model code. */
  runFilter?: RunFilter;
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
  /** Folders added under a program's folder once it reaches this stage, '/'-separated paths
   * (Candidate Package, DMPK/In vivo PK). Never compounds (FR-26). */
  folders?: string[];
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
  /** Unit-operation folders from each program's process flows, under the program's folder: the
   * flows and rules the file declares, and those a loader wrote to program_flow. None when absent. */
  flowFolders?: FlowFoldersConfig;
  /** Where the DMTA status (F-24) reads tested from, and functions that add events. Defaults when absent. */
  dmta?: DmtaConfig;
  /** Work units the portfolio uses (FR-06); none when absent. */
  workUnits?: WorkUnitKind[];
  /** Feeds that deliver work units from a source system, keyed by program code. */
  workUnitFeeds?: WorkUnitFeedConfig[];
  /** The ordered series maturity vocabulary, e.g. hit → validated → lead → optimized → candidate. */
  seriesMaturity?: string[];
  /** The registry feed of the compounds each program holds, with their roles. */
  compoundFeed?: CompoundFeedConfig;
  /** The program dataset query that data conditions are evaluated over. */
  programData?: ProgramDataConfig;
  /** Further datasets of every program, each a query with @program: a SAR table beside the
   * process dataset, assay results beside the compound table. The program page opens each. */
  programDatasets?: ProgramDatasetConfig[];
  /** Rules over the program dataset that flag compounds (or lots) needing attention: potent but
   * hERG not measured, a lead without solubility. Same condition language as the criteria. */
  attention?: AttentionRule[];
  /** When syncDuePortfolios syncs this portfolio: `every 6h`, or a cron line in UTC (`0 6 * * 1-5`).
   * Without it the portfolio syncs only on call. */
  schedule?: string;
  /** User logins told, in the platform's notifications, about a scheduled sync that changed something. */
  owners?: string[];
  /** Stage advance driven by data, off when absent: a program moves one stage when the data alone
   * meets the next gate, never past `through`. Portal authority only. */
  autoAdvance?: AutoAdvanceConfig;
  /** The portfolio's model catalogue: unit operations, run types, models and their cards. Checked
   * and loaded with the set, so a model may name a unit operation another portfolio declares. */
  catalogue?: ModelCatalogueConfig;
  /** Quick actions per persona (F-08): groups mapped to personas, a default persona, and the cards
   * each persona gets. Optional: without it every persona gets the core's cards. */
  personas?: PersonaConfig;
  /** The kinds of document the portfolio's programs hold (F-27). Optional: without it, one kind per
   * document criterion of the stage model, and `other`. */
  documentKinds?: DocumentKindConfig[];
  /** The optional approval workflow (F-30), off when absent: the gate decisions and artifacts that
   * a second person approves before they take effect, and how many approvers that takes. */
  approval?: ApprovalConfig;
  /** Where Ask SME questions go (F-51); the portfolio's approvers when absent. */
  askSme?: AskSmeConfig;
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
  errors.push(...validateWorkUnits(c.workUnits, c.workUnitFeeds, c.seriesMaturity));
  errors.push(...validateDmta(c.dmta, [...(c.programData ? [MAIN_DATASET] : []), ...(c.programDatasets ?? []).map((d) => d.code)]));
  errors.push(...validateFlowFolders(c.flowFolders, (c.stageModel?.stages ?? []).map((s) => s.code), c.workspaceRoot));
  errors.push(...validateCompoundFeed(c.compoundFeed));
  errors.push(...validateAskSme(c.askSme));
  if (c.programData !== undefined) {
    need(c.programData.connection, 'programData: connection');
    need(c.programData.sql, 'programData: sql');
    if (c.programData.sql && !c.programData.sql.includes('@program'))
      errors.push('programData: sql must select one program\'s rows with @program');
    if (c.programData.watermark !== undefined && (typeof c.programData.watermark !== 'string' || !c.programData.watermark.trim()))
      errors.push('programData: watermark must name a column');
  }
  errors.push(...validateProgramDatasets(c.programDatasets, c.programData?.viewers));
  const conditions = (c.stageModel?.stages ?? []).flatMap((s) => (s.criteria ?? []).filter((cr) => cr.condition));
  if (conditions.length > 0 && !c.programData)
    errors.push(`${conditions.length} criteria declare a condition but the portfolio has no programData to evaluate it over`);
  errors.push(...validateAttentionRules(c.attention));
  errors.push(...validateSchedule(c.schedule, c.owners));
  errors.push(...validatePersonas(c.personas));
  errors.push(...validateDocumentKinds(c.documentKinds));
  errors.push(...validateApproval(c.approval, c.stageAuthority));
  errors.push(...validateAutoAdvance(c.autoAdvance, c.stageAuthority, (c.stageModel?.stages ?? []).map((s) => s.code)));
  if ((c.attention ?? []).length > 0 && !c.programData)
    errors.push(`${c.attention!.length} attention rules are declared but the portfolio has no programData to evaluate them over`);
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
    errors.push(...validateStageFolders(s.code, s.folders ?? []).map((e) => `${where}: ${e}`));
    const criteria = new Set<string>();
    for (const cr of s.criteria ?? []) {
      if (criteria.has(cr.code))
        errors.push(`${where}: stage ${s.code} declares criterion "${cr.code}" twice`);
      criteria.add(cr.code);
      if (!['document', 'data_condition', 'decision'].includes(cr.kind))
        errors.push(`${where}: criterion ${s.code}/${cr.code} has kind "${cr.kind}"`);
      if (cr.kind === 'data_condition' && !cr.rule && !cr.condition)
        errors.push(`${where}: data condition ${s.code}/${cr.code} has no rule`);
      if (cr.condition && cr.kind !== 'data_condition')
        errors.push(`${where}: criterion ${s.code}/${cr.code} is a ${cr.kind}; only data conditions take a condition`);
      else if (cr.condition)
        errors.push(...validateCondition(cr.condition).map((e) => `${where}: condition ${s.code}/${cr.code}: ${e}`));
      errors.push(...validateRunCriterion(cr.artifact, cr.runFilter).map((e) => `${where}: criterion ${s.code}/${cr.code}: ${e}`));
    }
  }
  return errors;
}

/** Problems across the whole set of portfolios one instance runs. [alsoKnown] names portfolios the
 * instance already has, which a handoff may go to (a bundle imported beside them). */
export function validatePortfolioSet(portfolios: PortfolioConfig[], alsoKnown: string[] = []): string[] {
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
    if (p.handoff && !codes.has(p.handoff.to) && !alsoKnown.includes(p.handoff.to))
      errors.push(`${p.code}: hands off to "${p.handoff.to}", which is not configured`);
  }
  errors.push(...mergeCatalogues(portfolios.filter((p) => p.catalogue)
    .map((p) => ({portfolio: p.code, catalogue: p.catalogue!}))).errors);
  errors.push(...validateFlowSet(portfolios.map((p) => ({code: p.code, flowFolders: p.flowFolders,
    unitOperations: (p.catalogue?.unitOperations ?? []).map((u) => u.code)}))));
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
