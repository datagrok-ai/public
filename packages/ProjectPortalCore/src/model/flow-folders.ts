/** Unit-operation folders from a program's process flows: the leaf level of the workspace tree
 * (root → program → process type → `05 - ProA capture`). Optional per portfolio (`flowFolders`);
 * without it the builder makes no such folder.
 *
 * A program's flows come from two places. The `program_flow` rows a loader wrote (the steps in
 * `process_step`) win for their process type; for a process type they do not cover, the flows the
 * portfolio file declares are picked by its rules, by modality, highest priority first. The pick is
 * made in memory each build: nothing is written to `program_flow`.
 *
 * Folders, never compounds or runs (FR-26). Pure: the caller reads the rows. */
import {normalizeSpaceName} from './config';

export type FlowFolderPlacement = 'process_type' | 'program';
export const FLOW_FOLDER_PLACEMENTS: FlowFolderPlacement[] = ['process_type', 'program'];

export interface FlowStepConfig {
  step: number;
  /** A unit operation code of the catalogue. */
  unitOperation: string;
  /** The folder's name, numbered as people read the process: `05 - ProA capture`. */
  folderName: string;
  /** Steps sharing a value are parallel branches; folders all the same. */
  branchGroup?: string;
}

export interface FlowConfig {
  code: string;
  name: string;
  /** The split of flows (Drug Substance, Drug Product): a program runs one flow per process type. */
  processType: string;
  description?: string;
  steps: FlowStepConfig[];
}

export interface FlowRuleConfig {
  /** A program modality; absent means every modality. Case, spaces and punctuation do not matter. */
  modality?: string | null;
  flow: string;
  /** Higher wins within a process type. */
  priority?: number;
  note?: string;
}

export interface FlowFoldersConfig {
  /** Where a step's folder goes: under a folder named for the flow's process type (the default; a
   * workstream of that name is reused), or straight under the program's folder. */
  placement?: FlowFolderPlacement;
  /** The stage from which a program gets them; without it, every program that has a flow. */
  fromStage?: string;
  /** Flows for programs that have no `program_flow` row for a process type. */
  flows?: FlowConfig[];
  /** Which declared flow a program runs, by modality. */
  rules?: FlowRuleConfig[];
}

const CODE = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;

/** How a modality is compared: case, spaces and punctuation do not matter. */
export const normModality = (s: string | null | undefined): string => String(s ?? '').toLowerCase().replace(/[^a-z0-9]/g, '');

const spaceKey = (s: string) => normalizeSpaceName(s).toLowerCase();

/** Problems with a portfolio's `flowFolders`, as sentences. Unit operation codes are checked
 * against the catalogues of the set (validatePortfolioSet). */
export function validateFlowFolders(c: FlowFoldersConfig | undefined, stageCodes: string[], workspaceRoot?: string): string[] {
  if (c === undefined || c === null)
    return [];
  const errors: string[] = [];
  if (!workspaceRoot)
    errors.push('flowFolders need a workspaceRoot to hang under');
  if (c.placement !== undefined && !FLOW_FOLDER_PLACEMENTS.includes(c.placement))
    errors.push(`flowFolders.placement "${c.placement}" is not one of ${FLOW_FOLDER_PLACEMENTS.join(', ')}`);
  if (c.fromStage !== undefined && !stageCodes.includes(c.fromStage))
    errors.push(`flowFolders.fromStage "${c.fromStage}" is not a stage of the portfolio`);
  const flows = new Map<string, FlowConfig>();
  for (const f of c.flows ?? []) {
    const where = `flow ${f.code ?? '?'}`;
    if (!f.code || !CODE.test(f.code))
      errors.push(`flow code "${f.code ?? ''}" may hold only letters, digits, '_', '.', '-'`);
    else if (flows.has(f.code))
      errors.push(`flow ${f.code} is declared twice`);
    flows.set(f.code, f);
    if (!f.name)
      errors.push(`${where}: name is required`);
    if (!spaceKey(f.processType ?? ''))
      errors.push(`${where}: processType is required and must name a folder`);
    else if (f.processType.includes('/') || f.processType.includes('|'))
      errors.push(`${where}: processType "${f.processType}" may not contain '/' or '|'`);
    if (!f.steps?.length)
      errors.push(`${where}: has no steps`);
    const numbers = new Set<number>();
    const folders = new Set<string>();
    for (const s of f.steps ?? []) {
      if (!Number.isInteger(s.step) || s.step < 1 || numbers.has(s.step))
        errors.push(`${where}: step ${s.step} is not a positive number used once`);
      numbers.add(s.step);
      if (!s.unitOperation)
        errors.push(`${where} step ${s.step}: unitOperation is required`);
      const key = spaceKey(s.folderName ?? '');
      if (!key)
        errors.push(`${where} step ${s.step}: folderName is required and must hold a letter or digit`);
      else if (s.folderName.includes('/') || s.folderName.includes('|'))
        errors.push(`${where} step ${s.step}: folderName "${s.folderName}" may not contain '/' or '|'`);
      else if (folders.has(key))
        errors.push(`${where} step ${s.step}: folderName "${s.folderName}" is used twice (space names ignore case and punctuation)`);
      folders.add(key);
    }
  }
  const rules = new Set<string>();
  for (const r of c.rules ?? []) {
    const key = `${normModality(r.modality)}:${r.flow}`;
    if (!flows.has(r.flow))
      errors.push(`flow rule ${r.modality || 'any modality'} → ${r.flow}: ${r.flow} is not a declared flow`);
    if (r.priority !== undefined && !Number.isInteger(r.priority))
      errors.push(`flow rule ${r.modality || 'any modality'} → ${r.flow}: priority must be a whole number`);
    if (rules.has(key))
      errors.push(`flow rule ${r.modality || 'any modality'} → ${r.flow} is declared twice`);
    rules.add(key);
  }
  if ((c.flows ?? []).length > 0 && !(c.rules ?? []).length)
    errors.push('flowFolders declare flows but no rule picks them for a program');
  return errors;
}

/** A flow as one program runs it, from either source. */
export interface ProgramFlow {
  code: string;
  processType: string;
  /** `table`: a program_flow row; `portfolio`: picked by the portfolio file's rules. */
  source: 'table' | 'portfolio';
  steps: {step: number, folderName: string, unitOperation?: string | null}[];
}

/** A program_flow row joined to its flow and steps, as the caller read them. */
export interface TableFlow {
  code: string;
  processType: string;
  /** program_flow.source: a manual row overrides the rules for its process type. */
  rowSource: 'rule' | 'manual' | string;
  steps: {step: number, folderName: string, unitOperation?: string | null}[];
}

/** The flows a program of [modality] runs. Table rows first, a manual one over a rule one, per
 * process type; then, for each process type they leave open, the declared flow whose rule matches
 * with the highest priority (the first declared on a tie). */
export function flowsForProgram(modality: string | null | undefined, table: TableFlow[],
  config: FlowFoldersConfig | null | undefined): ProgramFlow[] {
  const out: ProgramFlow[] = [];
  const covered = new Set<string>();
  const byType = new Map<string, TableFlow[]>();
  for (const t of table) {
    const k = spaceKey(t.processType);
    byType.set(k, [...byType.get(k) ?? [], t]);
  }
  for (const [k, ts] of byType) {
    const manual = ts.filter((t) => t.rowSource === 'manual');
    for (const t of (manual.length ? manual : ts).sort((a, b) => a.code.localeCompare(b.code)))
      out.push({code: t.code, processType: t.processType, source: 'table', steps: t.steps});
    covered.add(k);
  }
  const flows = new Map((config?.flows ?? []).map((f) => [f.code, f]));
  const m = normModality(modality);
  const picked = new Map<string, {flow: FlowConfig, priority: number}>();
  for (const r of config?.rules ?? []) {
    const flow = flows.get(r.flow);
    if (!flow || (r.modality && normModality(r.modality) !== m))
      continue;
    const k = spaceKey(flow.processType);
    if (covered.has(k))
      continue;
    const priority = r.priority ?? 0;
    const was = picked.get(k);
    if (!was || priority > was.priority)
      picked.set(k, {flow, priority});
  }
  for (const {flow} of picked.values()) {
    out.push({code: flow.code, processType: flow.processType, source: 'portfolio',
      steps: flow.steps.map((s) => ({step: s.step, folderName: s.folderName, unitOperation: s.unitOperation}))});
  }
  return out.sort((a, b) => a.processType.localeCompare(b.processType) || a.code.localeCompare(b.code));
}

/** Whether a program at [ordinal] has reached [fromOrdinal]; no fromStage means every program. */
export function flowFoldersDue(fromOrdinal: number | null | undefined, ordinal: number | null): boolean {
  if (fromOrdinal === null || fromOrdinal === undefined)
    return true;
  return ordinal !== null && ordinal >= fromOrdinal;
}

/** The folder paths under a program's folder for [flows], in step order within each flow. A folder
 * two flows share is listed once. */
export function flowFolderPaths(flows: ProgramFlow[], placement: FlowFolderPlacement = 'process_type'): string[][] {
  const out: string[][] = [];
  const seen = new Set<string>();
  for (const f of flows) {
    for (const s of [...f.steps].sort((a, b) => a.step - b.step)) {
      const parts = placement === 'program' ? [s.folderName.trim()] : [f.processType.trim(), s.folderName.trim()];
      const key = parts.map(spaceKey).join('/');
      if (!seen.has(key) && parts.every((p) => spaceKey(p))) {
        seen.add(key);
        out.push(parts);
      }
    }
  }
  return out;
}

/** Problems across a set: a step naming a unit operation that no catalogue of the set declares. */
export function validateFlowSet(portfolios: {code: string, flowFolders?: FlowFoldersConfig,
  unitOperations: string[]}[]): string[] {
  const errors: string[] = [];
  const known = new Set(portfolios.flatMap((p) => p.unitOperations));
  for (const p of portfolios) {
    for (const f of p.flowFolders?.flows ?? []) {
      for (const s of f.steps ?? []) {
        if (s.unitOperation && !known.has(s.unitOperation))
          errors.push(`${p.code}: flow ${f.code} step ${s.step}: unit operation ${s.unitOperation} is not in the catalogue of any portfolio of the set`);
      }
    }
  }
  return errors;
}
