/** Choosing the run a gate criterion rests on (F-15): which runs a criterion takes, who may choose,
 * the side-by-side comparison of candidates, and the history of choices. Pure; the gate panel,
 * the dialog and selectRunForCriterion all ask the same rules. */
import {ProgramStatus} from './config';

/** What a criterion takes as evidence. Only model runs today. */
export type CriterionArtifactKind = 'model_run';
export const CRITERION_ARTIFACTS: CriterionArtifactKind[] = ['model_run'];

/** Which of a program's runs a criterion offers, by catalogue code. Each list narrows; an empty
 * or absent list does not. */
export interface RunFilter {
  /** Unit operation codes (`EX-PROD`). */
  unitOperations?: string[];
  /** Run type codes (`EX-SIM`). */
  runTypes?: string[];
  /** Model codes, the function's nqName (`ProjectPortalCore:FedBatchTiter`). */
  models?: string[];
}

const FILTER_KEYS = ['unitOperations', 'runTypes', 'models'] as const;

/** Problems with a criterion's `artifact` and `runFilter`, as sentences. */
export function validateRunCriterion(artifact: unknown, filter: RunFilter | undefined): string[] {
  const errors: string[] = [];
  if (artifact !== undefined && !(CRITERION_ARTIFACTS as unknown[]).includes(artifact))
    errors.push(`artifact "${artifact}" is not one of ${CRITERION_ARTIFACTS.join(', ')}`);
  if (filter === undefined)
    return errors;
  if (artifact === undefined)
    errors.push('runFilter needs artifact "model_run"');
  for (const [k, v] of Object.entries(filter)) {
    if (!(FILTER_KEYS as readonly string[]).includes(k))
      errors.push(`runFilter: "${k}" is not one of ${FILTER_KEYS.join(', ')}`);
    else if (!Array.isArray(v) || v.some((x) => typeof x !== 'string' || !x.trim()))
      errors.push(`runFilter: ${k} must be a list of codes`);
  }
  return errors;
}

/** A run as the filter sees it: its catalogue codes. */
export interface RunCodes {
  model: string;
  unitOperation: string | null;
  runType: string | null;
}

export function matchesRunFilter(run: RunCodes, filter: RunFilter | null | undefined): boolean {
  if (!filter)
    return true;
  const ok = (list: string[] | undefined, v: string | null) => !list?.length || (v !== null && list.includes(v));
  return ok(filter.unitOperations, run.unitOperation) && ok(filter.runTypes, run.runType) && ok(filter.models, run.model);
}

/** The filter in words, for the dialog: "Production bioreactor runs of type Simulated". */
export function describeRunFilter(filter: RunFilter | null | undefined, names: Record<string, string> = {}): string {
  const n = (codes?: string[]) => (codes ?? []).map((c) => names[c] ?? c).join(' or ');
  const parts: string[] = [];
  if (filter?.models?.length)
    parts.push(`of ${n(filter.models)}`);
  if (filter?.unitOperations?.length)
    parts.push(`on ${n(filter.unitOperations)}`);
  if (filter?.runTypes?.length)
    parts.push(`of type ${n(filter.runTypes)}`);
  return parts.length ? `Runs ${parts.join(', ')}` : 'Every run saved to the program';
}

export interface RunSelectionInput {
  isApprover: boolean;
  isContributor: boolean;
  programStatus: ProgramStatus;
  rationale?: string;
}

/** Why the current user may not select a run, empty when they may. Selecting a run is evidence,
 * not a stage decision: contributors may do it as well as approvers, also where the registry owns
 * the stage. A waived criterion stays waived; the run is still recorded. */
export function checkRunSelection(i: RunSelectionInput): string[] {
  const reasons: string[] = [];
  if (!i.isApprover && !i.isContributor)
    reasons.push('Only the portfolio\'s contributors and approvers select runs');
  if (['terminated', 'completed', 'transitioned'].includes(i.programStatus))
    reasons.push(`The program is ${i.programStatus}`);
  if (i.rationale !== undefined && !i.rationale.trim())
    reasons.push('Selecting a run needs a rationale: why this run over the others');
  return reasons;
}

/** A run's values for the comparison. */
export interface RunValues {
  title: string;
  inputs: Record<string, unknown>;
  outputs: Record<string, unknown>;
}

export interface ParamLabel {
  name: string;
  caption?: string | null;
  units?: string | null;
}

export interface ComparisonRow {
  section: 'input' | 'output';
  name: string;
  /** Caption with units: "Feed rate (% v/v per day)". */
  label: string;
  /** One per run, in the order given. */
  values: string[];
  /** Whether the runs disagree on it. */
  differs: boolean;
}

/** A value as the comparison shows it: numbers to 4 significant digits, a table by its size. */
export function formatRunValue(v: unknown): string {
  if (v === null || v === undefined || (typeof v === 'number' && Number.isNaN(v)))
    return '—';
  if (typeof v === 'number')
    return Number.isInteger(v) ? String(v) : String(Number(v.toPrecision(4)));
  const df = v as {rowCount?: number, columns?: {length: number}};
  if (typeof df === 'object' && typeof df.rowCount === 'number')
    return `table, ${df.rowCount} rows`;
  if (v instanceof Date)
    return v.toISOString().slice(0, 10);
  return String(v);
}

/** Inputs, then outputs, one row per parameter in the declared order; parameters a run has but
 * the declaration lacks follow. Labels come from [params] when given. */
export function compareRuns(runs: RunValues[], params: {inputs?: ParamLabel[], outputs?: ParamLabel[]} = {}):
  ComparisonRow[] {
  const rows: ComparisonRow[] = [];
  for (const section of ['input', 'output'] as const) {
    const declared = (section === 'input' ? params.inputs : params.outputs) ?? [];
    const names = [...declared.map((p) => p.name)];
    for (const r of runs) {
      for (const k of Object.keys(section === 'input' ? r.inputs : r.outputs))
        if (!names.includes(k))
          names.push(k);
    }
    for (const name of names) {
      const p = declared.find((x) => x.name === name);
      const values = runs.map((r) => formatRunValue((section === 'input' ? r.inputs : r.outputs)[name]));
      rows.push({section, name, label: `${p?.caption || name}${p?.units ? ` (${p.units})` : ''}`, values,
        differs: new Set(values).size > 1});
    }
  }
  return rows;
}

/** A recorded choice: an artifact row with its run's title and who chose it. */
export interface SelectionEntry {
  runTitle: string;
  role: 'candidate' | 'selected' | 'superseded';
  selectedOn: string | null;
  selectedBy: string | null;
  rationale: string | null;
}

/** The current choice and the ones it superseded, newest first. Candidates never chosen are left out. */
export function selectionHistory(entries: SelectionEntry[]): {selected: SelectionEntry | null, superseded: SelectionEntry[]} {
  const byDate = (a: SelectionEntry, b: SelectionEntry) => (b.selectedOn ?? '').localeCompare(a.selectedOn ?? '');
  return {
    selected: entries.find((e) => e.role === 'selected') ?? null,
    superseded: entries.filter((e) => e.role === 'superseded').sort(byDate),
  };
}

export function describeSelection(e: SelectionEntry): string {
  const when = e.selectedOn ? e.selectedOn.slice(0, 10) : '';
  return [e.runTitle, [e.selectedBy, when].filter((x) => x).join(', '), e.rationale].filter((x) => x).join(' · ');
}
