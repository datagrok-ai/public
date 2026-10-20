/** The save-to-program dialog's rules (F-14, FR-53): what it offers, what it fills in before the
 * user touches it, and what it refuses. Pure, so the dialog's logic is tested without a dialog. */
import {RUN_PURPOSES, RunAnnotations, RunPurpose} from './model-runs';

export interface Choice {
  code: string;
  name: string;
}

/** What the dialog can offer for one call. */
export interface SaveRunChoices {
  /** The portfolios the user can see. */
  portfolios: Choice[];
  /** The portfolios whose catalogue declares the call's model, in the catalogue's order. */
  modelPortfolios: string[];
  unitOperations: Choice[];
  /** The model's unit operations; the first is the default. */
  modelUnitOperations: string[];
  runTypes: Choice[];
  defaultRunType?: string | null;
  /** What the call already carries, when it was saved to a program before. */
  annotations?: RunAnnotations | null;
  /** Where the dialog was opened: a program page fixes the portfolio and the program. */
  context?: {portfolio?: string | null, program?: string | null} | null;
}

export interface SaveRunForm {
  portfolio: string;
  program: string;
  unitOperation: string;
  runType: string;
  purpose: RunPurpose | '';
  /** Empty: saveModelRun's default, `<model> · <program>`. */
  title: string;
  description: string;
}

/** The form before the user touches it. The portfolio comes from where the dialog was opened,
 * else from the call's own annotations, else from the model's catalogue (the first of its
 * portfolios the user can see), else from the only portfolio there is. A code the lists do not
 * hold is never filled in. */
export function saveRunDefaults(c: SaveRunChoices): SaveRunForm {
  const has = (list: Choice[], code?: string | null): code is string => !!code && list.some((x) => x.code === code);
  const a = c.annotations ?? null;
  const portfolio = [c.context?.portfolio, a?.portfolio, ...c.modelPortfolios].find((p) => has(c.portfolios, p)) ??
    (c.portfolios.length === 1 ? c.portfolios[0].code : '');
  const program = c.context?.program && c.context.portfolio === portfolio ? c.context.program :
    a?.program && a.portfolio === portfolio ? a.program : '';
  const unitOperation = [a?.unitOperation, ...c.modelUnitOperations].find((u) => has(c.unitOperations, u)) ?? '';
  const runType = [a?.runType, c.defaultRunType].find((r) => has(c.runTypes, r)) ?? '';
  const purpose = a?.purpose && (RUN_PURPOSES as readonly string[]).includes(a.purpose) ? a.purpose : '';
  return {portfolio, program, unitOperation, runType, purpose, title: a?.title ?? '', description: a?.description ?? ''};
}

/** The programs a portfolio offers, as the dialog lists them: the ones still running first, then
 * by code. Only what the user can see is passed in; row security decides that. */
export function programChoices(programs: {code: string, name: string, status?: string | null}[]): Choice[] {
  const ended = (s?: string | null) => ['terminated', 'completed', 'transitioned'].includes(s ?? '') ? 1 : 0;
  return [...programs].sort((a, b) => ended(a.status) - ended(b.status) || a.code.localeCompare(b.code))
    .map((p) => ({code: p.code, name: p.name === p.code ? p.code : `${p.code} · ${p.name}`}));
}

/** Why the form cannot be saved yet; empty when it can. */
export function validateSaveRun(f: SaveRunForm, programs: Choice[]): string[] {
  const errors: string[] = [];
  if (!f.portfolio)
    errors.push('Choose a portfolio');
  if (!f.program)
    errors.push('Choose a program');
  else if (!programs.some((p) => p.code === f.program))
    errors.push(`${f.program} is not a program of ${f.portfolio || 'the portfolio'} you can see`);
  if (f.purpose && !(RUN_PURPOSES as readonly string[]).includes(f.purpose))
    errors.push(`"${f.purpose}" is not a purpose`);
  return errors;
}

/** One line saying where the run went, for the dialog's result and the balloon. */
export function describeSavedRun(program: string, stage: string | null | undefined,
  folder: {path: string} | null, notes: string[]): string {
  const where = folder ? `linked into ${folder.path}` : 'not linked into a workspace';
  const already = notes.some((n) => n.startsWith('already recorded'));
  return already ? `Already saved to ${program}; nothing changed` :
    `Saved to ${program}${stage ? ` at ${stage}` : ''}, ${where}`;
}
