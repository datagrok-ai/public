// What a stage's form holds for each study table (the endpoint, the control group, the correction, alpha, the model,
// the plan…): kept as the user leaves the table (another page, Back, another study picked) and given back when the
// page shows that table again, so Back and Forward bring a page back exactly (QA 0.1.3 round 3 chaos 2). Kept in
// memory only: what is saved with the study is saved by the stage (its tags), not here.
import * as DG from 'datagrok-api/dg';

/** A form's values worth keeping: plain values, a column by its name. */
export type Kept = {[name: string]: string | number | boolean | null};

/** The values of `values` worth keeping: plain ones, a column by its name; tables (the study, an exposure table)
 *  and the names in `skip` left out. */
export function keptOf(values: {[name: string]: unknown}, skip: string[] = ['study']): Kept {
  const kept: Kept = {};
  for (const [name, v] of Object.entries(values)) {
    if (skip.includes(name))
      continue;
    if (v === null || v === undefined || ['string', 'number', 'boolean'].includes(typeof v))
      kept[name] = (v ?? null) as Kept[string];
    else if (v instanceof DG.Column)
      kept[name] = v.name;
  }
  return kept;
}

/** The call's values worth keeping. */
export const keptOfCall = (call: DG.FuncCall): Kept =>
  keptOf(Object.fromEntries(call.func.inputs.map((p) => [p.name, call.getParamValue(p.name)])));

/** Gives `kept` back to the form's fields (one it no longer has is skipped; one already so is not set again). */
export function applyKept(form: {getInput(name: string): {value: unknown} | null | undefined}, kept: Kept | null):
  void {
  for (const [name, v] of Object.entries(kept ?? {})) {
    const input = form.getInput(name);
    const now = input?.value instanceof DG.Column ? input.value.name : input?.value;
    if (input && (now ?? null) !== v)
      input.value = v;
  }
}

const memory = new WeakMap<DG.DataFrame, Map<string, Kept>>();
/** The page's own values (the endpoint, the evaluation day: not saved with the study) kept for the tab too, by the
 *  study's row (else its name), so a reload, or the study opened again, gives them back (QA 0.1.3 round 4 chaos N7).
 *  ponytail: sessionStorage, the tab's; a study setting if they ever need to follow it to another machine. */
const storeKey = (df: DG.DataFrame, stage: string): string | null => {
  let row = '';
  try {
    row = JSON.parse(df.getTag('pmx.row') || 'null')?.id ?? '';
  } catch {}
  // (a study not saved yet by its name: the demo opened again; a table that is no study's, not at all)
  return row ? `pmx.form:${stage}:row=${row}` : df.getTag('pmx.study') ? `pmx.form:${stage}:name=${df.name}` : null;
};
const STAGES = ['monitor', 'analyze'];
/** The values kept for the tab: the page's own only (what is saved with the study comes back with it). */
const PAGE_ONLY = ['endpoint', 'evaluationDay'];
const pageOnly = (kept: Kept): Kept => Object.fromEntries(Object.entries(kept).filter(([k]) => PAGE_ONLY.includes(k)));
/** Keeps `kept` as `stage`'s form values for `df`. */
export function rememberForm(df: DG.DataFrame, stage: string, kept: Kept): void {
  if (!memory.has(df))
    memory.set(df, new Map());
  memory.get(df)!.set(stage, kept);
  const key = storeKey(df, stage);
  try {
    if (key)
      sessionStorage.setItem(key, JSON.stringify(pageOnly(kept)));
  } catch {} // (storage off: kept for the page only)
}
/** `stage`'s form values kept for `df` (null: none). */
export function recallForm(df: DG.DataFrame, stage: string): Kept | null {
  const own = memory.get(df)?.get(stage);
  if (own)
    return own;
  const key = storeKey(df, stage);
  try {
    const s = key ? sessionStorage.getItem(key) : null;
    return s ? JSON.parse(s) : null;
  } catch {
    return null;
  }
}
/** Forgets what the stages kept for `df` (its changes discarded: it opens again as saved). */
export function forgetForm(df: DG.DataFrame): void {
  memory.delete(df);
  for (const stage of STAGES) {
    const key = storeKey(df, stage);
    try {
      if (key)
        sessionStorage.removeItem(key);
    } catch {}
  }
}
