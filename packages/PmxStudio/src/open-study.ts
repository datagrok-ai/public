/* What a file brought to the hub holds, and the stage it opens in: a saved plan (.pmxplan), a followed study
   (.pmxstudy) or a table of measurements (.csv, .tsv). Reads the dates only, never the groups, so a running study
   stays blind until it is analyzed. */
import * as DG from 'datagrok-api/dg';

import {PLAN_EXT, parsePlan} from './plan-file';

export const STUDY_EXT = 'pmxstudy';
export const OPENED_EXTS = [PLAN_EXT, STUDY_EXT, 'csv', 'tsv', 'txt'];

/** A table of measurements, as far as the hub reads it. */
export interface MeasuredTable {kind: 'table', file: string, df: DG.DataFrame, animals: number | null,
  groups: number | null, lastDay: number | null, found: string[], missing: string[], weight: boolean}
export interface PlanRead {kind: 'plan', file: string, text: string, name: string, line: string, end: number | null}
export interface StudyRead {kind: 'study', file: string, text: string}
export type FileRead = MeasuredTable | PlanRead | StudyRead | {kind: 'bad', file: string};

/** The columns Analyze needs, by the names a study table usually gives them. */
const ROLES: [string, RegExp][] = [
  ['Animal', /^(animal|mouse|subject|id\b|tumou?r[ _]?id)/i],
  ['Group', /(group|arm|treatment|cohort)/i],
  ['Day', /^(day|days|study[ _]?day|time)\b/i],
  ['Volume', /vol/i],
];
const WEIGHT = /(body[ _]?weight|^bw\b|^weight)/i;

const extOf = (name: string) => (name.split('.').pop() ?? '').toLowerCase();
const baseOf = (name: string) => name.replace(/\.[^.]+$/, '');

/** What one file holds. */
export function readFile(name: string, text: string): FileRead {
  const ext = extOf(name);
  try {
    if (ext === PLAN_EXT) {
      const p = parsePlan(text);
      return {kind: 'plan', file: name, text, name: p.name || baseOf(name), line: p.predicted.line ?? '',
        end: Number.isFinite(p.predicted.studyEnd) ? Math.round(p.predicted.studyEnd) : null};
    }
    if (ext === STUDY_EXT)
      return {kind: 'study', file: name, text};
    if (ext === 'csv' || ext === 'tsv' || ext === 'txt') {
      const df = DG.DataFrame.fromCsv(text);
      if (df.rowCount === 0)
        return {kind: 'bad', file: name};
      df.name = baseOf(name);
      const names = df.columns.names();
      const at = (re: RegExp) => names.find((c) => re.test(c));
      const roles = ROLES.map(([role, re]) => ({role, column: at(re)}));
      const count = (role: string) => {
        const c = roles.find((r) => r.role === role)!.column;
        return c ? df.col(c)!.stats.uniqueCount : null;
      };
      const day = roles.find((r) => r.role === 'Day')!.column;
      const dayCol = day ? df.col(day) : null;
      return {kind: 'table', file: name, df, animals: count('Animal'), groups: count('Group'),
        lastDay: dayCol?.isNumerical ? Math.round(dayCol.stats.max) : null,
        found: roles.filter((r) => r.column).map((r) => `${r.role} ← ${r.column}`),
        missing: roles.filter((r) => !r.column).map((r) => r.role), weight: !!at(WEIGHT)};
    }
  } catch {
    // a file that does not parse as its extension says is not a study file
  }
  return {kind: 'bad', file: name};
}

/** What the files brought so far hold together: one of each kind, the latest winning. */
export interface Brought {table?: MeasuredTable, plan?: PlanRead, study?: StudyRead}

export function bring(into: Brought, reads: FileRead[]): {brought: Brought, bad: string[]} {
  const brought = {...into};
  const bad: string[] = [];
  for (const r of reads) {
    if (r.kind === 'bad')
      bad.push(r.file);
    else if (r.kind === 'table')
      brought.table = r;
    else if (r.kind === 'plan')
      brought.plan = r;
    else
      brought.study = r;
  }
  return {brought, bad};
}

export type OpenStage = 'plan' | 'monitor' | 'analyze';

/** The stage a study is suggested to open in, from its dates alone: a table whose last measurement comes before
 *  its plan's end is still running. */
export function suggestedStage(b: Brought, monitorReady: boolean): OpenStage {
  // a followed study is Monitor's file, whether Monitor is here yet or not
  if (b.study)
    return 'monitor';
  if (b.table) {
    const t = b.table;
    const running = !!b.plan && b.plan.end !== null && t.lastDay !== null && t.lastDay < b.plan.end - 1;
    return running && monitorReady ? 'monitor' : 'analyze';
  }
  return 'plan';
}
