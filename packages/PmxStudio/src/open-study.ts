/* What a file brought to the hub holds, and the stage it opens in: a saved plan (.pmxplan), a followed study
   (.pmxstudy) or a table of measurements (.csv, .tsv), long or wide (wide-table.ts). Reads the dates only, never
   the groups, so a blinded study stays blind until it is analyzed; a wide table says how many animals it read in
   each group, as nothing may be lost unseen in turning it long. */
import * as DG from 'datagrok-api/dg';

import {PLAN_EXT, SavedPlan, parsePlan} from './plan-file';
import {STUDY_EXT, SavedStudy, parseStudy} from './study-file';
import {plainError} from './ui';
import {ROLES, WEIGHT, WideRead, WideRow, cellsOf, readWide, separatorOf} from './wide-table';

export const OPENED_EXTS = [PLAN_EXT, STUDY_EXT, 'csv', 'tsv', 'txt'];

/** A table of measurements, as far as the hub reads it. */
export interface MeasuredTable {kind: 'table', file: string, df: DG.DataFrame, animals: number | null,
  groups: number | null, lastDay: number | null, found: string[], missing: string[], weight: boolean,
  /** A wide table as read (wide-table.ts): what was read, or why nothing was. */
  wide: WideRead | null}
/** `saved`: the file as read, so opening it does not read it again. */
export interface PlanRead {kind: 'plan', file: string, text: string, saved: SavedPlan, name: string, line: string,
  end: number | null, blinded: boolean}
export interface StudyRead {kind: 'study', file: string, text: string, saved: SavedStudy}
/** `why`: what is wrong with it, in a sentence. */
export type FileRead = MeasuredTable | PlanRead | StudyRead | {kind: 'bad', file: string, why: string};

export {ROLES, WEIGHT};

const extOf = (name: string) => (name.split('.').pop() ?? '').toLowerCase();
const baseOf = (name: string) => name.replace(/\.[^.]+$/, '');

/** What one file holds. */
export function readFile(name: string, text: string): FileRead {
  const ext = extOf(name);
  try {
    if (ext === PLAN_EXT) {
      const p = parsePlan(text);
      return {kind: 'plan', file: name, text, saved: p, name: p.name || baseOf(name), line: p.predicted.line ?? '',
        end: Number.isFinite(p.predicted.studyEnd) ? Math.round(p.predicted.studyEnd) : null,
        blinded: !!p.study?.blinded};
    }
    if (ext === STUDY_EXT) {
      // read now, as a plan is: a damaged or newer study is said so here, not on opening it
      return {kind: 'study', file: name, text, saved: parseStudy(text)};
    }
    if (ext === 'csv' || ext === 'tsv' || ext === 'txt') {
      const read = readWide(cellsOf(text), separatorOf(text) === ';');
      const df = read ? longTable(read.error ? [] : read.rows) : DG.DataFrame.fromCsv(text);
      if (df.rowCount === 0 && !read)
        return {kind: 'bad', file: name, why: 'The table has no rows: open the file with the measurements.'};
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
        missing: roles.filter((r) => !r.column).map((r) => r.role), weight: !!at(WEIGHT), wide: read};
    }
  } catch (e: any) {
    // a plan or a study says what is wrong in a sentence; a table that does not parse, below
    if (ext === PLAN_EXT || ext === STUDY_EXT)
      return {kind: 'bad', file: name, why: plainError(e, 'The file is damaged: open another copy of it.')};
    console.warn(`PMX Studio: ${name}`, e);
    return {kind: 'bad', file: name, why: 'The table could not be read. Save it again as .csv or .tsv from your ' +
      'spreadsheet.'};
  }
  return {kind: 'bad', file: name, why: 'PMX Studio opens .pmxplan, .pmxstudy, .csv and .tsv files.'};
}

/** Long rows as the study table the stages take (emptyStudy's columns, without body weight; the source cells,
 *  only a check, are not kept). */
function longTable(rows: WideRow[]): DG.DataFrame {
  const cols = [DG.Column.fromStrings('animal', rows.map((r) => r.animal)),
    DG.Column.fromStrings('group', rows.map((r) => r.group)),
    DG.Column.fromFloat64Array('day', Float64Array.from(rows, (r) => r.day)),
    DG.Column.fromList(DG.TYPE.FLOAT, 'volume', rows.map((r) => r.volume))];
  return DG.DataFrame.fromColumns(cols);
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
export function suggestedStage(b: Brought): OpenStage {
  // a followed study is Monitor's file
  if (b.study)
    return 'monitor';
  if (b.table) {
    const t = b.table;
    const running = !!b.plan && b.plan.end !== null && t.lastDay !== null && t.lastDay < b.plan.end - 1;
    return running ? 'monitor' : 'analyze';
  }
  return 'plan';
}
