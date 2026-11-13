/* What a file brought to the hub holds, and the stage it opens in: a saved plan (.pmxplan), a followed study
   (.pmxstudy) or a table of measurements (.csv, .tsv), long or wide, read by study-reader.ts. Reads the dates only,
   never the groups, so a blinded study stays blind until it is analyzed; a table says how many animals it read in
   each group, and how, as nothing may be lost unseen in reading it. */
import * as DG from 'datagrok-api/dg';

import {PLAN_EXT, SavedPlan, parsePlan} from './plan-file';
import {STUDY_EXT, SavedStudy, VOLUME_TAG, WEIGHT_TAG, oneDecimal, parseStudy} from './study-file';
import {Answers, TableRead, readStudy} from './study-reader';
import {plainError} from './ui';
import {cellsOf, separatorOf} from './wide-table';

export const OPENED_EXTS = [PLAN_EXT, STUDY_EXT, 'csv', 'tsv', 'txt'];

/** A table of measurements, as far as the hub reads it. */
export interface MeasuredTable {kind: 'table', file: string, df: DG.DataFrame, animals: number | null,
  groups: number | null, lastDay: number | null, weight: boolean,
  /** The table as read (study-reader.ts), with the answers to its questions so far; `text`, the file, read again
   *  with each answer. */
  read: TableRead, text: string, answers: Answers}
/** `saved`: the file as read, so opening it does not read it again. */
export interface PlanRead {kind: 'plan', file: string, text: string, saved: SavedPlan, name: string, line: string,
  end: number | null, blinded: boolean}
export interface StudyRead {kind: 'study', file: string, text: string, saved: SavedStudy}
/** `why`: what is wrong with it, in a sentence. */
export type FileRead = MeasuredTable | PlanRead | StudyRead | {kind: 'bad', file: string, why: string};

const extOf = (name: string) => (name.split('.').pop() ?? '').toLowerCase();
const baseOf = (name: string) => name.replace(/\.[^.]+$/, '');

/** What one file holds; a table read with the `answers` given to its questions. */
export function readFile(name: string, text: string, answers: Answers = {}): FileRead {
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
      // an .xlsx plugs in here: a sheet's cells as text, read by readStudy as a CSV's are
      const cells = cellsOf(text);
      if (!cells.some((row) => row.some((c) => c !== '')))
        return {kind: 'bad', file: name, why: 'The table has no rows: open the file with the measurements.'};
      const read = readStudy(cells, answers, separatorOf(text) === ';');
      const df = studyTable(read);
      df.name = baseOf(name);
      const days = read.rows.map((r) => r.day);
      return {kind: 'table', file: name, df, animals: days.length ? new Set(read.rows.map((r) => r.animal)).size : null,
        groups: days.length ? read.groups.length : null, lastDay: days.length ? Math.round(Math.max(...days)) : null,
        weight: read.rows.some((r) => r.weight !== null), read, text, answers};
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

/** A read as the study table the stages take (emptyStudy's columns, body weight when read, the table's other
 *  columns as they were), with where its volumes come from. */
function studyTable(read: TableRead): DG.DataFrame {
  const rows = read.rows;
  const floats = (name: string, xs: (number | string | null)[]) => DG.Column.fromFloat64Array(name,
    Float64Array.from(xs, (x) => typeof x === 'number' ? x : DG.FLOAT_NULL));
  const cols = [DG.Column.fromStrings('animal', rows.map((r) => r.animal)),
    DG.Column.fromStrings('group', rows.map((r) => r.group)), floats('day', rows.map((r) => r.day)),
    floats('volume', rows.map((r) => r.volume))];
  const weight = rows.some((r) => r.weight !== null);
  if (weight)
    cols.push(floats('body_weight', rows.map((r) => r.weight)));
  const used = new Set(cols.map((c) => c.name));
  for (const [title, xs] of read.extra) {
    // two columns with one heading: the later one numbered
    let name = title;
    for (let k = 2; used.has(name); k++)
      name = `${title} (${k})`;
    used.add(name);
    cols.push(xs.every((x) => typeof x !== 'string') ? floats(name, xs) :
      DG.Column.fromStrings(name, xs.map((x) => x === null ? '' : String(x))));
  }
  const df = oneDecimal(DG.DataFrame.fromColumns(cols));
  if (weight)
    df.setTag(WEIGHT_TAG, 'body_weight');
  if (read.volume)
    df.setTag(VOLUME_TAG, JSON.stringify(read.volume));
  return df;
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
