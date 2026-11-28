/* The user's studies, as the Overview lists them: each saved plan and followed study in their Home, by its details
   (ID, line, what the groups receive, owner, tags) and how far it has gone, kept in a small index beside them (each
   file's row, written as it is saved): the list reads the index, not the studies, and a study is read when opened.
   The index is a cache: missing, damaged or of another format, it is rebuilt from the files. A study followed from
   a plan is listed once, at its furthest stage. Datagrok projects are not listed: a study shared as a project keeps its
   tags and marks in its table, and opens from Browse. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {parsePlan, PLAN_EXT, PLAN_FORMAT, plansFolder} from './plan-file';
import {parseStudy, savedGroupsShown, studiesFolder, STUDY_EXT, STUDY_FORMAT} from './study-file';
import {armParts, compoundsOf, parseMeta, treatmentOf} from './study-meta';
import {futureDays} from './monitor';
import {controlOf} from './study-effect';

export interface LibraryEntry {
  kind: 'plan' | 'study';
  /** 1 planned, 2 followed (Monitor), 3 finished (Analyze) */
  stage: 1 | 2 | 3;
  file: DG.FileInfo | null;
  text: string;
  id: string;
  name: string;
  line: string;
  compounds: string;
  /** the drugs its groups receive, each once (their details, else their names) */
  drugs: string;
  owner: string;
  /** who runs it, when not the owner (Study details' Run by; empty: the owner) */
  runBy?: string;
  tags: string[];
  /** ms since the epoch */
  updated: number;
  /** what each group receives, with its dose, route and schedule (else the compounds, else the arms) */
  treatment: string;
  /** day 0, YYYY-MM-DD (null: not fixed yet) */
  start: string | null;
  /** a followed study's last day measured and its next (from its usual gap); null for a plan */
  day: number | null;
  next: number | null;
  /** the study's last day as planned (null: no plan) */
  end: number | null;
  /** with domain storage on (store.ts): the row it comes from, read when opened, and the model as lists name it */
  row?: {kind: 'plan' | 'study', id: string, version: number};
  model?: string;
  /** a copy saved from a row shared to read: the original's row and whose it is */
  copyOf?: {id: string, owner: string};
  /** a row shared with the user to read: no Delete for them */
  readOnly?: boolean;
}

/** A plan or study file's entry (null: not one, or damaged; `strict`: what is wrong, thrown, for a document the app
 *  made itself, which must read). */
export function entryOf(text: string, updated: number, file: DG.FileInfo | null = null,
  strict = false): LibraryEntry | null {
  try {
    const x = JSON.parse(text);
    const study = x?.format === 'pmx-study' ? parseStudy(text) : null;
    const plan = study ? study.plan : parsePlan(text);
    const meta = parseMeta(study ? x.meta : plan?.study) ?? parseMeta(plan?.study);
    const days = (study?.rows ?? []).map((r) => r.day).filter((d): d is number => Number.isFinite(d));
    // a blinded study names no arm or compound outside itself until its groups are revealed (row, search, index)
    const shown = !study || savedGroupsShown(study);
    // what the groups receive: the details, else the plan's arms, else (a study opened from a table) its own groups
    const groups = study ? [...new Set(study.rows.map((r) => r.group).filter((g): g is string => !!g))] : [];
    const arms = (plan?.predicted.arms ?? []).map((a) => a.name);
    const compounds = !shown ? '' : compoundsOf(meta) || arms.join(', ') || groups.join(', ');
    const names = [...new Set([...(meta?.groups ?? []).map((g) => g.group), ...arms, ...groups])];
    // each drug once, however its arms write it ("Trametinib", "trametinib"): the first writing kept
    const byKey = new Map<string, string>();
    for (const p of names.flatMap((n) => armParts(n, meta)))
      byKey.set(p.key, byKey.get(p.key) ?? p.drug);
    const drugs = !shown ? '' : [...byKey.values()].join(', ');
    // what each treated group receives, the groups apart by "; " as treatmentOf writes them (a ", " or " + " between
    // groups read as one list of drugs, or as a combination): else the plan's arms, else the study's treated groups
    const control = controlOf(groups, meta);
    const treatment = !shown ? '' : treatmentOf(meta) || arms.join('; ') ||
      groups.filter((g) => g !== control).join('; ');
    // the model: the plan's, else the one chosen for the study
    const line = plan?.history?.name || plan?.predicted.line || study?.history?.name || study?.line || '';
    return {kind: study ? 'study' : 'plan', stage: study ? (study.finished === true ? 3 : 2) : 1, file, text,
      // (a followed study still named as a new table's default goes by its plan's name)
      id: meta?.id ?? '', name: (study && (!x.name || x.name === 'New study') && plan?.name) || (x.name ?? ''),
      line, compounds, drugs, owner: meta?.owner || plan?.author || '', runBy: meta?.runBy ?? '',
      tags: meta?.tags ?? [],
      updated, treatment, start: meta?.start ?? null,
      day: days.length ? days.reduce((a, b) => Math.max(a, b)) : null,
      next: days.length ? futureDays(days, plan, meta?.start)[0] : null,
      end: plan?.predicted.studyEnd ?? null};
  } catch (e) {
    if (strict)
      throw e;
    return null;
  }
}

/** One entry per study, the furthest stage kept (a plan and the study followed from it share an ID), newest first. */
export function libraryOf(entries: (LibraryEntry | null)[]): LibraryEntry[] {
  const kept = new Map<string, LibraryEntry>();
  entries.forEach((e, i) => {
    if (!e)
      return;
    const key = e.id || `#${i}`;
    const was = kept.get(key);
    if (!was || e.stage > was.stage || (e.stage === was.stage && e.updated > was.updated))
      kept.set(key, e);
  });
  return [...kept.values()].sort((a, b) => b.updated - a.updated);
}

/** Whether an entry has every word of `query` in its ID, name, line, drugs, treatment, owner, run by or tags. */
export function matches(e: LibraryEntry, query: string): boolean {
  const text = [e.id, e.name, e.line, e.drugs, e.treatment, e.owner, e.runBy ?? '', ...e.tags].join(' ').toLowerCase();
  return query.toLowerCase().split(/\s+/).filter((w) => w).every((w) => text.includes(w));
}

/** What the library needs of the files API (the package's tests give a fake one). */
export interface LibraryFiles {
  list(folder: string, recursive: boolean, ext: string): Promise<DG.FileInfo[]>;
  readAsText(file: DG.FileInfo | string): Promise<string>;
  writeAsText(path: string, text: string): Promise<void>;
}

/** Each file's row, by file name: its entry without the file and its text, and the file's date as the server lists
 *  it (a row whose date is not the file's is read again); a file read and found to be no PMX Studio file at all, its
 *  date and `unreadable` (one of ours this version cannot read gets no row: the next version may read it). */
type Index = {[file: string]: Omit<LibraryEntry, 'file' | 'text'> | {updated: number, unreadable: true}};

/** The index file's format: another (or none) is rebuilt. Any change to entryOf or LibraryEntry, or to what parsePlan
 *  or parseStudy accept, raises it. */
// (6: the rows of 0.1.3's first indexing kept a study's groups joined by " + ", as one combination: QA 0.1.3 A16)
export const INDEX_VERSION = 6;
/** Whether a file's text is a PMX Studio plan or study, of any version. */
const ours = (text: string): boolean => {
  try {
    return [PLAN_FORMAT, STUDY_FORMAT].includes(JSON.parse(text)?.format);
  } catch {
    return false;
  }
};
export const indexPath = () => `${grok.shell.user.project.name}:Home/PMX Studio/library-index.json`;
const readIndex = (fs: LibraryFiles): Promise<Index> => fs.readAsText(indexPath()).then((t) => {
  const x = JSON.parse(t);
  return x?.version === INDEX_VERSION && x.rows && typeof x.rows === 'object' && !Array.isArray(x.rows) ? x.rows : {};
}).catch(() => ({}));
/** A row as entryOf writes it (one that is not is read again from its file). */
const rowOk = (r: any): boolean => !!r && typeof r === 'object' && typeof r.name === 'string' &&
  typeof r.id === 'string' && Array.isArray(r.tags) && Number.isFinite(r.updated) &&
  (r.kind === 'plan' || r.kind === 'study') && [1, 2, 3].includes(r.stage);
const writeIndex = (fs: LibraryFiles, rows: Index) =>
  fs.writeAsText(indexPath(), JSON.stringify({version: INDEX_VERSION, rows}));
const rowOf = (e: LibraryEntry, updated: number): Index[string] => {
  const {file: _f, text: _t, ...row} = e;
  return {...row, updated};
};
const dateOf = (f: DG.FileInfo) => f.updatedOn?.valueOf() ?? NaN;

/** Reads the library for `show`, one reading at a time. Asked again before the reading has listed the folders, it
 *  joins it (the listing sees whatever changed); asked once they are listed, it reads again once that reading ends.
 *  A reading that fails keeps the list shown and says so (`warn`). */
export function libraryReader(show: (entries: LibraryEntry[]) => void, fs: LibraryFiles = grok.dapi.files,
  warn: (message: string) => unknown = (message) => grok.shell.warning(message)) {
  let [reading, listing, again] = [false, false, false];
  const read = (): void => {
    if (listing)
      return;
    if (reading) {
      again = true;
      return;
    }
    reading = listing = true;
    void myLibrary(fs, () => listing = false).then(show, () =>
      warn('Your studies could not be listed just now: open the Overview again in a moment.')).finally(() => {
      reading = listing = false;
      if (again) {
        again = false;
        read();
      }
    });
  };
  return read;
}

/** Runs the calls given to it at most `n` at a time, the rest queued in order (a slot is handed straight to the next
 *  call, so none slips in between). */
export function limited(n: number): <T>(fn: () => Promise<T>) => Promise<T> {
  let active = 0;
  const queue: (() => void)[] = [];
  return async (fn) => {
    if (active < n)
      active++;
    else
      await new Promise<void>((r) => queue.push(r));
    try {
      return await fn();
    } finally {
      const next = queue.shift();
      if (next)
        next();
      else
        active--;
    }
  };
}

/** How many files the library reads at once (the first listing may read hundreds). */
const READS_AT_ONCE = 6;

/** The user's plans and studies: listed (names and dates only) and each one's row from the index; a file the index
 *  has no row for, or an older one (saved before the index, or changed outside the app), is read once and its row
 *  kept. No study is read otherwise, and nothing is computed. A folder that fails to list fails the reading (a
 *  missing one lists nothing), and no row is dropped. */
export async function myLibrary(fs: LibraryFiles = grok.dapi.files, onListed = () => {}): Promise<LibraryEntry[]> {
  const [plans, studies, index] = await Promise.all([fs.list(plansFolder(), false, PLAN_EXT),
    fs.list(studiesFolder(), false, STUDY_EXT), readIndex(fs)]);
  onListed();
  const files = [...plans, ...studies];
  let changed = false;
  const read = limited(READS_AT_ONCE);
  const entries = await Promise.all(files.map(async (f) => {
    // (a file with no date is read every time: a change to it could not be seen; one of another date than its row,
    // a copy restored too, is read again)
    const updated = dateOf(f);
    const dated = Number.isFinite(updated);
    const row = index[f.fileName];
    if (dated && row?.updated === updated && (row as {unreadable?: unknown}).unreadable === true)
      return null;
    if (dated && rowOk(row) && row.updated === updated)
      return {...row, file: f, text: ''} as LibraryEntry;
    // (a file that fails to read, maybe in passing, is read again next time, and so is one of ours this version cannot
    // read; one that is no PMX Studio file is not, until it changes)
    const text = await read(() => fs.readAsText(f)).catch(() => null);
    const e = text === null ? null : entryOf(text, dated ? updated : 0, f);
    if (text !== null && dated && (e || !ours(text))) {
      index[f.fileName] = e ? rowOf(e, updated) : {updated, unreadable: true};
      changed = true;
    }
    return e && {...e, text: ''};
  }));
  // rows of files no longer there leave the index
  const names = new Set(files.map((f) => f.fileName));
  for (const name of Object.keys(index)) {
    if (!names.has(name)) {
      delete index[name];
      changed = true;
    }
  }
  if (changed) // (an index that cannot be written leaves the list as read)
    await writeIndex(fs, index).catch(() => {});
  return libraryOf(entries);
}
