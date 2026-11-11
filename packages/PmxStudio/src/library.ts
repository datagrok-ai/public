/* The user's studies, as the Overview lists them: each saved plan and followed study in their Home, by its details
   (ID, line, what the groups receive, owner, tags) and how far it has gone, kept in a small index beside them (each
   file's row, written as it is saved): the list reads the index, not the studies, and a study is read when opened.
   The index is a cache: missing, damaged or of another format, it is rebuilt from the files. A study followed from
   a plan is listed once, at its furthest stage. Datagrok projects are not listed: a study shared as a project keeps its
   tags and marks in its table, and opens from Browse. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {parsePlan, PLAN_EXT, PLAN_FORMAT, plansFolder} from './plan-file';
import {parseStudy, studiesFolder, STUDY_EXT, STUDY_FORMAT} from './study-file';
import {compoundsOf, parseMeta, treatmentOf} from './study-meta';
import {futureDays} from './monitor';

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
  owner: string;
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
}

/** A plan or study file's entry (null: not one, or damaged). */
export function entryOf(text: string, updated: number, file: DG.FileInfo | null = null): LibraryEntry | null {
  try {
    const x = JSON.parse(text);
    const study = x?.format === 'pmx-study' ? parseStudy(text) : null;
    const plan = study ? study.plan : parsePlan(text);
    const meta = parseMeta(study ? x.meta : plan?.study) ?? parseMeta(plan?.study);
    const days = (study?.rows ?? []).map((r) => r.day).filter((d): d is number => Number.isFinite(d));
    const compounds = compoundsOf(meta) || (plan?.predicted.arms ?? []).map((a) => a.name).join(', ');
    return {kind: study ? 'study' : 'plan', stage: study ? 2 : 1, file, text, id: meta?.id ?? '', name: x.name ?? '',
      line: plan?.predicted.line ?? '', compounds, owner: meta?.owner || plan?.author || '', tags: meta?.tags ?? [],
      updated, treatment: treatmentOf(meta) || compounds, start: meta?.start ?? null,
      day: days.length ? days.reduce((a, b) => Math.max(a, b)) : null,
      next: days.length ? futureDays(days, plan, meta?.start)[0] : null,
      end: plan?.predicted.studyEnd ?? null};
  } catch {
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

/** Whether an entry has every word of `query` in its ID, name, line, treatment, owner or tags. */
export function matches(e: LibraryEntry, query: string): boolean {
  const text = [e.id, e.name, e.line, e.treatment, e.owner, ...e.tags].join(' ').toLowerCase();
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
const INDEX_VERSION = 1;
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
    const text = await fs.readAsText(f).catch(() => null);
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
