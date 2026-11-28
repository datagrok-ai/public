/* Plans and studies as rows of the package's domain tables (databases/pmx/schema.json, EMS): the document whole, as
   the .pmxplan / .pmxstudy file holds it, in a json column, and beside it what the Overview lists and searches by,
   computed by entryOf at each save. The row's id is the record's identity; its version guards an update against a
   save made elsewhere since it was opened. Study IDs come from a lab-wide counter (pmx.code, autoNumber).
   Access (D691): everyone may read the table and create (View and Edit on the table; Edit is EMS's insert grant); a
   row is its author's alone until shared, to read, from the Overview (share: as NCA, view only until the platform
   honours row Edit grants). Blinding: the document holds the allocation, so plan and study keep no audit trail
   (schema "audit": false): the platform's History has no earlier versions to show.
   The `domainStorage` setting is on by default (storage-setting.ts); switched off, the app keeps its files and none
   of this is used. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {LibraryEntry, entryOf, libraryOf, myLibrary} from './library';
import {PLAN_EXT, SavedPlan, myPlans, parsePlan, planKey, plansFolder, sha256Hex} from './plan-file';
import {STUDY_EXT, SavedStudy, parseStudy, studiesFolder, tableOf} from './study-file';
import {DEMO_OWNER, StudyMeta, emptyMeta, withId} from './study-meta';
import {domainStorageOn} from './storage-setting';
import {pmxDb} from './generated/db';

export type Kind = 'plan' | 'study';

/** A saved row: its id and the version it was read or written at (an update passes it back). */
export interface Saved {id: string, version: number}

/** A row of pmx.plan or pmx.study (the columns the store writes; the server adds the system ones). `owner` and
 *  `copy_owner` are user ids (their names come as `~caption_owner` / `~caption_copy_owner` when a query asks);
 *  `document` is the document itself (json), `start` day 0 at midnight UTC (read back as a Dayjs). */
export interface DocRow {
  id: string; version: number; updated_on?: {valueOf(): number}; author_id?: string;
  name: string; code?: string; plan_id?: string | null; stage?: number; line?: string; model?: string;
  owner?: string | null; run_by?: string; tags?: string[]; compounds?: string; drugs?: string[]; treatment?: string;
  start?: string | {toISOString(): string} | null; last_day?: number | null; next_day?: number | null;
  end_day?: number | null; document: object | string; source?: string; copy_of?: string; copy_owner?: string | null;
  copy_file_owner?: string; '~caption_owner'?: string | null; '~caption_copy_owner'?: string | null;
}
export interface CodeRow {id: string, version: number, number?: number}

/** A row's document as text, for the file readers: json comes back as an object, keys reordered (a row written
 *  before schema 1.5.0 and not migrated held the text itself). */
export const docText = (d: object | string): string => typeof d === 'string' ? d : JSON.stringify(d);
/** Day 0 as the document writes it (YYYY-MM-DD), from the datetime column (or a 1.4.1 row's text). */
export const dayOf = (x: DocRow['start']): string | null => x == null ? null :
  (typeof x === 'string' ? x : x.toISOString()).slice(0, 10) || null;
/** A whole study day for an int column (null: none). */
const whole = (x: number | null | undefined): number | null => x == null || !Number.isFinite(x) ? null : Math.round(x);
/** The user a new row belongs to. */
const meId = (): string | null => grok.shell.user?.id ?? null;

/** What the store needs of a domain table: the part of grok.dapi.domains.table it uses (tests give a memory one). */
export interface StoreTable<TRow> {
  insert(row: Partial<TRow>): Promise<{id: string, version?: number}[]>;
  promote?(id: string): Promise<unknown>;
  update(id: string, values: Partial<TRow>, options: {version: number}): Promise<{id: string, version: number}>;
  get(id: string, options?: {withAccess?: boolean}): Promise<(TRow & {'~can_edit'?: boolean}) | null>;
  query(spec: {columns?: string[], sort?: string, search?: string, limit?: number, filter?: Where,
    deleted?: 'only', withAccess?: boolean, captions?: string[]}): Promise<(TRow & {'~can_edit'?: boolean})[]>;
  delete(id: string): Promise<unknown>;
  restore(id: string): Promise<unknown>;
}
/** Equalities AND-ed, in the condition-tree form grok.dapi.domains takes ([c1, 'and', c2]). */
export type Where = ({property: string, operator: '=', value: string} | 'and')[];
export interface Tables {plan: StoreTable<DocRow>, study: StoreTable<DocRow>, code: StoreTable<CodeRow>}

/** The package's tables on the server, through the client grok api generates from the schema (src/generated/db.ts;
 *  its rows are typed as the store reads them, so the cast is only of the narrower StoreTable view). */
export const domainTables = (): Tables =>
  ({plan: pmxDb.plans, study: pmxDb.studies, code: pmxDb.codes}) as unknown as Tables;

/** Whether plans and studies are saved as rows (package setting `domainStorage`); off, files as before. */
export {domainStorageOn};

/** The listed columns of a document: its entry, as entryOf computes it for the Overview. */
export function columnsOf(kind: Kind, doc: SavedPlan | SavedStudy): Omit<DocRow, 'id' | 'version'> {
  const document = JSON.stringify(doc);
  // (the app's own document: what keeps it from reading is said, not hidden as another kind of file)
  const e = entryOf(document, 0, null, true);
  if (!e || e.kind !== kind)
    throw new Error(`Not a PMX Studio ${kind}`);
  const meta = kind === 'study' ? (doc as SavedStudy).meta : (doc as SavedPlan).study;
  const plan = kind === 'study' ? (doc as SavedStudy).plan : doc as SavedPlan;
  // (the owner is the user, set when the row is made: columnsOf never changes whose a row is)
  return {name: e.name, code: e.id, line: e.line, model: plan?.history?.name || e.line, run_by: meta?.runBy ?? '',
    tags: e.tags, compounds: e.compounds, drugs: e.drugs ? e.drugs.split(', ') : [], treatment: e.treatment,
    start: e.start ? `${e.start}T00:00:00Z` : null, end_day: whole(e.end), ...kind === 'study' ?
      {stage: e.stage, last_day: whole(e.day), next_day: whole(e.next)} : {}, document: JSON.parse(document)};
}

/** Saves a plan or study: a new row, or (`at`) the row it was opened from, refused if saved elsewhere since. */
export async function save(kind: Kind, doc: SavedPlan | SavedStudy, at: Saved | null = null,
  extra: {planId?: string | null, fileOwner?: string} = {}, t: Tables = domainTables(), me?: string):
  Promise<Saved & {code?: string, copyOf?: string, meta?: StudyMeta}> {
  // only the row the view was opened from (or saved to) is updated; anything else is a new row (firstSave: a new ID
  // when its ID is taken, and the user's copy when it came from someone else's file, `fileOwner`, as read when the
  // file was opened; what the owner field says now never makes a copy): no row is ever replaced by another's data.
  // `meta`: the details as saved, when they changed (the view shows them)
  // a row shared to read is refused here, as the server would, but without its "Access Denied" (the view offers a copy)
  // one deleted, here or elsewhere since (another tab, a colleague): said so, and the view offers it as a new one (QA
  // 0.1.3 round 3 chaos 11: it said "shared with you to read")
  const row = at && !trashed.has(`${kind}:${at.id}`) ? await t[kind].get(at.id, {withAccess: true}) : null;
  if (at && !row)
    throw Object.assign(new Error(`This ${kind} was deleted (it is in the trash).`), {code: 'deleted', kind});
  if (row?.['~can_edit'] === false)
    throw Object.assign(new Error(`This ${kind} is shared with you to read`), {code: 'forbidden'});
  const first = at ? {doc} : await firstSave(kind, doc, t, me, extra.fileOwner);
  const now = first.doc;
  const values = {...columnsOf(kind, now), ...kind === 'study' && extra.planId !== undefined ?
    {plan_id: extra.planId} : {}, ...first.copyOf ? {copy_file_owner: first.copyOf} : {}};
  const saved = at ? await t[kind].update(at.id, values, {version: at.version}) :
    await t[kind].insert({...values, owner: meId()}).then(([r]) => ({id: r.id, version: r.version ?? 1}));
  if (kind === 'plan')
    planRows.set(planKey(now as SavedPlan)!, {id: saved.id, version: saved.version});
  return {...saved, ...first.code ? {code: first.code, meta: metaIn(kind, now)} : {},
    ...first.copyOf ? {copyOf: first.copyOf} : {}};
}

/** What the user is told when a first save gave the record a new ID (`was`: the ID it had). */
export const newRowNote = (kind: Kind, saved: {code?: string, copyOf?: string}, was: string): string =>
  saved.copyOf ? `Saved as your copy of ${saved.copyOf}'s ${kind}, ${saved.code}.` :
    `Saved as a new ${kind}, ${saved.code}: one with the ID ${was} already exists.`;

/** A study's details, or a plan's. */
const metaIn = (kind: Kind, doc: SavedPlan | SavedStudy): StudyMeta =>
  (kind === 'study' ? (doc as SavedStudy).meta : (doc as SavedPlan).study) ?? emptyMeta();

/** The names an owner field may hold for the user (`me` given: that one). */
const myNames = (me?: string): string[] => (me !== undefined ? [me] : [grok.shell.user?.login,
  grok.shell.user?.friendlyName]).filter((n): n is string => !!n);

/** Who a plan or study names as its owner. */
export const ownerOf = (kind: Kind, doc: SavedPlan | SavedStudy): string => kind === 'study' ?
  (doc as SavedStudy).meta?.owner ?? '' : (doc as SavedPlan).study?.owner || (doc as SavedPlan).author || '';
/** The owner of the file a plan (the view's start) was opened from, as read then; not saved anywhere. */
const fileOwners = new WeakMap<object, string>();
export const markFromFile = (start: object, owner: string) => fileOwners.set(start, owner);
export const fileOwnerOf = (start: object | null | undefined): string | undefined =>
  start ? fileOwners.get(start) : undefined;
/** A study table opened from a file: its owner as read then (a tag no save reads). */
export const FROM_FILE_TAG = 'pmx.fromFile';
export const markTableFromFile = (df: DG.DataFrame, owner: string) => df.setTag(FROM_FILE_TAG, owner);
export const tableFileOwner = (df: DG.DataFrame): string | undefined => df.getTag(FROM_FILE_TAG) || undefined;

/** A plan or study about to become a new row, from wherever it came (a file, a drop, the Open dialog, a view): one
 *  of someone else's (its owner is not the user) becomes the user's copy, with a new ID and whose it was in copyOf;
 *  the user's own keeps its ID unless a row of its kind has it already, when it gets a new one. */
async function firstSave<T extends SavedPlan | SavedStudy>(kind: Kind, doc: T, t: Tables, me?: string,
  fileOwner?: string): Promise<{doc: T, code?: string, copyOf?: string}> {
  // a copy only for a file someone else owned; a record made in the app is the user's, whatever its owner field says
  const owner = fileOwner ?? '';
  const mine = myNames(me);
  // (the demo's files are no one's: saved, the user's own with an ID of its own, never "copy of demo", nor the study
  // followed from it, QA 0.1.3 browser check)
  const demo = owner.toLowerCase() === DEMO_OWNER;
  const others = !demo && !!owner && mine.length > 0 && !mine.some((n) => n.toLowerCase() === owner.toLowerCase());
  const code = others || demo ? await nextCode(t) : await takenCode(kind, doc, t);
  if (!code)
    return {doc};
  return {doc: withMeta(kind, doc, {id: code, ...(others || demo) && mine[0] ? {owner: mine[0]} : {},
    ...others ? {copyOf: owner} : {}}), code, ...others ? {copyOf: owner} : {}};
}

/** `doc` with its study details changed by `patch` (a study's meta, a plan's study). */
function withMeta<T extends SavedPlan | SavedStudy>(kind: Kind, doc: T, patch: Partial<StudyMeta>): T {
  if (kind === 'study') {
    const s = doc as SavedStudy;
    return {...s, meta: {...s.meta ?? emptyMeta(), ...patch}} as T;
  }
  const p = doc as SavedPlan;
  return {...p, study: {...p.study ?? emptyMeta(), ...patch}} as T;
}

/** A new study ID when `doc`'s is already a visible row's of its kind (null: free, or none). */
async function takenCode(kind: Kind, doc: SavedPlan | SavedStudy, t: Tables): Promise<string | null> {
  const code = columnsOf(kind, doc).code;
  if (!code)
    return null;
  const [same] = await t[kind].query({columns: ['name'], limit: 1,
    filter: [{property: 'code', operator: '=', value: code}]});
  return same ? nextCode(t) : null;
}

/** A plan or study row, read through the file's own checks (a damaged or newer document is refused as a file is). */
export async function open(kind: 'plan', id: string, t?: Tables): Promise<{doc: SavedPlan, saved: Saved}>;
export async function open(kind: 'study', id: string, t?: Tables): Promise<{doc: SavedStudy, saved: Saved}>;
export async function open(kind: Kind, id: string, t: Tables = domainTables()):
  Promise<{doc: SavedPlan | SavedStudy, saved: Saved}> {
  const row = await t[kind].get(id);
  if (!row) {
    // (one in the trash said so: a link to it opened a blank stage with no word, QA 0.1.3 round 3 chaos 14)
    const trash = await t[kind].query({columns: ['name'], deleted: 'only', limit: 1,
      filter: [{property: 'id', operator: '=', value: id}]}).catch(() => []);
    // (with its code: a save that reads the row first offers "Save as new" as every save of a deleted one does)
    throw Object.assign(new Error(trash.length ? `This ${kind} was deleted (it is in the trash). Restore it from the ` +
      'Overview (Show deleted) to open it.' : `This ${kind} is not there any more, or it is not shared with you.`),
    {code: trash.length ? 'deleted' : 'not-found', kind});
  }
  const saved = {id: row.id, version: row.version};
  const read = kind === 'study' ? parseStudy(docText(row.document)) : parsePlan(docText(row.document));
  // (a 1.4.1 row, its document still text, is read as the 1.5.0 migration reads it)
  const doc = typeof row.document === 'string' ? await runBy1_4(kind, read) : read;
  if (kind === 'plan')
    planRows.set(planKey(doc as SavedPlan)!, saved);
  return {doc, saved};
}

/** A 1.4.1 document's Run by: its owner text when that is not the login of any user (a login is a user, never free
 *  text; a copy kept the original owner's login there), as the 1.5.0 migration moves it. */
async function runBy1_4<T extends SavedPlan | SavedStudy>(kind: Kind, doc: T): Promise<T> {
  const meta = kind === 'study' ? (doc as SavedStudy).meta : (doc as SavedPlan).study;
  const owner = meta?.owner ?? '';
  return !meta || meta.runBy || !owner || await isLogin(owner) ? doc : withMeta(kind, doc, {runBy: owner});
}
/** Whether a text is a user's login here (a quote is in no login). */
const isLogin = async (s: string): Promise<boolean> => !/["\\]/.test(s) &&
  !!(await grok.dapi.users.filter(`login = "${s}"`).first().catch(() => null));

/** The rows a deep link names: by its row id, that row; by its study ID (PMX-…, what users see), the rows of its kind
 *  the user can see with that ID, newest first (several: the link opens a chooser). */
export async function rowsOf(kind: Kind, idOrCode: string, t: Tables = domainTables()): Promise<LinkedRow[]> {
  if (/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i.test(idOrCode))
    return [{id: idOrCode, name: '', owner: '', stage: 1, updated: 0}];
  const rows = await t[kind].query({columns: kind === 'study' ? ['name', 'stage'] : ['name'], captions: ['owner'],
    sort: '!updated_on', limit: 20, filter: [{property: 'code', operator: '=', value: idOrCode}]});
  if (!rows.length)
    throw new Error(`No ${kind} ${idOrCode} is saved in PMX Studio, or it is not shared with you.`);
  return rows.map((r) => ({id: r.id, name: r.name ?? '', owner: r['~caption_owner'] ?? '', stage: kind === 'study' ?
    r.stage ?? 2 : 1, updated: r.updated_on?.valueOf() ?? 0}));
}
/** A row a deep link by study ID may mean: when several the user can see have it, they pick one. */
export type LinkedRow = {id: string, name: string, owner: string, stage: number, updated: number};

/** The record a deep link names (`?study=<id>` or `?plan=<id>`), with domain storage on; null otherwise. */
export function linkedRow(search = window.location.search): {kind: Kind, id: string} | null {
  if (!domainStorageOn())
    return null;
  const q = new URLSearchParams(search);
  const [kind, id] = q.get('study') ? ['study' as const, q.get('study')!] : ['plan' as const, q.get('plan')];
  return id ? {kind, id} : null;
}

/** The rows moved to the trash in this session: a view still open on one does not save into it (QA pass 4: the
 *  view of a deleted plan stayed open and editable, and its Save would write a row the Overview no longer lists). */
const trashed = new Set<string>();
/** Whether the table's study row was moved to the trash in this session (the Study field lists it no more). */
export const inTrash = (df: DG.DataFrame): boolean => {
  const at = rowOfTable(df);
  return !!at && trashed.has(`study:${at.id}`);
};
/** Moves a plan or study to the trash (soft delete: restorable). */
export const remove = (kind: Kind, id: string, t: Tables = domainTables()) =>
  t[kind].delete(id).then((x) => (trashed.add(`${kind}:${id}`), x));
/** Brings a plan or study back from the trash (its author, or whoever may delete it). */
export const restore = (kind: Kind, id: string, t: Tables = domainTables()) =>
  t[kind].restore(id).then((x) => (trashed.delete(`${kind}:${id}`), x));

/** An Overview entry with the row it comes from. */
export type StoreEntry = LibraryEntry & {row: Saved & {kind: Kind}, model: string};

/** The declared columns a list reads (the server refuses a system column here: id, version, updated_on and author_id
 *  always come with the rows). */
const LISTED = ['name', 'code', 'stage', 'line', 'model', 'run_by', 'tags', 'compounds',
  'drugs', 'treatment', 'start', 'last_day', 'next_day', 'end_day', 'copy_of', 'copy_owner', 'copy_file_owner'];
/** The users a list names, by the platform's display name of each (a ref caption, one join each). */
const NAMED = ['owner', 'copy_owner'];
const entryOfRow = (kind: Kind, r: DocRow & {'~can_edit'?: boolean}): StoreEntry => {
  const copyOwner = r['~caption_copy_owner'] || r.copy_file_owner || '';
  return {kind, stage: kind === 'study' ? (r.stage === 3 ? 3 : 2) : 1,
    file: null, text: '', id: r.code ?? '', name: r.name, line: r.line || r.model || '', compounds: r.compounds ?? '',
    drugs: (r.drugs ?? []).join(', '), owner: r['~caption_owner'] ?? '', runBy: r.run_by ?? '',
    tags: r.tags ?? [], updated: r.updated_on?.valueOf() ?? 0, treatment: r.treatment ?? '', start: dayOf(r.start),
    day: r.last_day ?? null, next: r.next_day ?? null, end: r.end_day ?? null, model: r.model || r.line || '',
    row: {kind, id: r.id, version: r.version},
    ...r.copy_of || copyOwner ? {copyOf: {id: r.copy_of ?? '', owner: copyOwner}} : {},
    ...r['~can_edit'] === false ? {readOnly: true} : {}};
};

/** The plans and studies the user can see, newest first, one per study ID (its furthest stage), as the Overview lists
 *  them; `search` narrows on the server (name, ID, model, model key, run by, treatment; not the owner, a user column).
 *  No document is read. */
export async function list(search = '', limit = 1000, t: Tables = domainTables(), deleted = false):
  Promise<StoreEntry[]> {
  const spec = {columns: LISTED, captions: NAMED, sort: '!updated_on', limit, withAccess: true,
    ...search ? {search} : {},
    ...deleted ? {deleted: 'only' as const} : {}};
  const [plans, studies] = await Promise.all([t.plan.query({...spec, columns: LISTED.filter((c) =>
    !['stage', 'last_day', 'next_day'].includes(c))}), t.study.query(spec)]);
  return libraryOf([...plans.map((r) => entryOfRow('plan', r)), ...studies.map((r) => entryOfRow('study', r))]) as
    StoreEntry[];
}

/** With domain storage off: how many plans and studies the user can see as rows (the Overview lists files then). */
export async function keptRows(t: Tables = domainTables()): Promise<number> {
  // (as the list counts them, a plan and its study one; up to 1001: past 1000, the line says "more than 1000")
  return (await list('', 1001, t)).length;
}
/** The Overview's line for them, with storage off (null: none). The name is the setting's (package.json). */
export const keptRowsLine = (n: number): string | null => n ? `${n > 1000 ? 'More than 1000' : n} ${n === 1 ?
  'study is' : 'studies are'} kept in Datagrok storage: switch "Save in Datagrok" back on to see or export ` +
  `${n === 1 ? 'it' : 'them'}.` : null;

/** Brings a .pmxplan / .pmxstudy file's text in as a new row, once: a file already brought in (the same text, by its
 *  sha256 in `source`, among the rows the user can see) answers its row and writes nothing. No row is ever replaced.
 *  A file of someone else's (its owner is not `me`) comes in as Save as my copy does: the user's, with a new ID, and
 *  whose it was beside it (`copyOf`). The user's own keeps its ID unless a row of its kind has it already, when it
 *  gets a new one (`renumbered`). A damaged file is refused as when it is opened. */
export async function importFile(kind: Kind, text: string, t: Tables = domainTables(),
  me?: string): Promise<{saved: Saved, imported: boolean, renumbered?: boolean, copyOf?: string}> {
  const source = await sha256Hex(text);
  const [was] = await t[kind].query({columns: ['name'], limit: 1,
    filter: [{property: 'source', operator: '=', value: source}]});
  if (was)
    return {saved: {id: was.id, version: was.version}, imported: false};
  const doc = kind === 'plan' ? parsePlan(text) : parseStudy(text);
  const first = await firstSave(kind, doc, t, me, ownerOf(kind, doc));
  // ponytail: two tabs importing the same file at the same instant can both pass the check above and make two rows;
  // schema "idempotency" would close it, but the server only applies it when a table is created (EMS 1.28.0)
  const [r] = await t[kind].insert({...columnsOf(kind, first.doc), source, owner: meId(),
    ...first.copyOf ? {copy_file_owner: first.copyOf} : {}});
  return {saved: {id: r.id, version: r.version ?? 1}, imported: true, ...first.copyOf ? {copyOf: first.copyOf} :
    first.code ? {renumbered: true} : {}};
}

/** Shares a plan or study of the user's with groups, to read: the row is promoted to an entity, then granted View
 *  (edit sharing waits for the platform to honour row Edit grants, as in NCA). */
export async function share(kind: Kind, id: string, groups: DG.Group[], t: Tables = domainTables()): Promise<void> {
  await t[kind].promote?.(id);
  const [entity] = await grok.dapi.getEntities([id]);
  for (const g of groups)
    await grok.dapi.permissions.grant(entity, g, false);
}

/** What bringing the files in did: rows written, files already in, files saved under a new ID (one with theirs was
 *  there already), damaged files. */
export interface ImportReport {imported: number, already: number, renumbered: string[], damaged: string[]}

/** Brings every plan and study file of the user's PMX Studio folders in as rows (importFile each: run again, it
 *  writes nothing new). The files are left as they are. */
export async function importFolders(t: Tables = domainTables(), fs: {list(folder: string, recursive: boolean,
  ext: string): Promise<DG.FileInfo[]>, readAsText(f: DG.FileInfo): Promise<string>} = grok.dapi.files):
  Promise<ImportReport> {
  const report: ImportReport = {imported: 0, already: 0, renumbered: [], damaged: []};
  const [plans, studies] = await Promise.all([fs.list(plansFolder(), false, PLAN_EXT),
    fs.list(studiesFolder(), false, STUDY_EXT)]);
  // plans first: a study followed from a plan finds the plan's row by its study ID once both are in
  for (const [kind, files] of [['plan', plans], ['study', studies]] as const) {
    for (const f of files) {
      try {
        const r = await importFile(kind, await fs.readAsText(f), t);
        if (!r.imported)
          report.already++;
        else {
          report.imported++;
          if (r.renumbered)
            report.renumbered.push(f.fileName);
        }
      } catch {
        report.damaged.push(f.fileName);
      }
    }
  }
  return report;
}

/** The Overview's first opening with rows on: the user's files brought in (importFolders: only read, idempotent),
 *  then the rows listed, so the studies they had as files are there without a click. */
export async function firstOpen(t: Tables = domainTables(), fs?: Parameters<typeof importFolders>[1]):
  Promise<{entries: StoreEntry[], report: ImportReport}> {
  const report = await importFolders(t, fs);
  return {entries: await list('', 1000, t), report};
}

/** What the first opening says: how many studies are now saved, and the files that could not be read (null: nothing
 *  was brought in and nothing failed). */
export function firstOpenNotice(r: ImportReport): {info: string | null, warning: string | null} {
  const n = r.imported + r.already;
  const renumbered = r.renumbered.length ? ` ${r.renumbered.join(', ')}: saved under a new ID, one with theirs ` +
    'being there already.' : '';
  return {info: n ? `Your ${n} ${n === 1 ? 'plan or study is' : 'plans and studies are'} now saved in Datagrok ` +
    '(the files are left as they were).' + renumbered : null,
  warning: r.damaged.length ? `Could not be read, so not brought in: ${r.damaged.join(', ')}` : null};
}

/** A new lab-wide study ID, PMX-<year>-<number>: the number from the shared counter, never given twice (the year is
 *  the year it is given; the number does not restart with it). */
export async function nextCode(t: Tables = domainTables()): Promise<string> {
  const [r] = await t.code.insert({});
  const n = (await t.code.get(r.id))?.number;
  if (!n)
    throw new Error('No study ID could be given just now: try again.');
  return `PMX-${new Date().getFullYear()}-${String(n).padStart(4, '0')}`;
}

/* What the views use, each behind domainStorageOn(): the row a plan or a study table was opened from or saved to, so
   that the next save updates it at its version. */

/** The tag of a study table with the row it was opened from or saved to ({id, version}). */
export const ROW_TAG = 'pmx.row';
/** The plan rows read in this session, by plan key (name and inputs): a plan opened from one, however it reaches
 *  Plan (the Overview, a picker, the Open strip), saves back to it. */
const planRows = new Map<string, Saved>();
export const rowOfPlan = (p: {name: string, inputs: unknown} | null | undefined): Saved | null =>
  planRows.get(planKey(p) ?? '') ?? null;
export function rowOfTable(df: DG.DataFrame | null | undefined): Saved | null {
  try {
    const x = JSON.parse(df?.getTag(ROW_TAG) || 'null');
    return typeof x?.id === 'string' && Number.isFinite(x.version) ? {id: x.id, version: x.version} : null;
  } catch {
    return null;
  }
}
export const setRowOn = (df: DG.DataFrame, s: Saved) =>
  df.setTag(ROW_TAG, JSON.stringify({id: s.id, version: s.version}));

/** A study row's table: the one already open for that row at that version (QA 0.1.3: each reopen added another copy
 *  to the workspace, "X (2)" … "X (6)"), else a new one, an older copy of the row closed when no table view shows it.
 *  `shown`: whether a page still shows a table (native.ts inUse): one no page shows is never taken up again, however
 *  it was left (a Discard, a page closed), and the row is read as saved; its copy is closed (QA 0.1.3 round 4 chaos
 *  N1: a discarded correction came back with the study, and the next Ctrl+S saved it). */
export function tableOfRow(doc: SavedStudy, row: Saved, tables: DG.DataFrame[] = grok.shell.tables,
  shown?: (df: DG.DataFrame) => boolean): DG.DataFrame {
  const copies = tables.filter((t) => rowOfTable(t)?.id === row.id);
  const same = copies.find((t) => rowOfTable(t)!.version === row.version && (!shown || shown(t)));
  if (same)
    return same;
  for (const t of copies) {
    if (shown ? !shown(t) : ![...grok.shell.tableViews].some((v) => v.dataFrame?.id === t.id))
      grok.shell.closeTable(t);
  }
  const df = tableOf(doc);
  setRowOn(df, row);
  return df;
}

/** `m` with an ID: its own, else a new one (the lab-wide counter with domain storage, else the user's own). */
export const withCode = async (m: StudyMeta): Promise<StudyMeta> =>
  m.id ? m : domainStorageOn() ? {...m, id: await nextCode()} : withId(m);

/** What a refused save says to the user. */
export function saveError(e: any, what: string): string {
  const code = e?.code;
  return code === 'version-conflict' ? `This ${what} was saved elsewhere since you opened it. To keep ` +
    'your changes, save them under another name with File › Save as…; to see that version, open it again from the ' +
    'Overview.' :
    code === 'deleted' ? `This ${what} was deleted (it is in the trash). Restore it from the Overview ` +
    '(Show deleted), or save your changes as a new one with File › Save as….' :
    code === 'forbidden' || code === 'not-found' ? `This ${what} is shared with you to read (or it is not there any ` +
    'more): open it and save a copy of your own.' :
    // (no answer from the server: said in words, not the browser's "ClientException: Failed to fetch, uri=…", QA
    // 0.1.3 round 3 chaos 21)
    offline(e) ? 'Not saved: no connection to the server. Your changes are kept; press Ctrl+S again when you are ' +
    'back online.' :
      `The ${what} could not be saved: ${e?.message ?? e}`;
}

/** Whether an error is the server not answering at all (offline, the server down). */
const offline = (e: any): boolean => /failed to fetch|networkerror|network request failed|load failed/i
  .test(String(e?.message ?? e ?? ''));

/** Whether a save was refused because the row is shared with the user to read only, or is gone (deleted). */
export const isReadOnly = (e: any): boolean => ['forbidden', 'not-found', 'deleted'].includes(e?.code);

/** Save as my copy: `doc`, the user's changes to a row shared with them to read, as a new row of theirs, marked as a
 *  copy wherever it is listed (" (copy)" in its name, its own study ID, the original's row and owner beside it). The
 *  original is not touched. */
export async function saveCopy<T extends SavedPlan | SavedStudy>(kind: Kind, doc: T, source: Saved,
  t: Tables = domainTables(), owner: string = grok.shell.user?.login ?? ''): Promise<{saved: Saved, doc: T}> {
  const [from] = await t[kind].query({columns: ['name', 'owner'], captions: ['owner'], limit: 1,
    filter: [{property: 'id', operator: '=', value: source.id}]}).catch(() => []);
  const code = await nextCode(t);
  // "X (copy)", then "X (copy 2)", ... among the rows the user can see (a copy of a copy copies the original's name)
  const base = doc.name.replace(/ \(copy( \d+)?\)$/, '');
  const taken = new Set((await t[kind].query({columns: ['name'], search: `${base} (copy`, limit: 1000}))
    .map((r) => r.name));
  let n = 1;
  while (taken.has(`${base} (copy${n === 1 ? '' : ` ${n}`})`))
    n++;
  // the copy is the user's: they own it; whose the original is goes beside it (copy_owner, the user) and in its
  // details (their name)
  const original = from?.owner || from?.author_id || null;
  const originalName = from?.['~caption_owner'] ?? '';
  const copy = withMeta(kind, {...doc, name: `${base} (copy${n === 1 ? '' : ` ${n}`})`} as T,
    {id: code, owner, ...originalName ? {copyOf: originalName} : {}});
  const [r] = await t[kind].insert({...columnsOf(kind, copy), owner: meId(), copy_of: source.id,
    copy_owner: original});
  const saved = {id: r.id, version: r.version ?? 1};
  if (kind === 'plan')
    planRows.set(planKey(copy as SavedPlan)!, saved);
  return {saved, doc: copy};
}

/** The user's plans and studies for the Overview: rows with domain storage, else the files. */
export const myEntries = (): Promise<LibraryEntry[]> => domainStorageOn() ? list() : myLibrary();

/** An entry's document text, read now: its row's, else its file's (the list holds none). */
export async function entryText(e: LibraryEntry, t?: Tables): Promise<string> {
  if (e.text)
    return e.text;
  if (e.row) {
    const row = await (t ?? domainTables())[e.row.kind].get(e.row.id);
    if (!row)
      throw new Error(`"${e.name}" is not there any more, or it is no longer shared with you.`);
    if (e.row.kind === 'plan') {
      try {
        planRows.set(planKey(parsePlan(docText(row.document)))!, {id: row.id, version: row.version});
      } catch {} // (a damaged document is refused where it is opened)
    }
    return docText(row.document);
  }
  return e.file ? grok.dapi.files.readAsText(e.file) : '';
}

/** Tables in memory, with the server's rules the store relies on: ids and versions stamped, an update at another
 *  version refused (code 'version-conflict'), a required name and document, delete out of every query, the code
 *  table numbered from 1, every row authored by `author`, a json document's keys reordered as jsonb does (here:
 *  sorted), a user column's caption its name in `names` (else the id). For the package's tests. */
export function memoryTables(author = 'me', names: {[id: string]: string} = {}): Tables {
  const jsonb = (x: unknown): unknown => Array.isArray(x) ? x.map(jsonb) : x && typeof x === 'object' ?
    Object.fromEntries(Object.keys(x).sort().reverse().map((k) => [k, jsonb((x as any)[k])])) : x;
  const stored = (r: any) => r.document && typeof r.document === 'object' ? {...r, document: jsonb(r.document)} : r;
  let clock = 0;
  const table = <TRow extends {id: string, version: number}>(numbered = false): StoreTable<TRow> => {
    const rows = new Map<string, TRow>();
    const trash = new Map<string, TRow>();
    let counter = 0;
    const err = (code: string, message: string) => Object.assign(new Error(message), {code});
    const check = (r: any) => {
      if (!numbered && (!r.name || !r.document))
        throw err('validation', 'Value can\'t be empty');
    };
    return {
      async insert(row) {
        const r = stored({...row, id: crypto.randomUUID(), version: 1, updated_on: ++clock, author_id: author,
          ...numbered ? {number: ++counter} : {}}) as unknown as TRow;
        check(r);
        rows.set(r.id, r);
        return [{id: r.id, version: 1}];
      },
      async update(id, values, {version}) {
        const was = rows.get(id);
        if (!was)
          throw err('not-found', `No row ${id}`);
        if (was.version !== version)
          throw err('version-conflict', `Expected version ${version}, the row is at ${was.version}`);
        const r = stored({...was, ...values, id, version: version + 1, updated_on: ++clock});
        check(r);
        rows.set(id, r);
        return {id, version: r.version};
      },
      async get(id, options) {
        const r = rows.get(id) as any;
        return r ? {...r, ...options?.withAccess ? {'~can_edit': r.author_id === author} : {}} : null;
      },
      async query({columns, sort, search, limit, filter, deleted, withAccess, captions}) {
        const words = (search ?? '').toLowerCase();
        const hit = (r: any) => (!words || ['name', 'code', 'line', 'model', 'run_by', 'treatment']
          .some((c) => String(r[c] ?? '').toLowerCase().includes(words))) &&
          (filter ?? []).every((c) => c === 'and' || r[c.property] === c.value);
        const desc = sort?.startsWith('!');
        const key = sort?.replace('!', '');
        const out = [...(deleted ? trash : rows).values()].filter(hit)
          .sort((a: any, b: any) => !key ? 0 : (a[key] < b[key] ? -1 : a[key] > b[key] ? 1 : 0) * (desc ? -1 : 1))
          .slice(0, limit ?? Infinity);
        // (as the server: a system column in `columns` is refused, and the system columns always come)
        const SYSTEM = ['id', 'version', 'created_on', 'updated_on', 'author_id'];
        if (columns?.some((c) => SYSTEM.includes(c)))
          throw err('validation', `Unknown or inaccessible column '${columns.find((c) => SYSTEM.includes(c))}'`);
        // (withAccess: the author edits; a row stamped by another is shared to read here)
        const access = (r: any) => withAccess ? {'~can_edit': r.author_id === author} : {};
        const picked = (r: any) => columns ? Object.fromEntries([...SYSTEM, ...columns].map((c) => [c, r[c]])) : r;
        const named = (r: any) => Object.fromEntries((captions ?? []).map((c) =>
          [`~caption_${c}`, r[c] == null ? null : names[r[c]] ?? r[c]]));
        return out.map((r: any) => ({...picked(r), ...access(r), ...named(r)}) as TRow);
      },
      async delete(id) {
        const r = rows.get(id);
        if (!r)
          throw err('not-found', `No row ${id}`);
        rows.delete(id);
        trash.set(id, r);
      },
      async restore(id) {
        const r = trash.get(id);
        if (!r)
          throw err('not-found', `No deleted row ${id}`);
        trash.delete(id);
        rows.set(id, r);
      },
    };
  };
  return {plan: table<DocRow>(), study: table<DocRow>(), code: table<CodeRow>(true)};
}

/** The user's saved plans (or `kind`, studies) as a picker lists them: by name, each read when picked (rows, or the
 *  files). `entry`: it as the Overview lists it, to open as its row opens (openEntryOnce). */
export async function savedPlans(t?: Tables, kind: Kind = 'plan'):
  Promise<{name: string, read: () => Promise<DG.FileInfo>, entry: LibraryEntry}[]> {
  if (!domainStorageOn()) {
    const files = kind === 'plan' ? await myPlans() : await grok.dapi.files.list(studiesFolder(), false, STUDY_EXT);
    return files.map((f) => ({name: f.name, read: async () => f,
      entry: {kind, name: f.name, file: f, text: ''} as unknown as LibraryEntry}));
  }
  const rows = await (t ?? domainTables())[kind].query({columns: ['name', 'code'], sort: '!updated_on', limit: 1000});
  return rows.map((r) => {
    const entry = {kind, row: {kind, id: r.id, version: r.version}, name: r.name, file: null,
      text: ''} as unknown as LibraryEntry;
    return {name: r.code ? `${r.name} (${r.code})` : r.name, entry, read: async () => DG.FileInfo.fromString(
      `${r.name}.${kind === 'plan' ? PLAN_EXT : STUDY_EXT}`, await entryText(entry, t))};
  });
}
