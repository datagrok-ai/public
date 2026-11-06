/* The study as a thing of its own, from its plan to its report: an ID (PMX-2026-014), a description, who owns it, its
   project and tags, the day it starts (day 0, randomization), and what each group receives (compound, dose, route,
   schedule). It is written once in Plan and travels with the study: in the .pmxplan (SavedPlan.study), on the
   measurements table (META_TAG) and in the .pmxstudy (SavedStudy.meta). The Overview lists studies by it. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** What one group receives. */
export interface GroupDetails {
  group: string, compound: string, dose: number | null, unit: string, route: string, schedule: string,
}

export interface StudyMeta {
  /** PMX-<year>-<number>, given when the study is first saved or followed (empty until then) */
  id: string;
  description: string;
  /** the login of who runs it */
  owner: string;
  project: string;
  tags: string[];
  /** day 0 (randomization), YYYY-MM-DD; null: not fixed yet */
  start: string | null;
  groups: GroupDetails[];
  /** the study's notebook: what happened on a day, by whom; it never changes the data */
  notes: DayNote[];
  /** whoever enters the measurements must not see the groups: Monitor hides them until someone reveals them */
  blinded: boolean;
}

export const NOTE_KINDS = ['Observation', 'Incident', 'Dose', 'Measurement'] as const;
export type NoteKind = typeof NOTE_KINDS[number];

/** A note on a study day: its kind, the animal it is about (empty: the study), what, who and when. */
export interface DayNote {day: number, kind: NoteKind, animal: string, text: string, who: string, when: string}

/** The tag of a measurements table that carries its study's details (JSON). */
export const META_TAG = 'pmx.meta';
/** The tag that carries the plan (plan-file.ts's PLAN_TAG; here to read its details without a cycle). */
const PLAN_TAG = 'pmx.plan';

export const emptyMeta = (): StudyMeta => ({id: '', description: '', owner: grok.shell.user?.login ?? '',
  project: '', tags: [], start: null, groups: [], notes: [], blinded: false});

/** Details read from a file or a tag someone may have edited: each field checked, the rest from emptyMeta. */
export function parseMeta(x: unknown): StudyMeta | null {
  if (!x || typeof x !== 'object')
    return null;
  const m = x as Partial<StudyMeta>;
  const text = (v: unknown) => typeof v === 'string' ? v : '';
  return {...emptyMeta(), id: text(m.id), description: text(m.description), owner: text(m.owner), project:
    text(m.project), tags: Array.isArray(m.tags) ? m.tags.filter((t) => typeof t === 'string') : [],
  start: typeof m.start === 'string' && /^\d{4}-\d{2}-\d{2}$/.test(m.start) ? m.start : null,
  groups: Array.isArray(m.groups) ? m.groups.filter((g) => g && typeof g.group === 'string').map((g) => ({
    group: g.group, compound: text(g.compound), dose: typeof g.dose === 'number' && Number.isFinite(g.dose) ?
      g.dose : null, unit: text(g.unit), route: text(g.route), schedule: text(g.schedule)})) : [],
  notes: Array.isArray(m.notes) ? m.notes.filter((n) => n && Number.isInteger(n.day) && typeof n.text === 'string')
    .map((n) => ({day: n.day, kind: NOTE_KINDS.includes(n.kind) ? n.kind : 'Observation', animal: text(n.animal),
      text: n.text, who: text(n.who), when: text(n.when)})) : [], blinded: m.blinded === true};
}

/** The date of a study day, from day 0's (null: no date yet); UTC, so a day is a day wherever it is read. */
export function dateOf(start: string | null, day: number): Date | null {
  if (!start)
    return null;
  const d = new Date(`${start}T00:00:00Z`);
  d.setUTCDate(d.getUTCDate() + day);
  return d;
}

/** A date as the calendar writes it: 12 Nov (and the year when `year`). */
export const dateText = (d: Date, year = false) => d.toLocaleDateString('en-GB',
  {day: 'numeric', month: 'short', ...year ? {year: 'numeric'} : {}, timeZone: 'UTC'});

/** The weekday of a date, Monday 0 (as Plan's weekdays). */
export const weekdayOf = (iso: string) => (new Date(`${iso}T00:00:00Z`).getUTCDay() + 6) % 7;

/** The study day 0 falls on when `day` is the date `iso`. */
export function startFrom(iso: string, day: number): string {
  const d = new Date(`${iso}T00:00:00Z`);
  d.setUTCDate(d.getUTCDate() - day);
  return d.toISOString().slice(0, 10);
}

/** The days each group is dosed, up to `last`, from its schedule as protocols write it: QD×14, Q3D×4, QW×3, Q2W×2
 *  (every k days or weeks, n doses, the first on day 0). A schedule written otherwise (BIW, "days 0-4") marks
 *  nothing. ponytail: dosing always starts on day 0; a "from day N" needs a start field per group. */
export function dosingDays(groups: GroupDetails[], last: number): Map<number, string[]> {
  const days = new Map<number, string[]>();
  for (const g of groups) {
    const m = /^\s*Q(\d*)\s*(D|W)\s*[×xX*]\s*(\d+)\s*$/i.exec(g.schedule);
    if (!m)
      continue;
    const every = (Number(m[1]) || 1) * (m[2].toUpperCase() === 'W' ? 7 : 1);
    for (let k = 0; k < Number(m[3]) && k * every <= last; k++)
      days.set(k * every, [...days.get(k * every) ?? [], g.group]);
  }
  return days;
}

/** A note on `day` (`label`: the day as the calendar names it), asked for in a dialog; `save` gets it. */
export function noteDialog(day: number, label: string, animals: string[], save: (n: DayNote) => void): void {
  const kind = ui.input.choice('Kind', {value: 'Observation', items: [...NOTE_KINDS], nullable: false});
  const animal = ui.input.choice('Animal', {value: '', items: ['', ...animals], nullable: true,
    tooltipText: 'Empty: the note is about the study'});
  const text = ui.input.textArea('Note', {value: ''});
  const d = ui.dialog(`Note on ${label}`).add(ui.inputs([kind, animal, text]));
  const ok = () => d.getButton('OK');
  text.onChanged.subscribe(() => ok().disabled = text.value.trim() === '');
  d.onOK(() => save({day, kind: kind.value as NoteKind, animal: animal.value ?? '', text: text.value.trim(),
    who: grok.shell.user?.login ?? '', when: new Date().toISOString()})).show();
  ok().disabled = true;
}

/** A table's study details: its own, else its plan's (null: none yet). */
export function metaOf(df: DG.DataFrame | null | undefined): StudyMeta | null {
  const read = (tag: string, pick: (x: any) => unknown) => {
    try {
      return parseMeta(pick(JSON.parse(df?.getTag(tag) || 'null')));
    } catch {
      return null;
    }
  };
  return read(META_TAG, (x) => x) ?? read(PLAN_TAG, (p) => p?.study);
}

export const setMetaOn = (df: DG.DataFrame, m: StudyMeta) => df.setTag(META_TAG, JSON.stringify(m));

/** The next study ID of this user this year: PMX-2026-001, -002, ... ponytail: the counter is the user's own
 *  (grok.userSettings); a lab-wide sequence needs a shared counter (package settings) once studies are shared. */
export function nextStudyId(): string {
  const year = new Date().getFullYear();
  const key = `lastStudyId-${year}`;
  const n = (Number(grok.userSettings.getValue('PmxStudio', key)) || 0) + 1;
  grok.userSettings.add('PmxStudio', key, String(n));
  return `PMX-${year}-${String(n).padStart(3, '0')}`;
}

/** `m` with an ID (a new one when it has none). */
export const withId = (m: StudyMeta): StudyMeta => m.id ? m : {...m, id: nextStudyId()};

/** The groups' details in the order of `groups` (the study's groups now), each kept from `m` by name. */
export function groupsOf(m: StudyMeta, groups: string[]): GroupDetails[] {
  return groups.map((group) => m.groups.find((g) => g.group === group) ??
    {group, compound: '', dose: null, unit: 'mg/kg', route: '', schedule: ''});
}

/** One line on what the groups receive, for a list: the compounds, each once. */
export const compoundsOf = (m: StudyMeta | null) =>
  [...new Set((m?.groups ?? []).map((g) => g.compound.trim()).filter((c) => c))].join(', ');

/** What the groups receive, each compound once with its dose, route and schedule ('anti-PD-1 10 mg/kg i.p. Q3D×4 +
 *  …'); empty when no group names a compound. */
export const treatmentOf = (m: StudyMeta | null) => [...new Set((m?.groups ?? []).filter((g) => g.compound.trim())
  .map((g) => [g.compound.trim(), g.dose === null ? '' : `${g.dose} ${g.unit}`.trim(), g.route, g.schedule]
    .filter((s) => s).join(' ')))].join(' + ');

/** The study's details to edit: its fields, and its groups as a table (one row per group, typed or pasted into like
 *  any Datagrok grid); `read` gives them as edited. An ID left empty is given when the study is first saved. */
export function detailsForm(meta: StudyMeta, groups: string[]): {root: HTMLElement, read: () => StudyMeta} {
  const id = ui.input.string('ID', {value: meta.id, placeholder: 'Given when the study is first saved'});
  const description = ui.input.textArea('Description', {value: meta.description});
  const owner = ui.input.string('Owner', {value: meta.owner});
  const project = ui.input.string('Project', {value: meta.project});
  const tags = ui.input.string('Tags', {value: meta.tags.join(', '),
    tooltipText: 'Comma-separated, e.g. syngeneic, IO, combination'});
  const start = ui.input.date('Day 0', {tooltipText: 'Randomization, the first day of treatment'});
  if (meta.start)
    start.stringValue = meta.start;
  const blinded = ui.input.bool('Blinded entry', {value: meta.blinded, tooltipText: 'Whoever enters the ' +
    'measurements must not see the groups: Monitor hides them until someone reveals them, and the log says so'});
  const rows = groupsOf(meta, groups);
  const table = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('group', rows.map((g) => g.group)),
    DG.Column.fromStrings('compound', rows.map((g) => g.compound)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'dose', rows.map((g) => g.dose)),
    DG.Column.fromStrings('unit', rows.map((g) => g.unit)),
    DG.Column.fromStrings('route', rows.map((g) => g.route)),
    DG.Column.fromStrings('schedule', rows.map((g) => g.schedule))]);
  table.col('route')!.meta.choices = ['p.o.', 'i.p.', 'i.v.', 's.c.', 'i.t.'];
  const grid = DG.Viewer.grid(table, {allowEdit: true, showAddNewRowIcon: false, showRemoveRowIcon: false});
  grid.root.style.width = '640px';
  grid.root.style.height = `${Math.min(8, rows.length + 1) * 30 + 8}px`;
  const root = ui.divV([ui.inputs([id, description, owner, project, tags, start, blinded], 'pmx-details-fields'),
    ui.divText('What each group receives (the schedule as written in the protocol, e.g. Q3D×4 or QD×14):'),
    grid.root], 'pmx-details');
  const read = (): StudyMeta => {
    const day = (start.value as unknown as {format(f: string): string} | null)?.format('YYYY-MM-DD') ?? null;
    const cell = (c: string, i: number) => table.col(c)!.isNone(i) ? '' : String(table.get(c, i)).trim();
    return {id: id.value.trim() || meta.id, description: description.value.trim(), owner: owner.value.trim(),
      project: project.value.trim(), tags: tags.value.split(',').map((t) => t.trim()).filter((t) => t),
      start: day, groups: rows.map((g, i) => ({group: g.group, compound: cell('compound', i),
        dose: table.col('dose')!.isNone(i) ? null : table.get('dose', i), unit: cell('unit', i),
        route: cell('route', i), schedule: cell('schedule', i)})), notes: meta.notes,
      blinded: blinded.value};
  };
  return {root, read};
}
