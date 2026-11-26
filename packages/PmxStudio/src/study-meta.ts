/* The study as a thing of its own, from its plan to its report: an ID (PMX-2026-014), a description, who owns it, its
   project and tags, the day it starts (day 0, randomization), and what each group receives (compound, dose, route,
   schedule). It is written once in Plan and travels with the study: in the .pmxplan (SavedPlan.study), on the
   measurements table (META_TAG) and in the .pmxstudy (SavedStudy.meta). The Overview lists studies by it. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {okNamed, setUserSetting, showCentred, userSetting} from './ui';
import {domainStorageOn} from './storage-setting';
import {ArmPart, controlByName, describeArm, drugKey, parseArmName} from './arm-name';
import {DoseException, DosingRules, Schedule, doseLine, dosedDays, doseTotal, excepted, readSchedule, validDay}
  from './schedule';

/** The animals, as ARRIVE (item 8) and PREPARE (topic 10) ask for them; free text: '8 weeks', '20 g'. */
export interface AnimalDetails {species: string, strain: string, sex: string, age: string, weight: string,
  source: string}

/** The animal fields, each also the word a checklist names it by. */
const ANIMAL_FIELDS: (keyof AnimalDetails)[] = ['species', 'strain', 'sex', 'age', 'weight', 'source'];

/** The animals in one line ('mouse, C57BL/6, female, 8 weeks, 20 g, Charles River') and the fields still empty. */
export function animalsLine(a: AnimalDetails | undefined): {filled: string, missing: string[]} {
  const has = ANIMAL_FIELDS.filter((k) => (a?.[k] ?? '').trim());
  return {filled: has.map((k) => a![k].trim()).join(', '), missing: ANIMAL_FIELDS.filter((k) => !has.includes(k))};
}

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
  /** no dosing at weekends (schedule.ts DosingRules: daily doses skip them, others move to the nearest weekday);
   *  needs day 0's date, without which it changes nothing */
  weekendsOff: boolean;
  /** doses skipped or added on a day, for a group (schedule.ts excepted): the schedule as written never changes */
  exceptions: DoseException[];
  /** the doses given, as Monitor logs them (a day and a group); a missed dose is an exception, not a missing log */
  given: DoseGiven[];
  /** the animals used (old files have none: a mouse, the rest empty) */
  animals: AnimalDetails;
  /** the colours the user gave groups in a chart (colors.ts); a group without one has its automatic colour */
  colors: {[group: string]: string};
  /** each treated group's place in the palette, given when the group first appears in Plan (colors.ts slotsAfter):
   *  adding or renaming an arm recolours no other; a group without one takes a free place */
  slots?: {[group: string]: number};
  /** a copy of someone else's study (Save as my copy, or their file brought in): whose the original is */
  copyOf?: string;
  /** the arms stopped on a review and the review's day: read from the study's log (stopsOf) by metaOf, never
   *  read back from the details (parseMeta), so the log stays the one record of the decision; no dose after
   *  that day */
  stopped?: {[group: string]: number};
}

/** A group's doses on a day logged as given, by whom and when. */
export interface DoseGiven {day: number, group: string, who: string, when: string}

export const NOTE_KINDS = ['Observation', 'Incident', 'Dose', 'Measurement'] as const;
export type NoteKind = typeof NOTE_KINDS[number];

/** A note on a study day: its kind, the animal it is about (empty: the study), what, who and when. */
export interface DayNote {day: number, kind: NoteKind, animal: string, text: string, who: string, when: string}

/** The tag of a measurements table that carries its study's details (JSON). */
export const META_TAG = 'pmx.meta';
/** The tag that carries the plan (plan-file.ts's PLAN_TAG; here to read its details without a cycle). */
const PLAN_TAG = 'pmx.plan';
/** The study's log (study-file.ts's LOG_TAG, likewise): the decisions on its reviews. */
const LOG_TAG = 'pmx.log';
/** A value read back from a 32-bit Datagrok column, to the 7 digits it holds (0.3, not 0.30000001192…). */
export const f32 = (x: number): number => Number(x.toPrecision(7));

export const emptyMeta = (): StudyMeta => ({id: '', description: '', owner: grok.shell.user?.login ?? '',
  project: '', tags: [], start: null, groups: [], notes: [], blinded: false, weekendsOff: false, exceptions: [],
  given: [], animals: {species: 'mouse', strain: '', sex: '', age: '', weight: '', source: ''}, colors: {}});

/** Details read from a file or a tag someone may have edited: each field checked, the rest from emptyMeta. */
export function parseMeta(x: unknown): StudyMeta | null {
  if (!x || typeof x !== 'object')
    return null;
  const m = x as Partial<StudyMeta>;
  const text = (v: unknown) => typeof v === 'string' ? v : '';
  return {...emptyMeta(), id: text(m.id), description: text(m.description), owner: text(m.owner), project:
    text(m.project), tags: Array.isArray(m.tags) ? m.tags.filter((t) => typeof t === 'string') : [],
  // a date that exists: 2026-13-45 has the shape and would empty the calendar
  start: typeof m.start === 'string' ? validDay(m.start) : null,
  groups: Array.isArray(m.groups) ? m.groups.filter((g) => g && typeof g.group === 'string').map((g) => ({
    group: g.group, compound: text(g.compound), dose: typeof g.dose === 'number' && Number.isFinite(g.dose) ?
      g.dose : null, unit: text(g.unit), route: text(g.route), schedule: text(g.schedule)})) : [],
  notes: Array.isArray(m.notes) ? m.notes.filter((n) => n && Number.isInteger(n.day) && typeof n.text === 'string')
    .map((n) => ({day: n.day, kind: NOTE_KINDS.includes(n.kind) ? n.kind : 'Observation', animal: text(n.animal),
      text: n.text, who: text(n.who), when: text(n.when)})) : [], blinded: m.blinded === true,
  weekendsOff: m.weekendsOff === true,
  exceptions: list(m.exceptions).filter((e) => e.change === 'skip' || e.change === 'add').map((e) => ({day: e.day,
    group: e.group, change: e.change, note: text(e.note), who: text(e.who), when: text(e.when)})),
  given: list(m.given).map((g) => ({day: g.day, group: g.group, who: text(g.who), when: text(g.when)})),
  animals: Object.fromEntries(ANIMAL_FIELDS.map((k) => [k, text((m.animals as any)?.[k]) ||
    (k === 'species' ? 'mouse' : '')])) as unknown as AnimalDetails,
  colors: m.colors && typeof m.colors === 'object' ? Object.fromEntries(Object.entries(m.colors)
    .filter(([, c]) => typeof c === 'string' && /^#[0-9a-f]{6}$/i.test(c))) : {},
  ...m.slots && typeof m.slots === 'object' ? {slots: Object.fromEntries(Object.entries(m.slots)
    .filter(([, k]) => Number.isInteger(k) && k >= 0))} : {},
  ...text(m.copyOf) ? {copyOf: text(m.copyOf)} : {}};
}

/** The entries of a dated list read from a file: an object with a whole day and a group. */
const list = (x: unknown): any[] => Array.isArray(x) ?
  x.filter((e) => e && Number.isInteger(e.day) && typeof e.group === 'string') : [];

/** Who writes an entry, and when: now. */
const stamp = () => ({who: grok.shell.user?.login ?? '', when: new Date().toISOString()});

/** `m` with `group`'s exception on `day` set to `change` (null: none); one a day and group, the latest wins. */
export function withException(m: StudyMeta, day: number, group: string, change: 'skip' | 'add' | null, note = ''):
  StudyMeta {
  const rest = m.exceptions.filter((e) => !(e.day === day && e.group === group));
  return {...m, exceptions: change === null ? rest : [...rest, {day, group, change, note, ...stamp()}]};
}

/** `m` with the doses of `day` logged as given for exactly `groups` (the other groups' logs of that day removed). */
export function withGiven(m: StudyMeta, day: number, groups: string[]): StudyMeta {
  const kept = m.given.filter((g) => g.day !== day || groups.includes(g.group));
  const fresh = groups.filter((g) => !kept.some((k) => k.day === day && k.group === g))
    .map((group) => ({day, group, ...stamp()}));
  return {...m, given: [...kept, ...fresh]};
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

/** A group's schedules as read (schedule.ts readSchedule): its details' schedule, else each of its name's (a
 *  combination may have two); continuous exposure from the schedule or the details' route, never from the name. */
export function schedulesOf(g: GroupDetails): Schedule[] {
  const texts = g.schedule.trim() ? [g.schedule] : parseArmName(g.group).map((p) => p.schedule ?? '');
  return (texts.length ? texts : ['']).map((t) => readSchedule(t, g.route)).filter((s): s is Schedule => s !== null);
}

/** The days each group is dosed, up to `last`, from its schedules as protocols write them (schedule.ts readSchedule:
 *  QD, Q3D×4, BIW, from day 3, d3–d16, 5 on / 2 off; continuous exposure every day), under the study's rules
 *  (`rules`: day 0's date and no dosing at weekends). A schedule it cannot read marks nothing. */
export function dosingDays(groups: GroupDetails[], last: number,
  rules: (DosingRules & Pick<StudyMeta, 'stopped'>) | null = null): Map<number, string[]> {
  // a group's days from all its rows (a combination has one per drug), each day once
  const own = new Map<string, Set<number>>();
  for (const g of groups) {
    const ds = own.get(g.group) ?? new Set<number>();
    schedulesOf(g).flatMap((s) => dosedDays(s, last, rules)).forEach((d) => ds.add(d));
    own.set(g.group, ds);
  }
  const days = new Map<number, string[]>();
  for (const [group, ds] of own) {
    // the study's exceptions for the group; one with no schedule read is not dosed, whatever exceptions it kept
    // (a stopped arm's doses end on its review's day)
    for (const d of ds.size ? excepted([...ds].sort((a, b) => a - b), group, rules, last) : []) {
      if (d <= (rules?.stopped?.[group] ?? Infinity))
        days.set(d, [...days.get(d) ?? [], group]);
    }
  }
  return days;
}

/** A note on `day` (`label`: the day as the calendar names it), asked for in a dialog; `save` gets it. */
export function noteDialog(day: number, label: string, animals: string[], save: (n: DayNote) => void): void {
  const kind = ui.input.choice('Kind', {value: 'Observation', items: [...NOTE_KINDS], nullable: false});
  // one empty choice, the study itself ('' with nullable gave two); a plan has no animals yet, so no field at all
  const animal = ui.input.choice<string>('Animal', {value: '', items: ['', ...animals], nullable: false,
    tooltipText: 'Empty: the note is about the study'});
  const text = ui.input.textArea('Note', {value: ''});
  // what is noted first (typed at once), then its kind and, if any, its animal
  const d = ui.dialog(`Note on ${label}`).add(ui.inputs([text, kind, ...(animals.length ? [animal] : [])]));
  d.onOK(() => save({day, kind: kind.value as NoteKind, animal: animal.value ?? '', text: text.value.trim(),
    ...stamp()}));
  const ok = okNamed(d, 'Add note');
  text.onChanged.subscribe(() => ok.disabled = text.value.trim() === '');
  ok.disabled = true;
  showCentred(d);
  text.input.focus();
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
  const m = read(META_TAG, (x) => x) ?? read(PLAN_TAG, (p) => p?.study);
  if (m) {
    try {
      const stopped = stopsOf(JSON.parse(df?.getTag(LOG_TAG) || '[]'));
      if (Object.keys(stopped).length)
        m.stopped = stopped;
    } catch {
      // a log edited by hand: no stops read from it
    }
  }
  return m;
}

/** The log's entry for a decision on a review: `arm` stopped, its stop undone, or continued. */
export const decisionEntry = (arm: string, review: number, what: 'stop' | 'stop undone' | 'continue') =>
  `${arm}: ${what} (review on day ${review})`;

/** The arms stopped as the log tells it, in its order (a stop undone takes it back), each with its review's day. */
export function stopsOf(log: {what: string}[]): {[group: string]: number} {
  const out: {[group: string]: number} = {};
  for (const l of Array.isArray(log) ? log : []) {
    const m = /^(.+): (stop|stop undone) \(review on day (-?[\d.]+)\)$/.exec(String(l?.what ?? ''));
    if (m && m[2] === 'stop')
      out[m[1]] = Number(m[3]);
    else if (m)
      delete out[m[1]];
  }
  return out;
}

export const setMetaOn = (df: DG.DataFrame, m: StudyMeta) => df.setTag(META_TAG, JSON.stringify(m));

/** The next study ID of this user this year: PMX-2026-001, -002, ... ponytail: the counter is the user's own
 *  (grok.userSettings); a lab-wide sequence needs a shared counter (package settings) once studies are shared. */
export function nextStudyId(): string {
  const year = new Date().getFullYear();
  const key = `lastStudyId-${year}`;
  const n = (Number(userSetting(key)) || 0) + 1;
  setUserSetting(key, String(n)); // sent at once: a reload before the 10 s sync would give the same ID again
  return `PMX-${year}-${String(n).padStart(3, '0')}`;
}

/** `m` with an ID (a new one when it has none). */
export const withId = (m: StudyMeta): StudyMeta => m.id ? m : {...m, id: nextStudyId()};

/** The groups' details in the order of `groups` (the study's groups now), each kept from `m` by name: one row per
 *  compound (rowsOf). */
export function groupsOf(m: StudyMeta, groups: string[]): GroupDetails[] {
  return groups.flatMap((group) => rowsOf(m, group));
}

/** A group's rows of details, one per compound: as many as its name gives drugs (one for a single drug or none):
 *  those kept (beyond that number left out: a combination renamed to one drug), then a row for each drug no kept row
 *  gives; each empty cell read from its name. */
export function rowsOf(m: StudyMeta | null | undefined, group: string): GroupDetails[] {
  const parts = parseArmName(group);
  const n = Math.max(1, parts.length);
  const kept = (m?.groups ?? []).filter((g) => g.group === group).slice(0, n);
  // the drug each kept row gives: its compound's, else the one at its place
  const given = new Set(kept.map((g, i) => g.compound.trim() ? parts.findIndex((p) => p.key === drugKey(g.compound)) :
    i));
  const free = parts.map((_, i) => i).filter((i) => !given.has(i));
  return [...kept.map((g, i) => fromName(g, i)), ...Array.from({length: n - kept.length}, (_, j) =>
    fromName({group, compound: '', dose: null, unit: '', route: '', schedule: ''}, free[j] ?? kept.length + j))];
}

/** A group's details with each empty cell read from its name: its drug, or for a combination the drug the row's
 *  compound names (none typed: the `i`th drug of the name). A cell the user typed is never changed. */
export function fromName(g: GroupDetails, i = 0): GroupDetails {
  const parts = parseArmName(g.group);
  const typed = g.compound.trim() ? drugKey(g.compound) : null;
  const p = parts.length === 1 ? parts[0] : typed !== null ? parts.find((q) => q.key === typed) ?? null :
    parts[i] ?? null;
  const dosed = g.dose === null && p?.dose;
  return {...g, compound: g.compound || (p?.drug ?? ''), dose: dosed ? p!.dose!.value : g.dose,
    unit: g.unit || (dosed ? p!.dose!.unit : 'mg/kg'), route: g.route || (p?.route ?? ''),
    schedule: g.schedule || (p?.schedule ?? '')};
}

/** The groups renamed between `before` and `now` (each group's animals, by ID), old name to new: a group gone and one
 *  new that hold exactly the same animals, and no other group does. A group removed and another added with other
 *  animals is no rename. */
export function renamedGroups(before: Map<string, string[]>, now: Map<string, string[]>): Map<string, string> {
  const key = (animals: string[]) => [...new Set(animals)].sort().join('\u0000');
  const fresh = [...now].filter(([g, a]) => !before.has(g) && a.length).map(([g, a]) => [g, key(a)]);
  const gone = [...before].filter(([g, a]) => !now.has(g) && a.length).map(([g, a]) => [g, key(a)]);
  const out = new Map<string, string>();
  for (const [from, k] of gone) {
    const to = fresh.filter(([, f]) => f === k);
    if (to.length === 1 && gone.filter(([, o]) => o === k).length === 1)
      out.set(from, to[0][0]);
  }
  return out;
}

/** A row's details without what was only read from its name (fromName; `i`, its place among the group's rows), so
 *  the new name fills those again. */
function typedOnly(g: GroupDetails, i: number): GroupDetails {
  const n = fromName({group: g.group, compound: '', dose: null, unit: '', route: '', schedule: ''}, i);
  const sameDose = g.dose === n.dose && g.unit === n.unit;
  return {...g, compound: g.compound === n.compound ? '' : g.compound, dose: sameDose ? null : g.dose,
    unit: sameDose ? '' : g.unit, route: g.route === n.route ? '' : g.route,
    schedule: g.schedule === n.schedule ? '' : g.schedule};
}

/** `m` with the details of each group in `renames` (old name to new: renamedGroups, or Plan's edit of an arm's name)
 *  under its new name, only what was typed kept (its new name fills the rest); one whose new name has details of its
 *  own stays as it is. */
export function renamedDetails(m: StudyMeta, renames: Map<string, string>): StudyMeta {
  if (!renames.size)
    return m;
  const moved = (g: {group: string}) => renames.has(g.group) && !m.groups.some((d) => d.group === renames.get(g.group));
  // a renamed group's exceptions and doses given go with it, one a day and group (a moved one wins)
  const to = <T extends {group: string}>(x: T): T => moved(x) ? {...x, group: renames.get(x.group)!} : x;
  const same = (a: {day: number, group: string}, b: {day: number, group: string}) =>
    a.day === b.day && a.group === b.group;
  const exceptions = m.exceptions.map(to).filter((e, i, all) => !all.slice(i + 1).some((x) => same(x, e)));
  const given = m.given.map(to).filter((g, i, all) => all.findIndex((x) => same(x, g)) === i);
  // each row's place among its group's rows
  const seen = new Map<string, number>();
  const at = (g: GroupDetails) => {
    const i = seen.get(g.group) ?? 0;
    seen.set(g.group, i + 1);
    return i;
  };
  // a renamed group keeps the colour the user gave it, unless its new name has one
  const colors = {...m.colors};
  for (const [from, to] of renames) {
    if (from in m.colors && !(to in m.colors))
      colors[to] = m.colors[from];
    delete colors[from];
  }
  // and its place in the palette
  const slots = m.slots && {...m.slots};
  for (const [from, to] of renames) {
    if (slots && from in slots && !(to in slots))
      slots[to] = slots[from];
    delete slots?.[from];
  }
  return {...m, ...slots ? {slots} : {}, groups: m.groups.map((g) => {
    const i = at(g);
    return moved(g) ? {...typedOnly(g, i), group: renames.get(g.group)!} : g;
  }), exceptions, given, colors};
}

/** `m` with a row of details for each of `groups` (none given: as it is), its empty cells read from the name; a group
 *  in `renames` (renamedDetails) keeps what was typed for it under its new name; the rows of other groups kept,
 *  though no group reads them (ofGroups leaves them out). */
export function withGroups(m: StudyMeta, groups: string[], renames = new Map<string, string>()): StudyMeta {
  if (!groups.length)
    return m;
  const base = renamedDetails(m, renames);
  return {...base, groups: [...groupsOf(base, groups), ...base.groups.filter((g) => !groups.includes(g.group))]};
}

/** The details of the study's groups only (all of them when the groups are not known). */
export const ofGroups = (m: StudyMeta, groups?: string[]): GroupDetails[] =>
  groups ? m.groups.filter((g) => groups.includes(g.group)) : m.groups;

/** The name armParts reads for a group: its details written out as a name, a combination's rows joined by ' + '
 *  (each empty cell from its name), when the group has a row (its compound empty after a rename) and they name a
 *  compound; else the name itself (also a combination's whose rows name fewer compounds than its name). */
export function detailsName(group: string, meta: StudyMeta | null | undefined): string {
  // a row's compound may be empty and still hold what was typed: a rename leaves it for the new name to fill
  if (!meta?.groups.some((g) => g.group === group))
    return group;
  const rows = rowsOf(meta, group).filter((f) => f.compound.trim());
  if (!rows.length || rows.length < parseArmName(group).length)
    return group;
  return rows.map((f) => [f.compound, f.dose === null ? '' : `${f.dose} ${f.unit}`, f.route, f.schedule].join(' '))
    .join(' + ');
}

/** What a group receives: its details (each empty cell from its name) when they name a compound, a combination's
 *  drug by drug, else its name. `from`: what the details gave that the name did not say, or said otherwise ('dose',
 *  'schedule'...). */
export function armParts(group: string, meta: StudyMeta | null | undefined): ArmPart[] & {from?: string[]} {
  const named = parseArmName(group);
  const text = detailsName(group, meta);
  if (text === group)
    return named;
  const parts: ArmPart[] & {from?: string[]} = parseArmName(text);
  const said = (a: string, b: string | null | undefined) => a.replace(/[.\s]/g, '').toLowerCase() ===
    (b ?? '').replace(/[.\s]/g, '').toLowerCase();
  const rows = rowsOf(meta, group).filter((f) => f.compound.trim());
  parts.from = [...new Set(rows.flatMap((d, i) => {
    // the drug of the name this row is about
    const n = named.length === 1 ? named[0] : named.find((q) => q.key === drugKey(d.compound)) ?? named[i];
    return [['compound', said(d.compound, n?.drug)], ['dose', d.dose === null ||
      (n?.dose?.value === d.dose && said(d.unit, n.dose.unit))], ['route', !d.route || said(d.route, n?.route)],
    ['schedule', !d.schedule || said(d.schedule, n?.schedule)]].filter(([, same]) => !same).map(([k]) => k as string);
  }))];
  return parts;
}

/** One line on what the groups receive, for a list: the compounds, each once. */
export const compoundsOf = (m: StudyMeta | null) =>
  [...new Set((m?.groups ?? []).map((g) => g.compound.trim()).filter((c) => c))].join(', ');

/** What the groups receive, each compound with its dose, route and schedule ('anti-PD-1 10 mg/kg i.p. Q3D×4; …'),
 *  each group once: a combination's drugs (its rows) joined by " + ", the groups by "; " (a " + " between groups
 *  would read as a combination); empty when no group names a compound. */
export function treatmentOf(m: StudyMeta | null): string {
  const byGroup = new Map<string, string[]>();
  for (const g of m?.groups ?? []) {
    if (g.compound.trim()) {
      byGroup.set(g.group, [...byGroup.get(g.group) ?? [], [g.compound.trim(),
        g.dose === null ? '' : `${g.dose} ${g.unit}`.trim(), g.route, g.schedule].filter((s) => s).join(' ')]);
    }
  }
  return [...new Set([...byGroup.values()].map((drugs) => drugs.join(' + ')))].join('; ');
}

/** The study's details to edit: its fields, and its groups as a table (one row per compound: a combination has one
 *  per drug, its group's name on each; typed or pasted into like
 *  any Datagrok grid); `read` gives them as edited. An ID left empty is given when the study is first saved.
 *  `hideGroups` (a blinded study not yet revealed): the groups table is not shown and they are kept as they are.
 *  `taken`: the IDs of the user's other plans and studies, refused for this one; `problem` says what stops the
 *  details from being saved (null: nothing), and `onProblem` is told each time it may have changed. */
export function detailsForm(meta: StudyMeta, groups: string[], hideGroups = false,
  taken: Promise<string[]> | null = null):
  {root: HTMLElement, read: () => StudyMeta, problem: () => string | null, onProblem: (f: () => void) => void} {
  const id = ui.input.string('ID', {value: meta.id, placeholder: 'Given on first save'});
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
  // no dosing at weekends: daily doses skip them, the others move to the nearest weekday; it needs day 0's date
  const weekends = ui.input.bool('No dosing at weekends', {value: meta.weekendsOff, tooltipText:
    'Daily doses skip Saturdays and Sundays; a dose every few days or weekly moves to the nearest weekday'});
  const needsDate = ui.divText('Set day 0 to use it.', 'pmx-details-hint');
  // (saved as rows: blinding hides the groups in PMX Studio only, and the user should know it as they turn it on)
  const blindNote = ui.divText('Anyone the study is shared with can still see the groups in Datagrok’s own table ' +
    'views.', 'pmx-details-hint');
  const showBlindNote = () => blindNote.style.display = blinded.value && domainStorageOn() ? '' : 'none';
  showBlindNote();
  blinded.onChanged.subscribe(showBlindNote);
  // day 0 as typed, when it is a date that exists
  const iso = () => validDay((start.value as unknown as {format(f: string): string} | null)?.format('YYYY-MM-DD'));
  const dated = () => {
    const on = iso() !== null;
    weekends.enabled = on;
    needsDate.style.display = on ? 'none' : '';
  };
  start.onChanged.subscribe(dated);
  dated();
  // what cannot be saved, said under the fields: an ID another study has, a day 0 that does not exist (the date
  // input rolls 2030-02-30 over to 2 March)
  let others = new Set<string>();
  const problemLine = ui.divText('', 'pmx-details-hint');
  problemLine.style.color = 'var(--dg-red, #d9534f)';
  const listeners: (() => void)[] = [];
  const problem = (): string | null => {
    const typed = id.value.trim();
    if (typed && typed !== meta.id && others.has(typed))
      return `${typed} is already the ID of another of your studies or plans: give this one another.`;
    return badDate === null ? null : `Day 0: "${badDate}" is not a date that exists.`;
  };
  // as typed, not as the input shows it once it has rolled the day over: kept until the next thing typed
  let badDate: string | null = null;
  start.root.addEventListener('input', (e) => {
    const text = ((e.target as HTMLInputElement | null)?.value ?? '').trim();
    badDate = impossibleDate(text) ? text : null;
    checked();
  });
  const checked = () => {
    problemLine.textContent = problem() ?? '';
    problemLine.style.display = problemLine.textContent ? '' : 'none';
    listeners.forEach((l) => l());
  };
  void taken?.then((ids) => {
    others = new Set(ids.filter((x) => x));
    checked();
  }, () => {}); // (the list unread: the ID is not checked, as before)
  id.onChanged.subscribe(checked);
  checked();
  // the animals, a row of their own
  const a = meta.animals;
  const animal = {species: ui.input.string('Species', {value: a.species}),
    // (examples, read as such: no strain named, since one example would contradict half the models' hosts)
    strain: ui.input.string('Strain', {value: a.strain, placeholder: 'the tumor model\'s host strain'}),
    sex: ui.input.choice('Sex', {value: a.sex || null, items: ['female', 'male', 'both'], nullable: true}),
    age: ui.input.string('Age', {value: a.age, placeholder: 'e.g. 8 weeks', tooltipText: 'At the start of the study'}),
    weight: ui.input.string('Body weight', {value: a.weight, placeholder: 'e.g. 20 g',
      tooltipText: 'At the start of the study'}),
    source: ui.input.string('Source', {value: a.source, placeholder: 'e.g. a supplier, or in house',
      tooltipText: 'Supplier, or bred in house'})};
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
  grid.root.style.width = '100%';
  grid.root.style.height = `${Math.min(8, rows.length + 1) * 30 + 8}px`;
  // the columns share the grid's width, each in proportion to its own (no blank space after the last)
  ui.onSizeChanged(grid.root).subscribe(() => {
    const cols = Array.from({length: grid.columns.length - 1}, (_, i) => grid.columns.byIndex(i + 1)!);
    const room = grid.root.clientWidth - (grid.columns.rowHeader?.width ?? 0) - 18;
    const used = cols.reduce((a, c) => a + c.width, 0);
    if (room > 0 && used > 0 && Math.abs(room - used) > 2)
      cols.forEach((c) => c.width = Math.floor(c.width * room / used));
  });
  const fields = ui.inputs([id, description, owner, project, tags, start, blinded, weekends], 'pmx-details-fields');
  const animals = [ui.divText('Animals', 'pmx-details-head'), ui.inputs(Object.values(animal), 'pmx-details-fields')];
  const root = ui.divV([fields, problemLine, blindNote, needsDate, ...animals, ...hideGroups ?
    [ui.divText('What each group receives is hidden until the groups are revealed.')] :
    [ui.divText('What each group receives (the schedule as written in the protocol, e.g. Q3D×4 or QD×14):'),
      grid.root]], 'pmx-details');
  const read = (): StudyMeta => {
    const day = iso();
    const cell = (c: string, i: number) => table.col(c)!.isNone(i) ? '' : String(table.get(c, i)).trim();
    return {id: id.value.trim() || meta.id, description: description.value.trim(), owner: owner.value.trim(),
      project: project.value.trim(), tags: tags.value.split(',').map((t) => t.trim()).filter((t) => t),
      start: day, groups: hideGroups ? meta.groups : rows.map((g, i) => ({group: g.group,
        compound: cell('compound', i), dose: table.col('dose')!.isNone(i) ? null : f32(table.get('dose', i)),
        unit: cell('unit', i), route: cell('route', i), schedule: cell('schedule', i)})), notes: meta.notes,
      blinded: blinded.value, weekendsOff: weekends.value, exceptions: meta.exceptions, given: meta.given,
      animals: Object.fromEntries(ANIMAL_FIELDS.map((k) => [k, (animal[k].value ?? '').trim()])) as unknown as
        AnimalDetails, colors: meta.colors, ...meta.slots ? {slots: meta.slots} : {},
      ...meta.copyOf ? {copyOf: meta.copyOf} : {}};
  };
  return {root, read, problem, onProblem: (f) => listeners.push(f)};
}

/** Whether a day typed as y-m-d, d/m/y or m/d/y names a day no calendar has (2030-02-30; 31/04/2026 read either
 *  way): the date input would roll it over to another day. Anything else is the input's to read. */
export function impossibleDate(text: string): boolean {
  const n = text.trim().match(/^(\d{1,4})[-/.](\d{1,2})[-/.](\d{1,4})$/);
  if (!n)
    return false;
  const iso = (y: string, m: string, d: string) => validDay(`${y.padStart(4, '0')}-${m.padStart(2, '0')}-` +
    `${d.padStart(2, '0')}`);
  const [a, b, c] = [n[1], n[2], n[3]];
  return a.length === 4 ? iso(a, b, c) === null : c.length === 4 ? iso(c, b, a) === null && iso(c, a, b) === null :
    false;
}

/** What Plan reads in a group's name and details, in a few words: 'Cpd-A · 30 mg/kg · QD · p.o.', a combination's
 *  drugs joined by ' + ', 'control' for a control, what came from the details said so; '' when it reads nothing. */
export function parseLine(group: string, meta: StudyMeta | null = null): string {
  const parts = armParts(group, meta);
  if (!parts.length)
    return controlByName([group]) ? 'control' : '';
  return describeArm(parts) + (parts.from?.length ? ` (${parts.from.join(', ')} from the group details)` : '');
}

/** What arm `arm` adds up to over days 0 to `last`, drug by drug: its doses under the study's rules, the dose in all
 *  and the compound to prepare for `animals` of `weightKg` (schedule.ts doseTotal; null: 25 g); with `givenTo`, the
 *  doses up to that day against those. Display only. null: no drug of the arm has a schedule that gives doses. */
export function armDoses(arm: string, meta: StudyMeta | null, last: number, animals: number | null,
  weightKg: number | null, givenTo?: number, withAmount = true): {text: string, tip: string} | null {
  const parts = armParts(arm, meta);
  const dosed: number[] = [];
  // the schedule's own days, before the exceptions (a skip counts only where it took a dose away)
  const planned: number[] = [];
  const lines = parts.flatMap((p) => {
    const s = readSchedule(p.schedule);
    if (!s || s.kind === 'continuous')
      return [];
    // (a stopped arm's doses end on its review's day)
    const until = (d: number) => d <= (meta?.stopped?.[arm] ?? Infinity);
    const base = dosedDays(s, last, meta).filter(until);
    const days = excepted(base, arm, meta, last).filter(until);
    planned.push(...base);
    dosed.push(...days);
    const all = doseTotal(days.length * s.perDay, p.dose, animals, weightKg);
    const given = givenTo === undefined ? undefined :
      doseTotal(days.filter((d) => d <= givenTo).length * s.perDay, p.dose, animals, weightKg);
    const l = doseLine(all, animals, given, withAmount);
    return [{text: parts.length > 1 ? `${p.drug}: ${l.text}` : l.text, tip: l.tip}];
  });
  if (!lines.length)
    return null;
  // the exceptions and, in Monitor once the study logs doses, the dosing days so far logged as given
  // counted from the days they changed (one on a day the schedule does not dose, or past the end, changes nothing)
  const own = (meta?.exceptions ?? []).filter((e) => e.group === arm);
  const skipped = own.filter((e) => e.change === 'skip' && planned.includes(e.day)).length;
  const added = own.filter((e) => e.change === 'add' && dosed.includes(e.day) && !planned.includes(e.day)).length;
  const changes = [skipped ? `${skipped} skipped` : '', added ? `${added} added` : ''].filter((x) => x).join(', ');
  let logged = '';
  if (givenTo !== undefined && meta?.given.length) {
    const due = new Set(dosed.filter((d) => d <= givenTo));
    const done = new Set(meta.given.filter((g) => g.group === arm && due.has(g.day)).map((g) => g.day));
    logged = due.size ? `${done.size} of ${due.size} dosing days so far logged as given` : '';
  }
  return {text: [lines.map((l) => l.text).join('; '), changes && `${changes} by exception`, logged]
    .filter((x) => x).join(' · '), tip: lines[0].tip};
}

/** The dose per administration that gives `total` (mg/kg in all) over arm `arm`'s doses to day `last`, under the
 *  study's rules and exceptions: Plan's one-time help to match a total dose across schedules. null: no doses. */
export function doseForTotal(arm: string, meta: StudyMeta | null, last: number, total: number): number | null {
  const n = armParts(arm, meta).map((p) => {
    const s = readSchedule(p.schedule);
    return s && s.kind !== 'continuous' ? excepted(dosedDays(s, last, meta), arm, meta, last).length * s.perDay : 0;
  });
  return n.length === 1 && n[0] > 0 && total > 0 ? total / n[0] : null;
}
