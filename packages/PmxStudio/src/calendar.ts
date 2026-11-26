/* The study's calendar, shared by Plan and Monitor: a row a week, Monday to Sunday, from day 0 to the last day, each
   day a cell with its date, its marks and a tooltip; a day's menu (right click, the Menu key or Shift+F10 on a
   focused day) and a double click to note what happened. What a day shows is the stage's (`day`); the study's own
   (its date, its doses, its notes) the calendar adds. Built from u2's primitives (u2 has no calendar); call it in a
   scope, which owns the tooltips. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {div, divV, span, Tooltip} from '@datagrok-libraries/u2';

import {groupColors} from './colors';
import {WEEKDAY_NAMES} from './frames';
import {okNamed, showCentred} from './ui';
import {StudyMeta, dateOf, dateText, dosingDays, noteDialog, schedulesOf, startFrom, withException, withGiven}
  from './study-meta';

/** Where a day's menu opens (DG.Menu.show's options): at the mouse, or next to the focused day. */
export type MenuAt = {causedBy?: MouseEvent, element?: HTMLElement, nextToElement?: boolean};

export const CAL = {blue: '#2083d5', slate: '#5b6b7f', red: '#c0504d', dark: '#4d5261', green: '#2e8b57',
  amber: '#b7791f'};
export const WEEKDAYS = WEEKDAY_NAMES.map((d) => d.slice(0, 3));

/** A day's mark: a short label in its kind's colour. */
export function mark(text: string, color: string, filled: boolean): HTMLElement {
  const el = span(text);
  Object.assign(el.style, {display: 'inline-block', padding: '0 5px', borderRadius: '3px', fontSize: '10.5px',
    lineHeight: '15px', whiteSpace: 'nowrap', fontWeight: '500', color: filled ? '#fff' : color,
    background: filled ? color : `color-mix(in srgb, ${color} 12%, transparent)`});
  return el;
}

/** One small dot per arm dosed that day, in the arm's colour; its tooltip and aria text name the arm. */
export function doseDots(arms: string[], colors: {[arm: string]: string}): HTMLElement {
  const row = div(arms.map((a) => {
    const dot = span('', 'pmx-calendar-dot');
    Object.assign(dot.style, {display: 'inline-block', width: '8px', height: '8px', borderRadius: '50%',
      background: colors[a] ?? CAL.green, flex: 'none'});
    dot.title = a;
    return dot;
  }), 'pmx-calendar-dots');
  Object.assign(row.style, {display: 'flex', flexWrap: 'wrap', gap: '3px', alignItems: 'center', minHeight: '15px'});
  row.setAttribute('role', 'img');
  row.setAttribute('aria-label', `Dose: ${arms.join(', ')}`);
  return row;
}

/** A summary line: a mark, then what it means. */
export function marked(m: HTMLElement, text: string): HTMLElement {
  const line = div([m, span(text)]);
  Object.assign(line.style, {display: 'flex', alignItems: 'baseline', gap: '6px'});
  return line;
}

/** What a day shows: its marks (`study`: the study's doses and notes, placed by the stage), what it is (the
 *  tooltip's first line), a background and a border. */
export interface DayLook {marks: HTMLElement[], detail?: string, tint?: string, border?: string}

/** The weeks from day 0 (on weekday `startWeekday`, Monday 0) to `last`; `header`: the seven weekday labels. */
export function weeksGrid(o: {meta: StudyMeta, startWeekday: number, last: number, header?: HTMLElement[],
  day: (day: number, study: HTMLElement[]) => DayLook, onMenu: (day: number, at: MenuAt) => void,
  onNote: (day: number) => void, control?: string}): HTMLElement {
  const {meta, startWeekday: start, last} = o;
  const doses = dosingDays(meta.groups, last, meta);
  // the arms' colours as the charts have them (read here: a colour change draws the calendar again)
  const names = [...new Set(meta.groups.map((g) => g.group))];
  const colors = groupColors(names, o.control ?? 'Control', meta.colors, meta.slots);
  // an arm whose schedule starts after day 0: its first dose day
  const starts = new Map<number, string[]>();
  for (const g of names) {
    const days = [...doses].filter(([, gs]) => gs.includes(g)).map(([d]) => d);
    const late = meta.groups.filter((x) => x.group === g).flatMap(schedulesOf).some((s) => s.first > 0);
    if (days.length && late && Math.min(...days) > 0)
      starts.set(Math.min(...days), [...starts.get(Math.min(...days)) ?? [], g]);
  }
  const by = <T extends {day: number}>(xs: T[]) => {
    const on = new Map<number, T[]>();
    for (const x of xs)
      on.set(x.day, [...on.get(x.day) ?? [], x]);
    return on;
  };
  const notesOn = by(meta.notes);
  // those of a group no longer in the study are not shown
  const known = (x: {group: string}) => meta.groups.some((g) => g.group === x.group);
  const changesOn = by((meta.exceptions ?? []).filter(known));
  const givenOn = by((meta.given ?? []).filter(known));
  const grid = div([], 'pmx-calendar');
  Object.assign(grid.style, {display: 'grid', gridTemplateColumns: '56px repeat(7, minmax(0, 1fr))', gap: '4px'});
  grid.append(span(''), ...o.header ?? WEEKDAYS.map((d) => {
    const h = span(d, 'pmx-calendar-weekday');
    Object.assign(h.style, {fontSize: '11px', textAlign: 'center', padding: '2px 0',
      color: 'var(--dg-text-color-muted, #68727f)'});
    return h;
  }));
  const weeks = Math.floor((start + last) / 7) + 1;
  for (let w = 0; w < weeks; w++) {
    const label = span(`Week ${w + 1}`);
    Object.assign(label.style, {fontSize: '11px', color: 'var(--dg-text-color-muted, #68727f)', alignSelf: 'center'});
    grid.append(label);
    for (let d = 0; d < 7; d++) {
      const day = 7 * w + d - start;
      // as the demo draws a day: a soft chip (the weekend's clear), the day bold and its date muted, small marks
      const cell = divV([], 'pmx-calendar-day');
      Object.assign(cell.style, {minHeight: '40px', padding: '3px 5px', gap: '2px', borderRadius: '4px',
        border: `1px solid ${d >= 5 ? 'transparent' : 'var(--dg-grey-2, #e4e7ec)'}`, boxSizing: 'border-box',
        alignItems: 'flex-start', background: d >= 5 ? 'transparent' : 'var(--dg-bg, #fff)'});
      if (day < 0 || day > last) {
        cell.style.background = 'transparent';
        grid.append(cell);
        continue;
      }
      const on = dateOf(meta.start, day);
      const num = span(`Day ${day}`);
      Object.assign(num.style, {fontSize: '11px', fontWeight: '500', color: 'var(--dg-text-color, #4d5261)'});
      if (on) {
        const date = span(` · ${dateText(on)}`);
        Object.assign(date.style, {fontWeight: '400', color: 'var(--dg-text-color-muted, #68727f)'});
        num.append(date);
      }
      const notes = notesOn.get(day) ?? [];
      const changes = changesOn.get(day) ?? [];
      const skipped = changes.filter((c) => c.change === 'skip');
      const added = new Set(changes.filter((c) => c.change === 'add').map((c) => c.group));
      const dosed = doses.get(day) ?? [];
      const given = new Set((givenOn.get(day) ?? []).map((g) => g.group));
      // a day whose doses are all logged as given is filled; a skipped dose is marked in red
      const allGiven = dosed.length > 0 && dosed.every((g) => given.has(g));
      const starting = starts.get(day) ?? [];
      const study = [...dosed.length ? [doseDots(dosed, colors)] : [],
        ...allGiven ? [mark('Given', CAL.green, true)] : [],
        ...starting.length ? [mark('Start', CAL.green, false)] : [],
        ...skipped.length ? [mark('Skipped', CAL.red, false)] : [],
        ...notes.length > 0 ? [mark(notes.length === 1 ? '1 note' : `${notes.length} notes`, CAL.amber, false)] : []];
      const why = (c: {note: string, who: string}) => `${c.note ? ` (${c.note})` : ''}${c.who ? `, ${c.who}` : ''}`;
      const look = o.day(day, study);
      cell.append(num, ...look.marks);
      if (look.tint)
        cell.style.background = look.tint;
      if (look.border)
        cell.style.border = look.border;
      const tip = [`Day ${day}${on ? ` (${dateText(on, true)})` : ''}${look.detail ? `: ${look.detail}` : ''}`,
        ...dosed.length ? [`Dose: ${dosed.map((g) => g + (added.has(g) ? ' (added)' : '')).join(', ')}`] : [],
        ...starting.length ? [`Start: ${starting.map((g) => `${g} (from day ${day})`).join(', ')}`] : [],
        ...given.size ? [`Given: ${[...given].join(', ')}`] : [],
        ...skipped.map((c) => `Skipped: ${c.group}${why(c)}`),
        ...notes.map((n) => `${n.kind}${n.animal ? ` · ${n.animal}` : ''}: ${n.text} (${n.who})`)].join('\n');
      Tooltip.bind(cell, tip);
      cell.setAttribute('aria-label', tip.replace(/\n/g, '; '));
      cell.tabIndex = 0;
      cell.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        o.onMenu(day, {causedBy: e});
      });
      cell.addEventListener('keydown', (e) => {
        if (e.key === 'ContextMenu' || (e.key === 'F10' && e.shiftKey)) {
          e.preventDefault();
          o.onMenu(day, {element: cell, nextToElement: true});
        }
      });
      cell.addEventListener('dblclick', () => o.onNote(day));
      grid.append(cell);
    }
  }
  return grid;
}

/** What a day offers on the study's details: add or remove a note, skip or add a group's dose (an exception: the
 *  schedule as written stays), set or change its date (day 0's, from it), copy the date. `animals`: those a note can
 *  be about; `logGiven` (Monitor): the day's doses can be logged as given; `name`: a group as the user may see it
 *  (blinded: "Group 2"). */
export function dayActions(get: () => StudyMeta, set: (m: StudyMeta) => void, animals: () => string[] = () => [],
  opts: {logGiven?: boolean, name?: (group: string) => string, last?: () => number} = {}):
  {addNote: (day: number) => void, menu: (day: number, at: MenuAt) => void} {
  const name = opts.name ?? ((g: string) => g);
  // a dialog to pick some of `groups` (by the names the user sees) and why; `ok` gets the groups and the note
  // one group to choose from is ticked already; the button says the action and waits for a group (a dose change;
  // logging none given is an answer)
  const pick = (title: string, verb: string, groups: string[], ticked: string[],
    ok: (gs: string[], note: string) => void, withNote = true) => {
    const items = groups.map(name);
    const chosen = ui.input.multiChoice<string>('Groups', {value: (groups.length === 1 ? groups : ticked).map(name),
      items});
    const note = ui.input.string('Why', {value: '', placeholder: 'e.g. animal unwell, holiday',
      tooltipText: 'Optional: kept with the change'});
    const d = ui.dialog(title).add(ui.inputs(withNote ? [chosen, note] : [chosen])).onOK(() =>
      ok(groups.filter((_, i) => (chosen.value ?? []).includes(items[i])), note.value.trim()));
    const button = okNamed(d, verb);
    const sync = () => button.disabled = withNote && (chosen.value ?? []).length === 0;
    chosen.onChanged.subscribe(sync);
    sync();
    showCentred(d);
  };
  const doseItems = (popup: DG.Menu, day: number) => {
    const m = get();
    const scheduled = [...new Set(m.groups.filter((g) => schedulesOf(g).length).map((g) => g.group))];
    // the day's doses as the calendar places them (to its last day: a weekend dose moves within it)
    const dosed = dosingDays(m.groups, Math.max(day, opts.last?.() ?? day), m).get(day) ?? [];
    const changed = m.exceptions.filter((e) => e.day === day && m.groups.some((g) => g.group === e.group));
    if (dosed.length) {
      popup.item('Skip a dose…', () => pick(`Skip a dose on ${label(day)}`, 'Skip', dosed, [], (gs, note) =>
        set(gs.reduce((x, g) => withException(x, day, g, 'skip', note), get()))));
    }
    const undosed = scheduled.filter((g) => !dosed.includes(g) && !changed.some((e) => e.group === g));
    if (undosed.length) {
      popup.item('Add a dose…', () => pick(`Add a dose on ${label(day)}`, 'Add', undosed, [], (gs, note) =>
        set(gs.reduce((x, g) => withException(x, day, g, 'add', note), get()))));
    }
    // items of their own, not a submenu: a click on a submenu's parent closed the whole menu
    for (const e of changed) {
      popup.item(`Undo: ${name(e.group)} ${e.change === 'skip' ? 'skipped' : 'added'}`,
        () => set(withException(get(), day, e.group, null)));
    }
    if (opts.logGiven && dosed.length) {
      const logged = m.given.filter((g) => g.day === day).map((g) => g.group);
      popup.item('Log doses given…', () => pick(`Doses given on ${label(day)}`, 'Log', dosed,
        logged.length ? logged : dosed, (gs) => set(withGiven(get(), day, gs)), false));
    }
  };
  const label = (day: number) => {
    const d = dateOf(get().start, day);
    return `day ${day}${d ? ` (${dateText(d)})` : ''}`;
  };
  const addNote = (day: number) => noteDialog(day, label(day), animals(),
    (n) => set({...get(), notes: [...get().notes, n]}));
  const dateDay = (day: number) => {
    // the date the study gives it already, to change; the others follow from it
    const now = dateOf(get().start, day);
    const input = ui.input.date('Date', {tooltipText: 'Every other day of the study follows from it'});
    if (now)
      input.stringValue = now.toISOString().slice(0, 10);
    const d = ui.dialog(`Date of day ${day}`).add(input).onOK(() => {
      const iso = (input.value as unknown as {format(f: string): string} | null)?.format('YYYY-MM-DD');
      if (iso)
        set({...get(), start: startFrom(iso, day)});
    });
    okNamed(d, 'Set date');
    showCentred(d);
  };
  const menu = (day: number, at: MenuAt) => {
    const m = get();
    const popup = DG.Menu.popup().item('Add a note…', () => addNote(day));
    const own = m.notes.filter((n) => n.day === day);
    for (const n of own) {
      popup.item(`Remove note: ${n.kind}: ${n.text.length > 40 ? `${n.text.slice(0, 40)}…` : n.text}`,
        () => set({...get(), notes: get().notes.filter((x) => x !== n)}));
    }
    popup.separator();
    doseItems(popup, day);
    popup.separator();
    popup.item(m.start ? `Change the date of day ${day}…` : `Set the date of day ${day}…`, () => dateDay(day));
    const d = dateOf(m.start, day);
    if (d)
      popup.item('Copy the date', () => void navigator.clipboard.writeText(d.toISOString().slice(0, 10)));
    popup.show(at);
  };
  return {addNote, menu};
}

/** The study's days as an iCalendar file of all-day events (`events`: day and what is done, with the doses its
 *  groups' schedules give), for Google, Outlook or Apple Calendar to import. Needs day 0's date. `end`: the expected
 *  end; a day past it is in the file (the end is a forecast, the study may run on) but its events say so, so no one
 *  books a dose or a measurement past the end as if it were certain. */
export function icsOf(name: string, meta: StudyMeta, last: number, events: [number, string][], end?: number): string {
  if (!meta.start)
    throw new Error('Set the date of day 0 first: right-click a day in the calendar');
  const ymd = (d: Date) => d.toISOString().slice(0, 10).replace(/-/g, '');
  const text = (s: string) => s.replace(/[\\;,]/g, (c) => `\\${c}`).replace(/\n/g, '\\n');
  const stamp = `${new Date().toISOString().replace(/[-:]/g, '').slice(0, 15)}Z`;
  const all = [...events, ...[...dosingDays(meta.groups, last, meta)].map(([d, g]) => [d, `Dose: ${g.join(', ')}`])]
    .sort((a, b) => (a[0] as number) - (b[0] as number)) as [number, string][];
  const lines = ['BEGIN:VCALENDAR', 'VERSION:2.0', 'PRODID:-//Datagrok//PMX Studio//EN', 'CALSCALE:GREGORIAN'];
  // an event keeps its UID from one export to the next (its day and what it is), so importing a later export
  // updates the calendar's events and does not double them (the UID leaves out the "if it runs on" a moving end adds)
  all.forEach(([day, what]) => lines.push('BEGIN:VEVENT',
    `UID:${meta.id || name}-${day}-${what.replace(/[^a-z0-9]+/gi, '-')}@pmxstudio`,
    `DTSTAMP:${stamp}`, `DTSTART;VALUE=DATE:${ymd(dateOf(meta.start, day)!)}`,
    `DTEND;VALUE=DATE:${ymd(dateOf(meta.start, day + 1)!)}`, `SUMMARY:${text(`${name}: ${what} (day ${day})` +
      (end !== undefined && day > end ? `, only if the study runs past its expected end (day ${end})` : ''))}`,
    'END:VEVENT'));
  lines.push('END:VCALENDAR');
  return lines.join('\r\n');
}

/** The study day of `now` (local calendar date) from day 0's date; null: day 0 has no date. */
export function dayToday(start: string | null, now = new Date()): number | null {
  if (!start)
    return null;
  const today = Date.UTC(now.getFullYear(), now.getMonth(), now.getDate());
  return Math.round((today - Date.parse(`${start}T00:00:00Z`)) / 86400000);
}

/** Day 0's date against the table's last measured day: a sentence when the table has measurements past the study
 *  day today falls on (the date is not this table's, so no stage prints it, or a day dated from it, as fact); null
 *  when it can stand (no date, no measurements, or none past today). */
export function startMismatch(start: string | null, lastDay: number | null | undefined,
  now = new Date()): string | null {
  const today = dayToday(start, now);
  return today === null || lastDay === null || lastDay === undefined || today >= lastDay ? null :
    `Day 0 is dated ${dateText(dateOf(start, 0)!, true)}, but the table has measurements up to day ${lastDay}: ` +
    'the date does not match.';
}
