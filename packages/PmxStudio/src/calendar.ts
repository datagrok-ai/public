/* The study's calendar, shared by Plan and Monitor: a row a week, Monday to Sunday, from day 0 to the last day, each
   day a cell with its date, its marks and a tooltip; a day's menu (right click, the Menu key or Shift+F10 on a
   focused day) and a double click to note what happened. What a day shows is the stage's (`day`); the study's own
   (its date, its doses, its notes) the calendar adds. Built from u2's primitives (u2 has no calendar); call it in a
   scope, which owns the tooltips. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {div, divV, span, Tooltip} from '@datagrok-libraries/u2';

import {WEEKDAY_NAMES} from './frames';
import {DayNote, StudyMeta, dateOf, dateText, dosingDays, noteDialog, startFrom} from './study-meta';

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
  onNote: (day: number) => void}): HTMLElement {
  const {meta, startWeekday: start, last} = o;
  const doses = dosingDays(meta.groups, last);
  const notesOn = new Map<number, DayNote[]>();
  for (const n of meta.notes)
    notesOn.set(n.day, [...notesOn.get(n.day) ?? [], n]);
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
      const study = [...doses.has(day) ? [mark('Dose', CAL.green, false)] : [],
        ...notes.length > 0 ? [mark(notes.length === 1 ? '1 note' : `${notes.length} notes`, CAL.amber, false)] : []];
      const look = o.day(day, study);
      cell.append(num, ...look.marks);
      if (look.tint)
        cell.style.background = look.tint;
      if (look.border)
        cell.style.border = look.border;
      Tooltip.bind(cell, [`Day ${day}${on ? ` (${dateText(on, true)})` : ''}${look.detail ? `: ${look.detail}` : ''}`,
        ...doses.has(day) ? [`Dose: ${doses.get(day)!.join(', ')}`] : [],
        ...notes.map((n) => `${n.kind}${n.animal ? ` · ${n.animal}` : ''}: ${n.text} (${n.who})`)].join('\n'));
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

/** What a day offers on the study's details: add or remove a note, set or change its date (day 0's, from it), copy
 *  the date. `animals`: those a note can be about. */
export function dayActions(get: () => StudyMeta, set: (m: StudyMeta) => void, animals: () => string[] = () => []):
  {addNote: (day: number) => void, menu: (day: number, at: MenuAt) => void} {
  const label = (day: number) => {
    const d = dateOf(get().start, day);
    return `day ${day}${d ? ` (${dateText(d)})` : ''}`;
  };
  const addNote = (day: number) => noteDialog(day, label(day), animals(),
    (n) => set({...get(), notes: [...get().notes, n]}));
  const dateDay = (day: number) => {
    const input = ui.input.date(`Date of day ${day}`);
    ui.dialog(`Date ${label(day)}`).add(input).onOK(() => {
      const iso = (input.value as unknown as {format(f: string): string} | null)?.format('YYYY-MM-DD');
      if (iso)
        set({...get(), start: startFrom(iso, day)});
    }).show();
  };
  const menu = (day: number, at: MenuAt) => {
    const m = get();
    const popup = DG.Menu.popup().item('Add a note…', () => addNote(day));
    const own = m.notes.filter((n) => n.day === day);
    if (own.length > 0) {
      const remove = popup.group('Remove a note');
      for (const n of own) {
        remove.item(`${n.kind}: ${n.text.length > 40 ? `${n.text.slice(0, 40)}…` : n.text}`,
          () => set({...get(), notes: get().notes.filter((x) => x !== n)}));
      }
    }
    popup.item(m.start ? `Change the date of day ${day}…` : `Set the date of day ${day}…`, () => dateDay(day));
    const d = dateOf(m.start, day);
    if (d)
      popup.item('Copy the date', () => void navigator.clipboard.writeText(d.toISOString().slice(0, 10)));
    popup.show(at);
  };
  return {addNote, menu};
}

/** The study's days as an iCalendar file of all-day events (`events`: day and what is done, with the doses its
 *  groups' schedules give), for Google, Outlook or Apple Calendar to import. Needs day 0's date. */
export function icsOf(name: string, meta: StudyMeta, last: number, events: [number, string][]): string {
  if (!meta.start)
    throw new Error('Set the date of day 0 first: right-click a day in the calendar');
  const ymd = (d: Date) => d.toISOString().slice(0, 10).replace(/-/g, '');
  const text = (s: string) => s.replace(/[\\;,]/g, (c) => `\\${c}`).replace(/\n/g, '\\n');
  const stamp = `${new Date().toISOString().replace(/[-:]/g, '').slice(0, 15)}Z`;
  const all = [...events, ...[...dosingDays(meta.groups, last)].map(([d, g]) => [d, `Dose: ${g.join(', ')}`])]
    .sort((a, b) => (a[0] as number) - (b[0] as number)) as [number, string][];
  const lines = ['BEGIN:VCALENDAR', 'VERSION:2.0', 'PRODID:-//Datagrok//PMX Studio//EN', 'CALSCALE:GREGORIAN'];
  // an event keeps its UID from one export to the next (its day and what it is), so importing a later export
  // updates the calendar's events and does not double them
  all.forEach(([day, what]) => lines.push('BEGIN:VEVENT',
    `UID:${meta.id || name}-${day}-${what.replace(/[^a-z0-9]+/gi, '-')}@pmxstudio`,
    `DTSTAMP:${stamp}`, `DTSTART;VALUE=DATE:${ymd(dateOf(meta.start, day)!)}`,
    `DTEND;VALUE=DATE:${ymd(dateOf(meta.start, day + 1)!)}`, `SUMMARY:${text(`${name}: ${what} (day ${day})`)}`,
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
