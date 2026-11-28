// "Enter today's values": a measuring day typed in a dialog and written to the study table only on Save day, all of
// it at once, so a half-entered day never exists in the study. A draft (per user, study and day) keeps what was typed.
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {dayZeroAsk, mm3Text, ReadingCheck} from './monitor';
import {checkTyped} from './reading-check';
import {ReadingColumns, ellipsoid} from './study-table';
import {pastedNumber, plural} from './wide-table';
import {addLog} from './study-file';
import {DayNote, metaOf, setMetaOn} from './study-meta';
import {CAL} from './calendar';
import {dialogLine, modal, modalOf, okNamed, showCentred} from './ui';

/** An animal on the day's list: its last values before the day, and the range its volume should fall in (null: no
 *  range for this day). `gone`: why it is no longer on study (listed greyed, not counted, nothing typed); `note`: a
 *  word beside it (it missed the last session: "Last measured day 9"). */
export interface EntryAnimal {animal: string, group: string, lastDay: number | null, lastVolume: number | null,
  lastWeight: number | null, lower: number | null, upper: number | null, expected: string, gone?: string,
  note?: string}

/** How far past the study's last day (or its planned end) a day is asked about before it is entered. ponytail: a fixed
 *  two weeks, the calendar's own margin past a study's span. */
export const FAR_DAYS = 14;

/** A day as typed in the dialog: a whole study day, or why it is none (`error`), or a question first (`ask`: a day
 *  far past the study's last day measured and its planned end). */
export function readDay(text: string, o: {last: number | null, end: number | null, today?: number | null}):
  {day?: number, error?: string, ask?: string} {
  const s = text.trim();
  if (!s)
    return {error: 'Type the study day (0 = randomization).'};
  if (/^[+-]?\d+([.,]\d+)?$/.test(s) && !/^[+-]?\d+$/.test(s))
    return {error: `Whole days only: "${s}" is not a study day (0 is randomization, then 1, 2, 3…).`};
  if (!/^[+-]?\d+$/.test(s))
    return {error: `"${s}" is not a day: type a whole number (0 is randomization).`};
  const day = Number(s);
  if (day < 0)
    return {error: 'A study day is 0 or later: day 0 is randomization.'};
  const from = Math.max(o.last ?? 0, o.end ?? 0);
  if (day > from + FAR_DAYS) {
    const known = [...o.last !== null ? [`the last day measured is day ${o.last}`] : [],
      ...o.end !== null ? [`the planned end is day ${o.end}`] : []];
    return {day, ask: `Day ${day} is far past the study${known.length ? ` (${known.join(', ')})` : ''}. ` +
      `Enter values for day ${day}?`};
  }
  // a day still to come (`today`: the study day of today's date) asks first too (QA 0.1.3 round 4 chaos N9)
  if (o.today != null && day > o.today) {
    return {day, ask: `Day ${day} is ${day === o.today + 1 ? 'tomorrow' : 'still to come'} (today is day ` +
      `${o.today}). Enter values for day ${day}?`};
  }
  return {day};
}

/** What was typed for an animal, as typed; `missed`: why it was not measured. */
export interface Typed {volume?: string, length?: string, width?: string, weight?: string, missed?: string}
export type Values = {[animal: string]: Typed};

/** Why an animal was not measured. */
export const MISSED = ['Dead', 'Removed', 'Missed'];

/** The columns the day is written to. */
export interface DayColumns {animal: DG.Column, group: DG.Column, day: DG.Column, readings: ReadingColumns}

/** A typed cell as a number: null when empty, NaN when it is not one. */
export const typedNumber = (s: string | undefined): number | null => s === undefined ? null : pastedNumber(s);

/** The volume a row gives: typed, or from its caliper readings (null: not complete). */
export function stagedVolume(t: Typed, caliper: boolean): number | null {
  if (!caliper)
    return typedNumber(t.volume);
  const [l, w] = [typedNumber(t.length), typedNumber(t.width)];
  return l === null || w === null ? null : ellipsoid(l, w);
}

/** A row's checks, Monitor's own (checkReading): its volume and its body weight. */
export function rowChecks(a: EntryAnimal, t: Typed, caliper: boolean):
  {volume: ReadingCheck | null, weight: ReadingCheck | null} {
  if (t.missed || a.gone)
    return {volume: null, weight: null};
  const v = stagedVolume(t, caliper);
  const w = typedNumber(t.weight);
  return {volume: checkTyped(v, {before: a.lastVolume, lower: a.lower, upper: a.upper}),
    weight: checkTyped(w, {weight: true, before: a.lastWeight})};
}

/** The footer's counts: animals entered (a volume, or not measured with a reason), values to check, impossible. */
export function entrySummary(animals: EntryAnimal[], values: Values, caliper: boolean):
  {entered: number, toCheck: number, impossible: number} {
  let [entered, toCheck, impossible] = [0, 0, 0];
  // (an animal no longer on study is not counted)
  for (const a of animals.filter((x) => !x.gone)) {
    const t = values[a.animal] ?? {};
    const c = rowChecks(a, t, caliper);
    const kinds = [c.volume?.kind, c.weight?.kind];
    entered += t.missed || c.volume ? 1 : 0;
    toCheck += kinds.some((k) => k === 'again' || k === 'typo') ? 1 : 0;
    impossible += kinds.includes('impossible') ? 1 : 0;
  }
  return {entered, toCheck, impossible};
}

/** The footer's line: "48 of 60 entered · 2 to check". */
export function summaryText(s: {entered: number, toCheck: number, impossible: number}, n: number): string {
  return [`${s.entered} of ${n} entered`, ...s.toCheck ? [`${s.toCheck} to check`] : [],
    ...s.impossible ? [`${s.impossible} impossible`] : []].join(' · ');
}

/** The rows of `day` (by animal, the first one each), as entered before; on day 0, an animal's row with no day and
 *  no volume yet too (the animals listed before randomization take their day 0 there). */
function rowsOfDay(df: DG.DataFrame, c: DayColumns, day: number): Map<string, number> {
  const out = new Map<string, number>();
  for (let i = 0; i < df.rowCount; i++) {
    const a = String(c.animal.get(i));
    const on = c.day.isNone(i) ? day === 0 && c.readings.volume.isNone(i) : c.day.get(i) === day;
    if (on && !out.has(a))
      out.set(a, i);
  }
  return out;
}

/** The day's values already in the table, as typed text: editing a day starts from them. */
export function valuesOfDay(df: DG.DataFrame, c: DayColumns, day: number): Values {
  const text = (col: DG.Column | null | undefined, i: number) => !col || col.isNone(i) ? undefined :
    String(+Number(col.get(i)).toFixed(3));
  const out: Values = {};
  const r = c.readings;
  for (const [a, i] of rowsOfDay(df, c, day)) {
    const t: Typed = r.caliper ? {length: text(r.caliper.length, i), width: text(r.caliper.width, i)} :
      {volume: text(r.volume, i)};
    if (r.weight)
      t.weight = text(r.weight, i);
    out[a] = Object.fromEntries(Object.entries(t).filter(([, v]) => v !== undefined));
  }
  return out;
}

/** A row's readings as the log keeps them: "a1 150 mm³, 20.5 g". */
function readingText(c: DayColumns, i: number): string {
  const r = c.readings;
  const x = (col: DG.Column | null, unit: string) => !col || col.isNone(i) ? [] :
    [`${+Number(col.get(i)).toFixed(2)} ${unit}`];
  return `${c.animal.get(i)} ${[...x(r.volume, 'mm³'), ...x(r.weight, 'g')].join(', ')}`.trim();
}

/** Takes rows off the day: on day 0, an animal's row stays as it was before its day 0 was entered (no day, no
 *  readings: the animal is still listed); on any other day the row goes. */
function emptyRows(df: DG.DataFrame, c: DayColumns, rows: number[], day: number): void {
  const r = c.readings;
  if (day === 0) {
    for (const i of rows) {
      for (const col of [r.volume, r.weight, r.caliper?.length, r.caliper?.width, c.day])
        col?.set(i, null, false);
    }
  } else {
    const drop = new Set(rows);
    df.rows.removeWhereIdx((i) => drop.has(i));
  }
}

/** Writes the day: every row typed (a new row per animal, or the day's own row changed: a value erased is erased, a
 *  row with every reading erased taken off the day), the animals not measured as notes on the day, and one line in the
 *  log. Checked whole before anything is written: an impossible value writes nothing (thrown). Returns the log
 *  line. */
export function writeDay(df: DG.DataFrame, c: DayColumns, day: number, animals: EntryAnimal[], values: Values):
  string {
  const caliper = !!c.readings.caliper;
  const live = animals.filter((a) => !a.gone);
  const bad = live.filter((a) => {
    const k = rowChecks(a, values[a.animal] ?? {}, caliper);
    return k.volume?.kind === 'impossible' || k.weight?.kind === 'impossible';
  });
  if (bad.length)
    throw new Error(`Nothing saved: correct the impossible values first (${bad.map((a) => a.animal).join(', ')})`);
  const existing = rowsOfDay(df, c, day);
  const r = c.readings;
  let [added, changed] = [0, 0];
  const missed: [string, string][] = [];
  const emptied: number[] = [];
  const erased: string[] = [];
  for (const a of live) {
    const t = values[a.animal] ?? {};
    if (t.missed) {
      missed.push([a.animal, t.missed]);
      continue;
    }
    // a caliper reading half typed is no reading: the row stays as it is
    if (caliper && (t.length === undefined) !== (t.width === undefined))
      continue;
    // every reading the dialog shows, as typed (null: empty)
    const v = stagedVolume(t, caliper);
    const cells: [DG.Column, number | null][] = [[r.volume, v],
      ...r.caliper ? [[r.caliper.length, typedNumber(t.length)], [r.caliper.width, typedNumber(t.width)]] as
        [DG.Column, number | null][] : [],
      ...r.weight ? [[r.weight, typedNumber(t.weight)]] as [DG.Column, number | null][] : []];
    const none = cells.every(([, x]) => x === null);
    let i = existing.get(a.animal);
    if (i === undefined || c.day.isNone(i)) {
      if (none)
        continue;
      if (i === undefined) {
        i = df.rows.addNew(null, false).idx;
        c.animal.set(i, a.animal, false);
        c.group.set(i, a.group, false);
      }
      // (an animal listed before day 0: its row is its day 0)
      c.day.set(i, day, false);
      added++;
    } else if (cells.some(([col, x]) => x === null ? !col.isNone(i!) : col.isNone(i!) || col.get(i!) !== x)) {
      if (none) {
        erased.push(readingText(c, i));
        emptied.push(i);
        continue;
      }
      changed++;
    }
    for (const [col, x] of cells)
      col.set(i, x, false);
  }
  emptyRows(df, c, emptied, day);
  df.fireValuesChanged();
  if (missed.length) {
    const m = metaOf(df);
    if (m) {
      const [who, when] = [grok.shell.user?.login ?? '', new Date().toISOString()];
      setMetaOn(df, {...m, notes: [...m.notes, ...missed.map(([animal, why]): DayNote => ({day, animal, who, when,
        kind: why === 'Missed' ? 'Measurement' : 'Incident', text: `Not measured: ${why.toLowerCase()}`}))]});
    }
  }
  // what was written, in values (an animal's row is one value, whatever it holds): "Day 14: 2 values entered"
  const done = [...added ? [`${plural(added, 'value')} entered`] : [], ...changed ? [`${plural(changed, 'value')} ` +
    'changed'] : [], ...erased.length ? [`${plural(erased.length, 'value')} removed (${erased.join('; ')})`] : []];
  const text = `Day ${day}: ${done.join(', ') || 'nothing new'}` + (missed.length ? `; not measured: ${missed.map(
    ([a, why]) => `${a} (${why.toLowerCase()})`).join(', ')}` : '');
  addLog(df, text);
  return text;
}

/** Where a log entry's list of old values starts: the log shows the line before it, the values on a click. */
export const LOG_VALUES = ' · values: ';

/** Takes every value of `day` off the study (a day saved by mistake), the old values in the log: "Day 11 removed (60
 *  values)", the values after LOG_VALUES (QA 0.1.3 round 3: technician N19, chaos 26). Returns the log line, null
 *  when the day has none. */
export function removeDay(df: DG.DataFrame, c: DayColumns, day: number): string | null {
  const rows = Array.from({length: df.rowCount}, (_, i) => i).filter((i) => !c.day.isNone(i) && c.day.get(i) === day);
  if (!rows.length)
    return null;
  const text = `Day ${day} removed (${plural(rows.length, 'value')})${LOG_VALUES}` +
    rows.map((i) => readingText(c, i)).join('; ');
  emptyRows(df, c, rows, day);
  df.fireValuesChanged();
  addLog(df, text);
  return text;
}

/** Cells pasted from a spreadsheet: a row a line, its cells as typed (a tab between them). */
export const pastedRows = (text: string): string[][] => text.replace(/\r/g, '').replace(/\n+$/, '').split('\n')
  .map((l) => l.split('\t').map((s) => s.trim()));

const SETTINGS = 'PmxStudio';
export const draftKey = (study: string, day: number) => `day-draft-${study}-${day}`;

/** The day's draft: what was typed, when, and the animals the day listed, in order (absent: kept before they were,
 *  or not given) (null: none). */
export function loadDraft(study: string, day: number): {values: Values, when: string, animals?: string[]} | null {
  try {
    const d = JSON.parse(grok.userSettings.getValue(SETTINGS, draftKey(study, day)) ?? 'null');
    if (!d || typeof d.values !== 'object')
      return null;
    const animals = Array.isArray(d.animals) && d.animals.every((a: unknown) => typeof a === 'string') ?
      d.animals as string[] : undefined;
    return {values: d.values, when: d.when, ...animals ? {animals} : {}};
  } catch {
    return null;
  }
}

/** Keeps the day's draft (none: forgets it), with the animals the day lists. */
export function saveDraft(study: string, day: number, values: Values | null, animals?: string[]): void {
  if (values && Object.keys(values).length) {
    grok.userSettings.add(SETTINGS, draftKey(study, day), JSON.stringify({values, when: new Date().toISOString(),
      ...animals?.length ? {animals} : {}}));
  } else
    grok.userSettings.delete(SETTINGS, draftKey(study, day));
  void grok.userSettings.flush();
}

/** What a dialog for the day needs from Monitor. */
export interface DayEntry {
  df: DG.DataFrame;
  cols: DayColumns;
  /** the study's key for its drafts */
  study: string;
  day: number;
  /** the animals on study on a day, in the order they are measured */
  animals: (day: number) => EntryAnimal[];
  /** a word on the day (not a planned day, already entered), or null */
  warning: (day: number) => string | null;
  /** the study's last day measured and its planned end: a day far past both is asked about first */
  limits?: {last: number | null, end: number | null, today?: number | null};
  /** the day's date, in words (null: day 0 has no date) */
  date?: (day: number) => string | null;
  showGroup: boolean;
  /** the day written: its log line (day 0 saved while the table's day 0 date is not today is asked about first,
   *  dayZeroQuestion, and today's answer moves that date in the table's details before this runs) */
  saved: (text: string) => void;
}

/** `values` with an animal's row replaced; an emptied row is no row (typed and erased is not a change). */
const withRow = (values: Values, animal: string, t: Typed): Values => {
  const {[animal]: _, ...rest} = values;
  return Object.keys(t).length ? {...rest, [animal]: t} : rest;
};

const TONE: {[k in ReadingCheck['kind']]: [string, string]} = {ok: ['', ''], again: [CAL.amber, 'Measure again'],
  typo: [CAL.amber, 'Likely a typo'], impossible: [CAL.red, 'Impossible']};

/** The one question of a day 0 saved for a date that is not today (`ask`: dayZeroAsk's): Day 0 is today (`go` with
 *  today's date, YYYY-MM-DD), Keep the planned date (`go` with none), or Cancel (nothing). */
export function dayZeroQuestion(ask: string, go: (start?: string) => void): DG.Dialog | null {
  const q = modal('day0-date', 'Enter day 0 values');
  if (!q)
    return null;
  const now = new Date();
  const today = new Date(Date.UTC(now.getFullYear(), now.getMonth(), now.getDate())).toISOString().slice(0, 10);
  const planned = /planned for ([^;]+);/.exec(ask)?.[1] ?? 'the planned date';
  q.add(dialogLine(ask)).onOK(() => go(today));
  okNamed(q, 'Day 0 is today');
  q.addButton(`Keep ${planned}`, () => {
    q.close();
    go();
  });
  showCentred(q);
  return q;
}

/** Shows the dialog. */
export function enterDayDialog(o: DayEntry): void {
  const d = modal('enter-day', 'Enter today\'s values');
  if (!d)
    return;
  const caliper = !!o.cols.readings.caliper;
  const weighs = !!o.cols.readings.weight;
  let day = o.day;
  let animals: EntryAnimal[] = [];
  let initial: Values = {};
  let values: Values = {};
  // typed as text: "2.5" or "abc" is said, not dropped
  const dayInput = ui.input.string('Day', {value: String(day)});
  const dateLabel = ui.span([], {style: {color: 'var(--grey-5, #68727f)', marginLeft: '8px', alignSelf: 'center'}});
  dayInput.root.append(dateLabel);
  const dayError = dialogLine('');
  dayError.style.color = CAL.red;
  dayError.hidden = true;
  const warning = dialogLine('');
  warning.style.color = CAL.amber;
  const draftLine = ui.div([]);
  // a day saved by mistake: taken off whole, asked first
  const removeLine = ui.div([]);
  const body = document.createElement('tbody');
  const footer = ui.divText('');
  const table = document.createElement('table');
  Object.assign(table.style, {borderCollapse: 'collapse', width: '100%'});
  const head = ['Animal', 'Last', 'Expected, mm³', ...caliper ? ['Length, mm', 'Width, mm', 'Volume, mm³'] :
    ['Volume, mm³'], ...weighs ? ['Body weight, g'] : [], '', ''];
  table.innerHTML = `<thead><tr>${head.map((h) => `<th style="text-align:left;padding:2px 6px;font-weight:500">` +
    `${h}</th>`).join('')}</tr></thead>`;
  table.append(body);
  const scroll = ui.div([table]);
  Object.assign(scroll.style, {maxHeight: '55vh', overflowY: 'auto', margin: '6px 0'});

  const dirty = () => [...new Set([...Object.keys(values), ...Object.keys(initial)])]
    .some((a) => JSON.stringify(values[a] ?? {}) !== JSON.stringify(initial[a] ?? {}));
  // the day's draft, with the animals it lists: a study with no animals in its table yet has them only here (their IDs
  // as typed, QA 0.1.3 r5: reopened, the default IDs came back and the values typed under M1…M9 were lost)
  const draft = (at: number, v: Values | null) => saveDraft(o.study, at, v, animals.map((a) => a.animal));
  let saving: ReturnType<typeof setTimeout> | undefined;
  const keep = () => {
    clearTimeout(saving);
    saving = setTimeout(() => draft(day, dirty() ? values : null), 400);
  };
  const canSave = () => dirty() && dayError.hidden && fieldIsDay();
  // what day 0 saved for a date that is not today asks (null: another day, or day 0 is today)
  const dateAsk = () => day === 0 ? dayZeroAsk(metaOf(o.df)?.start) : null;
  const moveStart = (start?: string) => {
    if (start)
      setMetaOn(o.df, {...metaOf(o.df)!, start});
  };
  const write = (at: number, list: EntryAnimal[], typed: Values) => {
    try {
      const text = writeDay(o.df, o.cols, at, list, typed);
      saveDraft(o.study, at, null);
      o.saved(text);
    } catch (e: any) {
      // (kept: the draft has every value)
      draft(at, typed);
      grok.shell.warning(String(e?.message ?? e));
    }
  };
  let okButton: HTMLButtonElement | null = null;
  const refresh = () => {
    const s = entrySummary(animals, values, caliper);
    footer.textContent = summaryText(s, animals.filter((a) => !a.gone).length);
    const b = okButton;
    if (b) {
      // Save day writes the day the field shows, and only that one: a day typed and not taken yet (a far day still
      // asked about) saves nothing (QA 0.1.3 round 3: chaos 3, 27)
      const bad = !dayError.hidden;
      const other = !bad && !fieldIsDay();
      b.disabled = bad || other || s.impossible > 0 || !dirty();
      b.title = bad ? 'Type a whole study day first' : other ? `Confirm day ${(dayInput.value ?? '').trim()} first` :
        s.impossible ? 'Correct the impossible values first' : !dirty() ? 'Nothing typed yet' : '';
    }
  };

  // the reading fields in each row's order: a paste fills its row's next fields from its next columns
  const keys: (keyof Typed)[] = [...caliper ? ['length', 'width'] as const : ['volume'] as const,
    ...weighs ? ['weight'] as const : []];
  const field = (a: EntryAnimal, key: keyof Typed, row: number, paint: () => void) => {
    const e = document.createElement('input');
    e.inputMode = 'decimal';
    e.value = values[a.animal]?.[key] ?? '';
    e.dataset.col = key;
    e.dataset.row = String(row);
    Object.assign(e.style, {width: '72px', padding: '2px 4px', border: '1px solid var(--grey-2, #dbdcdf)',
      borderRadius: '3px'});
    e.addEventListener('input', () => {
      const t = {...values[a.animal]};
      if (e.value.trim())
        t[key] = e.value.trim();
      else
        delete t[key];
      values = withRow(values, a.animal, t);
      paint();
      refresh();
      keep();
    });
    // cells pasted from a spreadsheet fill down from here, and across (volume, then weight) when they are several
    // columns
    e.addEventListener('paste', (ev) => {
      const cells = pastedRows(ev.clipboardData?.getData('text') ?? '');
      if (cells.length < 2 && (cells[0]?.length ?? 0) < 2)
        return;
      ev.preventDefault();
      let taken = 0;
      keys.slice(keys.indexOf(key)).forEach((k, j) => {
        const all = [...body.querySelectorAll<HTMLInputElement>(`input[data-col="${k}"]`)];
        const room = all.filter((x) => Number(x.dataset.row) >= row && !x.disabled).slice(0, cells.length);
        taken = Math.max(taken, room.length);
        room.forEach((x, n) => {
          if (cells[n][j] === undefined)
            return;
          x.value = cells[n][j];
          x.dispatchEvent(new Event('input'));
        });
      });
      // more rows than animals from here: the rest said, not dropped silently (QA 0.1.3 round 4 chaos 19)
      if (cells.length > taken) {
        grok.shell.warning(`${plural(cells.length - taken, 'pasted row')} left out: the list has ${plural(taken,
          'animal')} from here on.`);
      }
    });
    return e;
  };

  const render = () => {
    body.replaceChildren(...animals.map((a, row) => {
      const tr = document.createElement('tr');
      tr.style.borderTop = '1px solid var(--grey-1, #eceef0)';
      const name = ui.span([a.animal]);
      if (o.showGroup && a.group)
        name.append(ui.span([` · ${a.group}`], {style: {color: 'var(--grey-4, #9aa5b1)'}}));
      const status = ui.span([]);
      const live = ui.span([]);
      const inputs: HTMLInputElement[] = [];
      const paint = () => {
        const t = values[a.animal] ?? {};
        const c = rowChecks(a, t, caliper);
        const worst = [c.volume, c.weight].filter((x) => x && x.kind !== 'ok')
          .sort((x, y) => (y!.kind === 'impossible' ? 1 : 0) - (x!.kind === 'impossible' ? 1 : 0))[0];
        for (const e of inputs) {
          const k = e.dataset.col === 'weight' ? c.weight : c.volume;
          const [color] = TONE[k?.kind ?? 'ok'];
          e.style.borderColor = color || 'var(--grey-2, #dbdcdf)';
          e.title = a.gone ?? k?.why ?? '';
          e.disabled = !!t.missed || !!a.gone;
        }
        const v = stagedVolume(t, caliper);
        live.textContent = caliper && v !== null && Number.isFinite(v) ? mm3Text(v) : '';
        // (an animal off the study, or one that missed the last session, says so in grey)
        const quiet = a.gone ?? (t.missed ? `Not measured: ${t.missed.toLowerCase()}` : worst ? null : a.note ?? null);
        // (a value refused says why in the words Correct uses, not "Impossible": QA 0.1.3 round 4 chaos 8)
        status.textContent = quiet ?? (worst ? worst.kind === 'impossible' ? worst.why : TONE[worst.kind][1] : '');
        status.style.color = quiet ? 'var(--grey-4, #9aa5b1)' : worst ? TONE[worst.kind][0] : '';
        status.title = worst?.why ?? '';
        tr.style.opacity = t.missed || a.gone ? '0.6' : '';
      };
      const cells = keys.map((k) => {
        const e = field(a, k, row, paint);
        inputs.push(e);
        return e;
      });
      // not measured, with its reason; back to measured
      const more = ui.iconFA('ellipsis-h', (ev: MouseEvent) => {
        const menu = DG.Menu.popup();
        const set = (why?: string) => {
          const t = {...values[a.animal]};
          if (why)
            t.missed = why;
          else
            delete t.missed;
          values = withRow(values, a.animal, t);
          paint();
          refresh();
          keep();
        };
        if (a.gone)
          return;
        if (values[a.animal]?.missed)
          menu.item('Measured after all', () => set());
        else if (!(a.animal in initial)) {
          for (const why of MISSED)
            menu.item(`Not measured: ${why.toLowerCase()}`, () => set(why));
        }
        menu.show({causedBy: ev});
      }, 'Not measured (dead, removed, missed)');
      const td = (...x: (HTMLElement | string)[]) => {
        const c = document.createElement('td');
        c.style.padding = '2px 6px';
        c.append(...x);
        return c;
      };
      const last = a.lastVolume === null ? '—' : `${mm3Text(a.lastVolume)} mm³, day ${a.lastDay}`;
      tr.append(td(name), td(last), td(a.expected), ...cells.map((e) => td(e)),
        ...caliper ? [td(live)] : [], td(status), td(more));
      paint();
      return tr;
    }));
    refresh();
  };

  const load = (n: number) => {
    day = n;
    animals = o.animals(day);
    initial = valuesOfDay(o.df, o.cols, day);
    const draft = loadDraft(o.study, day);
    values = {...initial, ...draft?.values};
    dayLines();
    draftLine.replaceChildren(...draft ? [ui.divH([ui.span([`Your draft is back (typed ${new Date(draft.when)
      .toLocaleString()}). `]), ui.link('Discard it', () => {
      saveDraft(o.study, day, null);
      load(day);
    })])] : []);
    removeLine.replaceChildren(...Object.keys(initial).length ? [ui.link('Remove this day\'s values…',
      () => askRemove(), `Takes every value of day ${day} off the study; the log keeps them`)] : []);
    render();
  };
  const askRemove = () => {
    const q = modal('enter-day-remove', `Remove day ${day}'s values?`);
    if (!q)
      return;
    q.add(dialogLine(`Every value of day ${day} (${plural(Object.keys(initial).length, 'animal')}) is taken off the ` +
      'study now. The log keeps the old values.')).onOK(() => {
      clearTimeout(saving);
      const text = removeDay(o.df, o.cols, day);
      saveDraft(o.study, day, null);
      if (text)
        o.saved(text);
      d.close();
    });
    okNamed(q, 'Remove');
    showCentred(q);
  };
  // the day's date and its word (not planned, already entered), for the day shown; none under a Day error: they
  // would be the last day's, not the one typed (QA 0.1.3 round 3: technician N13, designer U1)
  const dayLines = () => {
    const w = dayError.hidden ? o.warning(day) : null;
    warning.textContent = w ?? '';
    warning.hidden = !w;
    dateLabel.textContent = dayError.hidden ? o.date?.(day) ?? '' : '';
  };
  // the field says the day the list is for (a far day asked about and not taken yet: not)
  const fieldIsDay = () => readDay(dayInput.value ?? '', o.limits ?? {last: null, end: null}).day === day;
  // the far day's question open: tied to the value it asks about, any change of the field closes it unanswered
  let farAsk: DG.Dialog | null = null;
  const toDay = (n: number) => {
    dayError.hidden = true;
    if (n === day) {
      dayLines();
      return refresh();
    }
    clearTimeout(saving);
    draft(day, dirty() ? values : null);
    load(n);
  };
  dayInput.input.addEventListener('input', () => {
    const q = farAsk;
    farAsk = null;
    q?.close();
    refresh();
  });
  dayInput.input.addEventListener('change', () => {
    const typed = dayInput.value ?? '';
    const r = readDay(typed, o.limits ?? {last: null, end: null});
    dayError.textContent = r.error ?? '';
    dayError.hidden = !r.error;
    if (r.error) {
      dayLines();
      return refresh();
    }
    if (!r.ask)
      return toDay(r.day!);
    // a day far from the study (9999 for 9): asked, else back to the day shown
    farAsk?.close();
    const q = modal('enter-day-far', 'Check the day');
    if (!q)
      return;
    farAsk = q;
    let took = false;
    q.add(dialogLine(r.ask)).onOK(() => {
      // only the value it asked about: the field changed since, nothing is taken
      if (farAsk !== q || (dayInput.value ?? '') !== typed)
        return;
      took = true;
      farAsk = null;
      toDay(r.day!);
    });
    // (closed any other way: the day shown stays; closed by a change of the field, what is typed stays)
    q.onClose.subscribe(() => {
      if (farAsk !== q)
        return;
      farAsk = null;
      if (!took)
        dayInput.value = String(day);
      refresh();
    });
    okNamed(q, `Enter day ${r.day}`);
    showCentred(q);
  });

  // Enter and Shift+Enter: the next and the previous field; the arrows: the same field of the next and previous animal
  scroll.addEventListener('keydown', (e) => {
    const at = e.target as HTMLInputElement;
    if (at.tagName !== 'INPUT')
      return;
    const all = [...body.querySelectorAll<HTMLInputElement>('input')].filter((x) => !x.disabled);
    let to: HTMLInputElement | undefined;
    if (e.key === 'Enter')
      to = all[all.indexOf(at) + (e.shiftKey ? -1 : 1)];
    else if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
      const same = all.filter((x) => x.dataset.col === at.dataset.col);
      to = same[same.indexOf(at) + (e.key === 'ArrowDown' ? 1 : -1)];
    } else
      return;
    e.preventDefault();
    e.stopPropagation();
    to?.focus();
    to?.select();
  });

  d.add(ui.inputs([dayInput])).add(dayError).add(warning).add(draftLine).add(removeLine)
    .add(dialogLine('Nothing reaches the study until Save day', 'Each value is checked as you type it, against the ' +
      'range Monitor expects: outside it, measure again. Enter moves on, Shift+Enter back, the arrows up and down. ' +
      'What you type is kept as a draft until the day is saved.'))
    .add(scroll).add(footer)
    .onOK(() => {
      clearTimeout(saving);
      // (Enter can reach OK while Save day is off: nothing typed is nothing to write, nor a day not read, nor one the
      // field no longer shows)
      if (!canSave())
        return;
      const [at, list, typed] = [day, animals, values];
      const ask = dateAsk();
      if (!ask)
        return write(at, list, typed);
      // (only if the platform's Enter got here first, the dialog closing: the draft keeps every value and its IDs)
      draft(at, typed);
      dayZeroQuestion(ask, (start) => {
        moveStart(start);
        write(at, list, typed);
      });
    });
  okButton = okNamed(d, 'Save day');
  // day 0 saved for a date that is not today: asked first, whichever way it was typed or pasted (QA 0.1.3 round 4
  // designer D5), the entry staying open under the question: Cancel there closes only the question, everything typed
  // still in place (QA 0.1.3 r5: it closed the entry too). The date is the table's, not the caller's: Monitor's cached
  // details of a followed plan had none (QA 0.1.3 r5)
  const askDate = (e: Event) => {
    const ask = canSave() ? dateAsk() : null;
    if (!ask)
      return;
    e.preventDefault();
    e.stopImmediatePropagation();
    clearTimeout(saving);
    const [at, list, typed] = [day, animals, values];
    draft(at, typed);
    dayZeroQuestion(ask, (start) => {
      moveStart(start);
      write(at, list, typed);
      d.close();
    });
  };
  d.root.addEventListener('click', (e) => {
    if (okButton && (e.target as HTMLElement).closest('button, .ui-btn') === okButton)
      askDate(e);
  }, true);
  // (Enter outside the table, where it moves on, saves as the button does: through the question)
  d.root.addEventListener('keydown', (e) => {
    if (e.key === 'Enter' && !scroll.contains(e.target as Node) && !modalOf('day0-date'))
      askDate(e);
  }, true);

  // Cancel and Esc with values typed: asked first (the draft keeps them either way)
  const cancel = d.getButton('CANCEL');
  const ask = (e: Event) => {
    if (!dirty())
      return;
    e.preventDefault();
    e.stopImmediatePropagation();
    clearTimeout(saving);
    draft(day, values);
    const q = modal('enter-day-close', 'Close without saving the day?');
    if (!q)
      return;
    q.add(dialogLine(`What you typed stays as a draft for day ${day}, offered when you open this again.`))
      .onOK(() => d.close());
    q.addButton('Discard it', () => {
      saveDraft(o.study, day, null);
      q.close();
      d.close();
    });
    okNamed(q, 'Keep the draft');
    showCentred(q);
  };
  d.root.addEventListener('click', (e) => {
    if (cancel && (e.target as HTMLElement).closest('button, .ui-btn') === cancel)
      ask(e);
  }, true);
  const esc = (e: KeyboardEvent) => {
    // (the question's own Esc closes the question, not this)
    if (!d.root.isConnected)
      return window.removeEventListener('keydown', esc, true);
    if (e.key === 'Escape' && !['enter-day-close', 'enter-day-far', 'enter-day-remove'].some((p) => modalOf(p)))
      ask(e);
  };
  window.addEventListener('keydown', esc, true);
  load(day);
  showCentred(d, {width: 760});
  (body.querySelector('input:not([disabled])') as HTMLInputElement | null)?.focus();
}
