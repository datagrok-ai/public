/* The Plan view: the declaration form is generated from pmxPlanStudy's own signature (u2 funcForm, grouped by
   the parameters' categories), so the function Grokky reads and the form people fill in are one definition; the
   treated arms and the multiplicity correction have editors of their own (arms-input.ts), and the measurement and
   randomization weekdays are set on the study calendar. The results are the function's outputs: the study's
   figures, the groups, the predicted volumes (band-chart.ts), the study calendar, then, folded, the power by
   group size, by end of study (platform line charts) and by measurement schedule, and the study protocol. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BasicTable, buttonWithIcon, ChoiceInput, computed, Dialog, div, divH, divV, dot, loader, Reason, signal, span,
  Splitter, tableFromMap, TextInput, Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal, StatValue} from '@datagrok-libraries/u2';
import {funcForm, ViewHost, viewers} from '@datagrok-libraries/u2/dg';
import type {ViewContent} from '@datagrok-libraries/u2/dg';

import {PlanFrames, PlanInputs, WEEKDAY_NAMES, armColors, declarationOf, planFrames, tgisOf, weekdaysOf,
  weekdaysText} from './frames';
import {ArmsInput, CorrectionInput, armsFrame} from './arms-input';
import {MIN_GROUP, finite} from './engine';
import {askedOf, autoRun, busyBar, busyOverlay, callInputs, cardGrid, copyRow, declarationPane, DECLARATION_MIN,
  infoNote, notes,
  notesBox, num, pct, RESULTS_MIN, resultsFrame, settleForm, sized, staleNote,
  statCard, swatch, volumeAxis, wrapCells} from './ui';
import {openStage} from './stages';
import {panelGrid, Slot} from './panels';
import {studyHeader} from './study-header';
import {followPlan} from './monitor-view';
import {SavedPlan, downloadPlan, parsePlan, planOf, planSaved, plansFolder, savePlan} from './plan-file';
import {DayNote, StudyMeta, dateOf, dateText, dosingDays, emptyMeta, noteDialog, startFrom, weekdayOf, withId}
  from './study-meta';
import {LineInput, lastLine} from './line-input';
import {BandChart} from './band-chart';

/** Where the plan's first stage is run: Monitor takes the implanted pool to randomization. */
/** Plan's panels, in order and width (twelfths): the groups, their growth and power, the calendar, the trade-offs,
 *  the scenarios, the protocol. */
const PLAN_LAYOUT: Slot[] = [{id: 'groups', span: 12}, {id: 'growth', span: 8}, {id: 'power', span: 4},
  {id: 'calendar', span: 12}, {id: 'end', span: 6}, {id: 'schedule', span: 6}, {id: 'scenarios', span: 12},
  {id: 'protocol', span: 12}];

const monitorNote = () => infoNote('<b>From randomization on, in Monitor.</b> Follow this study to form the groups, ' +
  'balanced on tumor volume, and to see each session who needs attention and how the study is going.');

const WEEKDAYS = WEEKDAY_NAMES.map((d) => d.slice(0, 3));

/** What a section with nothing to compare says: no arm reaches the target power with any group size. */
const UNREACHABLE = 'No group size up to 40 reaches the target power: lower it, or set the animals per arm.';

/** Where a study's total has no number: the limit is per group, and totals over 40 are common. */
const OVER_40 = 'over 40 per group';

/** Where a day's menu opens (DG.Menu.show's options): at the mouse, or next to the focused day. */
type MenuAt = {causedBy?: MouseEvent, element?: HTMLElement, nextToElement?: boolean};

/** A day's mark in the calendar: a short label in its kind's colour. */
function mark(text: string, color: string, filled: boolean): HTMLElement {
  const el = span(text);
  Object.assign(el.style, {display: 'inline-block', padding: '1px 6px', borderRadius: '3px', fontSize: '11px',
    lineHeight: '16px', whiteSpace: 'nowrap', fontWeight: '500', color: filled ? '#fff' : color,
    background: filled ? color : `color-mix(in srgb, ${color} 12%, transparent)`});
  return el;
}

/** The study as a wall calendar: a row a week, Monday to Sunday, day 0 on the randomization weekday (the date is
 *  known only when the tumors reach the volume), each day with what is done and its details in the tooltip; its
 *  header the lab's measurement weekdays, each one clickable (`toggle`); above it, the randomization weekday
 *  (`startDay`, its picker) and what the weeks add up to. The study's details add each day's date once day 0 has
 *  one, the doses its groups' schedules give and the notebook's notes; a day's menu (`onMenu`: right click, or the
 *  Menu key on a focused day) and a double click (`onNote`) are the study's. Built from u2's primitives (u2 has no
 *  calendar); call it in a scope, which owns the tooltips. */
function calendarGrid(p: DG.DataFrame, calendar: DG.DataFrame, startDay: HTMLElement,
  toggle: (weekday: number) => void, study: {meta: StudyMeta, onMenu: (day: number, at: MenuAt) => void,
    onNote: (day: number) => void}): HTMLElement {
  const v = (c: string) => p.get(c, 0);
  const end = v('study_end') as number;
  const endpointDay = v('endpoint_day') as number;
  const humane = v('humane_limit');
  const watch = v('humane_watch_day') as number | null;
  const details = new Map<number, string>(Array.from({length: calendar.rowCount},
    (_, i) => [calendar.get('day', i), calendar.get('activity', i)]));
  const measured = new Set(String(v('measurement_days')).split(', ').map(Number));
  const weekdays = new Set(weekdaysOf(v('measurement_weekdays')));
  const [BLUE, SLATE, RED, DARK, GREEN, AMBER] = ['#2083d5', '#5b6b7f', '#c0504d', '#4d5261', '#2e8b57', '#b7791f'];
  const meta = study.meta;
  const date = (day: number) => dateOf(meta.start, day);
  const doses = dosingDays(meta.groups, end);
  const notesOn = new Map<number, DayNote[]>();
  for (const n of meta.notes)
    notesOn.set(n.day, [...notesOn.get(n.day) ?? [], n]);

  // day 0 falls on the randomization weekday: the measurement weekdays are the lab's, seen from it
  const start = (v('start_weekday') ?? 0) as number;
  const weekday = (day: number) => (start + day) % 7;

  // above the weeks: the randomization weekday (picked in the sentence), the end and the closer watch
  const q90 = v('duration_q90');
  const range = typeof q90 === 'number' && Number.isFinite(q90) ? `day ${Math.round(v('duration_q10'))} to ` +
    `${Math.round(q90)}` : `from day ${Math.round(v('duration_q10'))}`;
  const first = divH([mark('Randomize', BLUE, true), span('Randomize on a'), startDay,
    span(`(day 0${meta.start ? `, ${dateText(date(0)!, true)}` : ''}) · ${measured.size} measurement days.`)]);
  Object.assign(first.style, {alignItems: 'center', gap: '6px', flexWrap: 'wrap'});
  const lines: HTMLElement[] = [first];
  // each line after the mark its day carries in the weeks below
  const marked = (m: HTMLElement, text: string) => {
    const line = divH([m, span(text)]);
    Object.assign(line.style, {alignItems: 'baseline', gap: '6px'});
    return line;
  };
  if (v('last_day_chosen')) {
    lines.push(marked(mark('End', DARK, true), `The study ends on day ${end}, as chosen` +
      (weekdays.has(weekday(end)) ? '.' : `: a ${WEEKDAY_NAMES[weekday(end)]}, not one of your measurement days.`)));
  } else {
    const reach = `the median control reaches ${v('endpoint')} mm³`;
    lines.push(marked(mark('End', DARK, true), end === endpointDay ? `The study ends on day ${end}, when ${reach} ` +
      `(in 8 of 10 studies, ${range}).` : `${reach[0].toUpperCase()}${reach.slice(1)} on day ${endpointDay} (in 8 ` +
      `of 10 studies, ${range}), a ${WEEKDAY_NAMES[weekday(endpointDay)]}: the study ends on the next measurement ` +
      `day, day ${end} (${WEEKDAY_NAMES[weekday(end)]}).`));
  }
  if (watch !== null) {
    lines.push(marked(mark('Watch controls', RED, false), `From day ${watch}, check the controls more often: an ` +
      `animal over ${humane} mm³ is removed.`));
  }
  if (doses.size > 0) {
    const dosed = meta.groups.filter((g) => [...doses.values()].some((gs) => gs.includes(g.group)));
    lines.push(marked(mark('Dose', GREEN, false), dosed.map((g) => `${g.group} ${g.schedule}` +
      (g.compound ? ` (${g.compound})` : '')).join('; ') + ', from day 0.'));
  }
  const summary = divV(lines, 'pmx-calendar-summary');
  Object.assign(summary.style, {gap: '4px', lineHeight: '1.45', marginBottom: 'var(--dg-space-m)'});

  // the weeks: a header of the lab's weekdays (click one to measure on it or not), then a cell a day up to the end
  const grid = div([], 'pmx-calendar');
  Object.assign(grid.style, {display: 'grid', gridTemplateColumns: '56px repeat(7, minmax(0, 1fr))', gap: '4px'});
  grid.append(span(''), ...WEEKDAYS.map((d, i) => {
    const on = weekdays.has(i);
    const h = span(d, 'pmx-calendar-weekday');
    Object.assign(h.style, {fontSize: '11px', textAlign: 'center', padding: '2px 0', borderRadius: '3px',
      cursor: 'pointer', userSelect: 'none', fontWeight: on ? '600' : '400',
      color: on ? SLATE : 'var(--dg-text-color-muted, #9aa3ad)',
      background: on ? `color-mix(in srgb, ${SLATE} 12%, transparent)` : 'transparent'});
    Tooltip.bind(h, on ? (weekdays.size > 1 ? `Measured on ${WEEKDAY_NAMES[i]}s: click to stop` :
      'The only measurement day') : `Not measured on ${WEEKDAY_NAMES[i]}s: click to measure`);
    if (!on || weekdays.size > 1)
      h.addEventListener('click', () => toggle(i));
    return h;
  }));
  // calendar weeks (Monday first): the days before randomization and after the end left blank
  const weeks = Math.floor((start + end) / 7) + 1;
  for (let w = 0; w < weeks; w++) {
    const label = span(`Week ${w + 1}`);
    Object.assign(label.style, {fontSize: '11px', color: 'var(--dg-text-color-muted, #68727f)', alignSelf: 'center'});
    grid.append(label);
    for (let d = 0; d < 7; d++) {
      const day = 7 * w + d - start;
      const cell = divV([], 'pmx-calendar-day');
      Object.assign(cell.style, {minHeight: '48px', padding: '4px 6px', gap: '3px', borderRadius: '4px',
        border: '1px solid var(--dg-border-color, #e1e6ec)', boxSizing: 'border-box', alignItems: 'flex-start'});
      if (day < 0 || day > end) {
        cell.style.borderStyle = 'dashed';
        cell.style.opacity = '0.4';
        grid.append(cell);
        continue;
      }
      const on = date(day);
      const num = span(`Day ${day}${on ? ` · ${dateText(on)}` : ''}`);
      Object.assign(num.style, {fontSize: '11px', color: 'var(--dg-text-color-muted, #68727f)'});
      cell.append(num);
      if (day === 0)
        cell.append(mark('Randomize', BLUE, true));
      if (measured.has(day))
        cell.append(mark('Measure', SLATE, false));
      if (doses.has(day))
        cell.append(mark('Dose', GREEN, false));
      const notes = notesOn.get(day) ?? [];
      if (notes.length > 0)
        cell.append(mark(notes.length === 1 ? '1 note' : `${notes.length} notes`, AMBER, false));
      if (day === watch)
        cell.append(mark('Watch controls', RED, false));
      if (day === end) {
        cell.append(mark('End', DARK, true));
        cell.style.borderColor = DARK;
      }
      // the days the controls need a closer watch, tinted
      if (watch !== null && day >= watch)
        cell.style.background = `color-mix(in srgb, ${RED} 6%, transparent)`;
      const detail = details.get(day);
      const tip = [`Day ${day}${on ? ` (${dateText(on, true)})` : ''}${detail ? `: ${detail}` : ''}`,
        ...doses.has(day) ? [`Dose: ${doses.get(day)!.join(', ')}`] : [],
        ...notes.map((n) => `${n.kind}${n.animal ? ` · ${n.animal}` : ''}: ${n.text} (${n.who})`)];
      Tooltip.bind(cell, tip.join('\n'));
      // the day's own menu, from the mouse or the keyboard; a double click notes what happened
      cell.tabIndex = 0;
      cell.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        study.onMenu(day, {causedBy: e});
      });
      cell.addEventListener('keydown', (e) => {
        if (e.key === 'ContextMenu' || (e.key === 'F10' && e.shiftKey)) {
          e.preventDefault();
          study.onMenu(day, {element: cell, nextToElement: true});
        }
      });
      cell.addEventListener('dblclick', () => study.onNote(day));
      grid.append(cell);
    }
  }
  return divV([summary, grid], 'pmx-calendar-box');
}

/** One row per measurement schedule: its measurement days, the animals it needs and each arm's power with the
 *  plan's animals (`own`: the plan's measurements a week). */
function schedulesTable(sc: DG.DataFrame, own: number | undefined): HTMLElement {
  type Schedule = [number, {measurements: number, nTotal: number | null, power: {[arm: string]: number}}];
  const byWeek = new Map<number, Schedule[1]>();
  for (let i = 0; i < sc.rowCount; i++) {
    const k = sc.get('per_week', i);
    const row: Schedule[1] = byWeek.get(k) ?? {measurements: sc.get('measurements', i),
      nTotal: finite(sc.get('n_total', i)), power: {}};
    row.power[sc.get('arm', i)] = sc.get('assurance', i);
    byWeek.set(k, row);
  }
  const armNames = [...new Set(sc.col('arm')!.toList() as string[])];
  const columns = [
    {header: 'Measurements a week', render: ([k]: Schedule) => k === own ? `${k} · this plan` : String(k)},
    {header: 'Measurement days', render: ([, s]: Schedule) => String(s.measurements)},
    {header: 'Analyzable animals', render: ([k, s]: Schedule) => animals(s.nTotal ?? OVER_40, k === own)},
    ...armNames.map((arm) => ({header: armNames.length > 1 ? `Power, ${arm}` : 'Power',
      render: ([, s]: Schedule) => pct(s.power[arm])})),
  ];
  const table = new BasicTable<Schedule>({items: [...byWeek], keyOf: ([k]) => String(k),
    columns: columns.map((c, j) => ({...c, align: j >= 1 ? 'center' : 'left'}))});
  table.root.style.width = '100%';
  return table.root;
}

/** A count of animals, set apart from the table's other figures. */
function animals(n: unknown, primary: boolean): HTMLElement {
  const el = span(n === null || n === undefined ? '—' : String(n));
  Object.assign(el.style, {fontVariantNumeric: 'tabular-nums', fontWeight: primary ? '600' : '500',
    fontSize: primary ? '15px' : '13px', color: primary ? 'var(--dg-accent, #2083d5)' : 'inherit'});
  return el;
}

/** One row per assumed TGI: the animals to implant and randomize on the plan's schedule (`own`), and the analyzable
 *  animals on each schedule. */
function scenariosTable(sc: DG.DataFrame, own: number | undefined): HTMLElement {
  type Row = {tgi: number, implant: number | null, randomise: number | null, byWeek: Map<number, number | null>};
  const rows = new Map<number, Row>();
  for (let i = 0; i < sc.rowCount; i++) {
    const tgi = sc.get('tgi', i);
    const row = rows.get(tgi) ?? {tgi, implant: finite(sc.get('implant_total', i)),
      randomise: finite(sc.get('randomise_total', i)), byWeek: new Map()};
    row.byWeek.set(sc.get('per_week', i), finite(sc.get('n_total', i)));
    rows.set(tgi, row);
  }
  const weeks = [...new Set(sc.col('per_week')!.toList() as number[])];
  const columns = [
    {header: 'Assumed TGI', render: (r: Row) => `${r.tgi}%`},
    {header: 'To implant', render: (r: Row) => animals(r.implant, true)},
    {header: 'To randomize', render: (r: Row) => animals(r.randomise, false)},
    ...weeks.map((k) => ({header: `Analyzable, ${k} a week${k === own ? ' · this plan' : ''}`,
      render: (r: Row) => animals(r.byWeek.get(k) ?? OVER_40, false)})),
  ];
  const table = new BasicTable<Row>({items: [...rows.values()], keyOf: (r) => String(r.tgi),
    columns: columns.map((c, j) => ({...c, align: j >= 1 ? 'center' : 'left'}))});
  table.root.style.width = '100%';
  wrapCells(table.root);
  return table.root;
}

/** A .pmxplan file opened from Files or dropped on the platform: Plan with it (no table to open). */
export async function openPlanFile(content: string): Promise<DG.DataFrame[]> {
  const view = grok.shell.addView(await planView(parsePlan(content)));
  // as stages.ts opens a stage: the app's root before the view's path
  view.path = '/apps/PmxStudio/plan';
  return [];
}

/** The Plan view; given a saved plan, its declaration in the form and planned again. */
export async function planView(initial?: SavedPlan): Promise<DG.ViewBase> {
  return ViewHost.create({name: 'PMX Studio · Plan', content: (self) => planContent(self, initial)}).view;
}

/** What the view shows, built under its host, which owns it all. */
function planContent(self: ViewHost, initial?: SavedPlan): ViewContent {
  const call = DG.Func.byName('PmxStudio:pmxPlanStudy').prepare();
  // the engine itself on the form's inputs (with `extra` over them), not the registered function's FuncCall (see
  // callInputs); a wrong field rejects, as the call did
  const designed = (extra: object = {}) => Promise.resolve().then(() => {
    const x = callInputs<PlanInputs & {tgiScenarios?: string | null}>(call, extra);
    return planFrames(declarationOf(x), tgisOf(x.tgiScenarios));
  });
  const result = signal<PlanFrames | undefined>(undefined);
  // the study's details (header's name, the calendar's dates, doses and notes); an ID once the plan leaves Plan
  // (saved, downloaded or followed)
  const meta = signal<StudyMeta>(initial?.study ?? emptyMeta());
  // the form's and the calendar's replans share the call, one run at a time: a change made during a run is planned
  // next, and the run it overtook is never shown
  let running: Promise<void> | null = null;
  let again = false;
  const error = signal<string | null>(null);
  const armsInput = new ArmsInput({label: 'Arms'});
  const correction = new CorrectionInput({label: 'Multiple comparisons', arms: armsInput.count});
  // a plan that no longer matches the form is dimmed until the next run
  const stale = signal(false);
  // captions above every field: the panel reads top to bottom, left-aligned, at any width
  // no Run: each change replans once the changes stop
  const auto = autoRun(() => form, () => replan());
  const form = funcForm(call, {layout: 'tall',
    overrides: {arms: {input: armsInput}, correction: {input: correction},
      line: {input: new LineInput({label: 'Tumor model'})}},
    onInputChanged: () => {
      stale.value = result.peek() !== undefined;
      auto.changed();
    }});
  settleForm(form, {n: 'Auto (target power)', lastDay: 'Auto (endpoint reached)',
    tgiScenarios: 'Off, or e.g. 30, 50, 70'});

  // busy while the view plans (the first plan, a change of the form, a change of randomization day)
  const planning = signal(true);
  // a change made on the calendar replans the study too (the measurement days set the power and the animals), but
  // under a loader of the calendar's own, the page left in view
  const refreshing = signal(false);
  const busy = computed(() => planning.value || refreshing.value);
  const replan = (fromCalendar = false): Promise<void> => {
    if (running) {
      again = true;
      return running;
    }
    const flag = fromCalendar ? refreshing : planning;
    flag.value = true;
    const once = (): Promise<void> => {
      const sent = askedOf(call);
      const next = () => {
        again = false;
        return once();
      };
      return designed().then((c) => {
        if (again)
          return next();
        error.value = null;
        stale.value = askedOf(call) !== sent;
        result.value = c;
      }, (e) => {
        if (again)
          return next();
        error.value = String(e?.message ?? e);
      });
    };
    running = once().finally(() => {
      running = null;
      flag.value = false;
    });
    return running;
  };

  const plan = computed(() => result.value?.plan);
  const cell = <T = number>(name: string) => computed(() => (plan.value?.get(name, 0) ?? undefined) as T | undefined);
  const groups = computed(() => (cell('treated_arms').value ?? 1) + 1);
  // the animals in the order the study uses them, then its end; what each count covers in its tooltip
  const card = (label: string, value: ReadonlySignal<StatValue | undefined>, tip: () => string,
    format?: (x: StatValue) => string) => statCard(form, label, value, tip, format);
  // a count no group size reaches has no number (a column with no value reads as an empty string)
  const count = (name: string) => computed(() => plan.value ? finite(plan.value.get(name, 0)) ?? '—' : undefined);
  // each card by the plan column it shows, so a change made on the calendar can point at what it moved
  const cardOf: {[column: string]: HTMLElement} = {
    implant_total: card('Animals to implant', count('implant_total'), () => `Covers ${pct(cell('pre_loss').value)} ` +
      'that never reach randomization: no take, or out of the volume window'),
    randomise_total: card('To randomize', count('randomise_total'), () => `Covers ${pct(cell('post_loss').value)} ` +
      'lost after randomization: ulceration, deaths, technical failures'),
    n_total: card('Analyzable', count('n_total'), () => cell('n_fixed').value ? `${groups.value} groups of the ` +
      'animals per arm you set' : `${groups.value} groups, ${pct(cell('power_target').value)} power per comparison`),
    study_end: card('End of study', cell('study_end'), () => cell('last_day_chosen').value ? 'Chosen; the median ' +
      `control reaches the endpoint on day ${cell('endpoint_day').value}` :
      cell('study_end').value === cell('endpoint_day').value ? 'When the median control reaches the endpoint' :
        `The first measurement day after the median control reaches the endpoint (day ${cell('endpoint_day').value})`,
    (x) => `Day ${x}`),
  };
  const cards = cardGrid(Object.values(cardOf), 'pmx-plan-cards', 140);

  // the engine's warnings, each under the section it is about, one line each with the dot on the text's first line
  // (u2 has no banner, so a warning dot and the text)
  const groupNotes = notesBox();
  const growthNotes = notesBox();
  growthNotes.style.margin = 'var(--dg-space-s) 0 var(--dg-space-m)';
  const analysis = div([]);
  let protocolText = '';
  const groupsTable = div([]);
  const scheduleBox = div([]);
  // the power by measurement schedule is the design's costliest scan: the plan takes it along only while its
  // section is open; opened later, it is computed then for the plan in view
  const scheduleOpen = signal(false);
  const schedules = signal<DG.DataFrame | undefined>(undefined);
  const calendarBox = div([]);
  // what a change on the calendar did: the counts and the end it moved (their cards flash), and a power that fell
  // below the target (with the animals per arm fixed; else the animals grow to keep it)
  const calendarLoader = loader('Updating the plan…');
  const calendarStatus = divV([], 'pmx-calendar-status');
  Object.assign(calendarStatus.style, {gap: '4px', margin: '0 0 var(--dg-space-s)', lineHeight: '1.45'});
  const flash = (el: HTMLElement) => {
    Object.assign(el.style, {transition: 'box-shadow 0.4s', boxShadow: '0 0 0 2px var(--dg-accent, #2083d5)'});
    setTimeout(() => el.style.boxShadow = '', 1800);
  };
  const reportChange = (before: PlanFrames, after: PlanFrames) => {
    const was = (c: string) => before.plan.get(c, 0);
    const now = (c: string) => after.plan.get(c, 0);
    const moved: string[] = [];
    for (const [c, what] of [['n_total', 'analyzable animals'], ['implant_total', 'animals to implant']]) {
      if (was(c) !== now(c))
        moved.push(`${finite(was(c)) ?? '—'} → ${finite(now(c)) ?? '—'} ${what}`);
    }
    if (was('study_end') !== now('study_end'))
      moved.push(`the study ends on day ${now('study_end')} instead of ${was('study_end')}`);
    for (const c of ['implant_total', 'randomise_total', 'n_total', 'study_end']) {
      if (was(c) !== now(c))
        flash(cardOf[c]);
    }
    // only what moved: the calendar itself shows the change was planned
    const lines: HTMLElement[] = moved.length > 0 ? [span(`With these days: ${moved.join(', ')}.`)] : [];
    const a = after.armResults;
    const target = now('power_target');
    for (let i = 0; i < a.rowCount; i++) {
      const power = a.get('assurance', i);
      if (typeof power === 'number' && power < target) {
        const line = divH([dot('warning'), span(`${a.get('arm', i)} now has ${pct(power)} power, below the ` +
          `${pct(target)} target: measure more often or add animals.`)]);
        Object.assign(line.style, {alignItems: 'center', gap: 'var(--dg-space-s)'});
        lines.push(line);
      }
    }
    calendarStatus.replaceChildren(...lines);
  };
  // the calendar sets two of the form's own fields, the randomization weekday and the measurement weekdays, and
  // the plan follows the field
  const setField = (name: string, value: string) => {
    const field = form.getInput(name);
    if (field)
      field.value = value;
    call.setParamValue(name, value);
    const before = result.peek();
    void replan(true).then(() => {
      const after = result.peek();
      if (before && after && after !== before)
        reportChange(before, after);
    });
  };
  const startDay = new ChoiceInput({items: WEEKDAY_NAMES, value: 'Monday', nullable: false});
  Object.assign(startDay.root.style, {width: '130px', flex: '0 0 auto'});
  form.effect(() => {
    const s = cell('start_weekday').value;
    if (s !== undefined)
      startDay.value = WEEKDAY_NAMES[s];
  });
  form.effect(() => {
    const chosen = startDay.value;
    // against what was last asked, not the plan on the page: a run may be on its way to another day
    if (chosen && chosen !== (call.getParamValue('startWeekday') ??
      WEEKDAY_NAMES[plan.peek()?.get('start_weekday', 0) ?? 0]))
      setField('startWeekday', chosen);
  });
  const toggleWeekday = (i: number) => {
    const days = new Set(weekdaysOf(call.getParamValue('measurementDays') ??
      plan.peek()?.get('measurement_weekdays', 0)));
    if (days.has(i))
      days.delete(i);
    else
      days.add(i);
    if (days.size > 0)
      setField('measurementDays', weekdaysText([...days].sort()));
  };
  // a dated day 0 sets the randomization weekday; the picker then follows the date
  form.effect(() => {
    const start = meta.value.start;
    startDay.disabled = start ? 'Day 0 has a date: change it from a day\'s menu in the calendar, or in the ' +
      'study\'s details' : null;
    if (start && startDay.peek() !== WEEKDAY_NAMES[weekdayOf(start)])
      startDay.value = WEEKDAY_NAMES[weekdayOf(start)];
  });
  const dayLabel = (day: number) => {
    const d = dateOf(meta.peek().start, day);
    return `day ${day}${d ? ` (${dateText(d)})` : ''}`;
  };
  const addNote = (day: number) => noteDialog(day, dayLabel(day), [],
    (n) => meta.value = {...meta.peek(), notes: [...meta.peek().notes, n]});
  const dateDay = (day: number) => {
    const input = ui.input.date(`Date of day ${day}`);
    ui.dialog(`Date ${dayLabel(day)}`).add(input).onOK(() => {
      const iso = (input.value as unknown as {format(f: string): string} | null)?.format('YYYY-MM-DD');
      if (iso)
        meta.value = {...meta.peek(), start: startFrom(iso, day)};
    }).show();
  };
  const dayMenu = (day: number, at: MenuAt) => {
    const m = meta.peek();
    const menu = DG.Menu.popup().item('Add a note…', () => addNote(day));
    const own = m.notes.filter((n) => n.day === day);
    if (own.length > 0) {
      const remove = menu.group('Remove a note');
      for (const n of own) {
        remove.item(`${n.kind}: ${n.text.length > 40 ? `${n.text.slice(0, 40)}…` : n.text}`,
          () => meta.value = {...meta.peek(), notes: meta.peek().notes.filter((x) => x !== n)});
      }
    }
    menu.item(m.start ? `Change the date of day ${day}…` : `Set the date of day ${day}…`, () => dateDay(day));
    const d = dateOf(m.start, day);
    if (d)
      menu.item('Copy the date', () => void navigator.clipboard.writeText(d.toISOString().slice(0, 10)));
    menu.show(at);
  };
  // the calendar: the plan's days with the study's dates, doses and notes; built in its own scope, so its tooltips
  // go with the next build
  form.effect(() => {
    const r = result.value;
    const m = meta.value;
    if (r) {
      calendarBox.replaceChildren(calendarGrid(r.plan, r.calendar, startDay.root, toggleWeekday,
        {meta: m, onMenu: dayMenu, onNote: addNote}));
    }
  });
  const detectable =signal<{of: PlanFrames, mde: Map<string, number | null>} | null>(null);
  const endSummary = div([], 'pmx-plan-end-summary');
  Object.assign(endSummary.style, {margin: '0 0 var(--dg-space-s)', lineHeight: '1.45'});
  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    const v = (c: string) => r.plan.get(c, 0);
    groupNotes.replaceChildren(...notes(v('notes_groups')));
    growthNotes.replaceChildren(...notes(v('notes_growth')));

    // the control first, then each arm: what it expects, the animals it gets and the power they give (a table built
    // in this run is disposed with it)
    const a = r.armResults;
    const colors = armColors(a.col('arm')!.toList());
    type Row = {name: string, cells: (string | HTMLElement)[]};
    const fixed = v('n_fixed');
    // the smallest detectable effects arrive after the plan (their search takes several times as long)
    const known = detectable.value?.of === r ? detectable.value.mde : null;
    const rows: Row[] = [{name: 'Control', cells: [swatch('Control', colors['Control']), '—', '—',
      animals(finite(v('n_control')) ?? 'over 40', true), animals(finite(v('randomise_control')), false), '—',
      ...(fixed ? ['—'] : [])]}];
    for (let i = 0; i < a.rowCount; i++) {
      const g = (c: string) => a.get(c, i);
      const metric = g('metric') === 'TC' ? 'T/C' : 'TGI';
      rows.push({name: g('arm'), cells: [swatch(g('arm'), colors[g('arm')]), g('therapy'),
        `${metric} ${num(g('effect'))}%`, animals(finite(g('n')) ?? 'over 40', true),
        animals(finite(g('randomise')), false),
        pct(g('assurance')),
        ...(fixed ? [!known ? (stale.value || error.value !== null ? '—' : loader()) :
          known.get(g('arm')) === null ? 'out of reach' : `${metric} ${Math.round(known.get(g('arm'))!)}%`] : [])]});
    }
    const headers = ['Group', 'Therapy', 'Expected effect', 'Analyzable', 'To randomize', 'Power',
      ...(fixed ? ['Smallest detectable effect'] : [])];
    const table = new BasicTable<Row>({items: rows, keyOf: (x) => x.name, columns: headers.map((header, j) => ({header,
      render: (x: Row) => x.cells[j], align: j >= 2 ? 'center' : 'left'}))});
    table.root.style.width = '100%';
    wrapCells(table.root);
    groupsTable.replaceChildren(table.root);
  });
  form.effect(() => {
    const sc = schedules.value;
    // (nothing comes for a plan out of date or a failed run: no loader then)
    const waiting = stale.value || error.value !== null ? span('Compared once the plan is up to date.') :
      loader('Comparing measurement schedules…');
    scheduleBox.replaceChildren(!sc ? waiting : sc.rowCount > 0 ?
      schedulesTable(sc, plan.peek()?.get('per_week', 0)) : span(UNREACHABLE));
  });

  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    const v = (c: string) => r.plan.get(c, 0);
    // what ending later than the plan's end buys and costs, read off the power-by-end curves
    const du = r.duration;
    endSummary.textContent = du.rowCount > 0 ? '' : UNREACHABLE;
    if (du.rowCount > 0 && du.col('assurance')!.stats.missingValueCount < du.rowCount) {
      const e = v('study_end');
      // (a day with no power left out)
      const rowsOf = Array.from({length: du.rowCount}, (_, i) => ({arm: du.get('arm', i) as string,
        day: du.get('day', i) as number, power: finite(du.get('assurance', i)), removed: du.get('past_humane', i)}))
        .filter((x): x is typeof x & {power: number} => x.power !== null);
      const gain = Math.max(...[...new Set(rowsOf.map((x) => x.arm))].map((arm) => {
        const own = rowsOf.filter((x) => x.arm === arm);
        const atE = own.find((x) => x.day === e)?.power ?? 0;
        return Math.max(...own.filter((x) => x.day >= e).map((x) => x.power)) - atE;
      }));
      const last = rowsOf[rowsOf.length - 1];
      const removedAt = (day: number) => Math.round(100 * (rowsOf.find((x) => x.day === day)?.removed ?? 0));
      const points = Math.round(100 * gain);
      endSummary.textContent = `Ending after day ${e} ${points < 1 ? 'adds no power' : `adds at most ${points} ` +
        `points of power`}, while the controls removed at ${v('humane_limit')} mm³ go from ${removedAt(e)}% on day ` +
        `${e} to ${removedAt(last.day)}% on day ${last.day}.`;
    }


    // the study protocol, in the order of the ARRIVE essentials: design, sample size and its basis, animals,
    // randomization, outcome and analysis, multiplicity, end, humane limit, measurements
    const k = v('treated_arms');
    const corr = v('correction');
    const ar = r.armResults;
    const armRows = Array.from({length: ar.rowCount}, (_, i) => ({name: ar.get('arm', i) as string,
      therapy: String(ar.get('therapy', i)).toLowerCase(), effect: `${ar.get('metric', i) === 'TC' ? 'T/C' : 'TGI'} ` +
        `${num(ar.get('effect', i))}%`, n: finite(ar.get('n', i)), power: ar.get('assurance', i)}));
    const plan: {[row: string]: string} = {
      'Design': `${v('line')} syngeneic tumors; ${k + 1} groups: a control and ` + armRows.map((x) =>
        `${x.name} (${x.therapy}, expected ${x.effect})`).join(', ') + '; each treated group compared with the control',
      'Sample size': armRows.map((x) => x.n === null ? `${x.name}: more than 40 analyzable animals per group ` +
        `would be needed for ${x.effect}` : `${x.name}: ${x.n} analyzable animals per group give ${pct(x.power)} ` +
        `power for ${x.effect}`).join('; ') + `, two-sided α = ${num(v('alpha'), 2)}. ` +
        (v('n_fixed') ? 'The animals per group were set; the power is what they give' :
          `Sized for ${pct(v('power_target'))} power`) + `, averaged over the growth of ${v('line')} controls in ` +
        `${v('atlas_papers')} published studies (historical-control prior); at least ${MIN_GROUP} per group`,
      'Animals': (finite(v('n_total')) === null ? 'Not sized: no group size up to 40 reaches the target power. ' :
        `${v('implant_total')} to implant, ${v('randomise_total')} to randomize, ${v('n_total')} analyzable. `) +
        `${pct(v('pre_loss'))} expected lost before randomization (no take, out of the volume ` +
        `window) and ${pct(v('post_loss'))} after it (ulceration, deaths, technical failures); ` +
        (v('loss_rule') === 'expected' ? 'the counts cover the average loss' :
          `each count is enough with ${pct(v('confidence'))} probability`) +
        (k < 2 ? '' : '. The control is as large as the largest treated group'),
      'Randomization': `When the tumors reach ${v('v_rand')} mm³, into groups balanced on tumor volume; the ` +
        'randomization day is day 0',
      'Outcome and analysis': 'The growth rate of each animal’s tumor (slope of log volume over its ' +
        'measurements), compared between each treated group and the control by a two-sample t test',
      'Multiple comparisons': k < 2 ? 'One treated group: none needed' : corr === 'dunnett' ?
        `Dunnett's many-to-one test over ${k} comparisons (each at α = ${num(v('alpha_per_comparison'), 3)} for ` +
        'the largest group)' :
        corr === 'bonferroni' ? `Bonferroni over ${k} comparisons (each at α = ${num(v('alpha_per_comparison'), 3)})` :
          `None: each of the ${k} comparisons at α = ${num(v('alpha'), 2)} (exploratory)`,
      'End of study': v('last_day_chosen') ? `Day ${v('study_end')}, chosen. The expected effects are the ones on ` +
        `day ${v('endpoint_day')}, when the median control reaches ${v('endpoint')} mm³` :
        `The first measurement day once the median control reaches ${v('endpoint')} mm³ (expected day ` +
        `${v('study_end')})`,
      'Humane limit': `${v('humane_limit')} mm³: an animal above it is removed, whatever the study day; its ` +
        'measurements up to removal stay in the analysis',
      'Measurements': `Tumor volume and body weight on ${v('measurement_weekdays')}: days ` +
        `${v('measurement_days')} from randomization`,
    };
    analysis.replaceChildren(tableFromMap(plan));
    // the text the Copy button gives: the same rows
    protocolText = [`Study protocol (${r.plan.name})`, ...Object.entries(plan).map(([k, t]) => `${k}: ${t}`)]
      .join('\n');
  });

  const {row: scaleRow, log} = volumeAxis();
  const growth = new BandChart(computed(() => result.value?.curves), {x: 'day', lower: 'q05', upper: 'q95',
    line: 'median', split: 'arm', xTitle: 'Day after randomization', yTitle: 'Tumor volume', log});


  // the assurance of every arm on one value column split by arm (several Y columns would each get an axis), each
  // arm's size in its legend entry, the target a bare formula line
  const power = viewers.lineChart(computed(() => {
    const r = result.value;
    if (!r)
      return undefined;
    // each arm's size in its legend entry: it reads level and hides with its curve
    const pc = r.powerCurve;
    const a = r.armResults;
    const label = new Map(Array.from({length: a.rowCount}, (_, i) => [a.get('arm', i) as string,
      `${a.get('arm', i)} · ${finite(a.get('n', i)) ?? 'over 40'} animals`]));
    // on a power axis the declared target reads by itself: a bare formula line (a one-space title, so that neither a
    // label nor the formula is drawn)
    const target = r.plan.get('power_target', 0);
    const df = DG.DataFrame.fromObjects(Array.from({length: pc.rowCount}, (_, i) => i)
      .filter((i) => finite(pc.get('assurance', i)) !== null).map((i) => ({'Analyzable animals': pc.get('n', i),
        'Power': 100 * pc.get('assurance', i), 'Arm': label.get(pc.get('arm', i))})))!;
    df.col('Power')!.meta.units = '%';
    df.col('Power')!.meta.format = '0.0';
    const colors = armColors(a.col('arm')!.toList());
    df.col('Arm')!.meta.colors.setCategorical(Object.fromEntries([...label].map(([arm, l]) => [l, colors[arm]])));
    df.meta.formulaLines.addLine({title: ' ', formula: `\${Power} = ${target * 100}`, color: '#9aa3ad'});
    return df;
  }), {xColumnName: 'Analyzable animals', yColumnNames: ['Power'], splitColumnNames: ['Arm'],
    legendVisibility: 'Always', showXAxis: true, showYAxis: true,
    // a result, not an exploration: the X and Y pickers stay (they are the line chart's axis titles), the split
    // and aggregation ones go
    showTopPanel: false, showAggrSelectors: false, showSplitSelector: false});

  // each arm's power and the controls' state by the day the study ends, on one percent axis (series by legend, the
  // target one of them); the endpoint day and a later end are the frame's formula lines
  const byEnd = viewers.lineChart(computed(() => {
    const r = result.value;
    if (!r || r.duration.rowCount === 0)
      return undefined;
    const t = r.duration;
    const colors = armColors(r.armResults.col('arm')!.toList());
    const [ENDPOINT, HUMANE] = [`Controls over the endpoint (${r.plan.get('endpoint', 0)} mm³)`,
      `Controls removed (over ${r.plan.get('humane_limit', 0)} mm³)`];
    // one curve after another (each arm's power, then the controls' two shares), one decimal in the tooltips
    const rows: {[k: string]: unknown}[] = [];
    const point = (day: number, value: number, curve: string) =>
      rows.push({'Last day of study': day, 'Percent': 100 * value, 'Curve': curve, '~marked': true});
    for (let i = 0; i < t.rowCount; i++) {
      if (finite(t.get('assurance', i)) !== null)
        point(t.get('day', i), t.get('assurance', i), `Power, ${t.get('arm', i)}`);
    }
    const powers = [...new Set(t.col('arm')!.toList() as string[])].map((arm) => `Power, ${arm}`);
    const first = t.get('arm', 0);
    for (const [col, curve] of [['past_endpoint', ENDPOINT], ['past_humane', HUMANE]]) {
      for (let i = 0; i < t.rowCount; i++) {
        if (t.get('arm', i) === first)
          point(t.get('day', i), t.get(col, i), curve);
      }
    }
    // the target as a series of its own, named in the legend (a formula line's label floats away from its line)
    const [e, end, target] = [r.plan.get('endpoint_day', 0), r.plan.get('study_end', 0), r.plan.get('power_target', 0)];
    const TARGET = `Target ${pct(target)}`;
    for (const day of new Set(rows.map((x) => x['Last day of study'] as number)))
      rows.push({'Last day of study': day, 'Percent': 100 * target, 'Curve': TARGET, '~marked': false});
    const df = DG.DataFrame.fromObjects(rows)!;
    df.col('Percent')!.meta.format = '0.0';
    // the legend in the chart's reading order: the arms' power on top, the target, then the controls
    df.col('Curve')!.setCategoryOrder([...powers, TARGET, ENDPOINT, HUMANE]);
    df.col('Curve')!.meta.colors.setCategorical({...Object.fromEntries(Object.entries(colors)
      .map(([arm, c]) => [`Power, ${arm}`, c])), [TARGET]: '#9aa3ad', [ENDPOINT]: colors['Control'],
    [HUMANE]: '#c0504d'});
    df.meta.formulaLines.addLine({title: `Endpoint, day ${e}`, formula: `\${Last day of study} = ${e}`});
    if (end !== e) {
      df.meta.formulaLines.addLine({title: `End of study, day ${end}`, formula: `\${Last day of study} = ${end}`,
        style: 'solid', color: '#4d5261'});
    }
    return df;
  }), {xColumnName: 'Last day of study', yColumnNames: ['Percent'], splitColumnNames: ['Curve'],
    markersVisibilityColumnName: '~marked',
    legendVisibility: 'Always', showXAxis: true, showYAxis: true,
    showTopPanel: false, showAggrSelectors: false, showSplitSelector: false});

  const scenarioEnd = div([], 'pmx-plan-scenario-end');
  const scenarioBox = div([]);
  const noScenarios = infoNote('<b>No TGIs to assume.</b> Enter them under Explore in the form.');
  form.effect(() => {
    const has = !!result.value?.scenarios.rowCount;
    noScenarios.hidden = has;
    scenarioEnd.hidden = !has;
    scenarioBox.hidden = !has;
  });
  // the figures answer the plan's question; each panel, one closer question
  const panels = panelGrid([
    {id: 'groups', title: 'Experimental groups', hint: 'The animals each group needs for its comparison with the ' +
      'control', content: [groupsTable, groupNotes, monitorNote()]},
    {id: 'growth', title: 'Predicted tumor growth', hint: 'Median and 90% range of each group\u2019s tumor volume, ' +
      'to the end of the study', content: [growthNotes, scaleRow, sized(growth, 380)]},
    {id: 'power', title: 'Power by group size', hint: 'With more or fewer animals', content: [sized(power, 380)]},
    {id: 'calendar', title: 'Study calendar', hint: 'Week by week from randomization. Click a weekday to measure ' +
      'on it or not; double-click a day to note what happened, right-click it for more', content: [calendarLoader,
      calendarStatus, calendarBox]},
    {id: 'end', title: 'Power by end of study', hint: 'Ending earlier or later', content: [endSummary,
      sized(byEnd, 280)]},
    {id: 'schedule', title: 'Measurement schedule', hint: 'Measuring more or less often', content: [scheduleBox]},
    {id: 'scenarios', title: 'Assumed TGI', hint: 'What the study would need if the treatment had each TGI: ' +
      'scenarios you read across, not predictions. The animals of the whole study with every arm at that TGI, ' +
      'sized for the target power whatever animals per arm is set', content: [noScenarios, scenarioEnd, scenarioBox]},
    {id: 'protocol', title: 'Study protocol', hint: 'Design, sample size and analysis, ready to paste',
      content: [copyRow(() => protocolText, 'Study protocol'), analysis]},
  ], PLAN_LAYOUT, 'pmx-plan-panels');
  // an empty viewer says what columns it lacks: the charts hidden until the first plan gives them their data
  form.effect(() => [growth, power].forEach((v) => v.root.style.visibility = result.value ? '' : 'hidden'));
  // (with no group size, no power by end: the chart would keep the last plan's)
  form.effect(() => byEnd.root.style.visibility = result.value?.duration.rowCount ? '' : 'hidden');
  // the schedules are the design's costliest scan: computed once their panel has been in view (and from then on,
  // so scrolling does not flip the form's input and dim the plan)
  const seen = new IntersectionObserver((entries) => {
    if (entries.some((e) => e.isIntersecting)) {
      scheduleOpen.value = true;
      seen.disconnect();
    }
  });
  seen.observe(scheduleBox);

  form.effect(() => {
    const r = result.value;
    if (!r?.scenarios.rowCount)
      return;
    scenarioEnd.textContent = `The study ends on day ${r.scenarios.get('study_end', 0)} in every scenario.`;
    scenarioBox.replaceChildren(scenariosTable(r.scenarios, r.plan.get('per_week', 0)));
  });
  const content = divV([cards, panels], 'pmx-plan-content');
  Object.assign(content.style, {gap: 'var(--dg-space-l)', minWidth: `${RESULTS_MIN}px`});

  // the first plan under u2's loader; every later one under the thin bar, what was shown left in view
  const overlay = busyOverlay('Designing the study…');
  const bar = busyBar();
  const failure = div([], 'pmx-plan-error');
  Object.assign(failure.style, {color: 'var(--dg-error-color, #c62828)'});
  form.effect(() => {
    overlay.hidden = !planning.value || result.value !== undefined;
    bar.hidden = !planning.value || result.value === undefined;
  });
  // a calendar change: its own loader, the calendar dimmed, the page in view; a plan from the form clears its note
  form.effect(() => {
    calendarLoader.hidden = !refreshing.value;
    calendarBox.style.opacity = refreshing.value ? '0.5' : '';
    if (planning.value)
      calendarStatus.replaceChildren();
  });
  form.effect(() => {
    failure.hidden = error.value === null;
    failure.replaceChildren(...(error.value === null ? [] : [divH([dot('error'), span(error.value)])]));
  });

  // the view follows the form by itself: what it shows is behind only while a field is missing or invalid
  const outdated = staleNote('A field is missing or invalid: the plan shown is for the last complete inputs');
  form.effect(() => {
    outdated.hidden = !stale.value || form.missingRequired.value.length + form.invalidFields.value.length === 0;
    content.style.opacity = stale.value ? '0.5' : '';
  });

  const results = resultsFrame([failure, outdated, content], overlay, bar);
  const declaration = declarationPane(form.root);

  // the plan carries the schedules while their section is open; opened on a plan without them, they are computed
  // for it (the form's own inputs, the plan unchanged)
  void form.settled.then(() => form.effect(() => call.setParamValue('compareSchedules', scheduleOpen.value)));
  form.effect(() => {
    const r = result.value;
    schedules.value = r && r.schedules.rowCount > 0 ? r.schedules : undefined;
  });
  form.effect(() => {
    const r = result.value;
    if (!scheduleOpen.value || schedules.value || !r || stale.value || busy.value)
      return;
    const fill = designed({compareSchedules: true, tgiScenarios: null});
    void fill.then((c) => {
      if (result.peek() === r)
        schedules.value = c.schedules;
    }).catch((e) => error.value = String(e?.message ?? e));
  });

  // with the animals per arm set, the plan shows first and each arm's smallest detectable effect fills in after it,
  // from the form's own inputs (the plan unchanged)
  let searching: PlanFrames | null = null;
  form.effect(() => {
    const r = result.value;
    if (!r || !r.plan.get('n_fixed', 0) || stale.value || busy.value || searching === r || detectable.peek()?.of === r)
      return;
    searching = r;
    const fill = designed({compareSchedules: false, detectableEffect: true, tgiScenarios: null});
    void fill.then((c) => {
      const a = c.armResults;
      detectable.value = {of: r, mde: new Map(Array.from({length: a.rowCount}, (_, i) =>
        [a.get('arm', i) as string, finite(a.get('mde', i))]))};
    }).catch((e) => error.value = String(e?.message ?? e)).finally(() => {
      if (searching === r)
        searching = null;
    });
  });

  // the view opens at once; the default declaration (or the saved plan's) is planned behind the loader, once the
  // form has written its defaults into the call
  void form.settled.then(() => {
    call.setParamValue('detectableEffect', false);
    // a new plan opens on the tumor model picked last
    const last = initial ? null : lastLine();
    for (const [name, value] of Object.entries(initial?.inputs ?? (last ? {line: last} : {}))) {
      const v = name === 'arms' ? armsFrame(value as SavedPlan['inputs']['arms']) : value;
      const field = form.getInput(name);
      if (field)
        field.value = v;
      call.setParamValue(name, v);
    }
    auto.cancel();
    return replan();
  });

  // the plan in view as a file: saved to the user's plans (shared from Files like any file) or downloaded, to send
  // (the plan's predictive checks are computed then, not on every replan: planOf)
  const ready = (): PlanFrames | null => {
    const r = result.peek();
    if (!r || stale.peek()) {
      grok.shell.warning(r ? 'The plan does not match the form yet: complete every field, or wait for it to update' :
        'The plan is not ready yet');
      return null;
    }
    return r;
  };
  const defaultName = (r: PlanFrames) => initial?.name ?? `${r.plan.get('line', 0)} plan`;
  const stamped = () => meta.value = withId(meta.peek());
  const planFor = (name: string, r: PlanFrames) => planOf(name, call, r, stamped());
  const saveDialog = () => {
    const r = ready();
    if (!r)
      return;
    const where = span(`Saved to ${plansFolder()}. Share it from Files, or pick it in Analyze's Plan field.`);
    Object.assign(where.style, {color: 'var(--dg-text-color-muted, #68727f)', maxWidth: '360px', lineHeight: '1.45'});
    // a plan of the same name is replaced only when Save is pressed again for it (the dialog says so and stays open)
    let replacing: string | null = null;
    self.show(self.build(() => {
      const name = new TextInput({label: 'Name', value: defaultName(r), nullable: false});
      return new Dialog('Save plan').add(name).add(where).onOK(async () => {
        const chosen = name.value.trim() || defaultName(r);
        if (chosen !== replacing && await planSaved(chosen)) {
          replacing = chosen;
          throw new Error(`A plan named "${chosen}" is already saved: press Save again to replace it, or rename it`);
        }
        const path = await savePlan(await planFor(chosen, r));
        grok.shell.info(`Plan saved: ${path}`);
      }, {text: 'Save'});
    }));
  };
  const saveButton = buttonWithIcon('Save plan', 'save', saveDialog);
  saveButton.title = 'Save the plan to your files, to compare the study with it in Analyze';
  const saving = signal(false);
  const downloadButton = buttonWithIcon('Download', 'download', () => {
    const r = ready();
    if (!r)
      return;
    saving.value = true;
    void planFor(defaultName(r), r).then(downloadPlan).catch((e) => error.value = String(e?.message ?? e))
      .finally(() => saving.value = false);
  });
  downloadButton.title = 'Download the plan as a file, to send it to someone';
  const followButton = buttonWithIcon('Follow this study', 'eye', () => {
    const r = ready();
    if (!r)
      return;
    saving.value = true;
    void planFor(defaultName(r), r).then(followPlan).catch((e) => error.value = String(e?.message ?? e))
      .finally(() => saving.value = false);
  });
  followButton.title = 'Start following this study in Monitor, with this plan';
  form.effect(() => {
    const off = !result.value ? 'The plan is not ready yet' : busy.value ? 'The plan is being updated' :
      saving.value ? 'The plan is being prepared' : null;
    [saveButton, downloadButton, followButton].forEach((b) => Reason.apply(b, off));
  });

  const panes = new Splitter([declaration, results], {direction: 'horizontal', sizes: [0.3, 0.7],
    minSize: DECLARATION_MIN});
  return {
    root: panes,
    ribbon: [[studyHeader(form, 'plan', computed(() => {
      const r = result.value;
      if (!r)
        return undefined;
      const arms = r.armResults.col('arm')!.toList() as string[];
      const colors = armColors(arms);
      return {name: defaultName(r), groups: ['Control', ...arms].map((name) => ({name, color: colors[name]})),
        status: `${r.plan.get('line', 0)} · planned`, meta: meta.value, setMeta: (m) => meta.value = m};
    }), (s) => void openStage(s))], [saveButton, downloadButton, followButton]],
    path: '/plan',
    status: computed(() => busy.value ? 'Designing the study…' : result.value ? result.value.plan.name : ''),
  };
}
