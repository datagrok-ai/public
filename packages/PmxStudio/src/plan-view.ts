/* The Plan view: the declaration form is generated from pmxPlanStudy's own signature (u2 funcForm, grouped by
   the parameters' categories), so the function Grokky reads and the form people fill in are one definition; the
   treated arms and the multiplicity correction have editors of their own (arms-input.ts), and the measurement and
   randomization weekdays are set on the study calendar. The results are the function's outputs: the study's
   figures, the groups, the predicted volumes (band-chart.ts), the study calendar, then, folded, the power by
   group size, by end of study (platform line charts) and by measurement schedule, and the study protocol. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {BasicTable, buttonWithIcon, ChoiceInput, computed, Dialog, div, divH, divV, dot, loader, signal, span,
  Splitter, tableFromMap, TextInput, Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal, StatValue} from '@datagrok-libraries/u2';
import {appView, funcForm, viewers} from '@datagrok-libraries/u2/src/dg/index.js';

import {PlanFrames, WEEKDAY_NAMES, armColors, weekdaysOf, weekdaysText} from './frames';
import {ArmsInput, CorrectionInput, armsFrame} from './arms-input';
import {MIN_GROUP, finite} from './engine';
import {block, busyOverlay, cardGrid, copyRow, declarationPane, DECLARATION_MIN, formRunning, infoNote, notes,
  notesBox, num, outputsOf, pct, RESULTS_MIN, resultsFrame, resultTabs, settleForm, sized, staleNote, stageStrip,
  statCard, swatch, volumeAxis} from './ui';
import {openStage} from './stages';
import {SavedPlan, downloadPlan, parsePlan, planOf, planSaved, plansFolder, savePlan} from './plan-file';
import {BandChart} from './band-chart';

/** Where the plan's first stage is run: Monitor takes the implanted pool to randomization. */
const monitorNote = () => infoNote('<b>From implantation to randomization, in Monitor.</b> Monitor follows every ' +
  'implanted animal until its tumor enters the randomization window, flags the ones that do not take, and forms ' +
  'the groups automatically, balanced on tumor volume, as soon as there are enough.');

const WEEKDAYS = WEEKDAY_NAMES.map((d) => d.slice(0, 3));

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
 *  (`startDay`, its picker) and what the weeks add up to. Built from u2's primitives (u2 has no calendar); call it
 *  in a scope, which owns the tooltips. */
function calendarGrid(p: DG.DataFrame, calendar: DG.DataFrame, startDay: HTMLElement,
  toggle: (weekday: number) => void): HTMLElement {
  const v = (c: string) => p.get(c, 0);
  const end = v('study_end') as number;
  const endpointDay = v('endpoint_day') as number;
  const humane = v('humane_limit');
  const watch = v('humane_watch_day') as number | null;
  const details = new Map<number, string>(Array.from({length: calendar.rowCount},
    (_, i) => [calendar.get('day', i), calendar.get('activity', i)]));
  const measured = new Set(String(v('measurement_days')).split(', ').map(Number));
  const weekdays = new Set(weekdaysOf(v('measurement_weekdays')));
  const [BLUE, SLATE, RED, DARK] = ['#2083d5', '#5b6b7f', '#c0504d', '#4d5261'];

  // day 0 falls on the randomization weekday: the measurement weekdays are the lab's, seen from it
  const start = (v('start_weekday') ?? 0) as number;
  const weekday = (day: number) => (start + day) % 7;

  // above the weeks: the randomization weekday (picked in the sentence), the end and the closer watch
  const q90 = v('duration_q90');
  const range = typeof q90 === 'number' && Number.isFinite(q90) ? `day ${Math.round(v('duration_q10'))} to ` +
    `${Math.round(q90)}` : `from day ${Math.round(v('duration_q10'))}`;
  const first = divH([mark('Randomize', BLUE, true), span('Randomize on a'), startDay,
    span(`(day 0) · ${measured.size} measurement days.`)]);
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
      const num = span(`Day ${day}`);
      Object.assign(num.style, {fontSize: '11px', color: 'var(--dg-text-color-muted, #68727f)'});
      cell.append(num);
      if (day === 0)
        cell.append(mark('Randomize', BLUE, true));
      if (measured.has(day))
        cell.append(mark('Measure', SLATE, false));
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
      if (detail)
        Tooltip.bind(cell, `Day ${day}: ${detail}`);
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
    const row: Schedule[1] = byWeek.get(k) ?? {measurements: sc.get('measurements', i), nTotal: sc.get('n_total', i),
      power: {}};
    row.power[sc.get('arm', i)] = sc.get('assurance', i);
    byWeek.set(k, row);
  }
  const armNames = [...new Set(sc.col('arm')!.toList() as string[])];
  const columns = [
    {header: 'Measurements a week', render: ([k]: Schedule) => k === own ? `${k} · this plan` : String(k)},
    {header: 'Measurement days', render: ([, s]: Schedule) => String(s.measurements)},
    {header: 'Analyzable animals', render: ([k, s]: Schedule) => animals(s.nTotal ?? 'over 40', k === own)},
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

/** A .pmxplan file opened from Files or dropped on the platform: Plan with it (no table to open). */
export async function openPlanFile(content: string): Promise<DG.DataFrame[]> {
  const view = grok.shell.addView(await planView(parsePlan(content)));
  // as stages.ts opens a stage: the app's root before the view's path
  view.path = '/apps/PmxStudio/plan';
  return [];
}

/** The Plan view; given a saved plan, its declaration in the form and planned again. */
export async function planView(initial?: SavedPlan): Promise<DG.ViewBase> {
  const call = DG.Func.byName('PmxStudio:pmxPlanStudy').prepare();
  const read = (c: DG.FuncCall) => outputsOf<PlanFrames>(c);
  const result = signal<PlanFrames | undefined>(undefined);
  // the newest run: the form's Run and the view's own replans share the call, and a run that finishes after a newer
  // one started is dropped
  let latest = 0;
  const error = signal<string | null>(null);
  const armsInput = new ArmsInput({label: 'Arms'});
  const correction = new CorrectionInput({label: 'Multiple comparisons', arms: armsInput.count});
  // a plan that no longer matches the form is dimmed until the next run
  const stale = signal(false);
  // captions above every field: the panel reads top to bottom, left-aligned, at any width
  const form = funcForm(call, {showRun: true, layout: 'tall',
    overrides: {arms: {input: armsInput}, correction: {input: correction}},
    onInputChanged: () => stale.value = result.peek() !== undefined,
    onRun: (c) => {
      latest++;
      error.value = null;
      stale.value = false;
      result.value = read(c as DG.FuncCall);
    }});
  settleForm(form, {n: 'Auto (target power)', lastDay: 'Auto (endpoint reached)'});

  // busy while the view plans (the first plan, a change of randomization day: not saved as runs) and while the
  // form runs one (funcForm says so in aria-busy)
  const planning = signal(true);
  const running = formRunning(form);
  // a change made on the calendar replans the study too (the measurement days set the power and the animals), but
  // under a loader of the calendar's own, the page left in view
  const refreshing = signal(false);
  const blocking = computed(() => planning.value || running.value);
  const busy = computed(() => blocking.value || refreshing.value);
  const replan = (fromCalendar = false) => {
    (fromCalendar ? refreshing : planning).value = true;
    const mine = ++latest;
    return call.call().then((c) => {
      if (mine !== latest)
        return;
      error.value = null;
      stale.value = false;
      result.value = read(c);
    }).catch((e) => error.value = String(e?.message ?? e))
      .finally(() => (fromCalendar ? refreshing : planning).value = false);
  };

  const plan = computed(() => result.value?.plan);
  const cell = <T = number>(name: string) => computed(() => (plan.value?.get(name, 0) ?? undefined) as T | undefined);
  const groups = computed(() => (cell('treated_arms').value ?? 1) + 1);
  // the animals in the order the study uses them, then its end; what each count covers in its tooltip
  const card = (label: string, value: ReadonlySignal<number | undefined>, tip: () => string,
    format?: (x: StatValue) => string) => statCard(form, label, value, tip, format);
  // each card by the plan column it shows, so a change made on the calendar can point at what it moved
  const cardOf: {[column: string]: HTMLElement} = {
    implant_total: card('Animals to implant', cell('implant_total'), () => `Covers ${pct(cell('pre_loss').value)} ` +
      'that never reach randomization: no take, or out of the volume window'),
    randomise_total: card('To randomize', cell('randomise_total'), () => `Covers ${pct(cell('post_loss').value)} ` +
      'lost after randomization: ulceration, deaths, technical failures'),
    n_total: card('Analyzable', cell('n_total'), () => cell('n_fixed').value ? `${groups.value} groups of the ` +
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
        moved.push(`${was(c) ?? '—'} → ${now(c) ?? '—'} ${what}`);
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
      field.value.value = value;
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
      startDay.value.value = WEEKDAY_NAMES[s];
  });
  form.effect(() => {
    const chosen = startDay.value.value;
    if (chosen && chosen !== WEEKDAY_NAMES[plan.peek()?.get('start_weekday', 0) ?? 0] && !busy.peek())
      setField('startWeekday', chosen);
  });
  const toggleWeekday = (i: number) => {
    if (busy.peek())
      return;
    const days = new Set(weekdaysOf(plan.peek()?.get('measurement_weekdays', 0)));
    if (days.has(i))
      days.delete(i);
    else
      days.add(i);
    if (days.size > 0)
      setField('measurementDays', weekdaysText([...days].sort()));
  };
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
    const rows: Row[] = [{name: 'Control', cells: [swatch('Control', colors['Control']), '—', '—',
      animals(v('n_control'), true), animals(v('randomise_control'), false), '—', ...(fixed ? ['—'] : [])]}];
    for (let i = 0; i < a.rowCount; i++) {
      const g = (c: string) => a.get(c, i);
      const metric = g('metric') === 'TC' ? 'T/C' : 'TGI';
      rows.push({name: g('arm'), cells: [swatch(g('arm'), colors[g('arm')]), g('therapy'),
        `${metric} ${num(g('effect'))}%`, animals(g('n') ?? 'over 40', true), animals(g('randomise'), false),
        pct(g('assurance')),
        ...(fixed ? [g('mde') === null ? 'out of reach' : `${metric} ${Math.round(g('mde'))}%`] : [])]});
    }
    const headers = ['Group', 'Therapy', 'Expected effect', 'Analyzable', 'To randomize', 'Power',
      ...(fixed ? ['Smallest detectable effect'] : [])];
    const table = new BasicTable<Row>({items: rows, keyOf: (x) => x.name, columns: headers.map((header, j) => ({header,
      render: (x: Row) => x.cells[j], align: j >= 2 ? 'center' : 'left'}))});
    table.root.style.width = '100%';
    groupsTable.replaceChildren(table.root);
  });
  form.effect(() => {
    const sc = schedules.value;
    scheduleBox.replaceChildren(sc ? schedulesTable(sc, plan.peek()?.get('per_week', 0)) :
      loader('Comparing measurement schedules…'));
  });

  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    const v = (c: string) => r.plan.get(c, 0);
    // what ending later than the plan's end buys and costs, read off the power-by-end curves
    const du = r.duration;
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

    // built in this run's own scope, so its tooltips go with the next run
    calendarBox.replaceChildren(calendarGrid(r.plan, r.calendar, startDay.root, toggleWeekday));

    // the study protocol, in the order of the ARRIVE essentials: design, sample size and its basis, animals,
    // randomization, outcome and analysis, multiplicity, end, humane limit, measurements
    const k = v('treated_arms');
    const corr = v('correction');
    const ar = r.armResults;
    const armRows = Array.from({length: ar.rowCount}, (_, i) => ({name: ar.get('arm', i) as string,
      therapy: String(ar.get('therapy', i)).toLowerCase(), effect: `${ar.get('metric', i) === 'TC' ? 'T/C' : 'TGI'} ` +
        `${num(ar.get('effect', i))}%`, n: ar.get('n', i), power: ar.get('assurance', i)}));
    const plan: {[row: string]: string} = {
      'Design': `${v('line')} syngeneic tumors; ${k + 1} groups: a control and ` + armRows.map((x) =>
        `${x.name} (${x.therapy}, expected ${x.effect})`).join(', ') + '; each treated group compared with the control',
      'Sample size': armRows.map((x) => x.n === null ? `${x.name}: more than 40 analyzable animals per group ` +
        `would be needed for ${x.effect}` : `${x.name}: ${x.n} analyzable animals per group give ${pct(x.power)} ` +
        `power for ${x.effect}`).join('; ') + `, two-sided α = ${num(v('alpha'), 2)}. ` +
        (v('n_fixed') ? 'The animals per group were set; the power is what they give' :
          `Sized for ${pct(v('power_target'))} power`) + `, averaged over the growth of ${v('line')} controls in ` +
        `${v('atlas_papers')} published studies (historical-control prior); at least ${MIN_GROUP} per group`,
      'Animals': `${v('implant_total')} to implant, ${v('randomise_total')} to randomize, ${v('n_total')} ` +
        `analyzable. ${pct(v('pre_loss'))} expected lost before randomization (no take, out of the volume ` +
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

  const {input: scale, row: scaleRow, log} = volumeAxis();
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
      `${a.get('arm', i)} · ${a.get('n', i) ?? 'over 40'} animals`]));
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

  // the figures answer the plan's question; each tab, one closer question
  const tabs = resultTabs('plan', [
    {id: 'design', label: 'Design', intro: 'The animals each group needs and how its tumors are expected to grow.',
      content: [
        block('Experimental groups', 'The animals each group needs for its comparison with the control', groupsTable,
          groupNotes, monitorNote()),
        block('Predicted tumor growth', 'Median and 90% range of each group’s tumor volume, to the end of the study',
          growthNotes, scaleRow, sized(growth, 380))]},
    {id: 'calendar', label: 'Calendar', intro: 'When to randomize, measure and end the study, week by week.', content: [
      block('Study calendar', 'Week by week from randomization. Click a weekday to measure on it or not',
        calendarLoader, calendarStatus, calendarBox)]},
    {id: 'alternatives', label: 'Trade-offs', intro: 'How the power changes with more or fewer animals, an earlier ' +
      'or later end, or measuring more or less often.', content: [
      block('Power by group size', 'With more or fewer animals', sized(power, 260)),
      block('Power by end of study', 'Ending earlier or later', endSummary, sized(byEnd, 280)),
      block('Measurement schedule', 'Measuring more or less often', scheduleBox)]},
    {id: 'protocol', label: 'Protocol', intro: 'The plan written out, ready to paste into your study protocol.',
      content: [
        block('Study protocol', 'Design, sample size and analysis, ready to paste',
          copyRow(() => protocolText, 'Study protocol'), analysis)]},
  ]);
  // an empty viewer says what columns it lacks: the charts hidden until the first plan gives them their data
  form.effect(() => [growth, power, byEnd].forEach((v) => v.root.style.visibility = result.value ? '' : 'hidden'));
  // the schedules are the design's costliest scan: computed while their tab is open
  form.effect(() => scheduleOpen.value = tabs.activeTab.value === 'alternatives');
  const content = divV([cards, tabs.root], 'pmx-plan-content');
  Object.assign(content.style, {gap: 'var(--dg-space-l)', minWidth: `${RESULTS_MIN}px`});

  // the first plan, and every later run, under u2's loader; what was shown stays, dimmed, until the new plan is in
  const overlay = busyOverlay('Designing the study…');
  const failure = div([], 'pmx-plan-error');
  Object.assign(failure.style, {color: 'var(--dg-error-color, #c62828)'});
  form.effect(() => overlay.hidden = !blocking.value);
  // a calendar change: its own loader, the calendar dimmed, the page in view; a plan from the form clears its note
  form.effect(() => {
    calendarLoader.hidden = !refreshing.value;
    calendarBox.style.opacity = refreshing.value ? '0.5' : '';
    if (blocking.value)
      calendarStatus.replaceChildren();
  });
  form.effect(() => {
    failure.hidden = error.value === null;
    failure.replaceChildren(...(error.value === null ? [] : [divH([dot('error'), span(error.value)])]));
  });

  const outdated = staleNote('The inputs changed: run to update the plan');
  form.effect(() => {
    outdated.hidden = !stale.value;
    content.style.opacity = stale.value ? '0.5' : '';
  });

  const results = resultsFrame([failure, outdated, content], overlay);
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
    const f = DG.Func.byName('PmxStudio:pmxPlanStudy');
    const fill = f.prepare({...Object.fromEntries(f.inputs.map((p) => [p.name, call.getParamValue(p.name)])),
      compareSchedules: true});
    void fill.call().then((c) => {
      if (result.peek() === r)
        schedules.value = c.outputs['schedules'];
    }).catch((e) => error.value = String(e?.message ?? e));
  });

  // the view opens at once; the default declaration (or the saved plan's) is planned behind the loader, once the
  // form has written its defaults into the call
  void form.settled.then(() => {
    for (const [name, value] of Object.entries(initial?.inputs ?? {})) {
      const v = name === 'arms' ? armsFrame(value as SavedPlan['inputs']['arms']) : value;
      const field = form.getInput(name);
      if (field)
        field.value.value = v;
      call.setParamValue(name, v);
    }
    return replan();
  });

  // the plan in view as a file: saved to the user's plans (shared from Files like any file) or downloaded, to send
  // (the plan's predictive checks are computed then, not on every replan: planOf)
  const ready = (): PlanFrames | null => {
    const r = result.peek();
    if (!r || stale.peek()) {
      grok.shell.warning(r ? 'The inputs changed: run to update the plan first' : 'The plan is not ready yet');
      return null;
    }
    return r;
  };
  const defaultName = (r: PlanFrames) => initial?.name ?? `${r.plan.get('line', 0)} plan`;
  const saveDialog = () => {
    const r = ready();
    if (!r)
      return;
    const name = new TextInput({label: 'Name', value: defaultName(r), nullable: false});
    const where = span(`Saved to ${plansFolder()}. Share it from Files, or pick it in Analyze's Plan field.`);
    Object.assign(where.style, {color: 'var(--dg-text-color-muted, #68727f)', maxWidth: '360px', lineHeight: '1.45'});
    // a plan of the same name is replaced only when Save is pressed again for it (the dialog says so and stays open)
    let replacing: string | null = null;
    new Dialog('Save plan').add(name).add(where).onOK(async () => {
      const chosen = name.value.value.trim() || defaultName(r);
      if (chosen !== replacing && await planSaved(chosen)) {
        replacing = chosen;
        throw new Error(`A plan named "${chosen}" is already saved: press Save again to replace it, or rename it`);
      }
      const path = await savePlan(await planOf(chosen, call, r));
      grok.shell.info(`Plan saved: ${path}`);
    }, {text: 'Save'}).show();
  };
  const saveButton = buttonWithIcon('Save plan', 'save', saveDialog);
  saveButton.title = 'Save the plan to your files, to compare the study with it in Analyze';
  const saving = signal(false);
  const downloadButton = buttonWithIcon('Download', 'download', () => {
    const r = ready();
    if (!r)
      return;
    saving.value = true;
    void planOf(defaultName(r), call, r).then(downloadPlan).catch((e) => error.value = String(e?.message ?? e))
      .finally(() => saving.value = false);
  });
  downloadButton.title = 'Download the plan as a file, to send it to someone';
  form.effect(() => [saveButton, downloadButton].forEach((b) => b.disabled = !result.value || busy.value ||
    saving.value));

  const panes = new Splitter([declaration, results], {direction: 'horizontal', sizes: [0.3, 0.7],
    minSize: DECLARATION_MIN});
  return appView({
    name: 'PMX Studio · Plan',
    content: panes,
    ribbon: [[stageStrip('plan', (s) => void openStage(s))], [saveButton, downloadButton]],
    path: '/plan',
    own: [form, armsInput, correction, scale, growth, startDay, tabs],
    status: computed(() => busy.value ? 'Designing the study…' : result.value ? result.value.plan.name : ''),
  });
}
