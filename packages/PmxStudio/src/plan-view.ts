/* The Plan view: the declaration form is generated from pmxPlanStudy's own signature (u2 funcForm, grouped by
   the parameters' categories), so the function Grokky reads and the form people fill in are one definition; the
   treated arms and the multiplicity correction have editors of their own (arms-input.ts), and the measurement and
   randomization weekdays are set on the study calendar. The results are the function's outputs: the study's
   figures, the groups, the predicted volumes (band-chart.ts), the study calendar, then, folded, the power by
   group size, by end of study (platform line charts) and by measurement schedule, and the study protocol. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BasicTable, buttonWithIcon, ChoiceInput, computed, Dialog, div, divH, divV, dot, effect, emptyState, icon,
  loader, Menu, NumberInput, Reason, Form, signal, span, tableFromMap,
  TextInput, Tooltip} from '@datagrok-libraries/u2';
import {Emitter} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';
import {funcForm, ViewHost, viewers} from '@datagrok-libraries/u2/dg';
import type {ViewContent} from '@datagrok-libraries/u2/dg';

import {PlanFrames, PlanInputs, WEEKDAY_NAMES, declarationOf, planFrames, tgisOf, weekdaysOf, weekdaysText}
  from './frames';
import {StudyLink, armColors, follow, hasOwnColors, linkStudy, paint, resetColors, slotsAfter} from './colors';
import {ARCHETYPES, Therapy} from './archetypes';
import {Arm, ArmsInput, CorrectionInput, KIND_WORDS, PlanVolumes, armsFrame, lackingImmunity} from './arms-input';
import {ENDPOINT_RISK, MIN_GROUP, StudyDeclaration, combinationFailed, combinationN, combinationNote, combinationsOf,
  designCombination,
  LAB_SENSITIVE, NOTE_TIP, diameterOf, finite, isFiniteNumber, labSpreadGroup, safeEndpoint,
  CORRECTION_NAME, correctionOf, correctionOver, synergyClause} from './engine';
import {ArmPart, drugKey, parseArmName, shortName, therapyOf} from './arm-name';
import {askedOf, autoRun, block, callInputs, cardGrid, closeNamed, copyRow, copyText, fixedValue, hintLink, notes,
  notesBox, num, pct, settleForm, showCentred, sized, markStage, statCard, swatch, volumeAxis, wrapCells, stagePanes,
  whenSettled, docSubtitle, metricWord, setUserSetting, userSetting, modal, showOnce, objectPanel, showDepth, inWords,
  noteLine}
  from './ui';
import {PmxGroup} from './pmx-objects';
import {checklistBlocks, checklistDialog} from './arrive';
import {PREPARE_CREDIT, prepareDialog, prepareItems} from './prepare';
import {ethicsParts, ethicsText} from './ethics-text';
import {DOCX_TYPE, docx} from './docx';
import {ARRIVE_STUDY_PLAN_URL, FormSection, NTS_CREDIT, NTS_URL, NtsLanguage, PlanDocFacts, STUDY_PLAN_CREDIT,
  formBlocks, formItems, ntsForm, ntsTitle, studyPlanForm} from './ethics-documents';
import {canReturn, keepOnLeaving, linkRow, openStage, showStage} from './stages';
import {customizePage, panelGrid, Slot} from './panels';
import {stagePage, studyHeader} from './study-header';
import {DEMO_PLAN_ID, demoFinishedTable, openStudyDialog, stageLinks, takesFiles} from './hub-view';
import type {OpenAt} from './hub-view';
import {historyDialog} from './history-dialog';
import {followPlan, followedFrom} from './monitor-view';
import {PLAN_ATLAS_NOTE, SavedPlan, SavedShown, atlasNote, downloadPlan, framesOf, parsePlan, planOf, planSaved,
  plansFolder, planKey, sameAsked, savePlan, shownOf, shownPath, engineFingerprint} from './plan-file';
import {sameJson, maxOf, minOf} from './format';
import {scheduleText} from './schedule';
import {labHumane} from './study-file';
import {Saved, domainStorageOn, fileOwnerOf, isReadOnly, newRowNote, open as openSaved, rowOfPlan,
  rowsOf as studyRowsOf, save, saveCopy, saveError, withCode} from './store';
import {dirtyViews, discardedViews, leaveAsks, offerCopy, saveGate} from './native';
import {StudyMeta, armDoses, dateOf, dateText, doseForTotal, dosingDays, emptyMeta, weekdayOf, renamedDetails,
  onlyGroups, withGroups}
  from './study-meta';
import {CAL, MenuAt, WEEKDAYS, dayActions, icsOf, mark, marked, weeksGrid} from './calendar';
import {Exports, download, reasonOf, relabel, stageMenu, reportHtml, stageToolbox} from './report';
import {entered, esc, printSheet, sig} from './sheet';
import {LineInput, lastLine, tumorsOf} from './line-input';
import {historyNotesOf} from './history-models';
import {historyModelName, hostClassOf, modelName} from './priors';
import {BandChart} from './band-chart';
import {OwnStudy, StudySource, findAnchor, findPatientsDose, findUndosed, previousOfTable} from './previous-studies';
import {effectFromStudy} from './study-effect';
import {drugPkInWorkspace, forgetNca} from './exposure';
import {mm3Text} from './monitor';

// where Plan's ethics documents sit in its Export menu
const ETHICS = 'Documents';

/** A group's swatch and short name, a combination wrapping only at its "+"; its full name (and `more`) on hover. */
function shortSwatch(name: string, color: string, more = ''): HTMLElement {
  const short = shortName(name);
  const el = swatch(short, color, name);
  el.lastElementChild!.replaceChildren(...short.split(' + ').flatMap((p, i) => {
    const part = span(p);
    part.style.whiteSpace = 'nowrap';
    return i ? [' + ', part] : [part];
  }));
  el.title = more ? `${name}\n${more}` : name;
  return el;
}

/** A combination's note in a few words: the animals it needs to beat each drug alone (the sentence on hover). */
function comboLine(name: string, row: {[k: string]: unknown}): string {
  const need = finite(row['n_vs_best']);
  return `To show ${shortName(name)} beats each drug alone: ${need === null ? 'over 40' : need} animals per group`;
}

/** Plan's panels, in order and width (twelfths): the growth and the groups (the plan's answer), then folded until
 *  asked the calendar, what a smaller effect would need, the trade-offs of the size, end and schedule, the protocol. */
export const PLAN_LAYOUT: Slot[] = [{id: 'growth', span: 12}, {id: 'groups', span: 12},
  {id: 'calendar', span: 12, folded: true}, {id: 'scenarios', span: 12, folded: true},
  {id: 'whatif', span: 12, folded: true}, {id: 'protocol', span: 12, folded: true}];

/** A sized group is a Monte Carlo estimate (another seed moves some by one animal): said in the protocol and in the
 *  group size's header, not as a line of its own */
const SIZED_BY_SIMULATION = 'estimated by simulation, accurate to about one animal';

/** An arm's power: at its group size; with none (no size up to 40 reaches the target), the power at the largest size
 *  simulated, said so; '—' when there is none. */
export function powerText(n: number | null, power: number | null, curve: DG.DataFrame, arm: string): string {
  if (n !== null)
    return pct(power);
  let at: {n: number, p: number} | null = null;
  for (let i = 0; i < curve.rowCount; i++) {
    const [m, p] = [finite(curve.get('n', i)), finite(curve.get('assurance', i))];
    if (curve.get('arm', i) === arm && m !== null && p !== null && (!at || m > at.n))
      at = {n: m, p};
  }
  return at ? `${pct(at.p)} with ${at.n}` : '—';
}

/** What a section with nothing to compare says: no arm reaches the target power with any group size. */
const UNREACHABLE = 'No group size up to 40 reaches the target power: lower it, or set the animals per arm.';

/** Where a study's total has no number: the limit is per group, and totals over 40 are common. */
const OVER_40 = 'over 40 per group';
/** The widest the form's column opens (the user may drag it wider): its fields and arm cards fill about this much. */
const DECLARATION_MAX = 720;

/** A count of the plan `p` (its column `name`): with some arm unsized, the control's and the sized arms' (the
 *  `_sized` column); null: no group size reaches the power. */
const countOf = (p: DG.DataFrame, name: string): number | null =>
  finite(p.get(name, 0)) ?? (p.col(`${name}_sized`) ? finite(p.get(`${name}_sized`, 0)) : null);

/** The analyzable animals per group, as the card says it: "13 per group", or "13–14 per group" when each arm is
 *  sized for its own effect (the control as the largest); null: a group with no size. */
export function perGroupText(sizes: (number | null)[]): string | null {
  if (!sizes.length || sizes.some((n) => n === null))
    return null;
  const [lo, hi] = [minOf(sizes as number[]), maxOf(sizes as number[])];
  return lo === hi ? `${lo} per group` : `${lo}–${hi} per group`;
}

/** The longest name a plan or a study is saved under (the Overview and the ribbon show it whole). */
export const NAME_MAX = 100;

/** Why the endpoint volume cannot be planned (null: it can): the controls start at the volume at randomization, so an
 *  endpoint at or below it is reached on day 0. */
export const endpointProblem = (endpoint: unknown, vRand: unknown): string | null =>
  typeof endpoint === 'number' && typeof vRand === 'number' && endpoint <= vRand ?
    `Must be above the volume at randomization (${vRand} mm³)` : null;

/** The planned review days of the Planned reviews field (empty: none), or why they are not days: each a whole study
 *  day after randomization, once, up to the latest end Plan allows (`max`). */
export function reviewDays(text: unknown, max = 120): {days: number[]} | {problem: string} {
  const parts = typeof text === 'string' ? text.split(/[,;\s]+/).filter((x) => x) : [];
  const days: number[] = [];
  for (const x of parts) {
    const d = Number(x);
    if (!/^\d+$/.test(x) || d < 1 || d > max)
      return {problem: `"${x}" is not a study day: write whole days from 1 to ${max}, such as 7, 14`};
    if (days.includes(d))
      return {problem: `Day ${d} is written twice`};
    days.push(d);
  }
  return {days: days.sort((a, b) => a - b)};
}

/** What the lab mark says: `n` the analyzable animals per group were labs as far apart as is plausible (null: over
 *  40; undefined: being worked out; 'failed': not known), `central` the card's own (null: over 40). A number only
 *  when it is above the central one; else, or not known, the line without one. */
export function labMarkText(n: number | null | undefined | 'failed', central: number | null | undefined): string {
  const lead = 'For this model the number depends more than usual on how much laboratories differ';
  const notAbove = n !== undefined && n !== 'failed' && central !== undefined &&
    (central === null || (n !== null && n <= central));
  return n === 'failed' || notAbove ? `${lead}: it may need more animals.` :
    n === undefined ? `${lead}. Working out how many animals it could need…` :
      `${lead}: if they differ as much as is plausible, it would be ` +
      (n === null ? 'over 40 per group.' : `${n} per group.`);
}

/** The discreet mark on the analyzable card of a plan on a line whose group size depends more than usual on how far
 *  apart labs are (engine LAB_SENSITIVE), hidden on the others. Hovered, focused or clicked, it says the plan's own
 *  group were labs as far apart as is plausible: the plan designed again at that spread (labSpreadGroup),
 *  asked then and kept for the plan in view; `central` the card's own group, which it must exceed to be said.
 *  `owner` builds its effect. */
export function labMark(owner: {build: <T>(f: () => T) => T}, plan: ReadonlySignal<StudyDeclaration | undefined>,
  central: () => number | null | undefined): HTMLElement {
  const mark = icon('info-circle', {cls: 'pmx-lab-mark'});
  Object.assign(mark, {tabIndex: 0, role: 'button'});
  mark.setAttribute('aria-label', 'How much the animals per group depend on how laboratories differ');
  let asked: {key: string, n?: number | null | 'failed'} | undefined;
  let over = false;
  const show = () => {
    const d = plan.peek();
    if (!d)
      return;
    const key = JSON.stringify(d);
    if (asked?.key !== key) {
      const a: typeof asked = {key};
      asked = a;
      labSpreadGroup(d).then((n) => a.n = n, (e) => {
        console.error('PMX Studio:', e);
        a.n = 'failed';
      }).then(() => {
        if (over && asked === a)
          Tooltip.show(labMarkText(a.n, central()), mark);
      });
    }
    over = true;
    Tooltip.show(labMarkText(asked.n, central()), mark);
  };
  const hide = () => {
    over = false;
    Tooltip.hide();
  };
  mark.addEventListener('pointerenter', show);
  mark.addEventListener('focus', show);
  mark.addEventListener('pointerleave', hide);
  mark.addEventListener('blur', hide);
  // the card's own click and press (set the animals yourself) are not the mark's
  for (const type of ['click', 'pointerdown', 'keydown'] as const) {
    mark.addEventListener(type, (e) => {
      if (type === 'keydown' && !['Enter', ' '].includes((e as KeyboardEvent).key))
        return;
      e.stopPropagation();
      e.preventDefault();
      if (type !== 'pointerdown')
        show();
    });
  }
  owner.build(() => effect(() => {
    const d = plan.value;
    mark.style.display = d && LAB_SENSITIVE.includes(d.line) ? '' : 'none';
  }));
  return mark;
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
    onNote: (day: number) => void, onSelect?: (day: number) => void, reviews?: number[]}): HTMLElement {
  const v = (c: string) => p.get(c, 0);
  const end = v('study_end') as number;
  const endpointDay = v('endpoint_day') as number;
  const humane = v('humane_limit');
  const watch = finite(v('humane_watch_day'));
  const details = new Map<number, string>(Array.from({length: calendar.rowCount},
    (_, i) => [calendar.get('day', i), calendar.get('activity', i)]));
  const measured = new Set(String(v('measurement_days')).split(', ').map(Number));
  const weekdays = new Set(weekdaysOf(v('measurement_weekdays')));
  const {blue: BLUE, slate: SLATE, red: RED, dark: DARK, green: GREEN} = CAL;
  const meta = study.meta;
  const date = (day: number) => dateOf(meta.start, day);
  const doses = dosingDays(meta.groups, end, meta);

  // day 0 falls on the randomization weekday: the measurement weekdays are the lab's, seen from it
  const start = (v('start_weekday') ?? 0) as number;
  const weekday = (day: number) => (start + day) % 7;

  // above the weeks: the randomization weekday (picked in the sentence), the end and the closer watch
  const q90 = v('duration_q90');
  // (in words, not "in 8 of 10 studies": designer wording 4)
  const range = isFiniteNumber(q90) ? `usually between day ${Math.round(v('duration_q10'))} and day ` +
    `${Math.round(q90)}` : `usually from day ${Math.round(v('duration_q10'))}`;
  // the weekday is day 0's date's once it has one: said beside the picker it locks (a disabled field shows no hover)
  const locked = icon('info-circle');
  locked.title = 'Day 0 has a date, so it sets the weekday: change the date from a day\'s menu in the calendar, or ' +
    'in the study\'s details';
  const first = divH([mark('Randomize', BLUE, true), span('Randomize on a'), startDay, ...meta.start ? [locked] : [],
    span(`(day 0${meta.start ? `, ${dateText(date(0)!, true)}` : ''}) · ${measured.size} measurement days.`)]);
  Object.assign(first.style, {alignItems: 'center', gap: '6px', flexWrap: 'wrap'});
  const lines: HTMLElement[] = [first];
  if (v('last_day_chosen')) {
    lines.push(marked(mark('End', DARK, true), `The study ends on day ${end}, as chosen` +
      (weekdays.has(weekday(end)) ? '.' : `: a ${WEEKDAY_NAMES[weekday(end)]}, not one of your measurement days.`)));
  } else {
    const reach = `the median control reaches ${v('endpoint')} mm³`;
    lines.push(marked(mark('End', DARK, true), end === endpointDay ? `The study ends on day ${end}, when ${reach} ` +
      `(${range}).` : `${reach[0].toUpperCase()}${reach.slice(1)} on day ${endpointDay} ` +
      `(${range}), a ${WEEKDAY_NAMES[weekday(endpointDay)]}: the study ends on the next measurement ` +
      `day, day ${end} (${WEEKDAY_NAMES[weekday(end)]}).`));
  }
  if (watch !== null) {
    lines.push(marked(mark('Watch controls', RED, false), `From day ${watch}, check the controls more often: an ` +
      `animal over ${humane} mm³ is removed.`));
  }
  if (doses.size > 0) {
    const dosed = meta.groups.filter((g) => [...doses.values()].some((gs) => gs.includes(g.group)));
    // each group from its first dose (a schedule may start later, "from day 3"); no dosing at weekends said once
    const firstOf = (g: string) => minOf([...doses].filter(([, gs]) => gs.includes(g)).map(([d]) => d));
    // one short entry per drug; the arms it covers go to the tooltip, so the line stays within the calendar
    // (one drug however its arms write it, "Cpd-A", "cpd-a": its first writing, each of its schedules once)
    const perDrug = new Map<string, {name: string, schedules: Map<string, string>, from: number}>();
    for (const g of dosed) {
      // (a group with no details: its drug and schedule as its name gives them, never the name's raw "q3dx4", QA 0.1.3
      // round 4 designer U5)
      const named = parseArmName(g.group);
      const name = g.compound || (named.length === 1 ? named[0].drug : '') || g.group;
      const e = perDrug.get(drugKey(name)) ?? {name, schedules: new Map(), from: Infinity};
      // (as typed in the details too, in the one writing every schedule is shown in: "Q3D×4")
      for (const sched of g.schedule.trim() ? [g.schedule] : named.map((p) => p.schedule ?? '').filter((x) => x)) {
        if (!e.schedules.has(sched.toLowerCase()))
          e.schedules.set(sched.toLowerCase(), scheduleText(sched));
      }
      e.from = Math.min(e.from, firstOf(g.group));
      perDrug.set(drugKey(name), e);
    }
    const dose = marked(mark('Dose', GREEN, false), [...perDrug.values()].map((e) =>
      `${e.name} ${[...e.schedules.values()].join(' / ')} from day ${e.from}`.replace(/\s+/g, ' ')).join(' · ') +
      (meta.weekendsOff && meta.start ? ' · no dosing at weekends' : ''));
    Tooltip.bind(dose, dosed.map((g) => g.group).join(', '));
    lines.push(dose);
  }
  // the planned reviews (Planned reviews): those within the study, and the ones after its end said not held
  const reviews = new Set((study.reviews ?? []).filter((d) => d <= end));
  const after = (study.reviews ?? []).filter((d) => d > end);
  if (reviews.size || after.length) {
    lines.push(marked(mark('Review', DARK, false), (reviews.size ? `Review each treated arm on day ` +
      `${[...reviews].join(', ')}.` : '') + (after.length ? ` Day ${after.join(', ')} ${after.length > 1 ? 'fall' :
      'falls'} after the end (day ${end}): no review then.` : '')));
  }
  const summary = divV(lines, 'pmx-calendar-summary');
  Object.assign(summary.style, {gap: '4px', lineHeight: '1.45', marginBottom: 'var(--dg-space-m)'});

  // the weeks: a header of the lab's weekdays (click one to measure on it or not), then a cell a day up to the end
  const header = WEEKDAYS.map((d, i) => {
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
  });
  const grid = weeksGrid({meta, startWeekday: start, last: end, header, onMenu: study.onMenu, onNote: study.onNote,
    onSelect: study.onSelect,
    day: (day, own) => ({detail: details.get(day),
      marks: [...day === 0 ? [mark('Randomize', BLUE, true)] : [], ...measured.has(day) ?
        [mark('Measure', SLATE, false)] : [], ...own, ...day === watch ? [mark('Watch controls', RED, false)] : [],
      ...reviews.has(day) ? [mark('Review', DARK, false)] : [],
      ...day === end ? [mark('End', DARK, true)] : []],
      border: day === end ? `1px solid ${DARK}` : undefined,
      // the days the controls need a closer watch, tinted
      tint: watch !== null && day >= watch ? `color-mix(in srgb, ${RED} 6%, transparent)` : undefined})});
  return divV([summary, grid], 'pmx-calendar-box');
}

/** One row per measurement schedule: its measurement days, the animals it needs and each arm's power with the
 *  plan's animals (`own`: the plan's measurements a week). */
function schedulesTable(sc: DG.DataFrame, own: number | undefined): HTMLElement {
  type Schedule = [number, {measurements: number, nTotal: number | null, power: {[arm: string]: number | null}}];
  const byWeek = new Map<number, Schedule[1]>();
  for (let i = 0; i < sc.rowCount; i++) {
    const k = sc.get('per_week', i);
    const row: Schedule[1] = byWeek.get(k) ?? {measurements: sc.get('measurements', i),
      nTotal: finite(sc.get('n_total', i)), power: {}};
    row.power[sc.get('arm', i)] = finite(sc.get('assurance', i));
    byWeek.set(k, row);
  }
  const armNames = [...new Set(sc.col('arm')!.toList() as string[])];
  const columns = [
    {header: 'Measurements a week', render: ([k]: Schedule) => k === own ? `${k} · this plan` : String(k)},
    {header: 'Measurement days', render: ([, s]: Schedule) => String(s.measurements)},
    // (no total: some arm needs over 40 at that schedule; its power column says which)
    {header: 'Analyzable animals', render: ([k, s]: Schedule) => animals(s.nTotal ?? 'an arm over 40', k === own)},
    ...armNames.map((arm) => ({header: armNames.length > 1 ? `Power, ${arm}` : 'Power',
      render: ([, s]: Schedule) => pct(s.power[arm])})),
  ];
  const table = new BasicTable<Schedule>({items: [...byWeek], keyOf: ([k]) => String(k),
    columns: columns.map((c, j) => ({...c, align: j >= 1 ? 'center' : 'left'}))});
  table.root.style.width = '100%';
  // the plan's own row is the design's simulated studies, the others a smaller set of their own (Crux design.rs, the
  // SCHEDULES stream at mde_animals): a point or two between rows is that, not the schedule (QA: 81, 79, 78%)
  const noise = div(['A difference of one or two points of power between rows is simulation noise, not the ' +
    'schedule.'], 'ui-label-details');
  noise.style.marginTop = '4px';
  return divV([table.root, noise]);
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

/** The demo plan's earlier studies: the demo's own finished study, alone. */
export const demoStudies: StudySource = () => demoFinishedTable()
  .then((df) => [previousOfTable(df)].filter((s): s is NonNullable<typeof s> => s !== null));

/** A combination's drugs as its name gives them, each once (by drugKey, its first writing) with its own dose, route
 *  and schedule: "cpd-a 30 mg/kg p.o. qd". */
export function combinationDrugs(parts: ArmPart[]): string[] {
  const seen = new Set<string>();
  return parts.filter((p) => !seen.has(p.key) && !!seen.add(p.key)).map((p) =>
    [p.drug, p.dose && `${entered(p.dose.value)} ${p.dose.unit}`, p.route, p.schedule].filter(Boolean).join(' '));
}

/** A .pmxplan file opened from Files or dropped on the platform: Plan with it (no table to open). */
export async function openPlanFile(content: string): Promise<DG.DataFrame[]> {
  await showStage('plan', await planView(parsePlan(content)));
  return [];
}

/** What Plan opens with: a saved plan, or only its name, details and form. */
export type PlanStart = Pick<SavedPlan, 'name' | 'study' | 'provenance'> & {inputs: {[param: string]: unknown},
  predicted?: {shown?: SavedShown}, row?: Saved | null};

/** The saved row a Plan view opens on: the one it carries (a plan come back to, edited since: its content no longer
 *  names the row, QA 0.1.3 round 4 designer D2), else the row read for that plan. */
const rowOfStart = (initial?: PlanStart): Saved | null => initial?.row !== undefined ? initial.row : rowOfPlan(initial);

// each Plan view's call (its form's inputs, hidden ones too), by its host: what the package's tests read
const planCalls = new WeakMap<ViewHost, DG.FuncCall>();
/** The call behind a Plan view: its inputs as the next plan takes them (for the package's tests). */
export const planCallOf = (view: DG.ViewBase): DG.FuncCall | undefined => planCalls.get(ViewHost.of(view));
// each Plan view's Save, which File › Save presses: what the package's tests press and read
const planSaves = new WeakMap<ViewHost, HTMLElement>();
/** The Save behind a Plan view's File › Save (for the package's tests). */
export const planSaveOf = (view: DG.ViewBase): HTMLElement | undefined => planSaves.get(ViewHost.of(view));
type PlanFileFn = (name: string) => Promise<SavedPlan>;
const planFiles = new WeakMap<ViewHost, PlanFileFn>();
/** The plan a Plan view shows, as its file would keep it (for the package's tests). */
// the plans Plan has made, in any view (a remade arm is the engine's cache's, so Crux's calls do not count them)
let plansMade = 0;
/** How many plans Plan has made since the page opened (for the package's tests). */
export const plansMadeCount = () => plansMade;

export const planFileOf = (view: DG.ViewBase, name: string): Promise<SavedPlan> =>
  planFiles.get(ViewHost.of(view))!(name);

export async function planView(initial?: PlanStart): Promise<DG.ViewBase> {
  const exports: Exports = [];
  const host = ViewHost.create({name: 'PMX Studio', content: (self) => planContent(self, initial, exports)});
  const view = markStage(host.view, 'plan');
  // (opened from a saved row: its URL names it)
  linkRow(view, 'plan', rowOfStart(initial)?.id);
  return view;
}

/** What the view shows, built under its host, which owns it all. */
function planContent(self: ViewHost, initial: PlanStart | undefined, exports: Exports): ViewContent {
  const call = DG.Func.byName('PmxStudio:pmxPlanStudy').prepare();
  planCalls.set(self, call);
  // the engine itself on the form's inputs (with `extra` over them), not the registered function's FuncCall (see
  // callInputs); a wrong field rejects, as the call did
  // each plan's declaration, for what is computed from it after it shows (the endpoint offered)
  const declared = new WeakMap<PlanFrames, StudyDeclaration>();
  const designed = (extra: object = {}) => Promise.resolve().then(async () => {
    plansMade++;
    const x = callInputs<PlanInputs & {tgiScenarios?: string | null}>(call, extra);
    // (the combinations' sentences asked after the plan shows: their Crux calls kept off the plan's own wait)
    const d: StudyDeclaration = {...declarationOf(x), combinationNotes: false};
    const frames = await planFrames(d, tgisOf(x.tgiScenarios, d.endpoint, d.vRand));
    declared.set(frames, d);
    return frames;
  });
  const result = signal<PlanFrames | undefined>(undefined);
  // the study's details (header's name, the calendar's dates, doses and notes); an ID once the plan leaves Plan
  // (saved, downloaded or followed)
  const meta = signal<StudyMeta>(initial?.study ?? emptyMeta());
  // a plan's arms by name, and the study's colours for its groups
  const armsIn = (r: PlanFrames) => r.armResults.col('arm')!.toList() as string[];
  const colorsIn = (r: PlanFrames, m = meta.peek()) => armColors(armsIn(r), m.colors, m.slots);
  // the form's and the calendar's replans share the call, one run at a time: a change made during a run is planned
  // next, and the run it overtook is never shown
  let running: Promise<void> | null = null;
  let again = false;
  const error = signal<string | null>(null);
  // the "Total dose for…" dialog open, if one is: it goes with the view, or when the study is followed
  let doseDialog: DG.Dialog | null = null;
  self.own(() => doseDialog?.close());
  // (the fields carry no caption of their own: their section's title says it)
  const lineInput = new LineInput({});
  // the arms read the kind of treatment from their names and the model's host; TGI or T/C is the study's (Advanced)
  const effectAs = signal<'TGI' | 'TC'>('TGI');
  // a table opened, closed or edited: a PK study may have come, gone or changed (its NCA run is read again)
  const tables = signal(0);
  const watched = new WeakSet<DG.DataFrame>();
  const watch = () => {
    for (const t of grok.shell.tables.filter((t) => !watched.has(t))) {
      watched.add(t);
      // ponytail: a closed table's subscription lives until the view closes; few tables, so not worth tracking
      self.subscribe(DG.debounce(t.onDataChanged, 500), () => {
        forgetNca(t);
        tables.value++;
      });
    }
  };
  watch();
  self.subscribe(grok.events.onTableAdded, () => {
    watch();
    tables.value++;
  });
  self.subscribe(grok.events.onTableRemoved, () => tables.value++);
  const own = (): OwnStudy => ({id: meta.peek().id, name: initial?.name ?? ''});
  // the demo plan scales its arms from the demo's own finished study, never the user's library or a copy of someone
  // else's (QA 0.1.3 final designer #9)
  const demoSource = (): StudySource | undefined => meta.peek().id === DEMO_PLAN_ID ? demoStudies : undefined;
  const planVolumes = signal<PlanVolumes>({v0: 100, endpoint: 1000});
  const armsInput = new ArmsInput({metric: effectAs, volumes: planVolumes,
    host: computed(() => hostClassOf(lineInput.value, lineInput.models.value)),
    // an arm whose drug the user studied before may take its effect from that study, read on the plan's own day
    previous: {line: computed(() => lineInput.value ?? null),
      evaluationDay: computed(() => finite(result.value?.plan.get('endpoint_day', 0))),
      hostOf: (l) => hostClassOf(l, lineInput.models.value),
      // never the study this plan is of (opened from it): its arms would be estimated from themselves
      // (a study of the plan's own model first)
      find: (n) => findAnchor(n, own(), lineInput.value ?? null, demoSource()),
      undosed: (n) => findUndosed(n, own(), demoSource()),
      meta: () => meta.value,
      // the patients' dose under this plan's rules (no dosing at weekends), not the earlier study's
      patientsDose: (n) => findPatientsDose(n, meta.value, own(), demoSource()),
      renamed: (from, to) => meta.value = renamedDetails(meta.peek(), new Map([[from, to]])), tables,
      // within the doses of the drug's PK, as NCA reads a concentration table open in the workspace (as Analyze)
      // (by the earlier study's route: Crux scales within one route)
      estimate: async (a, name, day) => effectFromStudy(a, name, await drugPkInWorkspace(a.part.key, a.part.route),
        day)}});
  const correction = new CorrectionInput({label: 'Multiple comparisons', arms: armsInput.count});
  // a plan that no longer matches the form is dimmed until the next run
  const stale = signal(false);
  // captions above every field: the panel reads top to bottom, left-aligned, at any width
  // no Run: each change replans once the changes stop
  // the inputs the plan in view was made with: a form that settles back to them (filled again as the view opens, a
  // list of models read late) is no change, so no plan is made again under the page's bar
  let shownAsked: string | null = null;
  const auto = autoRun(() => form, () => result.peek() && !error.peek() && askedOf(call) === shownAsked ?
    stale.value = false : replan());
  // a change of the form since the plan was opened or saved (domain storage: what closing the view asks about)
  let edited = false;
  // whether the user has typed or clicked yet with this page current (its form, toolbox, menus, dialogs): a change
  // before is the form settling itself
  let touched = false;
  const touch = (e: Event) => touched ||= e.isTrusted && grok.shell.v?.dart === self.view.dart;
  for (const type of ['pointerdown', 'keydown', 'paste', 'drop'])
    document.addEventListener(type, touch, true);
  self.own(() => ['pointerdown', 'keydown', 'paste', 'drop'].forEach((t) =>
    document.removeEventListener(t, touch, true)));
  // the form's values as last saved here (null: not saved yet): the app refilling the form after a save (the details
  // as saved, the weekday they date) changes none of them, so it is no change, whenever its events arrive
  let savedInputs: string | null = null;
  // (a change of the form, and a save from wherever it was asked: the unsaved marker checks then, native.ts saveGate)
  const formChanges = new Emitter();
  // (the study's details with them: Study details and a day's exceptions in the calendar are changes of the plan)
  const inputsNow = () => JSON.stringify([...call.func.inputs.map((p) => {
    const v = call.getParamValue(p.name);
    return v instanceof DG.DataFrame ? v.toCsv() : v;
  }), meta.peek()]);
  // the study's details changed by the user (Study details, a day of the calendar, a colour): a change of the plan
  // as a field's is, marked unsaved and asked about on leaving (QA 0.1.3 designer D1, D4)
  const userMeta = (m: StudyMeta) => {
    meta.value = m;
    edited = true;
    formChanges.fire();
  };
  // (no table fills the arms: a study table open in the workspace, as when coming from Analyze, is not a plan's)
  const form = funcForm(call, {layout: 'tall', skipTableAutoFill: true,
    overrides: {arms: {input: armsInput}, correction: {input: correction}, line: {input: lineInput}},
    onInputChanged: (name) => {
      // (the schedules' comparison switched on is no change of the plan: Plan makes it again itself, the remade effect,
      // and a plan reopened on its saved tables, schedules compared, is not out of date for it)
      if (name === 'compareSchedules')
        return;
      // (the form filling itself as the view opens is no change: only one the user made once a plan is shown; an
      // old saved plan's form settling a value late, a moment after its first plan, read "Unsaved changes" with
      // nothing edited, QA 0.1.3 round 3 designer W2)
      edited ||= touched && result.peek() !== undefined;
      stale.value = result.peek() !== undefined;
      auto.changed();
      formChanges.fire();
    }});
  settleForm(form, {tgiScenarios: 'Off, or e.g. 30, 50, 70'});
  // the inputs of the first plan shown are what a change is measured against until a save (a plan opened from a row
  // and left as it was asks nothing on leaving: its form refilling itself as it opens is no change)
  // (taken again with each plan until the user changes something: what the form settles by itself is the plan as
  // saved)
  form.effect(() => {
    if (result.value && !edited)
      savedInputs = inputsNow();
  });
  // checks across fields, before the engine is asked (an invalid field holds the plan back, and says why)
  const vRandInput = form.getInput('vRand');
  form.getInput('endpoint')?.addValidator((x) => endpointProblem(x, vRandInput?.value));
  form.getInput('looks')?.addValidator((x) => {
    const r = reviewDays(x);
    return 'problem' in r ? r.problem : null;
  });

  // busy while the view plans (the first plan, a change of the form, a change of randomization day)
  const planning = signal(true);
  // a change made on the calendar replans the study too (the measurement days set the power and the animals), but
  // under a loader of the calendar's own, the page left in view
  const refreshing = signal(false);
  // the plan made again for the schedules' comparison: under that section's own loader, not the page's bar
  const scheduling = signal(false);
  const busy = computed(() => planning.value || refreshing.value || scheduling.value);
  const replan = (under: 'plan' | 'calendar' | 'schedules' = 'plan'): Promise<void> => {
    const flag = under === 'calendar' ? refreshing : under === 'schedules' ? scheduling : planning;
    if (running) {
      again = true;
      // a change of the form meanwhile is planned next, under the bar it would have had
      flag.value = true;
      return running;
    }
    flag.value = true;
    const once = (): Promise<void> => {
      // no treated arm yet (a new plan): nothing to design, the page says what to add
      if (armsInput.count.peek() === 0) {
        error.value = null;
        stale.value = false;
        result.value = undefined;
        return Promise.resolve();
      }
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
        shownAsked = sent;
        result.value = c;
      }, (e) => {
        if (again)
          return next();
        error.value = inWords(String(e?.message ?? e));
      });
    };
    running = once().finally(() => {
      running = null;
      planning.value = refreshing.value = scheduling.value = false;
    });
    return running;
  };

  const plan = computed(() => result.value?.plan);
  const cell = <T = number>(name: string) => computed(() => (plan.value?.get(name, 0) ?? undefined) as T | undefined);
  const groups = computed(() => (cell('treated_arms').value ?? 1) + 1);
  // each group's analyzable animals: the control, then each arm (null: no size reaches the power)
  const groupSizes = computed((): [string, number | null][] => {
    const a = result.value?.armResults;
    const arm = (i: number): [string, number | null] => [a!.get('arm', i), finite(a!.get('n', i))];
    return a ? [['Control', finite(cell('n_control').value)], ...Array.from({length: a.rowCount}, (_, i) => arm(i))] :
      [];
  });
  // the animals in the order the study uses them, then its end; what each count covers in its tooltip
  // the change arrows start again after a save and once a new plan has its first arm: nothing to compare with before
  // (QA 0.1.3 final designer #22)
  const epoch = signal(0);
  let armsCounted = armsInput.count.peek();
  form.effect(() => {
    const k = armsInput.count.value;
    if (armsCounted === 0 && k > 0)
      epoch.value = epoch.peek() + 1;
    armsCounted = k;
  });
  const card = (label: Parameters<typeof statCard>[1], value: Parameters<typeof statCard>[2], tip: () => string,
    format?: Parameters<typeof statCard>[4]) => statCard(form, label, value, tip, format, true, () => epoch.value);
  // a count no group size reaches has no number (a column with no value reads as an empty string); with only some
  // arm unsized, the control's and the sized arms' (the unsized one says why in the groups table)
  const count = (name: string) => computed(() => plan.value ? countOf(plan.value, name) ?? '—' : undefined);
  // the arms no group size reaches, named in the cards' tooltips while the cards count the others
  const unsizedNote = () => {
    const out = groupSizes.peek().filter(([g, n]) => g !== 'Control' && n === null).map(([g]) => g);
    return out.length && finite(cell('n_total').peek()) === null && finite(cell('n_total_sized').peek()) !== null ?
      `Without ${out.join(', ')}: no group size up to 40 reaches the power. ` : '';
  };
  // the animals per group and the end of study are computed, or set from their cards (a click asks the value); set,
  // the card says so and a click offers to change it or go back to the computed one (the function's `n` and
  // `lastDay`, not in the form)
  // a field of the form and its parameter of the call, set together (`n` and `lastDay`: the call's alone)
  const setInput = (name: string, value: unknown) => {
    const field = form.getInput(name);
    if (field)
      field.value = value;
    call.setParamValue(name, value);
  };
  const apply = (name: string, value: number | null) => {
    setInput(name, value);
    // one plan, now, and only with every field of the form valid (the deferred run's own gate)
    auto.changed(0);
  };
  // a value (checked again here, whatever called it, since it is the number of animals a study uses; empty is an
  // error here, going back is back's)
  const setParam = (name: string, raw: number, min: number, max: number) => {
    const r = fixedValue(raw, min, max);
    if ('error' in r || r.value === null)
      throw new Error('error' in r ? r.error : 'Enter a number');
    apply(name, r.value);
  };
  const settable = (el: HTMLElement, name: 'n' | 'lastDay', set: () => boolean, ask: {title: string, label: string,
    now: () => number | undefined, min: number, max: number, back: string}) => {
    // a button for the keyboard too
    Object.assign(el, {tabIndex: 0, role: 'button'});
    el.style.cursor = 'pointer';
    // its label underlined in dots, the app's one mark of what can be changed (the study's name, the arms' hints)
    el.classList.add('pmx-stat-settable');
    const open = () => showOnce(self, `set-${name}`, () => {
      // the computed value as the start, within the field's range (a computed end past it would block OK)
      const start = Math.min(ask.max, Math.max(ask.min, ask.now() ?? ask.min));
      const input = new NumberInput({label: ask.label, value: start, min: ask.min, max: ask.max,
        mode: 'int', nullable: false});
      // in a form: OK waits while the number is not valid; the check again on OK, empty is not "back" here
      return new Dialog(ask.title).add(new Form().add(input)).onOK(() => {
        if (input.value === null)
          throw new Error('Enter a number');
        setParam(name, input.value, ask.min, ask.max);
      }, {text: 'Set'});
    });
    const act = (anchor: HTMLElement) => {
      if (!result.peek())
        return;
      if (!set())
        return open();
      Menu.popup().item('Change…', open).item(ask.back, () => apply(name, null)).show({anchor});
    };
    el.addEventListener('click', () => act(el));
    el.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        act(el);
      }
    });
  };
  const nSet = () => !!cell('n_fixed').peek();
  const endSet = () => !!cell('last_day_chosen').peek();
  // each card by the plan column it shows, so a change made on the calendar can point at what it moved
  const cardOf: {[column: string]: HTMLElement} = {
    // (the animals to randomize, what the implanted come to, on its second line)
    implant_total: card('Animals to implant', count('implant_total'), () =>
      `${unsizedNote()}Covers ${pct(cell('pre_loss').value)} ` +
      'that never reach randomization: no take, or out of the volume window. Of them, ' +
      `${count('randomise_total').value ?? '—'} to randomize: covers ${pct(cell('post_loss').value)} lost after ` +
      'randomization (ulceration, deaths, technical failures)'),
    n_total: card(computed(() => {
      const per = perGroupText(groupSizes.value.map(([, n]) => n));
      return (cell('n_fixed').value ? 'Analyzable · set' : 'Analyzable') + (per ? ` · ${per}` : '');
    }),
    count('n_total'),
    () => {
      // groups of different sizes: each arm sized for its own effect, the control as the largest
      const sizes = groupSizes.value;
      const each = perGroupText(sizes.map(([, n]) => n))?.includes('–') ?
        ` (${sizes.map(([g, n]) => `${g} ${n}`).join(', ')}: each arm sized for its own effect)` : '';
      return unsizedNote() + (cell('n_fixed').value ? `${groups.value} groups of the animals per arm you set. ` +
        'Click to change ' +
        'it, or go back to what the power needs' : `${groups.value} groups${each}, ` +
        `${pct(cell('power_target').value)} power per comparison. Click to set the animals per group yourself`);
    }),
    study_end: card(computed(() => cell('last_day_chosen').value ? 'End of study · set' : 'End of study'),
      cell('study_end'), () => cell('last_day_chosen').value ? 'Set by you; the median control reaches the ' +
        `endpoint on day ${cell('endpoint_day').value}. Click to change it, or go back to the endpoint` :
        (cell('study_end').value === cell('endpoint_day').value ? 'When the median control reaches the endpoint' :
          'The first measurement day after the median control reaches the endpoint (day ' +
          `${cell('endpoint_day').value})`) + '. Click to end on a day of your own',
      (x) => `Day ${x}`),
  };
  settable(cardOf.n_total, 'n', nSet, {title: 'Animals per group', label: 'Analyzable animals per group',
    // the largest group (the control's): what every group then has; with some arm unsized, the control as the cards
    // count it (the largest sized arm), not the field's minimum
    now: () => (plan.peek() && countOf(plan.peek()!, 'n_control')) ?? undefined, min: 2, max: 100,
    back: 'Back to what the power needs'});
  settable(cardOf.study_end, 'lastDay', endSet, {title: 'End of study', label: 'Last day of study',
    now: () => cell('study_end').peek(), min: 5, max: 120, back: 'Back to when the controls reach the endpoint'});
  const toRandomize = div([], 'pmx-stat-sub');
  form.effect(() => toRandomize.textContent = plan.value ? `${count('randomise_total').value} to randomize` : '');
  cardOf.implant_total.append(toRandomize);
  cardOf.n_total.append(labMark(form, computed(() => {
    const r = result.value;
    return r ? declared.get(r) : undefined;
  }), () => plan.peek() ? finite(cell('n_control').peek()) : undefined));
  const cards = cardGrid(Object.values(cardOf), 'pmx-plan-cards', 140);

  // the engine's warnings, each under the section it is about, one line each with the dot on the text's first line
  // (u2 has no banner, so a warning dot and the text)
  const groupNotes = notesBox();
  const growthNotes = notesBox();
  // the last drawing that asked for a history model's notes
  let historyAsked = 0;
  // each plan's combination sentences, once in the notes (the synergy completes them after)
  const sentencesOf = new WeakMap<PlanFrames, Promise<HTMLElement[]>>();
  // the plan whose combination sentences and their synergy are in (or failed): the plan's other work after it shows
  // (the schedules made again, the detectable effects) waits for them, so they are not queued behind its calls
  const sentenced = signal<PlanFrames | null>(null);
  growthNotes.style.margin = 'var(--dg-space-s) 0 var(--dg-space-m)';
  const analysis = div([], 'pmx-protocol');
  let protocolText = '';
  let ethicsApplication = '';
  // the plan's facts for its ethics documents (the name and details are read when one is opened)
  let docFacts: Omit<PlanDocFacts, 'title' | 'meta'> | null = null;
  const groupsTable = div([]);
  const scheduleBox = div([]);
  // the power by measurement schedule is the design's costliest scan: the plan takes it along only while its
  // section is open; opened later, the plan in view is made again with it. A plan made without it has none to show
  // (its table is empty either way: a plan made with it and no arm reaching the target power has none either)
  const scheduleOpen = signal(false);
  const schedules = computed(() => {
    const r = result.value;
    return r && declared.get(r)?.schedules ? r.schedules : undefined;
  });
  const calendarBox = div([]);
  // what a change on the calendar did: the counts and the end it moved (their cards flash), and a power that fell
  // below the target (with the animals per arm fixed; else the animals grow to keep it)
  const calendarLoader = loader('Updating the plan…');
  const calendarStatus = divV([], 'pmx-calendar-status');
  // how the calendar is edited, said while the pointer is on it, the first three plans the user opens: its line
  // keeps its place hidden, so nothing under the pointer moves when it shows
  const HINT_KEY = 'calendarHintViews';
  let views = 0;
  try {
    views = Number(userSetting(HINT_KEY));
  } catch {} // the settings unreadable: the hint, as for a new user
  const hintViews = Number.isFinite(views) && views > 0 ? views : 0;
  const calendarHint = div(['Click a weekday to measure on it, or skip it; double-click a day to note what ' +
    'happened, right-click it for more'], 'ui-label-details pmx-context-hint');
  calendarHint.style.visibility = 'hidden';
  if (hintViews >= 3)
    calendarHint.style.display = 'none';
  else {
    try {
      setUserSetting(HINT_KEY, String(hintViews + 1));
    } catch {} // not counted: the hint once more next time
  }
  calendarBox.addEventListener('mouseenter', () => calendarHint.style.visibility = 'visible');
  calendarBox.addEventListener('mouseleave', () => calendarHint.style.visibility = 'hidden');
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
        flash(cardOf[c === 'randomise_total' ? 'implant_total' : c]);
    }
    // only what moved: the calendar itself shows the change was planned
    const lines: HTMLElement[] = moved.length > 0 ? [span(`With these days: ${moved.join(', ')}.`)] : [];
    const a = after.armResults;
    const target = now('power_target');
    for (let i = 0; i < a.rowCount; i++) {
      // (an arm with no size has no power: its empty cell is the platform's FLOAT_NULL, a tiny number that read 0%;
      // the groups table says what it reaches with 40)
      const power = finite(a.get('assurance', i));
      if (power !== null && power < target) {
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
    setInput(name, value);
    const before = result.peek();
    void replan('calendar').then(() => {
      const after = result.peek();
      if (before && after && after !== before)
        reportChange(before, after);
    });
  };
  // each arm's place in the palette, kept for life: one added (a combination too), removed, moved or renamed
  // recolours no other
  let armsBefore: string[] = [];
  form.effect(() => {
    const names = armsInput.names.value.filter((n) => n);
    const m = meta.peek();
    const slots = slotsAfter(m.slots ?? {}, armsBefore, names);
    armsBefore = names;
    if (!sameJson(slots, m.slots ?? {}))
      meta.value = {...m, slots};
  });
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
  // the plan's arms as groups, before a save stamps them: the calendar shows their doses and offers their exceptions
  const armed = (m: StudyMeta) => {
    const r = result.peek();
    return r ? withGroups(m, ['Control', ...armsIn(r)]) : m;
  };
  const {addNote, menu: dayMenu, select: daySelect} = dayActions(() => armed(meta.peek()), userMeta,
    undefined,
    {last: () => (result.peek()?.plan.get('study_end', 0) as number | undefined) ?? 0,
      // (the days the calendar's cells mark as measured, day 0 too)
      measured: () => {
        const r = result.peek();
        return r ? String(r.plan.get('measurement_days', 0)).split(', ').map(Number) : [];
      }});
  // the calendar: the plan's days with the study's dates, doses and notes; built in its own scope, so its tooltips
  // go with the next build
  // the calendar's events, for the .ics export: randomization, the measurement days and the end
  let planEvents: {last: number, events: [number, string][]} = {last: 0, events: []};
  // the Planned reviews the plan in view was made with (none while the field is not a list of days)
  const plannedReviews = () => {
    const r = reviewDays(call.getParamValue('looks'));
    return 'days' in r ? r.days : [];
  };
  // a review day after the study's end, marked under the field itself (the plan holds: the calendar and the protocol
  // say no review is held then)
  const looksAfter = div([], 'pmx-field-note');
  whenSettled(form, () => form.getInput('looks')?.root.append(looksAfter));
  form.effect(() => {
    const r = result.value;
    const end = r ? finite(r.plan.get('study_end', 0)) : null;
    const after = end === null ? [] : plannedReviews().filter((d) => d > end);
    looksAfter.replaceChildren(...after.length ? [noteLine(`Day ${after.join(', ')} ${after.length > 1 ? 'fall' :
      'falls'} after the end (day ${end}): no review then`)] : []);
  });
  // an arm removed takes its details and dose exceptions with it: the doses always the plan's arms' (calendar,
  // protocol, the overview, the study followed)
  form.effect(() => {
    const r = result.value;
    // (not while an arm's name is typed: its details move to the new name when it is left, armsInput's renamed)
    if (!r || document.activeElement?.closest('.pmx-arm-title'))
      return;
    const m = meta.peek();
    // (the arms as they stand too: a plan computed before a rename, shown after it, still has the old name)
    const kept = onlyGroups(m, ['Control', ...armsIn(r), ...armsInput.names.peek()]);
    if (kept !== m)
      meta.value = kept;
  });
  form.effect(() => {
    const r = result.value;
    const m = meta.value;
    if (r) {
      const end = r.plan.get('study_end', 0) as number;
      planEvents = {last: end, events: [[0, 'Randomize'], ...String(r.plan.get('measurement_days', 0)).split(', ')
        .map(Number).filter((d) => d > 0 && d < end).map((d): [number, string] => [d, 'Measure tumors']),
      [end, 'End of study (measure)']]};
      calendarBox.replaceChildren(calendarGrid(r.plan, r.calendar, startDay.root, toggleWeekday,
        {meta: armed(m), onMenu: dayMenu, onNote: addNote, onSelect: daySelect, reviews: plannedReviews()}));
    }
  });
  // an offered endpoint is for the plan in view: an edit hides it (its line keeps its place, the notes their height)
  form.effect(() => {
    const hide = stale.value;
    growthNotes.querySelectorAll<HTMLElement>('.pmx-endpoint-offer')
      .forEach((e) => e.style.visibility = hide ? 'hidden' : '');
  });
  // the plan's combinations, from its declaration, with its arms' kill rates
  const combinationsIn = (r: PlanFrames) => {
    const dc = declared.get(r);
    const rates = (dc?.arms ?? []).map((x, i) => ({archetype: ARCHETYPES[x.therapy],
      killRate: finite(r.armResults.get('kill_rate_per_day', i)) ?? NaN, n: finite(r.armResults.get('n', i))}));
    return {dc, rates, found: dc ? combinationsOf(dc, rates) : []};
  };
  const detectable =signal<{of: PlanFrames, mde: Map<string, number | null>} | null>(null);
  const endSummary = div([], 'pmx-plan-end-summary');
  Object.assign(endSummary.style, {margin: '0 0 var(--dg-space-s)', lineHeight: '1.45'});
  // the endpoint offer of the plan shown: its line is held while its search runs, so the notes keep their height
  let offerOf: {r: PlanFrames, line: HTMLElement} | null = null;
  // each combination against its drugs alone, asked once the plan shows (cached: one plan asks once), first among
  // the notes as the engine wrote it, marked (the synergy finds its line by the mark, not by the text, which an arm's
  // name may echo); one the engine could not work out says so, unmarked (no synergy). The notes and the synergy wait
  // on the same promise, which puts the lines in first
  const placeSentences = (r: PlanFrames, lines: HTMLElement[]) => {
    if (result.peek() === r && lines.length && !lines[0].isConnected)
      groupNotes.prepend(...lines);
  };
  const sentencesFor = (r: PlanFrames): Promise<HTMLElement[]> => {
    let p = sentencesOf.get(r);
    if (p)
      return p;
    const {dc, rates, found} = combinationsIn(r);
    p = (!dc || !found.length ? Promise.resolve([]) : Promise.all(found.map((c) =>
      designCombination({...dc, n: combinationN(dc, rates, c) ?? null}, rates[c.a], rates[c.b], dc.arms[c.ab].metric,
        c.interaction, false).then((row) => row ? comboLine(dc.arms[c.ab].name, row) + NOTE_TIP +
        combinationNote(dc.arms[c.ab].name, row, dc.power, `${dc.arms[c.a].name} and ${dc.arms[c.b].name}`) : null,
      (e) => {
        console.error('PMX Studio:', e);
        return null;
      })))
      .then((texts) => found.map((c, i) => {
        const [line] = notes(texts[i] ?? combinationFailed(dc!.arms[c.ab].name)) as HTMLElement[];
        if (texts[i])
          line.dataset.pmxCombination = String(c.ab);
        return line;
      })))
      .then((lines) => {
        placeSentences(r, lines);
        return lines;
      });
    sentencesOf.set(r, p);
    return p;
  };
  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    const v = (c: string) => r.plan.get(c, 0);
    groupNotes.replaceChildren(...notes(v('notes_groups')));
    // the combinations' sentences, first among the notes once asked (drawn again: put back)
    void sentencesFor(r).then((lines) => placeSentences(r, lines));
    const {dc, found} = combinationsIn(r);
    growthNotes.replaceChildren(...notes(v('notes_growth')));
    // a model not in the atlas: what its history gave and what it assumes, where the user decides
    // (only the last drawing's: two drawings of one plan before the notes arrive would each add them)
    const asked = ++historyAsked;
    void historyNotesOf(String(v('line')), finite(v('study_end'))).then((ns) => {
      if (asked === historyAsked && result.peek() === r && ns.length)
        growthNotes.append(...notes(ns.join(' | ')));
    }, () => {});
    // the controls of some studies may level off below the endpoint: under the warning (the engine's, same threshold)
    // the endpoint that avoids it, found behind the plan for the plan's own declaration, applied on a click; when none
    // does (lines whose controls level off at any volume), what to plan for instead: the study's end day, reached or
    // not. Not with a last day fixed: there the controls are not followed to the endpoint
    const d = declared.get(r);
    if (d && !d.lastDay && v('p_endpoint_not_reached') > ENDPOINT_RISK) {
      // (a plan's effect runs again when its smallest detectable effects arrive: the same line, searched once)
      if (offerOf?.r !== r) {
        const line = div([' '], 'ui-label-details pmx-endpoint-offer');
        line.style.marginLeft = '16px';
        offerOf = {r, line};
        // ponytail: the calls of an offer a newer plan overtook are not cancelled, only ignored (cached by inputs)
        void safeEndpoint(d).then((x) => {
          if (offerOf?.line !== line || result.peek() !== r)
            return;
          line.replaceChildren(...(x === null ?
            [`Plan to stop on day ${v('study_end')} even if some controls are still small`] :
            [`A ${x} mm³ endpoint avoids it: `,
              hintLink('use it', () => apply('endpoint', x), `Sets the endpoint volume to ${x} mm³`)]));
        }, () => {}); // a failed search leaves the line empty: it costs that line only
      }
      offerOf.line.style.visibility = stale.peek() ? 'hidden' : '';
      growthNotes.append(offerOf.line);
    }

    // the control first, then each arm: what it expects, the animals it gets and the power they give (a table built
    // in this run is disposed with it)
    const a = r.armResults;
    const colors = colorsIn(r);
    type Row = {name: string, cells: (string | HTMLElement)[]};
    const fixed = v('n_fixed');
    // the smallest detectable effects arrive after the plan (their search takes several times as long)
    const known = detectable.value?.of === r ? detectable.value.mde : null;
    // with some arm unsized, the control as the cards count it: as large as the largest sized arm
    const rows: Row[] = [{name: 'Control', cells: [swatch('Control', colors['Control']), '—',
      animals(countOf(r.plan, 'n_control') ?? 'over 40', true),
      animals(countOf(r.plan, 'randomise_control'), false), '—',
      ...(fixed ? ['—'] : [])]}];
    // a combination by its two arms, not by the kind its growth is modelled with (the more variable of the two)
    const comboOf = new Map(dc ? found.map((c) =>
      [dc.arms[c.ab].name, `${dc.arms[c.a].name} + ${dc.arms[c.b].name}`]) : []);
    for (let i = 0; i < a.rowCount; i++) {
      const g = (c: string) => a.get(c, i);
      const metric = metricWord(g('metric'));
      // the therapy on the group's hover, under its full name: the table stays compact
      const pair = comboOf.get(g('arm'));
      const therapy = pair ? `Combination of ${pair}; growth modelled as ` +
        `${String(g('therapy')).split(' (')[0].toLowerCase()}, the more variable of the two` : String(g('therapy'));
      rows.push({name: g('arm'), cells: [shortSwatch(g('arm'), colors[g('arm')], therapy),
        `${metric} ${num(g('effect'))}%`, animals(finite(g('n')) ?? 'over 40', true),
        animals(finite(g('randomise')), false),
        powerText(finite(g('n')), finite(g('assurance')), r.powerCurve, g('arm')),
        ...(fixed ? [!known ? (stale.value || error.value !== null ? '—' : loader()) :
          known.get(g('arm')) === null ? 'out of reach' : `${metric} ${Math.round(known.get(g('arm'))!)}%`] : [])]});
    }
    // (what it is, said in the header: designer wording 2)
    const metrics = new Set(Array.from({length: a.rowCount}, (_, i) => metricWord(a.get('metric', i))));
    const detectableHeader = `Smallest ${metrics.size === 1 ? [...metrics][0] : 'effect'} detectable ` +
      `(${pct(v('power_target'))} power)`;
    const headers = ['Group', 'Expected', 'Analyzable', 'To randomize', 'Power',
      ...(fixed ? [detectableHeader] : [])];
    const headerTips: {[h: string]: string} = {Expected: 'The effect each arm is expected to have',
      Analyzable: fixed ? 'Animals per group that reach the analysis' : `Animals per group, ${SIZED_BY_SIMULATION}`,
      [detectableHeader]: 'The smallest effect these groups would detect at the declared power'};
    // what each dosed arm adds up to over the study, for its randomized animals: in the group's panel (a click on its
    // row), with "Match a total…"
    const end = v('study_end') as number;
    const doseLine = (arm: string, i: number): HTMLElement | null => {
      const d = armDoses(arm, meta.peek(), end, finite(a.get('randomise', i)), null);
      if (!d)
        return null;
      const text = span(d.text);
      text.title = d.tip;
      // a one-time help: the dose per administration that gives a total over this arm's doses (it is never
      // redistributed later: a skipped dose lowers the total)
      const match = doseForTotal(arm, meta.peek(), end, 1) === null ? [] : [div([
        hintLink('Match a total…', () => {
          const answer = ui.divText('', {style: {fontWeight: '600', marginTop: '4px'}});
          ui.tooltip.bind(answer, 'Per administration over this arm\'s doses. Write it in the arm\'s name or its ' +
            'group details to use it.');
          const total = ui.input.float('Total dose (mg/kg)', {value: undefined, nullable: true,
            onValueChanged: (x) => {
              const per = x ? doseForTotal(arm, meta.peek(), end, x) : null;
              answer.textContent = per === null ? '' : `= ${sig(per)} mg/kg per dose`;
            }});
          doseDialog?.close();
          const d = modal('total-dose', `Total dose for ${arm}`);
          if (!d)
            return;
          doseDialog = showCentred(closeNamed(d).add(ui.inputs([total])).add(answer));
          total.input.focus();
        }, 'The dose per administration that adds up to a total you give')], 'pmx-plan-match')];
      return div([text, ...match], 'pmx-plan-doses');
    };
    // a group's panel: what it expects, the animals and the power, its doses
    const groupPanel = (x: Row, i: number) => {
      if (x.name === 'Control') {
        return objectPanel({title: 'Control', color: colors['Control'],
          subtitle: 'Vehicle, on the treated groups\' days',
          facts: [{label: 'Analyzable', value: String(countOf(r.plan, 'n_control') ?? 'over 40')},
            {label: 'To randomize', value: String(countOf(r.plan, 'randomise_control') ?? '—')}]});
      }
      const at = i - 1;
      const g = (c: string) => a.get(c, at);
      const power = finite(g('assurance'));
      const target = finite(cell('power_target').peek());
      const doses = doseLine(x.name, at);
      return objectPanel({title: x.name, color: colors[x.name], subtitle: String(g('therapy')),
        verdict: {text: `Expected ${metricWord(g('metric'))} ${num(g('effect'))}%` +
          (power === null ? '' : `, power ${pct(power)}`), tone: power !== null && target !== null && power < target ?
          'warn' : 'ok', sub: power !== null && target !== null && power < target ? `Below the ${pct(target)} target` :
          undefined},
        facts: [{label: 'Analyzable', value: String(finite(g('n')) ?? 'over 40')},
          {label: 'To randomize', value: String(finite(g('randomise')) ?? '—')}],
        sections: doses ? [{title: 'Doses', content: () => doses, open: true}] : []});
    };
    const table = new BasicTable<Row>({items: rows, keyOf: (x) => x.name, columns: headers.map((header, j) => ({header,
      render: (x: Row) => x.cells[j], align: j >= 1 ? 'center' : 'left',
      ...(headerTips[header] ? {tooltip: headerTips[header]} : {})})),
    onRowClick: (x, i) => showDepth(new PmxGroup(null, x.name, colors[x.name], {details: () => groupPanel(x, i)})),
    row: (_x, _i, tr) => {
      tr.style.cursor = 'pointer';
      tr.title = 'Its doses and power in the context panel';
    }});
    table.root.style.width = '100%';
    wrapCells(table.root);
    groupsTable.replaceChildren(table.root);
    // a narrow panel (a third of the row) scrolls the table sideways, no column cut (QA round 3 designer B7)
    Object.assign(groupsTable.style, {overflowX: 'auto', maxWidth: '100%', minWidth: '0'});
  });
  // each combination's synergy, searched behind the plan (several times its cost): its sentence completed in place.
  // Its own effect: it runs again when an out-of-date plan is current again (the same result, its answers dropped)
  form.effect(() => {
    const r = result.value;
    if (!r || stale.value)
      return;
    const {dc, rates, found} = combinationsIn(r);
    // the search for each combination's smallest synergy, then that synergy's own row (the effect it asks of the
    // combination), started with its sentence, not after it: the clause completes the sentence once both are in
    const searches = (dc ? found : []).map((c) => (async () => {
      const at = {...dc!, n: combinationN(dc!, rates, c) ?? null};
      const metric = dc!.arms[c.ab].metric;
      const row = await designCombination(at, rates[c.a], rates[c.b], metric, c.interaction, true);
      const m = finite(row?.['mde_synergy']);
      const needed = m !== null && m > 0 ?
        await designCombination(at, rates[c.a], rates[c.b], metric, m, false) : null;
      await sentencesFor(r);
      if (!row || result.peek() !== r || stale.peek())
        return;
      // the engine's sentence stays as it is: the synergy is added to its end, once (no sentence, as when the engine
      // failed on it: nothing to complete)
      const line = [...groupNotes.children].find((e): e is HTMLElement => e instanceof HTMLElement &&
        e.dataset.pmxCombination === String(c.ab) && e.dataset.pmxSynergy === undefined);
      const said = line?.lastElementChild as HTMLElement | null | undefined;
      if (said)
        said.title += synergyClause(row, needed, metric);
      if (line)
        line.dataset.pmxSynergy = '1';
    })().catch((e) => console.error('PMX Studio:', e))); // a failed search leaves the sentence without its synergy
    void Promise.all([sentencesFor(r), ...searches]).finally(() => sentenced.value = r);
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
      const gain = maxOf([...new Set(rowsOf.map((x) => x.arm))].map((arm) => {
        const own = rowsOf.filter((x) => x.arm === arm);
        const atE = own.find((x) => x.day === e)?.power ?? 0;
        return Math.max(...own.filter((x) => x.day >= e).map((x) => x.power)) - atE;
      }));
      const last = rowsOf[rowsOf.length - 1];
      const removedAt = (day: number) => Math.round(100 * (rowsOf.find((x) => x.day === day)?.removed ?? 0));
      const points = Math.round(100 * gain);
      const adds = points < 1 ? 'adds no power' : `adds at most ${points} point${points === 1 ? '' : 's'} of power`;
      endSummary.textContent = `Ending after day ${e} ${adds}, while the controls removed at ` +
        `${v('humane_limit')} mm³ go from ${removedAt(e)}% on day ${e} to ${removedAt(last.day)}% on day ${last.day}.`;
    }


    // the study protocol, in the order of the ARRIVE essentials: design, sample size and its basis, animals,
    // randomization, outcome and analysis, multiplicity, end, humane limit, measurements
    const k = v('treated_arms');
    const corr = v('correction');
    const corrName = CORRECTION_NAME[correctionOf(corr)];
    // each combination is also compared with its two drugs alone
    const {dc, found} = combinationsIn(r);
    const combined = !dc ? [] : found.map((c) => ({name: dc.arms[c.ab].name, a: dc.arms[c.a].name,
      b: dc.arms[c.b].name}));
    const ar = r.armResults;
    // a combination by its own drugs, each with the dose, route and schedule its name gives (not those of the arms of
    // the drug alone, QA 0.1.3 final designer #5); its growth is modelled as the more variable of their kinds
    // (its name already gives them: "a combination", not its name again in the parenthesis, designer wording 6)
    const drugsOf = (arm: string) => combinationDrugs(parseArmName(arm)).length > 1 ? 'a combination' : null;
    // a drug the name does not place, sized as acting through the immune system (the cautious n), is still said a
    // small molecule unless the user chose that (QA 0.1.3 final designer #25)
    const kindWords = (i: number) => {
      const t = ar.get('therapy_key', i) as Therapy;
      // (the declaration keeps each arm as the form gave it, armsOf, its therapy_source with it)
      return t === 'immune_agonist' && (dc?.arms[i] as Partial<Arm> | undefined)?.therapySource !== 'user' &&
        therapyOf(ar.get('arm', i)).placed === null ? KIND_WORDS.small_molecule : KIND_WORDS[t];
    };
    const armRows = Array.from({length: ar.rowCount}, (_, i) => ({name: ar.get('arm', i) as string,
      therapy: drugsOf(ar.get('arm', i)) ?? kindWords(i),
      effect: `${metricWord(ar.get('metric', i))} ${num(ar.get('effect', i))}%`, n: finite(ar.get('n', i)),
      power: ar.get('assurance', i)}));
    const unsized = armRows.filter((x) => x.n === null).map((x) => x.name);
    // a model not in the atlas, by its readable name (its line is its key)
    const own = historyModelName(v('line'));
    // the drugs that act through an immune system this model's mice lack (the arms' cards say it too)
    const host = hostClassOf(v('line'), lineInput.models.value);
    const lacking = [...new Set(Array.from({length: ar.rowCount}, (_, i) =>
      lackingImmunity(ar.get('arm', i), ar.get('therapy_key', i) as Therapy, host)).flat())];
    const plan: {[row: string]: string} = {
      'Design': `${tumorsOf(v('line'), lineInput.models.value?.[v('line')])}; ${k + 1} groups: a control and ` +
        armRows.map((x) => `${x.name} (${x.therapy}, expected ${x.effect})`).join(', ') +
        // (a combination's comparison with its two drugs, "each at α, not corrected", returns here and in the ethics
        // text when Analyze computes it: promise and calculation travel together)
        '; each treated group compared with the control' + (!lacking.length ? '' :
        `. Warning: ${lacking.join(' and ')} ${lacking.length > 1 ? 'act' : 'acts'} through the immune system, ` +
          'which the immunodeficient mice of this model lack; little or no effect is to be expected from ' +
          (lacking.length > 1 ? 'them' : 'it')),
      'Sample size': armRows.map((x) => x.n === null ? `${x.name}: more than 40 analyzable animals per group ` +
        `would be needed for ${x.effect}` : `${x.name}: ${x.n} analyzable animals per group give ${pct(x.power)} ` +
        `power for ${x.effect}`).join('; ') + `, two-sided α = ${num(v('alpha'), 2)}` +
        // one α, said once: the family's, and with a correction the level each group is sized at
        (k < 2 || corr === 'none' ? '' : `, ${corrName}-corrected over ${k} comparisons: each group sized at α = ` +
          `${num(v('alpha_per_comparison'), 3)}`) + '. ' +
        (v('n_fixed') ? 'The animals per group were set; the power is what they give' :
          `Sized for ${pct(v('power_target'))} power (${SIZED_BY_SIMULATION})`) + `, averaged over the growth of ` +
        (own ? `${own} controls in its own history, weighed against PMX Studio's PDX reference` :
          `${v('line')} controls in ${v('atlas_papers')} published studies (historical-control prior)`) +
        `; at least ${MIN_GROUP} per group`,
      'Animals': (finite(v('n_total')) !== null ?
        `${v('implant_total')} to implant, ${v('randomise_total')} to randomize, ${v('n_total')} analyzable. ` :
        // some arm unsized: the counts of the control and the sized arms, as the cards give them
        r.plan.col('n_total_sized') && finite(v('n_total_sized')) !== null ?
          `${v('implant_total_sized')} to implant, ${v('randomise_total_sized')} to randomize, ` +
          `${v('n_total_sized')} analyzable, without ${unsized.join(', ')}: no group size up to 40 reaches the ` +
          `target power for ${unsized.length > 1 ? 'them' : 'it'}. ` :
          'Not sized: no group size up to 40 reaches the target power. ') +
        `${pct(v('pre_loss'))} expected lost before randomization (no take, out of the volume ` +
        `window) and ${pct(v('post_loss'))} after it (ulceration, deaths, technical failures); ` +
        (v('loss_rule') === 'expected' ? 'the counts cover the average loss' :
          `these numbers leave a ${pct(v('confidence'))} chance that each group keeps enough animals`) +
        (finite(v('confidence_study')) === null ? '' :
          `; in ${pct(v('confidence_study'))} of studies every group does`) +
        (k < 2 ? '' : '. The control is as large as the largest treated group'),
      'Randomization': `When the tumors reach ${v('v_rand')} mm³, into groups balanced on tumor volume; the ` +
        'randomization day is day 0',
      'Outcome and analysis': 'The growth rate of each animal’s tumor (slope of log volume over its ' +
        'measurements), compared between each treated group and the control by Welch’s t test (each group its ' +
        'own variance)',
      'Multiple comparisons': k < 2 ? 'One treated group: none needed' : corr === 'none' ?
        `None: each of the ${k} comparisons at α = ${num(v('alpha'), 2)} (exploratory)` :
        corr === 'bh' ? `${correctionOver(corr, k)}: each group is sized as for Bonferroni, at ` +
          `${num(v('alpha_per_comparison'), 3)} (the conservative bound)` :
          `${corrName} over ${k} comparisons: the chance of any false positive among them stays at α = ` +
          `${num(v('alpha'), 2)}; each group is sized at ${num(v('alpha_per_comparison'), 3)}` +
          (corr === 'holm' ? ', the level of Holm’s first step' : ''),
      'End of study': v('last_day_chosen') ? `Day ${v('study_end')}, chosen. The expected effects are the ones on ` +
        `day ${v('endpoint_day')}, when the median control reaches ${v('endpoint')} mm³` :
        `The first measurement day once the median control reaches ${v('endpoint')} mm³ (expected day ` +
        `${v('study_end')})`,
      'Humane limit': `${v('humane_limit')} mm³ (a mean diameter of about ` +
        `${diameterOf(v('humane_limit')).toFixed(1)} cm): an animal above it is removed, whatever the study day; its ` +
        'measurements up to removal stay in the analysis',
      'Measurements': `Tumor volume and body weight on ${v('measurement_weekdays')}: days ` +
        `${v('measurement_days')} from randomization`,
    };
    const reviews = plannedReviews().filter((d) => d <= v('study_end'));
    if (reviews.length) {
      plan['Planned reviews'] = `Each treated group reviewed on day ${reviews.join(', ')}, against its chance of ` +
        'success; the days are fixed before the study starts';
    }
    analysis.replaceChildren(tableFromMap(plan));
    // the text the Copy button gives: the same rows
    protocolText = [`Study protocol (${r.plan.name})`, ...Object.entries(plan).map(([k, t]) => `${k}: ${t}`)]
      .join('\n');
    const ethicsFacts = {arms: armRows.map((x) => ({name: x.name, therapy: x.therapy, effect: x.effect,
      n: x.n, power: pct(x.power)})), combinations: combined, line: own ?? v('line'),
    papers: v('atlas_papers'), history: !!own,
    target: pct(v('power_target')),
    alpha: String(num(v('alpha'), 2)), alphaPer: String(num(v('alpha_per_comparison'), 3)), nFixed: !!v('n_fixed'),
    correction: correctionOf(corr),
    minGroup: MIN_GROUP, preLoss: pct(v('pre_loss')), postLoss: pct(v('post_loss')), confidence: pct(v('confidence')),
    lossRule: String(v('loss_rule')), implant: v('implant_total'), randomise: v('randomise_total'),
    nTotal: finite(v('n_total')), humane: v('humane_limit'), diameter: diameterOf(v('humane_limit')).toFixed(1),
    lastDayChosen: !!v('last_day_chosen'), end: v('study_end'), endpoint: v('endpoint')};
    ethicsApplication = ethicsText(ethicsFacts);
    // the control is as large as the largest treated group
    const nMax = maxOf(armRows.map((x) => x.n ?? 0));
    docFacts = {line: own ?? v('line'), protocol: plan, ethics: ethicsParts(ethicsFacts), end: v('study_end'),
      vRand: v('v_rand'), humane: v('humane_limit'), diameter: ethicsFacts.diameter,
      perGroup: [{name: 'Control', n: armRows.some((x) => x.n === null) ? null : nMax},
        ...armRows.map((x) => ({name: x.name, n: x.n}))],
      // (some arm unsized: the animals the cards count, the control's and the sized arms')
      implant: countOf(r.plan, 'implant_total'), measurements: String(v('measurement_days')).split(',').length,
      therapies: armRows.map((x) => x.therapy)};
  });

  // the protocol's Copy button, and beside it the same facts in the wording of an ethics application
  const protocolCopy = copyRow(() => protocolText, 'Study protocol');
  protocolCopy.append(buttonWithIcon('Copy for an ethics application', 'copy', () => !ethicsApplication ?
    grok.shell.warning('Plan the study first') : void navigator.clipboard.writeText(ethicsApplication)
      .then(() => grok.shell.info('Text for an ethics application copied'),
        () => grok.shell.warning('Text for an ethics application: the clipboard is not available'))));

  const {row: scaleRow, log} = volumeAxis();
  // the arms in the study's colours: a colour changed in any stage of the study redraws them
  const growthColors = computed(() => {
    const r = result.value;
    if (!r?.curves.col('arm'))
      return null;
    const colors = colorsIn(r, meta.value);
    paint(r.curves, 'arm', colors);
    return colors;
  });
  const growth = new BandChart(computed(() => result.value?.curves), {x: 'day', lower: 'q05', upper: 'q95',
    line: 'median', split: 'arm', xTitle: 'Day after randomization', yTitle: 'Tumor volume', log, fit: true,
    redrawOn: growthColors});


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
    const colors = colorsIn(r);
    paint(df, 'Arm', Object.fromEntries([...label].map(([arm, l]) => [l, colors[arm]])),
      Object.fromEntries([...label].map(([arm, l]) => [l, arm])));
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
    const colors = colorsIn(r);
    const [ENDPOINT, HUMANE] = [`Controls over the endpoint (${mm3Text(r.plan.get('endpoint', 0))} mm³)`,
      `Controls removed (over ${mm3Text(r.plan.get('humane_limit', 0))} mm³)`];
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
    paint(df, 'Curve', {...Object.fromEntries(Object.entries(colors)
      .map(([arm, c]) => [`Power, ${arm}`, c])), [TARGET]: '#9aa3ad', [ENDPOINT]: colors['Control'],
    [HUMANE]: CAL.amber}, Object.fromEntries(Object.keys(colors).map((arm) => [`Power, ${arm}`, arm])));
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

  // the animals the study needs as the assumed TGI goes, always drawn: after the plan, a design per TGI (the form's
  // own inputs, the plan unchanged); a plan given TGIs through the API draws its own
  const SWEEP = '30, 40, 50, 60, 70, 80, 90';
  const sensitivity = signal<{of: PlanFrames, sc: DG.DataFrame} | undefined>(undefined);
  let sweeping: PlanFrames | null = null;
  const sweepFailed = signal<PlanFrames | null>(null);
  form.effect(() => {
    const r = result.value;
    if (!r || stale.value || busy.value || sweeping === r || sensitivity.peek()?.of === r)
      return;
    if (r.scenarios.rowCount > 0) {
      sensitivity.value = {of: r, sc: r.scenarios};
      return;
    }
    sweeping = r;
    // (with the plan's own schedules setting: the sweep's first design is then the plan's, already computed)
    void designed({tgiScenarios: SWEEP}).then((c) => {
      if (result.peek() === r)
        sensitivity.value = {of: r, sc: c.scenarios};
    }).catch((e) => {
      // the plan itself stands: the chart says it has no sizes, the run's error line stays the run's
      if (result.peek() === r)
        sweepFailed.value = r;
      grok.shell.error(`The study could not be sized for each TGI: ${e?.message ?? e}`);
    }).finally(() => {
      if (sweeping === r)
        sweeping = null;
    });
  });
  const sensitivityData = computed(() => {
    // (never undefined: the chart is built over its first frame, and a column-less placeholder fails its split)
    const s = sensitivity.value;
    const sc = s && s.of === result.value ? s : undefined;
    const own = sc?.of.plan.get('per_week', 0);
    const rows = !sc ? [] : Array.from({length: sc.sc.rowCount}, (_, i) => i)
      .filter((i) => sc.sc.get('per_week', i) === own)
      .flatMap((i) => [['To implant', finite(sc.sc.get('implant_total', i))], ['Analyzable',
        finite(sc.sc.get('n_total', i))]].filter(([, n]) => n !== null).map(([what, n]) => ({
        'Assumed TGI': sc.sc.get('tgi', i), 'Animals': n, 'Animals of the study': what})));
    const df = rows.length === 0 ? DG.DataFrame.fromColumns([DG.Column.fromList('double', 'Assumed TGI', []),
      DG.Column.fromList('double', 'Animals', []), DG.Column.fromList('string', 'Animals of the study', [])]) :
      DG.DataFrame.fromObjects(rows)!;
    df.col('Assumed TGI')!.meta.units = '%';
    df.col('Animals of the study')!.meta.colors.setCategorical({'To implant': '#2083d5', 'Analyzable': '#5b6b7f'});
    return df;
  });
  const sensitivityChart = viewers.lineChart(sensitivityData, {xColumnName: 'Assumed TGI', yColumnNames: ['Animals'],
    splitColumnNames: ['Animals of the study'], legendVisibility: 'Always', showXAxis: true, showYAxis: true,
    markerType: 'circle', showMarkers: 'Always', showTopPanel: false, showAggrSelectors: false,
    showSplitSelector: false} as any);
  const sweepNote = div([], 'pmx-plan-sweep-note');
  form.effect(() => {
    const s = sensitivity.value;
    const r = result.value;
    sweepNote.textContent = !r ? '' : sweepFailed.value === r ? 'The sizes for each TGI could not be computed.' :
      !s || s.of !== r ? 'Sizing the study for each TGI…' : 'The whole study with ' +
      'every arm at the TGI on the axis, sized for the target power: a smaller effect costs animals fast. Scenarios ' +
      'to read across, not predictions.';
  });
  const scenarioEnd = div([], 'pmx-plan-scenario-end');
  const scenarioBox = div([]);
  form.effect(() => {
    const has = !!result.value?.scenarios.rowCount;
    scenarioEnd.hidden = !has;
    scenarioBox.hidden = !has;
  });
  // what each folded panel holds, in one line
  const calendarSummary = computed(() => {
    const r = result.value;
    if (!r)
      return '';
    const days = form.getInput('measurementDays')?.value;
    return [days ? `Measured ${days}` : '', `${r.plan.get('per_week', 0)} measurements a week`,
      `ends on day ${r.plan.get('study_end', 0)}`].filter((x) => x).join(' · ');
  });
  const scenariosSummary = computed(() => {
    const df = sensitivityData.value;
    const rows = Array.from({length: df.rowCount}, (_, i) => i)
      .filter((i) => df.get('Animals of the study', i) === 'To implant');
    if (rows.length < 2)
      return result.value ? 'Sizing the study for each TGI…' : '';
    const at = (i: number) => `${df.get('Animals', i)} to implant at TGI ${df.get('Assumed TGI', i)}%`;
    return `${at(rows[0])}, ${at(rows[rows.length - 1])}`;
  });
  // the figures answer the plan's question; each panel, one closer question
  const panels = panelGrid([
    {id: 'groups', title: 'Experimental groups', hint: 'The animals each group needs for its comparison with the ' +
      'control, and the power that gives. Click a group for its doses', content: [groupsTable, groupNotes]},
    {id: 'growth', title: 'Predicted tumor growth', hint: 'Median and 90% range of each group\u2019s tumor volume, ' +
      'to the end of the study', content: [growthNotes, scaleRow, sized(growth, 380)]},
    {id: 'calendar', title: 'Study calendar', hint: 'Week by week from randomization. Click a weekday to measure on ' +
      'it, or skip it; double-click a day to note what happened, right-click it for more', content: [calendarLoader,
      calendarStatus, calendarHint, calendarBox], summary: calendarSummary},
    {id: 'scenarios', title: 'Effect scenarios', hint: 'What the study would need if the treatment had each TGI: ' +
      'scenarios you read across, not predictions. The animals of the whole study with every arm at that TGI, ' +
      'sized for the target power whatever animals per arm is set', content: [
      sized(sensitivityChart, 280), sweepNote, scenarioEnd, scenarioBox], summary: scenariosSummary},
    {id: 'whatif', title: 'Group size, end of study and schedule', hint: 'More or fewer animals, an earlier or later ' +
      'end, other measuring days', content: [
      block('Power by group size', 'With more or fewer animals', sized(power, 280)),
      block('Power by end of study', 'Ending earlier or later', endSummary, sized(byEnd, 280)),
      block('Measurement schedule', 'Measuring more or less often', scheduleBox)],
    summary: computed(() => result.value ? 'power with more or fewer animals, an earlier or later end, other ' +
      'measuring days' : '')},
    {id: 'protocol', title: 'Study protocol', hint: 'Design, sample size and analysis, ready to paste; also in ' +
      'File › Documents', content: [protocolCopy, analysis],
    summary: computed(() => result.value ? 'design, sample size and analysis, ready to paste' : '')},
  ], PLAN_LAYOUT, 'pmx-plan-panels', 'plan', form);
  exports.push(['Customize this page', () => customizePage(panels), 'View']);
  // the study's colours: a colour changed in a chart is the plan's, and its frames follow it
  const link: StudyLink = {meta: () => meta.peek(), setMeta: userMeta, colors: () => {
    const r = result.peek();
    return armColors(r ? armsIn(r) : [], meta.peek().colors, meta.peek().slots);
  }, plan: () => planKey(initial)};
  self.own(linkStudy(link));
  form.effect(() => {
    const r = result.value;
    if (r)
      follow(link, Object.values(r).filter((x): x is DG.DataFrame => x instanceof DG.DataFrame));
  });
  exports.push(['Reset colors', () => resetColors(link), 'View',
    () => hasOwnColors(link) ? null : 'The colors are the automatic ones']);
  // an empty viewer says what columns it lacks: the charts hidden until the first plan gives them their data
  form.effect(() => [growth, power].forEach((v) => v.root.style.visibility = result.value ? '' : 'hidden'));
  // (with no group size, no power by end: the chart would keep the last plan's)
  form.effect(() => byEnd.root.style.visibility = result.value?.duration.rowCount ? '' : 'hidden');
  // (nor the sizes by TGI before their sweep: the chart would say its columns are missing)
  form.effect(() => sensitivityChart.root.style.visibility = sensitivityData.value.rowCount > 0 ? '' : 'hidden');
  // the schedules are the design's costliest scan: computed once the first plan shows, wherever the page is scrolled,
  // and from then on (QA: computed only when scrolled to, the section said "Comparing…" while nothing ran)
  form.effect(() => {
    if (result.value && !scheduleOpen.peek())
      scheduleOpen.value = true;
  });

  form.effect(() => {
    const r = result.value;
    if (!r?.scenarios.rowCount)
      return;
    scenarioEnd.textContent = `The study ends on day ${r.scenarios.get('study_end', 0)} in every scenario.`;
    scenarioBox.replaceChildren(scenariosTable(r.scenarios, r.plan.get('per_week', 0)));
  });
  // a plan made with another atlas says so, once, above its numbers
  const atlasLine = div([], 'ui-label-details');
  atlasLine.hidden = true;
  void atlasNote(initial?.provenance, PLAN_ATLAS_NOTE)
    .then((t) => t && (atlasLine.hidden = !(atlasLine.textContent = t)));
  const content = divV([atlasLine, cards, panels], 'pmx-plan-content');
  Object.assign(content.style, {gap: 'var(--dg-space-l)'});

  // a calendar change: its own loader, the calendar dimmed, the page in view; a plan from the form clears its note
  form.effect(() => {
    calendarLoader.hidden = !refreshing.value;
    calendarBox.style.opacity = refreshing.value ? '0.5' : '';
    if (planning.value)
      calendarStatus.replaceChildren();
  });

  // no treated arm yet: what to add, in place of a plan
  const noArms = emptyState('Add a treatment arm: the animals, the end of the study, the calendar and the power ' +
    'follow from it.', {icon: 'flask', actions: [buttonWithIcon('Add the first arm', 'plus', () => armsInput.add())]});
  noArms.style.marginTop = '96px';

  // the first plan under u2's loader, every later one under the thin bar
  const panes = stagePanes(form, content, {busy: planning, result, error, stale, before: [noArms],
    busyText: 'Designing the study…', staleText: 'the plan shown is for the last complete inputs'});
  // the form's column opens no wider than its fields use (at 2560 px a fifth left ~870 px of it empty: QA 0.1.3 final
  // designer #24): its default width only, once the page is first laid out; a sash the user drags stays where put
  const widths = new ResizeObserver(() => {
    const w = panes.root.clientWidth;
    if (w <= 0)
      return;
    widths.disconnect();
    if (panes.sizes.peek()[0] * w > DECLARATION_MAX)
      panes.sizes.value = [DECLARATION_MAX / w, 1 - DECLARATION_MAX / w];
  });
  widths.observe(panes.root);
  form.scope.own(() => widths.disconnect());
  // no arm: the empty state alone, not empty cards and a white chart under it. A box of its own hides the plan
  // (stagePanes shows `content` again whenever the error clears, as a new plan's first replan does: QA 0.1.3 final
  // designer #27); `display: contents` leaves the layout as it was
  const withArms = div([]);
  withArms.style.display = 'contents';
  content.replaceWith(withArms);
  withArms.append(content);
  form.effect(() => {
    const none = armsInput.count.value === 0;
    noArms.style.display = none ? '' : 'none';
    withArms.style.display = none ? 'none' : 'contents';
  });

  // the plan carries the schedules while their section is open; opened on a plan without them, the plan is made
  // again with them, by the same run as every other plan: a newer plan, a change made meanwhile or a failure are
  // then handled as for any plan, and the section waits only while that run is on (it once had a computation of its
  // own beside the plan's, whose answer could be dropped or never come, and its loader stayed)
  whenSettled(form, () => form.effect(() => call.setParamValue('compareSchedules', scheduleOpen.value)));
  // the arms typed under the other metric converted to the same treated volume (TGI 60 is T/C 46 from 100 to
  // 1000 mm³) before the cards read it: the plan stays the same
  whenSettled(form, () => form.effect(() => {
    const m = form.getInput('effectAs')?.value === 'T/C' ? 'TC' : 'TGI';
    const num = (name: string, x: number) => Number(form.getInput(name)?.peek()) || x;
    // (it sets effectAs, the cards' metric, with the arms: ArmsInput.sayAs)
    armsInput.sayAs(m, num('vRand', 100), num('endpoint', 1000));
  }));
  // the volumes a TGI is bounded by (complete regression, tgiMax), as the form holds them
  whenSettled(form, () => form.effect(() => {
    const num = (name: string, x: number) => Number(form.getInput(name)?.value) || x;
    planVolumes.value = {v0: num('vRand', 100), endpoint: num('endpoint', 1000)};
  }));
  let remade: PlanFrames | null = null;
  form.effect(() => {
    const r = result.value;
    if (!scheduleOpen.value || schedules.value || !r || stale.value || busy.value || remade === r ||
      sentenced.value !== r)
      return;
    remade = r; // once per plan: a run that fails says so, and is not tried again in a loop
    void replan('schedules');
  });

  // with the animals per arm set, the plan shows first and each arm's smallest detectable effect fills in after it,
  // from the form's own inputs (the plan unchanged)
  let searching: PlanFrames | null = null;
  form.effect(() => {
    const r = result.value;
    if (!r || !r.plan.get('n_fixed', 0) || stale.value || busy.value || searching === r ||
      detectable.peek()?.of === r || sentenced.value !== r)
      return;
    searching = r;
    const fill = designed({compareSchedules: false, detectableEffect: true, tgiScenarios: null});
    void fill.then((c) => {
      const a = c.armResults;
      detectable.value = {of: r, mde: new Map(Array.from({length: a.rowCount}, (_, i) =>
        [a.get('arm', i) as string, finite(a.get('mde', i))]))};
    }).catch((e) => error.value = inWords(String(e?.message ?? e))).finally(() => {
      if (searching === r)
        searching = null;
    });
  });

  // the view opens at once; the default declaration (or the saved plan's) is planned behind the loader, once the
  // form has written its defaults into the call
  whenSettled(form, () => {
    call.setParamValue('detectableEffect', false);
    // a new plan opens on the tumor model picked last, with no arm: not the table in view, which the platform gives
    // a table parameter by default (Monitor's study would read as 267 arms)
    const last = initial ? null : lastLine();
    const lab = labHumane();
    // a plan saved before TGI / T/C was the study's: its arms' own
    const saved = initial?.inputs;
    const savedArms = saved?.arms as Arm[] | null | undefined;
    const effects = saved && saved.effectAs === undefined && savedArms?.length ?
      {effectAs: metricWord(savedArms[0].metric)} : {};
    // only what the function takes: setParamValue throws on any other name, and the inputs after it went unset (a
    // plan saved by a build whose pmxPlanStudy had "looks", opened where the registered function lacks it)
    const params = new Set(call.func.inputs.map((p) => p.name));
    for (const [name, value] of Object.entries(saved ? {...saved, ...effects} :
      {arms: null, effectAs: 'TGI', ...last ? {line: last} : {}, ...lab ? {humane: lab} : {}})) {
      if (!params.has(name)) {
        console.warn(`PMX Studio: the plan's "${name}" is not an input of ${call.func.name}: not set`);
        continue;
      }
      setInput(name, name === 'arms' && value ? armsFrame(value as SavedPlan['inputs']['arms']) : value);
    }
    auto.cancel();
    // a plan saved with its tables, opened on the same inputs and computed by this app, atlas and engine: shown at
    // once as saved; otherwise (or edited meanwhile) computed, the busy overlay up, as any plan
    const shown = initial?.predicted?.shown;
    const kept = shown && sameAsked(shown.asked, initial!.inputs) ? framesOf(shown) : null;
    if (!kept)
      return replan();
    return engineFingerprint().then((e) => e === shown!.engine, () => false).then((same) => {
      if (!same || stale.peek() || result.peek() !== undefined)
        return replan();
      declared.set(kept, shown!.declared);
      error.value = null;
      shownAsked = askedOf(call);
      result.value = kept;
      planning.value = false;
    });
  });

  // the plan in view as a file: saved to the user's plans (shared from Files like any file) or downloaded, to send
  // (the plan's predictive checks are computed then, not on every replan: planOf)
  // what in the form stops it from planning: the plan in view then never catches up, so nothing waits for it (two arms
  // of one name, a field out of its range or empty)
  const formProblem = computed(() => {
    const [missing, invalid] = [form.missingRequired.value, form.invalidFields.value];
    return armsInput.validity.value ?? (missing.length ? `Fill in ${missing.join(', ')}` :
      invalid.length ? `Correct ${invalid.join(', ')}` : null);
  });
  const ready = (): PlanFrames | null => {
    const r = result.peek();
    const problem = formProblem.peek();
    // (nor while it is being computed: a save then would keep tables the page is replacing)
    if (!r || stale.peek() || problem || busy.peek()) {
      grok.shell.warning(problem ?? (!r ? 'The plan is not ready yet' : busy.peek() ?
        'The plan is being updated: save it in a moment' :
        'The plan does not match the form yet: wait for it to update'));
      return null;
    }
    return r;
  };
  // the plan's name once saved or opened from a file: what the header, the next Save and Monitor call it
  const named = signal<string | null>(initial?.name ?? null);
  const defaultName = (r: PlanFrames) => named.value ?? `${modelName(r.plan.get('line', 0))} plan`;
  // an ID, and each group's details with its empty cells read from its name
  // (with domain storage the ID is the lab's counter's, taken once, when the plan is first saved or downloaded)
  const stamped = async (r: PlanFrames) => meta.value =
    await withCode(withGroups(meta.peek(), ['Control', ...armsIn(r)]));
  // (with the plan's tables as shown: the file reopens on them, not recomputed)
  const planFor = async (name: string, r: PlanFrames): Promise<SavedPlan> => {
    const [p, d] = [await planOf(name, call, r, await stamped(r)), declared.get(r)];
    // (the inputs as the file keeps them: what it reopens with, so the same text then)
    const engine = d ? await engineFingerprint().catch(() => null) : null;
    return d && engine ? {...p, predicted: {...p.predicted, shown: shownOf(r, JSON.stringify(p.inputs), d, engine)}} :
      p;
  };
  planFiles.set(self, (name) => planFor(name, ready()!));
  // with domain storage: the plan's row (opened from it, or saved to it); a save updates it at its version
  let planRow = rowOfStart(initial);
  // opened from someone else's file (Open dialog, drop, a file in the Overview): its first save is the user's copy
  let fromFile = fileOwnerOf(initial);
  const saveRow = async (name: string, r: PlanFrames) => {
    const plan = await planFor(name, r);
    // (one balloon a save: a new row's note says it was saved)
    let said = false;
    try {
      const saved = await save('plan', plan, planRow, {fileOwner: fromFile});
      planRow = saved;
      linkRow(self.view, 'plan', saved.id);
      fromFile = undefined;
      if (saved.meta) {
        grok.shell.info(newRowNote('plan', saved, meta.peek().id));
        said = true;
        meta.value = saved.meta;
      }
    } catch (e) {
      // a plan shared with the user to read: their version goes to a copy of their own, if they say so
      if (!planRow || !isReadOnly(e))
        throw e;
      const from = planRow;
      return offerCopy(async () => {
        const {saved, doc} = await saveCopy('plan', plan, from);
        planRow = saved;
        linkRow(self.view, 'plan', saved.id);
        meta.value = doc.study ?? meta.peek();
        named.value = doc.name;
        edited = false;
        savedInputs = inputsNow();
        epoch.value++;
        formChanges.fire();
        grok.shell.info(`Saved as your copy: ${esc(doc.name)}`);
      }, e);
    }
    named.value = name;
    edited = false;
    savedInputs = inputsNow();
    epoch.value++;
    formChanges.fire();
    if (!said)
      grok.shell.info('Plan saved');
  };
  // a save asked while the plan is being made waits for it (the button says so), then saves; a plan that fails says
  // it was not saved
  const pendingSave = signal<((r: PlanFrames) => void) | null>(null);
  const whenReady = (act: (r: PlanFrames) => void) => {
    if (busy.peek() || stale.peek() && !formProblem.peek()) {
      pendingSave.value = act;
      return;
    }
    const r = ready();
    if (r)
      act(r);
  };
  form.effect(() => {
    const act = pendingSave.value;
    const failed = error.value;
    if (!act || busy.value || stale.value && !failed)
      return;
    pendingSave.value = null;
    // (outside the effect: the save reads signals of its own)
    queueMicrotask(() => {
      if (failed)
        return void grok.shell.error(`The plan could not be computed, so it was not saved: ${failed}`);
      const r = ready();
      if (r)
        act(r);
    });
  });
  const saveDialog = () => whenReady((r) => {
    // a row saved before keeps its name: Save writes it at once (Save as… is the dialog, for another name)
    if (domainStorageOn() && planRow && named.peek() !== null)
      return void saveRow(named.peek()!, r).catch((e) => grok.shell.error(saveError(e, 'plan'), {autoHide: false}));
    saveAsDialog(r);
  });
  const saveAsDialog = (r: PlanFrames) => {
    if (domainStorageOn()) {
      showOnce(self, 'save-plan', () => {
        const name = new TextInput({label: 'Name', value: defaultName(r), nullable: false});
        const note = span('Saved in PMX Studio, yours until you share it; the Overview lists it.');
        Object.assign(note.style, {color: 'var(--dg-text-color-muted, #68727f)', maxWidth: '360px'});
        // Save as… under the plan's own name replaces it only when Save is pressed again (the dialog says so)
        let replacing: string | null = null;
        return new Dialog('Save plan').add(name).add(note).onOK(async () => {
          const chosen = name.value.trim() || defaultName(r);
          if (chosen.length > NAME_MAX)
            throw new Error(`A name of at most ${NAME_MAX} characters`);
          if (planRow && chosen === named.peek() && chosen !== replacing) {
            replacing = chosen;
            throw new Error(`This plan is saved as "${chosen}": press Save again to replace it, or give another name ` +
              'for a new plan');
          }
          // under another name it is another study, with an ID and a row of its own (names may repeat)
          if (named.peek() !== null && chosen !== named.peek()) {
            meta.value = {...meta.peek(), id: ''};
            planRow = null;
          }
          await saveRow(chosen, r).catch((e) => {
            throw new Error(saveError(e, 'plan'));
          });
        }, {text: 'Save'});
      });
      return;
    }
    const where = span(`Saved to ${shownPath(plansFolder())}. Share it from Files; the Overview lists it.`);
    Object.assign(where.style, {color: 'var(--dg-text-color-muted, #68727f)', maxWidth: '360px', lineHeight: '1.45'});
    // a plan of the same name is replaced only when Save is pressed again for it (the dialog says so and stays open)
    let replacing: string | null = null;
    showOnce(self, 'save-plan', () => {
      const name = new TextInput({label: 'Name', value: defaultName(r), nullable: false});
      return new Dialog('Save plan').add(name).add(where).onOK(async () => {
        const chosen = name.value.trim() || defaultName(r);
        if (chosen.length > NAME_MAX)
          throw new Error(`A name of at most ${NAME_MAX} characters`);
        if (chosen !== replacing && await planSaved(chosen)) {
          replacing = chosen;
          throw new Error(`A plan named "${chosen}" is already saved: press Save again to replace it, or rename it`);
        }
        // saved under another name, it is another study, with an ID of its own: the Overview then lists both
        if (named.peek() !== null && chosen !== named.peek())
          meta.value = {...meta.peek(), id: ''};
        const path = await savePlan(await planFor(chosen, r));
        named.value = chosen;
        grok.shell.info(`Plan saved: ${shownPath(path)}`);
      }, {text: 'Save'});
    });
  };
  const saveButton = self.build(() => buttonWithIcon(computed(() => pendingSave.value ?
    'Saving when the plan is ready…' : 'Save plan'), 'save', saveDialog));
  planSaves.set(self, saveButton);
  // Ctrl+S, and closing with the form changed since the last save asks first (a plan never saved: Save opens its
  // name dialog, and the view stays)
  if (domainStorageOn()) {
    queueMicrotask(() => saveGate(self, 'plan', async () => {
      if (busy.peek() || stale.peek() && !formProblem.peek()) {
        saveDialog();
        throw new Error('saved once the plan is ready');
      }
      const r = ready();
      if (!r)
        throw new Error('not ready');
      if (!planRow || named.peek() === null) {
        saveAsDialog(r);
        throw new Error('named first');
      }
      await saveRow(named.peek()!, r).catch((e) => {
        grok.shell.error(saveError(e, 'plan'), {autoHide: false});
        throw e;
      });
    }, () => edited && inputsNow() !== savedInputs, undefined, formChanges));
  }
  saveButton.title = domainStorageOn() ? 'Save the plan in PMX Studio (Ctrl+S), to compare the study with it in ' +
    'Analyze' : 'Save the plan to your files, to compare the study with it in Analyze';
  if (domainStorageOn())
    exports.unshift(['Save plan as…', () => whenReady(saveAsDialog), 'Study']);
  const saving = signal(false);
  // the plan in view as a saved plan, handed to `to` (Download, Follow); the ribbon's buttons wait meanwhile
  const withPlan = (r: PlanFrames, to: (p: SavedPlan) => unknown) => {
    saving.value = true;
    void planFor(defaultName(r), r).then(to).catch((e) => grok.shell.error(String(e?.message ?? e)))
      .finally(() => saving.value = false);
  };
  // the plan as a file, from Export (as Monitor's study file)
  const downloadFile = () => {
    const r = ready();
    if (r)
      withPlan(r, downloadPlan);
  };
  const copyProtocol = () => {
    if (!ready())
      return;
    return copyText(protocolText, 'Study protocol');
  };
  // under a document's title: the study's ID, the plan it comes from and today's date
  const docSub = (plan: string) => docSubtitle(meta.peek().id, `From the plan ${plan}`);
  exports.push(['Report (PDF)…', async () => {
    const r = ready();
    // (a folded panel is left out until opened: reportHtml)
    if (r) {
      printSheet(await reportHtml({title: defaultName(r), stage: 'Plan', meta: meta.peek(), content,
        groups: ['Control', ...armsIn(r)]}));
    }
  }], ['Plan file (.pmxplan)', downloadFile],
  ['Calendar (.ics)', () => {
    const r = ready();
    if (!r)
      return;
    const name = defaultName(r);
    download(`${name}.ics`, icsOf(name, meta.peek(), planEvents.last, planEvents.events), 'text/calendar');
  }],
  // the documents of an ethics application: each opens its checklist, what the plan fills and what is left, with
  // Download (.docx) beside Copy; a checklist reopens, as it stands now, once the study's details it opened close
  ['Study plan (ARRIVE)…', function studyPlan() {
    const r = ready();
    if (r && docFacts) {
      ethicsDocument(defaultName(r), 'Study plan', studyPlanForm({...docFacts, title: defaultName(r),
        meta: meta.peek()}), STUDY_PLAN_CREDIT, {guide: 'ARRIVE study plan', url: ARRIVE_STUDY_PLAN_URL,
        asker: 'The study plan', reopen: studyPlan});
    }
  }, ETHICS],
  // (in English, as the app is: its items are English whatever the labels; Spanish is the Labels choice)
  ['Non-technical summary…', function summary(lang: NtsLanguage = 'en') {
    const r = ready();
    if (!r || !docFacts)
      return;
    // the template's labels in English or Spanish (the official Spanish wording); the plan's text stays English
    const labels = ui.input.choice('Labels', {value: lang === 'es' ? 'Español' : 'English',
      items: ['English', 'Español'], nullable: false, onValueChanged: (x) => {
        d.close();
        summary(x === 'Español' ? 'es' : 'en');
      }});
    const d = ethicsDocument(defaultName(r), ntsTitle(lang), ntsForm({...docFacts, title: defaultName(r),
      meta: meta.peek()}, lang), NTS_CREDIT, {guide: 'EU template (Decision 2020/569)', url: NTS_URL,
      asker: 'The template', reopen: () => summary(lang), options: labels.root});
  }, ETHICS],
  // what the plan settles before the study; the rest is the facility's and the ethics committee's
  ['PREPARE checklist…', function prepare() {
    const r = ready();
    if (r) {
      const items = prepareItems({line: historyModelName(r.plan.get('line', 0)) ?? r.plan.get('line', 0),
        protocol: protocolText,
        humaneLimit: r.plan.get('humane_limit', 0), lastDay: r.plan.get('study_end', 0),
        toImplant: finite(r.plan.get('implant_total', 0)), meta: meta.peek()});
      const plan = defaultName(r);
      const sub = docSub(plan);
      prepareDialog(items, {details: () => head.details(prepare), protocol: copyProtocol},
        async () => download(`${plan} - PREPARE checklist.docx`, await docx(checklistBlocks({title: 'PREPARE ' +
          'checklist', sub, items, credit: PREPARE_CREDIT, last: ['Notes', 'notes']}),
        {title: 'PREPARE checklist', landscape: true}), DOCX_TYPE));
    }
  }, ETHICS],
  ['Copy the study protocol', copyProtocol, ETHICS]);
  /** A document's checklist dialog, Download (.docx) beside Copy; `reopen` opens it again once the study's details
   *  it opened close. */
  const ethicsDocument = (plan: string, title: string, sections: FormSection[], credit: string,
    o: {guide: string, url: string, asker: string, reopen: () => void, options?: HTMLElement}) => {
    const sub = docSub(plan);
    return checklistDialog({title, guide: o.guide, url: o.url, asker: o.asker, options: o.options,
      download: async () => download(`${plan} - ${title}.docx`,
        await docx(formBlocks({title, sub, sections, credit}), {title}), DOCX_TYPE)},
    formItems(sections), {details: () => head.details(o.reopen)});
  };
  // the way on: the page's one primary action
  // (a plan already followed opens its study: QA 0.1.3 round 4 designer D8)
  const followed = signal(false);
  // (one name, the button's and the Study menu's and the toolbox's)
  const followText = computed(() => followed.value ? 'Open the followed study' : 'Follow this study');
  const followButton = buttonWithIcon(followText, 'arrow-right', () => void followStudy(), {iconRight: true});
  // changes not saved asked about first, before the plan is taken (D2: the plan handed to Monitor was the one before
  // the save, with the old ID): saved, the plan as saved goes on; discarded, the plan as last saved (D7: the edited
  // one went on, and Plan showed it again)
  const followStudy = async () => {
    const r = ready();
    if (!r)
      return;
    doseDialog?.close();
    doseDialog = null;
    const ask = leaveAsks.get(self.view.dart);
    if (dirtyViews.get(self.view.dart)?.() && ask && !await ask(false, 'before following it'))
      return;
    const saved = discardedViews.has(self.view.dart) ? await lastSaved().catch(() => null) : null;
    if (saved)
      return void followPlan(saved).catch((e) => grok.shell.error(String(e?.message ?? e)));
    withPlan(r, followPlan);
  };
  // the plan as last saved: its row, else the file it was opened from (null: never saved)
  const lastSaved = async (): Promise<SavedPlan | null> => planRow ? (await openSaved('plan', planRow.id)).doc :
    initial && 'format' in initial ? initial as SavedPlan : null;
  followButton.classList.add('pmx-go');
  // (whole up to its longest name, 250 px in studio.css: at 172 px it read "PEN THE FOLLOWED STUDY", cut at its start;
  // narrower, it ends in an ellipsis; inline: the package CSS is capped, QA 0.1.3 browser check)
  // (its label, before the icon: iconRight)
  const label = followButton.firstElementChild as HTMLElement | null;
  if (label)
    Object.assign(label.style, {minWidth: '0', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap'});
  // (in the Study menu too: on a narrow window the button gives way first, QA 0.1.3 designer D15)
  exports.push([() => followText.peek(), () => followButton.click(), 'Study', () => reasonOf(followButton)]);
  // a study followed from this plan, open or saved (by the plan's ID)
  // (looked up again as the plan gets its ID and after each save; a saved study looked up once per ID)
  let looked = '';
  form.effect(() => {
    void epoch.value;
    const id = meta.value.id;
    const plan = {name: '', inputs: {arms: []}, study: meta.value} as unknown as SavedPlan;
    followed.value = !!id && grok.shell.tables.some((t) => followedFrom(t, plan));
    if (id && !followed.peek() && domainStorageOn() && looked !== id) {
      looked = id;
      void studyRowsOf('study', id).then((rows) => followed.value ||= rows.length > 0 && meta.peek().id === id,
        () => {});
    }
  });
  followButton.title = 'Start following this study from randomization: Monitor opens with the plan\'s groups and ' +
    'calendar';
  form.effect(() => {
    const off = armsInput.count.value === 0 ? 'Add a treatment arm first' : formProblem.value ??
      (!result.value ? 'The plan is not ready yet' : busy.value ? 'The plan is being updated' :
        saving.value ? 'The plan is being prepared' : null);
    Reason.apply(followButton, off);
    // (Save is taken while the plan is updated: it waits for it, and says so)
    const waiting = pendingSave.value !== null;
    Reason.apply(saveButton, waiting ? 'Saving when the plan is ready…' :
      off === 'The plan is being updated' ? null : off);
  });

  // another stage shown in this view's place: the plan as it stands (its inputs, edited or not), to come back to
  keepOnLeaving(self, async () => {
    const r = result.peek();
    if (!r || stale.peek() || busy.peek())
      return initial ? () => planView(initial) : null;
    const p = await planOf(named.peek() ?? defaultName(r), call, r, meta.peek());
    const d = declared.get(r);
    const engine = d ? await engineFingerprint().catch(() => null) : null;
    const start: PlanStart = d && engine ? {...p, predicted: {...p.predicted, shown: shownOf(r,
      JSON.stringify(p.inputs), d, engine)}} : p;
    const row = planRow;
    return () => planView({...start, row});
  });

  const head = studyHeader(form, 'plan', computed(() => {
    const r = result.value;
    if (!r)
      return undefined;
    const colors = colorsIn(r);
    return {name: defaultName(r), groups: ['Control', ...armsIn(r)].map((name) => ({name, color: colors[name]})),
      status: `${modelName(r.plan.get('line', 0))} · planned · ${armsIn(r).length + 1} groups · ends on day ` +
        `${r.plan.get('study_end', 0)}`, meta: () => meta.value, setMeta: userMeta};
  }), (s) => void openStage(s), computed(() => {
    // a new plan with no arm yet: nothing to follow or analyze (unless a stage left behind has a study)
    const none = armsInput.count.value === 0;
    return {monitor: none && !canReturn('monitor') ? 'Follow a plan, or open a study with measurements, first' : null,
      analyze: none && !canReturn('analyze') ? 'Needs a study with measurements: follow the plan, then measure' :
        null};
  }));
  // Open… here opens a plan, or a table as a model's growth history (the model then picked here, as + Add your own
  // model does)
  const at: OpenAt = {stage: 'plan', history: (t) => historyDialog((key) => lineInput.value = key,
    {rows: t.read.rows, from: t.file})};
  const links = stageLinks(at, () => head.details());
  // Save plan as… is File's Save as…
  const saveAsAt = exports.findIndex((x) => x[0] === 'Save plan as…');
  if (saveAsAt >= 0)
    exports[saveAsAt] = ['Save as…', exports[saveAsAt][1], 'File', () => reasonOf(saveButton)];
  let rebuildMenu: (() => void) | null = null;
  queueMicrotask(() => rebuildMenu = stageMenu(self.view, exports, {newPlan: () => void planView().then((v) =>
    showStage('plan', v)), open: () => void openStudyDialog(undefined, at),
  save: {text: () => saveButton.textContent?.trim() || 'Save plan', run: () => saveButton.click(),
    why: () => reasonOf(saveButton)}}, links, (d) => self.own(d)));
  const toolbox = stageToolbox(exports, links.study, links.help, panels);
  // (the Study menu and the toolbox named as the button is, as it changes)
  form.effect(() => {
    void followText.value;
    rebuildMenu?.();
    relabel(toolbox);
  });
  return {
    root: takesFiles(stagePage(head.row, panes), at),
    ribbon: [[head.steps], [followButton]],
    toolbox,
    path: '/plan',
    status: computed(() => busy.value ? 'Designing the study…' : result.value ? result.value.plan.name : ''),
  };
}
