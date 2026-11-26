/* The Monitor view, built as Analyze is: the form generated from pmxMonitorStudy's signature, the results from its
   outputs. The groups are shown, unless the study is blinded (its details): then they stay hidden until someone
   reveals them, and that goes in the study's log. A study's measurements
   are an ordinary table (paste into its grid): Monitor only reads it, and saving the study is a gesture of its own.
   Tabs by question: who needs attention today, how the groups will be formed, how the study is growing, and whether to
   go on, only on the reviews fixed at the start (MONITOR_DESIGN.md). */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {badge, BasicTable, button, buttonWithIcon, computed, div, divH, divV, signal,
  span} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';
import {funcForm, ViewHost} from '@datagrok-libraries/u2/dg';
import type {TableInput, ViewContent} from '@datagrok-libraries/u2/dg';

import {analyzeView} from './analyze-view';
import {BandChart} from './band-chart';
import {weekdaysOf} from './frames';
import {isFiniteNumber} from './engine';
import {StudyLink, groupColors, hasOwnColors, linkStudy, paint, resetColors} from './colors';
import {allocateGroups, balanceNote, day0Weekday, MonitorFrames, monitorFrames, MonitorInputs, latestReadings,
  outsideBand, naturalCompare, promptsByChance, weightMissing} from './monitor';
import {FROZEN_ATLAS_NOTE, PLAN_TAG, SavedPlan, atlasNote, parsePlan, planKeyOf, shownPath} from './plan-file';
import {planInput} from './plan-vs-observed';
import {openStage, showStage} from './stages';
import {LIMITS_TAG, LINE_TAG, LOG_TAG, LOOKS_TAG, REVEALED, VOLUME_TAG, WEIGHT_TAG, addLog, groupsShown, logOf,
  newStudyTable, parseStudy, EXCLUDED, EXPOSURE_TAG, exposureOfTag, stamped, lineOf, setHistoryOn,
  restoreFields, roleColumns, saveStudy, studyNameOf, studyOf, tableOf, ticked} from './study-file';
import {TWIN_NOTE, savedHistory} from './history-models';
import {markLine, modal} from './ui';
import {WeightRow, groupWeights, readSchedule} from './schedule';
import {editLayout, panelTabs, Tab} from './panels';
import {stagePage, studyHeader} from './study-header';
import {openButton, takesFiles} from './hub-view';
import {animalMarks, animalPanel as panelOfAnimal, blindColumn, emptyStudy, openStudyTable, studyDataLine, studyMenu,
  studyOf as studyNamesOf,
  keepStudyField, tagStudy, withoutExcluded} from './study-table';
import {EmptyState, NoControl} from './crux';
import {META_TAG, StudyMeta, armDoses, dateOf, dateText, decisionEntry, emptyMeta, metaOf, renamedGroups, setMetaOn,
  stopsOf, withGroups, withId} from './study-meta';
import {renamedExposure} from './exposure';
import {Saved, domainStorageOn, isReadOnly, markTableFromFile, newRowNote, ownerOf, rowOfTable, save, saveCopy,
  saveError, setRowOn, tableFileOwner, withCode} from './store';
import {offerCopy, saveGate, tableSavers} from './native';
import {CAL, DayLook, dayActions, dayToday, icsOf, mark, marked, startMismatch, weeksGrid} from './calendar';
import {LineInput} from './line-input';
import {esc, printSheet, sheetHtml} from './sheet';
import {Exports, download, exportMenu, reportHtml} from './report';
import {studyInputsOf} from './explore';
import {rowsOf, tolerabilityTab} from './tolerability';
import {findCombinations, parseArmName} from './arm-name';
import {linePrior, modelName} from './priors';
import {IndependenceBand, bandWhy, bandWords, independenceBand, notCombined} from './synergy';
import {askedOf, autoRun, block, callInputs, cardGrid, dayText, exampleLine, infoNote, onTableEdits, TABLE_EDITS,
  dialogLine, isBodyWeight, okNamed, showCentred,
  keepNamedColumns, noteLine, notes, notesBox, pct, plainError, settleForm, sized, markStage, signed, stageOf,
  statCard, swatch, volumeAxis, wrapCells, stagePanes, whenSettled, showDepth, percentRange} from './ui';

/** Monitor's tabs and their panels, in order and width (twelfths): today's work and how the groups grow; the
 *  measurements and the groups formed; how the animals bear it; whether each group is worth continuing. */
const MONITOR_TABS: Tab[] = [
  {name: 'Today', layout: [{id: 'growth', span: 12}, {id: 'today', span: 12}, {id: 'calendar', span: 12}]},
  {name: 'Measurement sheet', layout: [{id: 'sheet', span: 12}, {id: 'enrollment', span: 12}]},
  {name: 'Welfare', layout: [{id: 'welfare', span: 12}]},
  {name: 'Decisions', layout: [{id: 'decisions', span: 12}]}];

/** Below this chance of success, on a planned review, Monitor suggests stopping the arm. A person decides. */
export const FUTILITY = 0.1;
/** An arm whose TGI interval starts at this (the effect well clear of none) has a clear result already in sight. */
export const CLEAR_TGI = 50;

/** The "chance of a clear result" in words. Only the low end is a figure: on the atlas's real arms a low chance held,
 *  a high one was too hopeful (MONITOR_DESIGN.md §8), so no high figure is shown. Above the low end it reads "not
 *  low", except for an arm whose TGI interval already starts far from no effect (`tgiLower` at CLEAR_TGI or more),
 *  where "not low" would read as doubt about a clearly effective arm: "High", a word from its interval, not the
 *  model's figure. */
export function chanceWords(success: number, tgiLower: number): string {
  if (success < FUTILITY)
    return `Low (under ${pct(FUTILITY)})`;
  return tgiLower >= CLEAR_TGI ? `High (TGI at least ${Math.round(tgiLower)}%)` : `Not low (${pct(FUTILITY)} or more)`;
}

interface AnimalRow {[c: string]: any}

/** The details as a blinded study shows them before its groups are revealed: each group's doses under "Group 1",
 *  "Group 2"… (the calendar and its .ics say when to dose, not what) */
const unnamed = (m: StudyMeta, revealed: boolean): StudyMeta => {
  if (revealed)
    return m;
  const name = blindName(m);
  const known = (x: {group: string}) => m.groups.some((g) => g.group === x.group);
  // the schedules a name gives (no schedule typed) read before the name goes, one per distinct schedule; the
  // exceptions and doses given under the same "Group k" (those of a group no longer in the study left out)
  return {...m, groups: m.groups.flatMap((g) => [...new Set<string>(g.schedule.trim() ? [g.schedule] :
    parseArmName(g.group).map((p) => p.schedule ?? ''))].map((schedule) => ({...g, group: name(g.group),
    schedule, compound: '', dose: null}))),
  exceptions: m.exceptions.filter(known).map((e) => ({...e, group: name(e.group)})),
  given: m.given.filter(known).map((g) => ({...g, group: name(g.group)})),
  ...m.stopped ? {stopped: Object.fromEntries(Object.entries(m.stopped).map(([g, d]) => [name(g), d]))} : {}};
};
/** A group as a blinded study names it before it is revealed: "Group k", numbered by group, not by row (a
 *  combination has a row per drug). */
const blindName = (m: StudyMeta) => (group: string) =>
  `Group ${[...new Set(m.groups.map((g) => g.group))].indexOf(group) + 1}`;
const grams = (x: unknown) => isFiniteNumber(x) ? x.toFixed(1) : '—';
const mm3 = (x: unknown) => isFiniteNumber(x) ? `${Math.round(x)}` : '—';
function fill<R>(host: HTMLElement, items: R[], key: (r: R) => string,
  columns: [string, (r: R) => string | HTMLElement][], onRowClick?: (r: R) => void): void {
  const t = new BasicTable<R>({items, keyOf: key, columns: columns.map(([header, render], j) =>
    ({header, render, align: j === 0 || header === 'Group' ? 'left' : 'center'})), onRowClick});
  if (onRowClick)
    t.root.style.cursor = 'pointer';
  t.root.style.width = '100%';
  wrapCells(t.root);
  host.replaceChildren(t.root);
}

/** The Monitor view; `initial`: the table of a study to follow. */
/** What an exported file of the study is named after: the study (its plan's name), else its ID, else its table; with
 *  no character a file name cannot hold. */
export function exportName(df: DG.DataFrame): string {
  const named = df.getTag(PLAN_TAG) ? studyNameOf(df) : '';
  return (named || metaOf(df)?.id || df.name || 'study').replace(/[\\/:*?"<>|]/g, '-').trim();
}

/** The tags a study is saved from (studyOf), besides its measurements. */
const SAVED_TAGS = [PLAN_TAG, META_TAG, LIMITS_TAG, LINE_TAG, VOLUME_TAG, EXPOSURE_TAG, LOOKS_TAG, LOG_TAG];

/** A file opened in place of a study just followed from Plan (still empty) is that study: it takes the plan, the
 *  details and the reviews the followed table carried, unless it carries a plan of its own. True when it took them. */
export function takeOverFollowed(followed: DG.DataFrame, opened: DG.DataFrame): boolean {
  if (opened === followed || followed.rowCount > 0 || !followed.getTag(PLAN_TAG) || opened.getTag(PLAN_TAG))
    return false;
  for (const t of SAVED_TAGS) {
    const v = followed.getTag(t);
    if (v && !opened.getTag(t))
      opened.setTag(t, v);
  }
  return true;
}

export async function monitorView(initial?: DG.DataFrame): Promise<DG.ViewBase> {
  const exports: Exports = [];
  const host = ViewHost.create({name: 'PMX Studio', content: (self) => monitorContent(self, initial,
    exports)});
  exportMenu(host.view, exports);
  return markStage(host.view, 'monitor');
}

/** Each study table's groups and their animals as they stood when Monitor first read it or last saved it: a group
 *  renamed since is found by its animals (study-meta.ts renamedGroups). */
const groupsBefore = new WeakMap<DG.DataFrame, Map<string, string[]>>();

/** What the view shows, built under its host, which owns it all (the form, the charts, the panes). */
function monitorContent(self: ViewHost, initial: DG.DataFrame | undefined, exports: Exports): ViewContent {
  const call = DG.Func.byName('PmxStudio:pmxMonitorStudy').prepare();
  const result = signal<MonitorFrames | undefined>(undefined);
  const error = signal<string | null>(null);
  // the error is an empty state (nothing measured yet), not a failure
  const empty = signal(false);
  const stale = signal(false);
  const busy = signal(false);
  // the log's version: what it records (a reveal, the blinding switched) changes what is shown
  const logVersion = signal(0);
  // the study table's values: a flag or an exclusion changes what the lists show without a new forecast
  const dataVersion = signal(0);
  let latest = 0;
  // no Run: each change forecasts again once the changes stop
  let preparing = false;
  const auto = autoRun(() => form, () => preparing ? undefined : run());
  const form = funcForm(call, {layout: 'tall', overrides: {line: {input: new LineInput({label: 'Tumor model'})}},
    onInputChanged: () => {
      stale.value = result.peek() !== undefined;
      auto.changed();
    }});
  settleForm(form, {control: 'Auto (the vehicle group)', evaluationDay: 'Auto (endpoint reached)'});
  keepNamedColumns(form, {bodyWeight: isBodyWeight});

  const studyTable = (): DG.DataFrame | null => {
    const v = call.getParamValue('study');
    return v instanceof DG.DataFrame ? v : null;
  };
  // the plan's text: the one picked in the form, else the study table's own
  const planText = () => call.getParamValue('plan') || studyTable()?.getTag(PLAN_TAG);
  // the study as it is saved: its table with the form's limits
  const savedStudy = (df: DG.DataFrame) =>
    studyOf(df, call.getParamValue('humane') ?? null, call.getParamValue('weightLimit') ?? null);
  // the groups, unless the study is blinded and no one has revealed them
  const revealed = computed(() => {
    void result.value;
    void logVersion.value;
    const df = studyTable();
    return !df || groupsShown(df);
  });
  // the Group input is not the user's to change while the groups are hidden
  form.effect(() => {
    const field = form.getInput('group');
    if (field)
      field.disabled = revealed.value ? null : 'Hidden until the groups are revealed';
  });
  const lockLooks = () => {
    const field = form.getInput('looks');
    if (field)
      field.disabled = studyTable()?.getTag(LOOKS_TAG) ? 'The reviews were fixed when this study started' : null;
  };
  // what goes with the table: its body weight column and the protocol's limits, for the study file and for Analyze
  const remember = () => {
    const df = studyTable();
    if (!df)
      return;
    const weight = call.getParamValue('bodyWeight');
    df.setTag(WEIGHT_TAG, weight instanceof DG.Column ? weight.name : typeof weight === 'string' ? weight : '');
    df.setTag(LIMITS_TAG, JSON.stringify({humane: call.getParamValue('humane') ?? null,
      weightLimit: call.getParamValue('weightLimit') ?? null}));
    // an empty field leaves the model the study was saved with
    const line = call.getParamValue('line');
    if (typeof line === 'string' && line) {
      df.setTag(LINE_TAG, line);
      setHistoryOn(df, line);
    }
  };
  const planState = {stale, error, result: computed(() => result.value && {analysis: result.value.status}),
    reading: Promise.resolve() as Promise<unknown>, changed: () => auto.changed()};
  const run = async (prepare?: () => Promise<void>) => {
    busy.value = true;
    const mine = ++latest;
    try {
      // a study the view fills in changes the form as it goes: none of those changes runs by itself
      preparing = prepare !== undefined;
      await prepare?.();
      preparing = false;
      auto.cancel();
      await planState.reading;
      // a measurements file opened in place of the study just followed: it is that study, with its plan
      const opened = studyTable();
      if (initial && opened && opened !== initial && takeOverFollowed(initial, opened))
        lockLooks();
      // the reviews fixed at the start hold, whatever the field says
      const fixed = studyTable()?.getTag(LOOKS_TAG);
      if (fixed)
        call.setParamValue('looks', fixed);
      // a blinded study's group column is hidden from the form's list, which would drop it: given directly
      const df0 = studyTable();
      if (df0 && !groupsShown(df0))
        call.setParamValue('group', roleColumns(df0).group);
      const sent = askedOf(call);
      tagStudy(studyInputsOf(call));
      const frames = await monitorFrames(withoutExcluded(callInputs<MonitorInputs>(call)));
      if (mine !== latest)
        return;
      error.value = null;
      empty.value = false;
      stale.value = askedOf(call) !== sent;
      // the reviews stay with the study from its first run
      const df = studyTable();
      if (df && !df.getTag(LOOKS_TAG))
        df.setTag(LOOKS_TAG, String(call.getParamValue('looks') ?? ''));
      lockLooks();
      remember();
      result.value = frames;
    } catch (e: any) {
      // only the last run speaks: an older one failing late would write over a newer result's place
      if (mine === latest) {
        empty.value = e instanceof EmptyState;
        error.value = e instanceof EmptyState ? e.message : e instanceof NoControl ? (revealed.peek() ? e.message :
          `There is no group "${e.control}" in the study table: check its spelling.`) :
          plainError(e, 'The forecast could not be computed. ' +
          'Try again; if it fails again, reload the page.');
        // what is shown is no longer what the table or the form says
        stale.value = result.peek() !== undefined;
      }
    } finally {
      if (mine === latest) {
        preparing = false;
        busy.value = false;
      }
    }
  };

  // the form reads the group column by name; after a rename it follows
  const followGroup = (df: DG.DataFrame) => {
    const input = form.getInput('group');
    let name: string | null = null;
    try {
      name = roleColumns(df).group.name;
    } catch {
      // no group column: nothing to follow
    }
    if (input && name && input.peek() !== name)
      input.value = name;
  };

  const setStudy = async (df: DG.DataFrame) => {
    // a blinded study opens with its group column hidden (`~group`), before any view of the table exists
    blindColumn(df);
    grok.shell.addTableView(df);
    // back to Monitor once the table is open; a view not yet shown is left to showStage: making it current would add
    // it, and showStage would add it a second time
    if ([...grok.shell.views].some((v) => v.dart === self.view.dart))
      grok.shell.v = self.view;
    const study = form.getInput('study');
    if (study)
      study.value = df.name;
    const looks = df.getTag(LOOKS_TAG);
    const field = form.getInput('looks');
    if (looks && field)
      field.value = looks;
    await form.settled;
    // the form's table change resets the call's other fields as it settles: the reviews go in after it
    if (looks) {
      if (field)
        field.value = looks;
      call.setParamValue('looks', looks);
    }
    restoreFields(form, df);
    followGroup(df);
    await form.settled;
    lockLooks();
  };
  const newStudy = () => run(async () => {
    const df = emptyStudy();
    await setStudy(df);
    tagStudy({study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: 'body_weight'});
    grok.shell.v = [...grok.shell.tableViews].find((v) => v.dataFrame.id === df.id) ?? grok.shell.v;
  });
  whenSettled(form, () => {
    keepStudyField(form.getInput('study') as unknown as TableInput, () => call.getParamValue('study') ?? null);
    const study = form.getInput('study') as unknown as {actions?: {id: string, icon: string, tooltip: string,
      run: () => void}[]} | null;
    if (study?.actions) {
      study.actions = [...study.actions, {id: 'new', icon: 'plus',
        tooltip: 'Start an empty study table, to type or paste the measurements into', run: () => void newStudy()}];
    }
    if (initial)
      void run(() => setStudy(initial));
  });
  planInput(form, call, planState);

  const start = divV([
    span('Open the study\'s measurements (one row per animal and day, day 0 the randomization), check the columns ' +
      'and the forecasts follow. Paste or type each session\'s values into the table: each is checked as you enter ' +
      'it, and the forecasts update.'),
    exampleLine(),
  ], 'pmx-monitor-start');
  Object.assign(start.style, {gap: 'var(--dg-space-m)', maxWidth: '760px'});

  const status = computed(() => result.value?.status);
  const cell = <T = number>(name: string) => computed(() => {
    const s = status.value;
    return (s && !s.col(name)?.isNone(0) ? s.get(name, 0) ?? undefined : undefined) as T | undefined;
  });
  const card = statCard.bind(null, form);
  const metaVersion = signal(0);
  const meta = computed(() => {
    void result.value;
    void metaVersion.value;
    void logVersion.value; // the arms stopped are read from the log
    return metaOf(studyTable()) ?? emptyMeta();
  });
  // a day's card says its date once day 0 has one, unless the table is measured past today's study day (the date the
  // calendar says does not match: a date from it would be one the data contradicts); the expected end, how it
  // stands against the plan's
  const dated = (label: string, day: ReadonlySignal<number | undefined>, also?: () => string) => computed(() => {
    const d = day.value;
    const today = dayToday(meta.value.start);
    const latest = status.value?.get('day', 0) as number | undefined;
    const date = d === undefined || (today !== null && latest !== undefined && today < latest) ? null :
      dateOf(meta.value.start, d);
    return [label, date ? dateText(date) : '', also?.() ?? ''].filter((x) => x).join(' · ');
  });
  const [needing, nextDay, endMedian, plannedEnd] = [cell('attention'), cell('next_day'), cell('end_median'),
    cell('planned_end')];
  const versusPlan = () => {
    const [e, p] = [endMedian.value, plannedEnd.value];
    if (e === undefined || p === undefined)
      return '';
    const d = Math.round(e - p);
    return d === 0 ? 'as planned' : `${Math.abs(d)} ${Math.abs(d) === 1 ? 'day' : 'days'} ` +
      `${d > 0 ? 'later' : 'earlier'} than planned`;
  };
  const attentionCard = card('Animals needing attention', needing,
    () => 'Animals listed under Today, each with its reason');
  form.effect(() => attentionCard.classList.toggle('pmx-stat-attn', (needing.value ?? 0) > 0));
  const cards = cardGrid([
    attentionCard,
    card(dated('Next measurement', nextDay), nextDay, () => 'The next day at the study\'s usual interval', dayText),
    card(dated('Expected end', endMedian, versusPlan), endMedian, () => {
      const [a, b] = [cell('end_lower').value, cell('end_upper').value];
      // not "8 in 10": halfway through, the atlas's MC38 studies fell in this range about half the time (§8)
      return `The day the median control reaches the endpoint, likely between ${dayText(a)} and ` +
        `${b === undefined ? 'later than the days forecast' : dayText(b)}`;
    }, dayText),
    card(dated('End in the plan', plannedEnd), plannedEnd,
      () => 'The last day the plan predicted (empty: no plan)', dayText),
  ], 'pmx-monitor-cards', 150);

  // Today
  const attention = div([], 'pmx-monitor-attention');
  const sheet = div([], 'pmx-monitor-sheet');
  // an animal's menu on its row (its name, the first cell), in the attention list and the measurement sheet alike
  for (const list of [attention, sheet]) {
    list.addEventListener('contextmenu', (e) => {
      const row = (e.target as HTMLElement).closest('tbody tr');
      const name = (row?.querySelector('[data-animal]') as HTMLElement | null)?.dataset.animal ??
        row?.querySelector('td')?.textContent?.trim();
      if (!name)
        return;
      e.preventDefault();
      studyMenu(studyNamesOf(studyTable()), 'animal', name, {causedBy: e});
    });
  }
  const noteBox = notesBox();
  // each value as it goes in, against what the last run expected: no run needed to hear "measure again"
  const entered = div([], 'pmx-monitor-entered');
  const checkEntered = () => {
    const [r, df] = [result.peek(), studyTable()];
    // the columns as the last run read them (the table's own names, whatever they are)
    let cols: ReturnType<typeof roleColumns> | null = null;
    try {
      cols = df ? roleColumns(df) : null;
    } catch {
      // not a study yet
    }
    if (!r || !df || !cols) {
      entered.replaceChildren();
      return;
    }
    const {animal, day, volume} = cols;
    const nextDay = r.status.get('next_day', 0);
    const off = outsideBand(Array.from({length: df.rowCount}, (_, i) => ({animal: String(animal.get(i)),
      day: day.get(i), volume: volume.isNone(i) ? NaN : volume.get(i)})), nextDay,
    new Map(rowsOf(r.animals).map((a) => [a.animal, a as any])));
    entered.replaceChildren(...off.length === 0 ? [] : [noteLine(`Measure again (${dayText(nextDay)}): ` +
      off.map((o) => `${o.animal} at ${mm3(o.volume)} mm³, expected ${mm3(o.lower)}–${mm3(o.upper)}`).join('; '),
    'warning')]);
  };
  let watching: {unsubscribe(): void}[] = [];
  // edits come in bursts (a paste, values typed one after another): the checks and the lists redone once after them
  let edits: ReturnType<typeof setTimeout> | undefined;
  const edited = () => {
    clearTimeout(edits);
    edits = setTimeout(() => {
      checkEntered();
      dataVersion.value++;
    }, 200);
  };
  const stopWatching = () => {
    clearTimeout(edits);
    watching.forEach((w) => w.unsubscribe());
    watching = [];
  };
  form.effect(() => {
    void result.value;
    stopWatching();
    const df = studyTable();
    if (df)
      watching = [df.onDataChanged.subscribe(edited), df.onRowsAdded.subscribe(edited)];
    checkEntered();
  });
  const groupCols = (): [string, (r: AnimalRow) => string][] => revealed.value ? [['Group', (a) => a.group]] : [];
  // an animal clicked: its curve, readings, weight and marks in the context panel (its group and colour only once
  // the groups are revealed), cleared with the view
  let shownAnimal: HTMLElement | null = null;
  self.own(() => {
    if (shownAnimal && grok.shell.o === shownAnimal)
      grok.shell.o = null;
  });
  const showAnimal = (a: AnimalRow) => {
    const n = studyNamesOf(studyTable());
    const r = result.peek();
    if (!n || !r)
      return;
    const shown = revealed.peek();
    shownAnimal = showDepth(panelOfAnimal(n, String(a.animal), {showGroup: shown,
      color: shown ? colorsOf(r)[String(a.group)] : undefined}));
  };
  // only for a study with body weights
  const weightCols = (r: MonitorFrames): [string, (a: AnimalRow) => string | HTMLElement][] =>
    r.animals.col('weight_change') ?
      [['Weight change', (a) => {
        // a gain says so (+4%), as a loss does (−4%)
        const el = span(typeof a.weight_change === 'number' ?
          (a.weight_change >= 0.5 ? '+' : '') + signed(a.weight_change) : '—');
        if (typeof a.last_weight === 'number')
          ui.tooltip.bind(el, `Last weighed: ${grams(a.last_weight)} g`);
        return el;
      }]] : [];
  form.effect(() => {
    const r = result.value;
    void dataVersion.value;
    if (!r)
      return;
    // an animal flagged or excluded says so by its name, as the study table colours its rows
    const n = studyNamesOf(studyTable());
    const marks = n ? animalMarks(n) : new Map<string, string>();
    const animalCell = (a: AnimalRow) => {
      const m = marks.get(String(a.animal));
      const el = span(String(a.animal));
      if (m) {
        el.append(' ', badge(m === 'excluded' ? 'Excluded' : m === 'partly' ? 'Some values excluded' : 'Flagged',
          {variant: m === 'flagged' ? 'warning' : 'default'}));
      }
      el.dataset.animal = String(a.animal);
      return el;
    };
    // by animal: the engine's order, group by group, would show the groups
    const all = rowsOf(r.animals).sort((a, b) => naturalCompare(String(a.animal), String(b.animal)));
    const listed = all.filter((a) => a.attention);
    if (listed.length === 0)
      attention.replaceChildren(noteLine('No animal needs attention today.', 'success'));
    else {
      fill(attention, listed, (a) => a.animal, [['Animal', animalCell], ...groupCols(),
        ['Why', (a) => a.attention], ['Last volume', (a) => `${mm3(a.last_volume)} mm³, day ${a.last_day}`],
        ...weightCols(r),
        ['Next expected', (a) => a.next_median === null ? '—' : `${mm3(a.next_lower)}–${mm3(a.next_upper)} mm³`]],
      showAnimal);
    }
    fill(sheet, all, (a) => a.animal, [['Animal', animalCell], ...groupCols(),
      ['Last volume, mm³', (a) => mm3(a.last_volume)], ...weightCols(r),
      ['Next expected, mm³', (a) => a.next_median === null ? '—' : `${mm3(a.next_lower)}–${mm3(a.next_upper)}`],
      ['Measure again outside, mm³',
        (a) => a.next_median === null ? '—' : `${mm3(a.check_lower)}–${mm3(a.check_upper)}`]], showAnimal);
    // a prompt to measure again, not a verdict: said with what it costs by chance (no correction over animals and
    // sessions, which would let a mistyped value through)
    const chance: number | null = r.status.col('check_chance')?.isNone(0) === false ?
      r.status.get('check_chance', 0) : null;
    const above: number = chance && r.status.col('check_above')?.isNone(0) === false ?
      r.status.get('check_above', 0) : 0.5;
    const side = above > 2 / 3 ? ', mostly above it' : above < 1 / 3 ? ', mostly below it' : '';
    // one line; what it costs by chance on hover
    const again = infoNote('A value outside its range: measure again; a second reading that agrees stands.');
    ui.tooltip.bind(again, chance ? `In the atlas's published controls a correct value fell outside its range about ` +
      `1 time in ${Math.round(1 / chance)} (the average over studies)${side}: with ${all.length} animals, ` +
      `${promptsByChance(all.length, chance)} by chance.` :
      `How often a correct value falls outside its range is not measured for ${r.status.get('line', 0)} yet.`);
    sheet.append(again);
    noteBox.replaceChildren(...notes(r.status.get('notes', 0)));
  });

  // Enrollment
  const allocationTable = div([], 'pmx-monitor-allocation');
  const propose = async () => {
    const df = studyTable();
    if (!df)
      throw new Error('Open the study\'s measurements first');
    const text = df.getTag(PLAN_TAG);
    if (!text)
      throw new Error('Propose groups needs the study\'s plan, for its groups and animals: add it in the Plan field.');
    const plan: SavedPlan = parsePlan(text);
    const [animal, day, volume] = ['animal', 'day', 'volume'].map((role) => roleCol(df, role));
    if (!animal || !day || !volume)
      throw new Error('Propose groups needs the animals, their days and volumes: pick those columns in the form');
    // each animal's latest volume
    const wName = df.getTag(WEIGHT_TAG);
    const weightCol = roleCol(df, 'bodyWeight') ?? (wName ? df.col(wName) : null);
    const latest = latestReadings(Array.from({length: df.rowCount}, (_, i) => ({animal: String(animal.get(i)),
      day: day.get(i), volume: volume.isNone(i) ? NaN : volume.get(i),
      weight: weightCol && !weightCol.isNone(i) ? Number(weightCol.get(i)) : null})));
    const groups = [{name: 'Vehicle', n: plan.predicted.nControl ?? 0},
      ...plan.predicted.arms.map((a) => ({name: a.name, n: a.n ?? 0}))];
    // what the engine would refuse, said as what to do (a retry or a reload changes neither)
    const unsized = groups.filter((g) => !(g.n > 0)).map((g) => g.name);
    if (unsized.length) {
      throw new Error(`Every group needs at least one animal, and the plan gives none to ${unsized.join(', ')}: size ` +
        `${unsized.length > 1 ? 'them' : 'it'} in Plan, save the plan and follow it again.`);
    }
    if (latest.animal.length === 0)
      throw new Error('No animal has a tumor volume yet: enter each animal\'s volume, then propose the groups.');
    // a new seed per study: it decides who gets what
    const seed = crypto.getRandomValues(new Uint32Array(1))[0];
    const {allocation, balance, weightUsed} = await allocateGroups(latest.animal, latest.volume, groups, seed,
      weightCol ? latest.weight : undefined);
    const note = balanceNote(!!weightCol, weightMissing(latest.volume, latest.weight), weightUsed);
    const group = new Map(rowsOf(allocation).map((r) => [r.animal, r.group]));
    const col = roleCol(df, 'group') ?? df.columns.addNewString(groupsShown(df) ? 'group' : '~group');
    for (let i = 0; i < df.rowCount; i++)
      col.set(i, group.get(String(animal.get(i))) ?? '');
    record(`Groups allocated (seed ${seed})${weightUsed ? '; balanced on body weight too' : ''}`);
    fill(allocationTable, rowsOf(balance), (b) => b.group, [['Group', (b) => b.group], ['Animals', (b) => String(b.n)],
      ['Geometric mean volume', (b) => `${mm3(b.geo_mean)} mm³`],
      ['Volume spread (CV)', (b) => `${Math.round(100 * Math.sqrt(Math.exp(b.sd_log * b.sd_log) - 1))}%`],
      ...weightUsed ? [['Body weight, g (mean ± sd)',
        (b: AnimalRow) => `${grams(b.weight_mean)} ± ${grams(b.weight_sd)}`] as [string, (b: AnimalRow) => string]] :
        []]);
    if (note)
      allocationTable.append(ui.divText(note));
  };
  // an action's failure is said where it happened and goes away (the run's error line stays the run's)
  const failed = (what: string) => (e: any) =>
    grok.shell.error(plainError(e, `${what}. Try again; if it fails again, reload the page.`));
  // groups already in the table are someone's allocation, perhaps the one the animals are being dosed by, so they
  // are replaced only when asked twice
  const proposeAsked = () => {
    const df = studyTable();
    const groups = df ? roleCol(df, 'group') : null;
    const allocated = !!df && !!groups &&
      Array.from({length: df.rowCount}, (_, i) => i).some((i) => !groups.isNone(i) && String(groups.get(i)).trim());
    if (!allocated)
      return void propose().catch(failed('The groups could not be proposed'));
    const d = modal('replace-groups', 'Replace the groups?');
    if (!d)
      return;
    d.add(dialogLine('Every animal gets a new random group', 'From its latest volume, in place of the one it has. ' +
        'The study\'s log keeps the new seed, not the old groups.'))
      .onOK(() => void propose().catch(failed('The groups could not be proposed')));
    okNamed(d, 'Replace');
    showCentred(d);
  };
  const proposeButton = buttonWithIcon('Propose groups', 'random', proposeAsked);
  proposeButton.title = 'Fill the table\'s group column with random groups alike in their tumors';
  // the animal room's sheet, from the groups in the table and the seed the log kept
  const sheetButton = buttonWithIcon('Print randomization sheet', 'print', () => {
    const df = studyTable();
    if (!df || !roleCol(df, 'group')) {
      grok.shell.warning('The sheet needs the study\'s groups: propose them, or open a table with a group column.');
      return;
    }
    const seed = [...logOf(df)].reverse().map((e) => /seed (\d+)/.exec(e.what)?.[1]).find((s) => s) ?? null;
    try {
      printSheet(sheetHtml(df.name, metaOf(df), df, seed));
      addLog(df, 'Randomization sheet printed (it shows the groups)');
      logVersion.value++;
    } catch (e: any) {
      failed('The sheet could not be printed')(e);
    }
  });
  sheetButton.title = 'A page for the animal room: each animal, its group and its volume at randomization, ' +
    'what each group receives, and lines to sign. It shows the groups.';

  // Growth: by group, once revealed
  const growthAxis = volumeAxis();
  const colorsOf = (r: MonitorFrames) => {
    const names = [...new Set(r.days.col('group')!.toList() as string[])];
    const m = metaOf(studyTable());
    const colors = groupColors(names, r.status.get('control', 0), m?.colors, m?.slots);
    return Object.fromEntries(names.map((n) => [n, colors[n]]));
  };
  // each combination found in the groups' names (unless the plan unticked it), its drugs started together: where A+B
  // would run if its two drugs acted independently, dashed beside it, and whether it slows growth more or less on the
  // last day measured (synergy.ts independenceBand); computed after each reading, shown when it arrives
  const independent = signal(new Map<string, {band: IndependenceBand, a: string, b: string}>());
  let bandedFor: MonitorFrames | undefined;
  form.effect(() => {
    const r = result.value;
    const s = studyInputsOf(call);
    if (!revealed.value || r === bandedFor)
      return;
    const mine = bandedFor = r;
    if (independent.peek().size)
      independent.value = new Map();
    if (!r || !s)
      return;
    const control = r.status.get('control', 0) as string;
    const treated = [...new Set(r.days.col('group')!.toList() as string[])].filter((g) => g !== control);
    const off = notCombined(planText());
    // the band reads A and B alone as they act inside A+B: not for drugs started on other days (a staggered design)
    const starts = (g: string) => parseArmName(g).map((p) => readSchedule(p.schedule)?.first ?? 0);
    const combos = findCombinations(treated).filter((c) => !off.includes(treated[c.ab]) &&
      new Set([c.a, c.b, c.ab].flatMap((i) => starts(treated[i]))).size === 1);
    if (!combos.length)
      return;
    void linePrior(r.status.get('line', 0)).then(async (lp) => {
      const logBeta = lp.study_prior?.mean_log_beta;
      if (logBeta === undefined)
        return;
      const found = new Map<string, {band: IndependenceBand, a: string, b: string}>();
      for (const c of combos) {
        const [a, b, ab] = [treated[c.a], treated[c.b], treated[c.ab]];
        const band = await independenceBand(s, {control, a, b, ab}, Math.exp(logBeta), r.status.get('day', 0))
          .catch(() => null);
        if (band)
          found.set(ab, {band, a, b});
      }
      if (bandedFor === mine && found.size)
        independent.value = found;
    }).catch(() => {});
  });
  const independentName = (ab: string) => `${ab} if independent`;
  // the dashed categories of the growth chart: the bands of independence (BandChart reads them each time it draws)
  const dashed: string[] = [];
  const growthData = computed(() => {
    const r = result.value;
    if (!r || !revealed.value)
      return undefined;
    // the bands of independence as more groups, in their combination's colour (a copy: the frame stays Monitor's)
    const bands = independent.value;
    const t = bands.size ? r.days.clone() : r.days;
    dashed.splice(0, dashed.length, ...[...bands.keys()].map(independentName));
    for (const [ab, {band: x}] of bands) {
      x.days.forEach((d, j) => {
        const i = t.rowCount;
        t.rows.addNew();
        t.set('group', i, independentName(ab));
        t.set('day', i, d);
        t.set('lower', i, Math.exp(x.base + x.lower[j]));
        t.set('median', i, Math.exp(x.base + x.mid[j]));
        t.set('upper', i, Math.exp(x.base + x.upper[j]));
      });
    }
    const colors = colorsOf(r);
    for (const ab of bands.keys())
      colors[independentName(ab)] = colors[ab];
    paint(t, 'group', colors);
    // what the study is read against: the endpoint, today and the expected end (the chart's scale follows them)
    const lines = t.meta.formulaLines;
    lines.clear();
    const endpoint = Number(call.getParamValue('endpoint'));
    if (endpoint > 0)
      lines.addLine({title: `Endpoint ${endpoint} mm³`, formula: `\${median} = ${endpoint}`});
    // today as the calendar has it: from day 0's date, else the last day measured; none before day 0
    const today = dayToday(meta.value.start) ?? r.status.get('day', 0) as number;
    if (today >= 0)
      lines.addLine({title: `Today, day ${today}`, formula: `\${day} = ${today}`, color: '#2083d5'});
    const end = r.status.get('end_median', 0);
    if (Number.isFinite(end) && end !== null)
      lines.addLine({title: `Expected end, day ${end}`, formula: `\${day} = ${end}`, color: '#9aa5b1'});
    // the notebook on the days it speaks of: what happened stays next to how the groups grew
    const notes = new Map<number, string[]>();
    for (const n of meta.value.notes)
      notes.set(n.day, [...notes.get(n.day) ?? [], n.animal ? `${n.animal}: ${n.text}` : n.text]);
    for (const [day, texts] of notes) {
      const text = texts.join('; ');
      lines.addLine({title: text.length > 32 ? `${text.slice(0, 31)}…` : text, formula: `\${day} = ${day}`,
        color: '#8b5cc4'});
    }
    return t;
  });
  const growth = new BandChart(growthData, {x: 'day', lower: 'lower', upper: 'upper', line: 'median', split: 'group',
    xTitle: 'Day after randomization', yTitle: 'Median volume (mm³)', log: growthAxis.log, markers: true, fit: true,
    redrawOn: meta, dashed});
  growth.root.dataset.pmx = 'monitor-growth';
  // the reveal where it is asked for, as in Decisions and Welfare (`reveal` is defined with the log, below)
  const growthHidden = divV([infoNote('<b>The groups are hidden.</b> Reveal them to see each group\'s growth.'),
    button('Reveal the groups', () => reveal())]);
  growthHidden.style.gap = 'var(--dg-space-m)';
  // under the chart, a line per combination: where it stands against its band on the last day measured
  const bandLines = div([], 'pmx-band-lines');
  bandLines.dataset.pmx = 'monitor-independence';
  form.effect(() => {
    bandLines.replaceChildren(...[...independent.value].map(([ab, x]) => {
      const line = ui.divText(bandWords(ab, x.band));
      line.title = bandWhy(x.a, x.b);
      return line;
    }));
  });
  form.effect(() => {
    growth.root.hidden = !revealed.value;
    growthAxis.row.hidden = !revealed.value;
    growthHidden.hidden = revealed.value;
    bandLines.hidden = !revealed.value;
  });

  // Calendar: the study by its days, today among them once day 0 has a date
  const setMeta = (m: StudyMeta) => {
    const df = studyTable();
    if (df) {
      setMetaOn(df, m);
      metaVersion.value++;
    }
  };
  // a column the form maps a role to, else the study table's own of that name
  const roleCol = (df: DG.DataFrame, role: string): DG.Column | null => {
    const c = call.getParamValue(role);
    return c instanceof DG.Column ? c : df.col(role) ?? df.col('~' + role);
  };
  // each group's animals, by ID
  const animalsByGroup = (df: DG.DataFrame): Map<string, string[]> => {
    const [animal, group] = [roleCol(df, 'animal'), roleCol(df, 'group')];
    const out = new Map<string, string[]>();
    for (let i = 0; animal && group && i < df.rowCount; i++) {
      if (!group.isNone(i) && !animal.isNone(i))
        out.set(String(group.get(i)), [...out.get(String(group.get(i))) ?? [], String(animal.get(i))]);
    }
    return out;
  };
  form.effect(() => {
    void result.value;
    const df = studyTable();
    if (df && !groupsBefore.has(df))
      groupsBefore.set(df, animalsByGroup(df));
  });
  const studyPlan = (): SavedPlan | null => {
    const text = planText();
    try {
      return text ? parsePlan(text) : null;
    } catch {
      return null;
    }
  };
  const calendarDays = dayActions(() => meta.peek(), setMeta, () => {
    const df = studyTable();
    const a = df ? roleCol(df, 'animal') : null;
    return a ? [...new Set(a.toList().map(String))].sort(naturalCompare) : [];
  }, {logGiven: true, last: () => monitorEvents.last,
    name: (g) => revealed.peek() ? g : blindName(meta.peek())(g)});
  const calendarSummary = divV([], 'pmx-calendar-summary');
  Object.assign(calendarSummary.style, {gap: '4px', lineHeight: '1.45', marginBottom: 'var(--dg-space-m)'});
  const calendarBox = div([]);
  let monitorEvents: {last: number, events: [number, string][], end?: number} = {last: 0, events: []};
  form.effect(() => {
    const r = result.value;
    const m = meta.value;
    const df = studyTable();
    if (!r || !df) {
      calendarSummary.replaceChildren();
      calendarBox.replaceChildren();
      return;
    }
    const dayCol = roleCol(df, 'day');
    const volCol = roleCol(df, 'volume');
    const measured = new Set<number>();
    for (let i = 0; i < df.rowCount; i++) {
      if (dayCol && volCol && !dayCol.isNone(i) && !volCol.isNone(i))
        measured.add(dayCol.get(i));
    }
    const latest = r.status.get('day', 0) as number;
    const plan = studyPlan();
    const startWeekday = day0Weekday(m.start, plan);
    const weekdays = new Set(weekdaysOf(plan?.inputs.measurementDays as string | undefined));
    // the plan's measurement days: day 0 and its weekdays (no plan: only the days measured count)
    const planned = (d: number) => plan !== null && (d === 0 || weekdays.has((startWeekday + d) % 7));
    const today = dayToday(m.start);
    const round = (x: number | undefined) => x === undefined ? undefined : Math.round(x);
    const [end, lo, hi] = [round(cell('end_median').value), round(cell('end_lower').value),
      round(cell('end_upper').value)];
    const planEnd = plan?.predicted.studyEnd ?? null;
    // the calendar's last day: the study's own span, and today while it is no more than two weeks past it (a study
    // opened again a year on is not a year of empty weeks, each a measurement missed); the likely end's range to a
    // week past the expected end at most (halfway through it can run weeks on, each an empty week to scroll)
    const span = Math.max(latest, end === undefined ? hi ?? 0 : Math.min(hi ?? end, end + 7), planEnd ?? 0);
    const last = Math.max(span, Math.min(today ?? 0, span + 14));
    // the day the study is on: today's, or without a date the last day measured
    const now = Math.min(today ?? latest, last);
    const looks = new Set(String(df.getTag(LOOKS_TAG) ?? '').split(/[,;\s]+/).filter((x) => x).map(Number));
    const missing = Array.from({length: Math.max(0, now)}, (_, d) => d).filter((d) => planned(d) && !measured.has(d));

    const lines: HTMLElement[] = [];
    const todayMark = mark('Today', CAL.blue, true);
    if (today === null) {
      lines.push(marked(todayMark, `Day 0 has no date yet: today is taken as the last day measured, day ` +
        `${latest}. Set the date from any day's menu (right click) and it follows today.`));
    } else if (measured.size && startMismatch(m.start, latest)) {
      // measurements past today's study day: the date given to day 0 is not this table's
      lines.push(marked(todayMark, `${startMismatch(m.start, latest)} Set day 0's date from any day's menu (right ` +
        'click).'));
    } else if (today < 0)
      lines.push(marked(todayMark, `Day 0 is in ${-today} days, on ${dateText(dateOf(m.start, 0)!, true)}.`));
    else {
      lines.push(marked(todayMark, `Day ${today}, ${dateText(dateOf(m.start, today)!, true)}: ${today} days in` +
        (end === undefined ? '.' : today <= end ? `, about ${end - today} to the expected end (day ${end}` +
          (lo !== undefined && hi !== undefined ? `, likely ${Math.max(lo, today)} to ${hi}).` : ').') :
          `, past the expected end (day ${end}).`)));
    }
    if (today !== null && planned(today) && !measured.has(today))
      lines.push(marked(mark('Measure today', CAL.blue, true), 'Today is a measurement day and has no values yet.'));
    lines.push(marked(mark('Measured', CAL.slate, true), `${measured.size} days with measurements, the last ` +
      `day ${latest}.`));
    if (missing.length > 0) {
      lines.push(marked(mark('No data', CAL.red, false), `The plan measures on day${missing.length > 1 ? 's' : ''} ` +
        `${missing.join(', ')}, with no values in the table.`));
    }
    // each dosed group: the doses its schedule gives (with the study's exceptions) up to today against the study's,
    // for the animals it has and their last body weights; once doses are logged, how many dosing days so far are
    if (revealed.value) {
      const groupCol = roleCol(df, 'group');
      const animalCol = roleCol(df, 'animal');
      const weightCol = df.col(df.getTag(WEIGHT_TAG) ?? '') ?? null;
      const excluded = df.col(EXCLUDED);
      const rows: WeightRow[] = [];
      for (let i = 0; groupCol && animalCol && i < df.rowCount; i++) {
        if (groupCol.isNone(i) || animalCol.isNone(i))
          continue;
        rows.push({group: String(groupCol.get(i)), animal: String(animalCol.get(i)),
          day: dayCol && !dayCol.isNone(i) ? Number(dayCol.get(i)) : -Infinity,
          weight: weightCol && !weightCol.isNone(i) ? Number(weightCol.get(i)) : null, excluded: ticked(excluded, i)});
      }
      for (const [group, {animals, kg}] of groupWeights(rows, weightCol?.meta.units ?? null)) {
        const d = armDoses(group, m, planEnd ?? span, animals, kg, Math.max(-1, now), false);
        if (d) {
          const line = marked(mark('Dose', CAL.green, false), `${group}: ${d.text}.`);
          line.title = `${d.tip} The doses so far are the schedule's up to today, less the skipped and with the ` +
            'added; log the doses given, skip or add one from a day\'s menu in the calendar.';
          lines.push(line);
        }
      }
    }
    // the days after the expected end stay on the calendar (the end is a forecast: the study may run on to its
    // likely range), but what falls on them is said to be conditional, here, on each day and in the .ics
    const runsOn = (d: number) => end !== undefined && d > end && d >= now && !measured.has(d);
    if (end !== undefined && last > end && Math.max(now, end + 1) <= last) {
      lines.push(marked(mark('If it runs on', CAL.amber, false), `Days ${Math.max(now, end + 1)}–${last} are past ` +
        `the expected end (day ${end}): their doses and measurements happen only if the study runs on.`));
    }
    calendarSummary.replaceChildren(...lines);

    const look = (d: number, own: HTMLElement[]): DayLook => {
      const marks: HTMLElement[] = [];
      if (d === today)
        marks.push(mark('Today', CAL.blue, true));
      if (d === 0)
        marks.push(mark('Randomize', CAL.blue, false));
      if (measured.has(d))
        marks.push(mark('Measured', CAL.slate, true));
      else if (planned(d)) {
        marks.push(d < now ? mark('No data', CAL.red, false) : d === today ? mark('Measure today', CAL.blue, true) :
          mark(runsOn(d) ? 'Measure if running' : 'Measure', CAL.slate, false));
      }
      if (looks.has(d))
        marks.push(mark('Review', CAL.dark, false));
      marks.push(...own);
      if (d === end)
        marks.push(mark('Expected end', CAL.dark, true));
      if (planEnd !== null && d === planEnd && d !== end)
        marks.push(mark('Planned end', CAL.dark, false));
      const detail = runsOn(d) && d !== today ? `past the expected end (day ${end}): its doses and measurements ` +
        'only if the study runs on' : d === today ? 'today' : measured.has(d) ? 'measured' : planned(d) ? (d < now ?
          'a measurement day with no values' : 'a measurement day') : undefined;
      // the days gone by shaded, so where the study stands reads at a glance; the likely end's days in amber
      return {marks, detail, border: d === today ? `2px solid ${CAL.blue}` : undefined,
        tint: d < now ? `color-mix(in srgb, ${CAL.slate} 16%, transparent)` :
          lo !== undefined && hi !== undefined && d >= lo && d <= hi ?
            `color-mix(in srgb, ${CAL.amber} 10%, transparent)` : undefined};
    };
    // from today (or the last day measured) on: what is still to do, for the .ics export
    monitorEvents = {last, events: Array.from({length: last - Math.max(0, now) + 1}, (_, i) => Math.max(0, now) + i)
      .flatMap((d): [number, string][] => [...planned(d) && !measured.has(d) ? [[d, 'Measure tumors'] as [number,
        string]] : [], ...looks.has(d) ? [[d, 'Planned review'] as [number, string]] : [],
      ...d === end ? [[d, 'Expected end of study'] as [number, string]] : []]), end};
    calendarBox.replaceChildren(weeksGrid({meta: unnamed(m, revealed.value), startWeekday, last, day: look,
      control: String(r.status.get('control', 0)), onMenu: calendarDays.menu,
      onNote: calendarDays.addNote}));
  });

  // Decisions
  const decisionBody = div([], 'pmx-monitor-decisions');
  const logBox = div([], 'pmx-monitor-log');
  const record = (what: string) => {
    const df = studyTable();
    if (df) {
      addLog(df, what);
      blindColumn(df);
      followGroup(df);
      logVersion.value++;
    }
  };
  const reveal = () => record(REVEALED);
  form.effect(() => {
    const r = result.value;
    void logVersion.value;
    if (!r)
      return;
    if (!revealed.value) {
      decisionBody.replaceChildren(infoNote('<b>The groups are hidden.</b> The chance of success tells how the ' +
        'groups ' +
        'differ: reveal them only when someone must decide.'), button('Reveal the groups', reveal));
      return;
    }
    const d = rowsOf(r.decisions).filter((a) => a.planned);
    if (!r.status.get('at_look', 0) || d.length === 0) {
      const next = r.status.get('next_look', 0);
      decisionBody.replaceChildren(infoNote(r.status.col('next_look')?.isNone(0) || !next ?
        'No planned review ahead. The chance of success is shown only on the reviews fixed at the start.' :
        `The next planned review is on ${dayText(next)}. The chance of success is shown only on the reviews fixed at ` +
        'the start: looking between them would make it easier to stop on chance.'));
      return;
    }
    const review = r.status.get('day', 0) as number;
    const log = studyTable() ? logOf(studyTable()!) : [];
    // a stop holds from the review it was taken on, whichever review this is, until undone; a continue is this
    // review's
    const stopped = stopsOf(log);
    const continued = (arm: string) => log.some((l) => l.what === decisionEntry(arm, review, 'continue'));
    fill(decisionBody, d, (a) => a.arm, [['Group', (a) => swatch(a.arm, colorsOf(r)[a.arm])],
      // only the low end as a figure (MONITOR_DESIGN.md §8); a clearly effective arm reads "High" (chanceWords)
      ['Chance of a clear result', (a) => chanceWords(a.success, a.tgi_lower)],
      ['TGI if it continues', (a) =>
        `${Math.round(a.tgi_median)}% (${percentRange(a.tgi_lower, a.tgi_upper)})`],
      ['Suggests', (a) => badge(a.success < FUTILITY ? 'Consider stopping' : 'Continue',
        {variant: a.success < FUTILITY ? 'warning' : 'success'})],
      // once decided, the decision in place of the buttons (the log keeps it, with the review's day); a stop ends
      // the arm's doses after its review (the calendar's too), so it can be undone
      ['You decide', (a) => {
        const stop = stopped[a.arm];
        if (stop !== undefined) {
          return divH([span(`Stopped on day ${stop}, no more doses`),
            button('Undo', () => record(decisionEntry(a.arm, stop, 'stop undone')))]);
        }
        return continued(a.arm) ? span('Continued') : divH([
          button('Continue', () => record(decisionEntry(a.arm, review, 'continue'))),
          button('Stop', () => record(decisionEntry(a.arm, review, 'stop')))]);
      }]]);
  });
  form.effect(() => {
    void logVersion.value;
    void result.value;
    const df = studyTable();
    const log = df ? logOf(df) : [];
    // while the groups are hidden, an entry names them as the calendar does ("Group k": a decision names its arm);
    // the longest names first, so a combination is not taken for one of its drugs
    const m = meta.value;
    const order = [...new Set([...m.groups.map((g) => g.group), ...result.value?.days.col('group')?.toList() ?? []])];
    const named = revealed.value ? [] : [...order].sort((a, b) => b.length - a.length);
    const what = (text: string) => named.reduce((t, g) => t.split(g).join(`Group ${order.indexOf(g) + 1}`), text);
    if (log.length === 0)
      logBox.replaceChildren(span('Nothing recorded yet.'));
    else {
      fill(logBox, log, (l) => l.when + l.what, [['When', (l) => l.when.slice(0, 16).replace('T', ' ')],
        ['Who', (l) => l.who], ['What', (l) => what(l.what)]]);
    }
  });
  // Welfare: Analyze's Tolerability on the weights so far; by group, so once revealed
  const welfareFrames = computed(() => {
    const r = result.value;
    if (!r)
      return undefined;
    const control = r.status.get('control', 0) as string;
    const today = r.status.get('day', 0) as number;
    const arms = r.tolerability.rowCount === 0 ? [] :
      [...new Set(r.tolerability.col('group')!.toList() as string[])].filter((g) => g !== control);
    return {...r, analysis: DG.DataFrame.fromColumns([DG.Column.fromStrings('control', [control])]),
      comparisons: DG.DataFrame.fromColumns([DG.Column.fromStrings('arm', arms)]),
      volumes: r.days.clone(DG.BitSet.create(r.days.rowCount, (i) => r.days.get('day', i) <= today))};
  });
  const welfare = tolerabilityTab(form, welfareFrames, colorsOf);
  const welfareBody = divV(welfare.content);
  const welfareHidden = divV([infoNote('<b>The groups are hidden.</b> Body weight by group tells how the groups ' +
    'differ. The animals to look at are under Today, without their group.'), button('Reveal the groups', reveal)]);
  welfareHidden.style.gap = 'var(--dg-space-m)';
  form.effect(() => {
    const has = (result.value?.tolerability.rowCount ?? 0) > 0;
    welfareHidden.hidden = !has || revealed.value;
    welfareBody.hidden = has && !revealed.value;
  });

  // what a study row had when opened or last saved: its measurements and the tags it is saved from (domain storage:
  // closing the view with anything else asks first)
  const savedAs = new WeakMap<DG.DataFrame, string>();
  const snapshot = (df: DG.DataFrame) => [df.toCsv(), ...SAVED_TAGS.map((t) => df.getTag(t) ?? '')].join('|');
  const saveStudyNow = (): Promise<void> => {
    const df = studyTable();
    if (!df)
      return Promise.reject(new Error('No study table'));
    remember();
    // through setMeta: the header shows the ID the study has just been given
    const groups = roleCol(df, 'group');
    const names = groups ? groups.categories.filter((g) => g) : [];
    // a group renamed in the table since it was read or saved (the same animals under a new name) takes its details,
    // kept exposure and analysis's TGI with it
    const now = animalsByGroup(df);
    const renames = renamedGroups(groupsBefore.get(df) ?? now, now);
    groupsBefore.set(df, now);
    if (renames.size && df.getTag(EXPOSURE_TAG))
      df.setTag(EXPOSURE_TAG, JSON.stringify(renamedExposure(exposureOfTag(df), renames) ?? {}));
    const details = withGroups(metaOf(df) ?? emptyMeta(), names, renames);
    if (!domainStorageOn()) {
      setMeta(withId(details));
      const said = (p: string) => grok.shell.info(`Study saved: ${shownPath(p)}`);
      void Promise.resolve().then(() => saveStudy(savedStudy(df), df.name))
        .then(said, failed('The study could not be saved'));
      return Promise.resolve();
    }
    // a row: the one the table was opened from or saved to, at its version (another's save since is refused)
    return withCode(details).then(async (m) => {
      setMeta(m);
      const s = await stamped(savedStudy(df));
      const at = rowOfTable(df);
      try {
        const saved = await save('study', s, at, {fileOwner: tableFileOwner(df)});
        setRowOn(df, saved);
        head.refresh(); // its ID is now its own: no longer "Not saved yet"
        markTableFromFile(df, '');
        // (its ID was another study's: saved as a new one under a new ID)
        if (saved.meta) {
          setMeta(saved.meta);
          grok.shell.info(newRowNote('study', saved, m.id));
        }
      } catch (e) {
        // a study shared with the user to read: their version goes to a copy of their own, if they say so
        if (!at || !isReadOnly(e))
          throw e;
        return offerCopy(async () => {
          const {saved, doc} = await saveCopy('study', s, at);
          setRowOn(df, saved);
          df.name = doc.name;
          setMeta(doc.meta!);
          savedAs.set(df, snapshot(df));
          grok.shell.info(`Saved as your copy: ${esc(doc.name)}`);
        });
      }
      savedAs.set(df, snapshot(df));
      grok.shell.info('Study saved');
    }).catch((e) => {
      grok.shell.error(saveError(e, 'study'));
      throw e;
    });
  };
  const saveButton = buttonWithIcon('Save study', 'save', () => void saveStudyNow().catch(() => {}));
  if (domainStorageOn()) {
    // the table as first shown, once the view has read it (what it writes on its first reading is not a change)
    form.effect(() => {
      const df = result.value ? studyTable() : null;
      if (df && rowOfTable(df) && !savedAs.has(df))
        savedAs.set(df, snapshot(df));
    });
    queueMicrotask(() => {
      saveGate(self, 'study', saveStudyNow, () => {
        const now = studyTable();
        return !!now && savedAs.has(now) && savedAs.get(now) !== snapshot(now);
      });
    });
  }
  saveButton.title = domainStorageOn() ? 'Save the plan, the measurements and the log in PMX Studio (yours until ' +
    'you share it)' : 'Save the plan, the measurements and the log to your files';
  let toAnalyze = false;
  const analyzeThis = () => {
    const df = studyTable();
    // (a second click while Analyze is being built would open it twice)
    if (!df || toAnalyze)
      return;
    toAnalyze = true;
    // Analyze shows the groups: that is revealing them
    if (!revealed.peek())
      reveal();
    // in place of an Analyze view already open
    [...grok.shell.views].find((v) => stageOf(v) === 'analyze')?.close();
    // Analyze's Save (and Ctrl+S) saves this table as this view does, while it is this view's table
    tableSavers.set(df, () => studyTable() === df ? saveStudyNow() :
      Promise.reject(new Error('Open the study in Monitor to save it')));
    void analyzeView({study: df}).then((v) => showStage('analyze', v)).catch(failed('Analyze could not be opened'))
      .finally(() => toAnalyze = false);
  };
  const analyzeButton = buttonWithIcon('Analyze this study', 'chart-line', analyzeThis);
  analyzeButton.title = 'Analyze this table as planned; in a blinded study this reveals the groups, and the log ' +
    'says so';

  // the three limits an animal is checked against, together, as the last run read them
  const limits = div([], 'pmx-monitor-limits');
  form.effect(() => {
    if (!result.value)
      return;
    const v = (k: string) => call.getParamValue(k);
    limits.replaceChildren(noteLine(`Endpoint ${v('endpoint')} mm³ (median control) · humane limit ` +
      `${v('humane')} mm³ · weight loss ${v('weightLimit') ?? 20}%`));
  });
  // the study's data in Datagrok's own table view; here, one line on it and the way there
  const dataLine = studyDataLine(() => studyInputsOf(call));
  self.own(dataLine.dispose);
  // after a run that failed too: the way to the table is what a study with no forecast yet needs
  form.effect(() => {
    void result.value;
    void error.value;
    dataLine.refresh();
  });
  // a study with animals and volumes but no groups yet has no forecast, so no panels: forming the groups is here
  const firstGroups = divH([buttonWithIcon('Propose groups', 'random', proposeAsked),
    span('Day-0 volumes in the table and no groups yet? Random groups alike in their tumors, from the plan.')]);
  Object.assign(firstGroups.style, {gap: 'var(--dg-space-m)', alignItems: 'center'});
  form.effect(() => firstGroups.hidden = !!result.value || !error.value || !studyTable()?.getTag(PLAN_TAG));
  const panels = panelTabs([
    {id: 'today', title: 'Needs attention', hint: 'Who needs attention today; a person decides what to do',
      content: [limits, entered, attention, noteBox]},
    {id: 'growth', title: 'Growth', hint: 'Each group so far and what is expected up to the end',
      content: [growthHidden, growthAxis.row, sized(growth, 340), bandLines]},
    {id: 'calendar', title: 'Study calendar', hint: 'The study by its days: today, the days measured and to measure, ' +
      'doses, reviews and the expected end. Double-click a day to note what happened, right-click it for more',
    content: [calendarSummary, calendarBox]},
    {id: 'sheet', title: 'Measurement sheet', hint: 'The range the next value of each animal should fall in: check ' +
      'each value as you measure it', content: [sheet]},
    {id: 'decisions', title: 'Decisions', hint: 'Whether each treated group is still worth continuing, and what was ' +
      'decided', content: [decisionBody, block('Log', 'Who did what, and when', logBox)]},
    {id: 'welfare', title: 'Welfare', hint: 'Whether the treatments are being tolerated, from the body weights so far',
      content: [welfareHidden, welfareBody]},
    {id: 'enrollment', title: 'Enrollment', hint: 'Random groups alike in their starting tumors',
      content: [divH([proposeButton, sheetButton]), allocationTable]},
  ], MONITOR_TABS, 'pmx-monitor-panels', 'monitor');
  exports.push(['Layout…', () => editLayout(panels.root), 'View']);
  // the study's colours: a colour changed in a chart is the study's, in every stage open on it
  const link: StudyLink = {meta: () => metaOf(studyTable()), setMeta, colors: () => {
    const r = result.peek();
    return r ? colorsOf(r) : {};
  }, table: studyTable, plan: () => planKeyOf(studyTable()),
  control: () => result.peek()?.status.get('control', 0) ?? null};
  self.own(linkStudy(link));
  exports.push(['Reset colors', () => resetColors(link), 'View',
    () => hasOwnColors(link) ? null : 'The colors are the automatic ones']);

  // a study whose plan was made with another atlas says so above its numbers (the plan the form gives, else the
  // table's); the latest answer stands
  const atlasLine = div([], 'ui-label-details');
  atlasLine.hidden = true;
  let atlasAsked = 0;
  form.effect(() => {
    void result.value;
    const asked = ++atlasAsked;
    void atlasNote(studyPlan()?.provenance, FROZEN_ATLAS_NOTE).then((t) => {
      if (asked === atlasAsked)
        atlasLine.hidden = !(atlasLine.textContent = t ?? '');
    });
  });
  // a model planned from its history: its twin is not validated for it (D694)
  const twinLine = markLine(TWIN_NOTE);
  form.effect(() => {
    void result.value;
    const df = studyTable();
    twinLine.hidden = !savedHistory(String(call.getParamValue('line') || (df ? lineOf(df) : '') || ''));
  });
  const content = divV([atlasLine, twinLine, cards, panels.root], 'pmx-monitor-content');
  const needStudy = () => studyTable() ?? (() => {
    throw new Error('Open a study first');
  })();
  exports.push(
    ['Report (PDF)…', async () => {
      const df = needStudy();
      // the report is what is on screen: not while that is behind the table or the form
      if (busy.peek() || stale.peek() || !result.peek())
        throw new Error('The forecasts are not up to date yet: export the report once they are');
      printSheet(await panels.all(() => reportHtml({title: df.name, stage: 'Monitor', meta: metaOf(df), content,
        groups: roleCol(df, 'group')?.categories.filter((g) => g),
        planAtlas: studyPlan() ? studyPlan()!.provenance?.atlas ?? null : undefined, log: logOf(df),
        lastDay: result.peek()!.status.get('day', 0) as number})));
    }],
    ['Study file (.pmxstudy)', () => {
      const df = needStudy();
      remember();
      void stamped(savedStudy(df))
        .then((s) => download(`${exportName(df)}.pmxstudy`, JSON.stringify(s), 'application/json'),
          failed('The study file could not be made'));
    }],
    ['Measurements (CSV)', () => download(`${needStudy().name}.csv`, csvOf(needStudy()), 'text/csv')],
    ['Measurements', () => openStudyTable(needStudy()), 'Data > Open as table'],
    ['Calendar (.ics)', () => {
      const df = needStudy();
      download(`${exportName(df)}.ics`, icsOf(exportName(df), unnamed(meta.peek(), revealed.peek()), monitorEvents.last,
        monitorEvents.events, monitorEvents.end),
      'text/calendar');
    }]);
  Object.assign(content.style, {gap: 'var(--dg-space-l)'});
  form.effect(() => {
    content.hidden = !result.value;
    // a study open with nothing measured yet: its empty state says what to do, not this too
    start.hidden = !!result.value || empty.value;
  });
  const panes = stagePanes(form, content, {busy, result, error, stale, empty,
    before: [start, firstGroups, dataLine.root],
    busyText: 'Forecasting the study…', staleText: 'the forecasts shown are for the last complete inputs'});
  // a session's values pasted or typed into the table: the forecasts follow, after a pause (one value, then the next)
  self.own(onTableEdits(form, computed(() => [result.value, error.value]),
    studyTable, () => auto.changed(TABLE_EDITS)));
  // the table's own subscriptions are not built under the host: dropped with it
  self.own(stopWatching);
  const head = studyHeader(form, 'monitor', computed(() => {
    const r = result.value;
    // the details changed (an ID given on saving, a note, day 0's date): the header shows them as they are
    void metaVersion.value;
    if (!r)
      return undefined;
    const colors = colorsOf(r);
    const day = cell('day').value;
    const df = studyTable();
    return {name: df ? studyNameOf(df) : '', groups: revealed.value ?
      Object.entries(colors).map(([name, color]) => ({name, color})) : [],
    status: [modelName(r.status.get('line', 0)), day === undefined ? '' : dayText(day),
      revealed.value ? '' : 'groups hidden'].filter((x) => x).join(' · '),
    groupsHidden: !revealed.value, meta: () => metaOf(df),
    // with the store, an ID is a study's once saved (a plan's or a file's may be taken: the store gives another)
    unsaved: () => domainStorageOn() && !!df && !rowOfTable(df), setMeta: df ? (m) => {
      const was = !!metaOf(df)?.blinded;
      setMeta(m);
      if (m.blinded !== was)
        record(m.blinded ? 'Blinded entry on' : 'Blinded entry off: the groups are shown');
    } : undefined, groupMenu: (g, e) => studyMenu(studyNamesOf(df), 'group', g, {causedBy: e})};
    // the header's Analyze step is this study's analysis (unless that would reveal a blinded study's groups: the
    // button does that, and says so); Plan opens with the study's plan, but not on a blinded study's hidden groups
  }), (s) => s === 'analyze' && studyTable() && revealed.peek() ? analyzeThis() : void openStage(s,
    revealed.peek() ? planText() : null));
  // saving is always in reach, as in Plan: a study's measurements, log and notes live in the table until saved
  form.effect(() => [saveButton, analyzeButton].forEach((b) => b.hidden = !result.value));
  return {
    root: takesFiles(stagePage(head.row, panes), {stage: 'monitor'}),
    ribbon: [[head.steps], [openButton({stage: 'monitor'}), saveButton, analyzeButton]],
    path: '/monitor',
    status: computed(() => busy.value ? 'Forecasting the study…' :
      result.value ? modelName(String(status.value?.get('line', 0))) : ''),
  };
}

/** Plan's "Follow this study": a new measurements table carrying the plan, opened in Monitor. */
export async function followPlan(plan: SavedPlan): Promise<void> {
  await showStage('monitor', await monitorView(newStudyTable(plan)));
}

/** A .pmxstudy file opened from Files or a drop. */
export async function openStudyFile(content: string, row: Saved | null = null): Promise<DG.DataFrame[]> {
  const s = parseStudy(content);
  const df = tableOf(s);
  if (row) // (a study row's document: its next save updates that row)
    setRowOn(df, row);
  else // (a file's: its first save is the user's copy when its owner is someone else)
    markTableFromFile(df, ownerOf('study', s));
  // a finished study (saved from Analyze, the Overview's stage 3) opens in Analyze: in Monitor its last animals would
  // read as needing attention and its next measurement day as due
  if (s.finished === true)
    await showStage('analyze', await analyzeView({study: df}));
  else
    await showStage('monitor', await monitorView(df));
  return [];
}

/** The measurements as CSV, without the group column while the groups are hidden. */
const csvOf = (df: DG.DataFrame): string => df.toCsv(groupsShown(df) ? undefined :
  {columns: df.columns.names().filter((c) => !c.startsWith('~'))});
