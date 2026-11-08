/* The Monitor view, built as Analyze is: the form generated from pmxMonitorStudy's signature, the results from its
   outputs. The groups are shown, unless the study is blinded (its details): then they stay hidden until someone
   reveals them, and that goes in the study's log. A study's measurements
   are an ordinary table (paste into its grid): Monitor only reads it, and saving the study is a gesture of its own.
   Tabs by question: who needs attention today, how the groups will be formed, how the study is growing, and whether to
   go on, only on the reviews fixed at the start (MONITOR_DESIGN.md). */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {badge, BasicTable, button, buttonWithIcon, computed, div, divH, divV, signal, span,
  Splitter} from '@datagrok-libraries/u2';
import type {ReadonlySignal, StatValue} from '@datagrok-libraries/u2';
import {funcForm, ViewHost} from '@datagrok-libraries/u2/dg';
import type {ViewContent} from '@datagrok-libraries/u2/dg';

import {exampleStudy} from './analyze';
import {analyzeView} from './analyze-view';
import {BandChart} from './band-chart';
import {WEEKDAY_NAMES, groupColors, weekdaysOf} from './frames';
import {allocateGroups, MonitorFrames, monitorFrames, MonitorInputs, outsideBand} from './monitor';
import {PLAN_TAG, SavedPlan, parsePlan} from './plan-file';
import {planInput} from './plan-vs-observed';
import {openStage, showStage} from './stages';
import {LIMITS_TAG, LOOKS_TAG, REVEALED, WEIGHT_TAG, addLog, groupsShown, logOf, newStudyTable, parseStudy,
  restoreFields, saveStudy, studyOf, tableOf} from './study-file';
import {editLayout, panelTabs, Tab} from './panels';
import {stagePage, studyHeader} from './study-header';
import {emptyStudy, studyDataLine, studyMenu, studyOf as studyNamesOf, tagStudy, withoutExcluded} from './study-table';
import {StudyMeta, dateOf, dateText, emptyMeta, metaOf, setMetaOn, weekdayOf, withId} from './study-meta';
import {CAL, DayLook, dayActions, dayToday, icsOf, mark, marked, weeksGrid} from './calendar';
import {LineInput} from './line-input';
import {printSheet, sheetHtml} from './sheet';
import {Exports, download, exportMenu, reportHtml} from './report';
import {studyInputsOf} from './explore';
import {tolerabilityTab} from './tolerability';
import {askedOf, autoRun, block, busyBar, busyOverlay, callInputs, cardGrid, dayText, declarationPane, DECLARATION_MIN,
  infoNote,
  onTableEdits, TABLE_EDITS,
  isBodyWeight, keepNamedColumns, noteLine, notes,
  notesBox, num, pct, resultsFrame, settleForm, sized, staleNote,
  markStage, signed, stageOf, statCard, swatch, volumeAxis, wrapCells} from './ui';

/** Monitor's tabs and their panels, in order and width (twelfths): today's work and how the groups grow; the
 *  measurements and the groups formed; how the animals bear it; whether each group is worth continuing. */
const MONITOR_TABS: Tab[] = [
  {name: 'Today', layout: [{id: 'growth', span: 12}, {id: 'today', span: 12}, {id: 'calendar', span: 12}]},
  {name: 'Measurements', layout: [{id: 'sheet', span: 12}, {id: 'enrollment', span: 12}]},
  {name: 'Welfare', layout: [{id: 'welfare', span: 12}]},
  {name: 'Decisions', layout: [{id: 'decisions', span: 12}]}];

/** Below this chance of success, on a planned review, Monitor suggests stopping the arm. A person decides. */
export const FUTILITY = 0.1;

interface AnimalRow {[c: string]: any}

const rowsOf = (df: DG.DataFrame): AnimalRow[] => Array.from({length: df.rowCount}, (_, i) =>
  Object.fromEntries(df.columns.names().map((c) => {
    const col = df.col(c)!;
    return [c, col.isNone(i) ? null : col.get(i)];
  })));

const mm3 = (x: unknown) => typeof x === 'number' && Number.isFinite(x) ? `${Math.round(x)}` : '—';
function fill<R>(host: HTMLElement, items: R[], key: (r: R) => string,
  columns: [string, (r: R) => string | HTMLElement][]): void {
  const t = new BasicTable<R>({items, keyOf: key, columns: columns.map(([header, render], j) =>
    ({header, render, align: j === 0 ? 'left' : 'center'}))});
  t.root.style.width = '100%';
  wrapCells(t.root);
  host.replaceChildren(t.root);
}

/** The Monitor view; `initial`: the table of a study to follow. */
export async function monitorView(initial?: DG.DataFrame): Promise<DG.ViewBase> {
  const exports: Exports = [];
  const host = ViewHost.create({name: 'PMX Studio', content: (self) => monitorContent(self, initial,
    exports)});
  exportMenu(host.view, exports);
  return markStage(host.view, 'monitor');
}

/** What the view shows, built under its host, which owns it all (the form, the charts, the panes). */
function monitorContent(self: ViewHost, initial: DG.DataFrame | undefined, exports: Exports): ViewContent {
  const call = DG.Func.byName('PmxStudio:pmxMonitorStudy').prepare();
  const result = signal<MonitorFrames | undefined>(undefined);
  const error = signal<string | null>(null);
  const stale = signal(false);
  const busy = signal(false);
  // the log's version: what it records (a reveal, the blinding switched) changes what is shown
  const logVersion = signal(0);
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
  // the groups, unless the study is blinded and no one has revealed them
  const revealed = computed(() => {
    void result.value;
    void logVersion.value;
    const df = studyTable();
    return !df || groupsShown(df);
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
  };
  const planState = {stale, error, result: computed(() => result.value && {analysis: result.value.status}),
    reading: Promise.resolve() as Promise<unknown>};
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
      // the reviews fixed at the start hold, whatever the field says
      const fixed = studyTable()?.getTag(LOOKS_TAG);
      if (fixed)
        call.setParamValue('looks', fixed);
      const sent = askedOf(call);
      tagStudy(studyInputsOf(call));
      const frames = await monitorFrames(withoutExcluded(callInputs<MonitorInputs>(call)));
      if (mine !== latest)
        return;
      error.value = null;
      stale.value = askedOf(call) !== sent;
      // the reviews stay with the study from its first run
      const df = studyTable();
      if (df && !df.getTag(LOOKS_TAG))
        df.setTag(LOOKS_TAG, String(call.getParamValue('looks') ?? ''));
      lockLooks();
      remember();
      result.value = frames;
    } catch (e: any) {
      error.value = String(e?.message ?? e);
    } finally {
      preparing = false;
      if (mine === latest)
        busy.value = false;
    }
  };

  const setStudy = async (df: DG.DataFrame) => {
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
    await form.settled;
    lockLooks();
  };
  // the example: Plan's example study halfway, today a planned review
  const openExample = () => run(async () => {
    const whole = await exampleStudy();
    const plan = parsePlan(whole.getTag(PLAN_TAG)!);
    const days = [...new Set(whole.col('day')!.toList() as number[])].sort((a, b) => a - b);
    const today = days.filter((d) => d <= plan.predicted.studyEnd / 2).pop()!;
    const df = whole.clone(DG.BitSet.create(whole.rowCount, (i) => whole.get('day', i) <= today));
    df.name = `Example ${plan.predicted.line} study, day ${today}`;
    df.setTag(PLAN_TAG, whole.getTag(PLAN_TAG)!);
    df.setTag(LOOKS_TAG, String(today));
    df.setTag(WEIGHT_TAG, 'body_weight');
    await setStudy(df);
  });
  const newStudy = () => run(async () => {
    const df = emptyStudy();
    await setStudy(df);
    tagStudy({study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: 'body_weight'});
    grok.shell.v = [...grok.shell.tableViews].find((v) => v.dataFrame.id === df.id) ?? grok.shell.v;
  });
  void form.settled.then(() => {
    const study = form.getInput('study') as unknown as {actions?: {id: string, icon: string, tooltip: string,
      run: () => void}[]} | null;
    if (study?.actions) {
      study.actions = [...study.actions, {id: 'example', icon: 'flask', tooltip: 'Follow an example study: MC38 as ' +
        'Plan predicts it, halfway through', run: () => void openExample()}, {id: 'new', icon: 'plus',
        tooltip: 'Start an empty study table, to type or paste the measurements into', run: () => void newStudy()}];
    }
    if (initial)
      void setStudy(initial);
  });
  planInput(form, call, planState);

  const start = divV([
    span('Open the study\'s measurements (one row per animal and day, day 0 the randomization), check the columns ' +
      'and the forecasts follow. Paste or type each session\'s values into the table: each is checked as you enter ' +
      'it, and the forecasts update.'),
    divH([button('Follow an example study', () => void openExample()),
      span('An MC38 study as Plan predicts it, halfway through.')]),
    infoNote('<b>Blinded entry.</b> When whoever enters the measurements must not know the groups, turn it on in ' +
      'the study\'s details (click its name): Monitor then shows each animal without its group until someone ' +
      'reveals them, and the log says so.'),
  ], 'pmx-monitor-start');
  Object.assign(start.style, {gap: 'var(--dg-space-m)', maxWidth: '760px'});
  (start.children[1] as HTMLElement).style.gap = 'var(--dg-space-m)';
  (start.children[1] as HTMLElement).style.alignItems = 'center';

  const status = computed(() => result.value?.status);
  const cell = <T = number>(name: string) => computed(() => {
    const s = status.value;
    return (s && !s.col(name)?.isNone(0) ? s.get(name, 0) ?? undefined : undefined) as T | undefined;
  });
  const card = (label: string | ReadonlySignal<string>, value: ReadonlySignal<StatValue | undefined>, tip: () => string,
    format?: (x: StatValue) => string) => statCard(form, label, value, tip, format);
  const metaVersion = signal(0);
  const meta = computed(() => {
    void result.value;
    void metaVersion.value;
    return metaOf(studyTable()) ?? emptyMeta();
  });
  // a day's card says its date once day 0 has one; the expected end, how it stands against the plan's
  const dated = (label: string, day: ReadonlySignal<number | undefined>, also?: () => string) => computed(() => {
    const d = day.value;
    const date = d === undefined ? null : dateOf(meta.value.start, d);
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
      const name = (e.target as HTMLElement).closest('.u2-table-row')?.querySelector('td')?.textContent?.trim();
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
    if (!r || !df || !['animal', 'day', 'volume'].every((c) => df.col(c))) {
      entered.replaceChildren();
      return;
    }
    const nextDay = r.status.get('next_day', 0);
    const off = outsideBand(Array.from({length: df.rowCount}, (_, i) => ({animal: String(df.get('animal', i)),
      day: df.get('day', i), volume: df.col('volume')!.isNone(i) ? NaN : df.get('volume', i)})), nextDay,
    new Map(rowsOf(r.animals).map((a) => [a.animal, a as any])));
    entered.replaceChildren(...off.length === 0 ? [] : [noteLine(`Measure again (${dayText(nextDay)}): ` +
      off.map((o) => `${o.animal} at ${mm3(o.volume)} mm³, expected ${mm3(o.lower)}–${mm3(o.upper)}`).join('; '),
    'warning')]);
  };
  let watching: {unsubscribe(): void}[] = [];
  const stopWatching = () => {
    watching.forEach((w) => w.unsubscribe());
    watching = [];
  };
  form.effect(() => {
    void result.value;
    stopWatching();
    const df = studyTable();
    if (df)
      watching = [df.onDataChanged.subscribe(checkEntered), df.onRowsAdded.subscribe(checkEntered)];
    checkEntered();
  });
  const groupCols = (): [string, (r: AnimalRow) => string][] => revealed.value ? [['Group', (a) => a.group]] : [];
  // only for a study with body weights
  const weightCols = (r: MonitorFrames): [string, (a: AnimalRow) => string][] => r.animals.col('weight_change') ?
    [['Body weight', (a) => typeof a.weight_change === 'number' ? signed(a.weight_change) : '—']] : [];
  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    // by animal: the engine's order, group by group, would show the groups
    const all = rowsOf(r.animals).sort((a, b) => String(a.animal).localeCompare(String(b.animal), undefined,
      {numeric: true}));
    const listed = all.filter((a) => a.attention);
    if (listed.length === 0)
      attention.replaceChildren(noteLine('No animal needs attention today.', 'success'));
    else {
      fill(attention, listed, (a) => a.animal, [['Animal', (a) => a.animal], ...groupCols(),
        ['Why', (a) => a.attention], ['Last volume', (a) => `${mm3(a.last_volume)} mm³, day ${a.last_day}`],
        ...weightCols(r),
        ['Next expected', (a) => a.next_median === null ? '—' : `${mm3(a.next_lower)}–${mm3(a.next_upper)} mm³`]]);
    }
    fill(sheet, all, (a) => a.animal, [['Animal', (a) => a.animal], ...groupCols(),
      ['Last volume', (a) => mm3(a.last_volume)], ...weightCols(r),
      ['Next expected', (a) => a.next_median === null ? '—' : `${mm3(a.next_lower)}–${mm3(a.next_upper)}`],
      ['Measure again outside', (a) => a.next_median === null ? '—' : `${mm3(a.check_lower)}–${mm3(a.check_upper)}`]]);
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
      throw new Error('Propose groups needs the study\'s plan, for its groups and animals');
    const plan: SavedPlan = parsePlan(text);
    // each animal's latest volume
    const last = new Map<string, {day: number, v: number}>();
    for (let i = 0; i < df.rowCount; i++) {
      const [a, d, v] = [String(df.get('animal', i)), df.get('day', i), df.get('volume', i)];
      if (!(last.get(a)?.day! >= d) && Number.isFinite(v))
        last.set(a, {day: d, v});
    }
    const groups = [{name: 'Vehicle', n: plan.predicted.nControl ?? 0},
      ...plan.predicted.arms.map((a) => ({name: a.name, n: a.n ?? 0}))];
    // a new seed per study: it decides who gets what
    const seed = crypto.getRandomValues(new Uint32Array(1))[0];
    const {allocation, balance} = await allocateGroups([...last.keys()], [...last.values()].map((x) => x.v), groups,
      seed);
    const group = new Map(rowsOf(allocation).map((r) => [r.animal, r.group]));
    const col = df.col('group') ?? df.columns.addNewString('group');
    for (let i = 0; i < df.rowCount; i++)
      col.set(i, group.get(String(df.get('animal', i))) ?? '');
    addLog(df, `Groups allocated (seed ${seed})`);
    fill(allocationTable, rowsOf(balance), (b) => b.group, [['Group', (b) => b.group], ['Animals', (b) => String(b.n)],
      ['Geometric mean volume', (b) => `${mm3(b.geo_mean)} mm³`],
      ['Spread (sd of log volume)', (b) => String(num(b.sd_log, 2))]]);
  };
  const proposeButton = buttonWithIcon('Propose groups', 'random', () => void propose().catch((e) =>
    error.value = String(e?.message ?? e)));
  proposeButton.title = 'Fill the table\'s group column with random groups alike in their tumors';
  // the animal room's sheet, from the groups in the table and the seed the log kept
  const sheetButton = buttonWithIcon('Print randomization sheet', 'print', () => {
    const df = studyTable();
    if (!df?.col('group')) {
      error.value = 'The sheet needs the study\'s groups: propose them, or open a table with a group column';
      return;
    }
    const seed = [...logOf(df)].reverse().map((e) => /seed (\d+)/.exec(e.what)?.[1]).find((s) => s) ?? null;
    try {
      printSheet(sheetHtml(df.name, metaOf(df), df, seed));
      addLog(df, 'Randomization sheet printed (it shows the groups)');
    } catch (e: any) {
      error.value = String(e?.message ?? e);
    }
  });
  sheetButton.title = 'A page for the animal room: each animal, its group and its volume at randomization, ' +
    'what each group receives, and lines to sign. It shows the groups.';

  // Growth: by group, once revealed
  const growthAxis = volumeAxis();
  const colorsOf = (r: MonitorFrames) => {
    const names = [...new Set(r.days.col('group')!.toList() as string[])];
    const colors = groupColors(names, r.status.get('control', 0));
    return Object.fromEntries(names.map((n) => [n, colors[n]]));
  };
  const growthData = computed(() => {
    const r = result.value;
    if (!r || !revealed.value)
      return undefined;
    r.days.col('group')!.meta.colors.setCategorical(colorsOf(r));
    // what the study is read against: the endpoint, today and the expected end (the chart's scale follows them)
    const lines = r.days.meta.formulaLines;
    lines.clear();
    const endpoint = Number(call.getParamValue('endpoint'));
    if (endpoint > 0)
      lines.addLine({title: `Endpoint ${endpoint} mm³`, formula: `\${median} = ${endpoint}`});
    lines.addLine({title: `Today, day ${r.status.get('day', 0)}`, formula: `\${day} = ${r.status.get('day', 0)}`,
      color: '#2083d5'});
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
    return r.days;
  });
  const growth = new BandChart(growthData, {x: 'day', lower: 'lower', upper: 'upper', line: 'median', split: 'group',
    xTitle: 'Study day', yTitle: 'Median volume (mm³)', log: growthAxis.log, markers: true, fit: true,
    redrawOn: meta});
  const growthHidden = infoNote('<b>The groups are hidden.</b> Reveal them to see each group\'s growth.');
  form.effect(() => {
    growth.root.hidden = !revealed.value;
    growthAxis.row.hidden = !revealed.value;
    growthHidden.hidden = revealed.value;
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
    return c instanceof DG.Column ? c : df.col(role);
  };
  const studyPlan = (): SavedPlan | null => {
    const text = call.getParamValue('plan') || studyTable()?.getTag(PLAN_TAG);
    try {
      return text ? parsePlan(text) : null;
    } catch {
      return null;
    }
  };
  const calendarDays = dayActions(() => meta.peek(), setMeta, () => {
    const df = studyTable();
    const a = df ? roleCol(df, 'animal') : null;
    return a ? [...new Set(a.toList().map(String))].sort((x, y) => x.localeCompare(y, undefined, {numeric: true})) : [];
  });
  const calendarSummary = divV([], 'pmx-calendar-summary');
  Object.assign(calendarSummary.style, {gap: '4px', lineHeight: '1.45', marginBottom: 'var(--dg-space-m)'});
  const calendarBox = div([]);
  let monitorEvents: {last: number, events: [number, string][]} = {last: 0, events: []};
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
    const startWeekday = m.start ? weekdayOf(m.start) :
      Math.max(0, WEEKDAY_NAMES.indexOf(String(plan?.inputs.startWeekday ?? 'Monday')));
    const weekdays = new Set(weekdaysOf(plan?.inputs.measurementDays as string | undefined));
    // the plan's measurement days: day 0 and its weekdays (no plan: only the days measured count)
    const planned = (d: number) => plan !== null && (d === 0 || weekdays.has((startWeekday + d) % 7));
    const today = dayToday(m.start);
    // the day the study is on: today's, or without a date the last day measured
    const now = today ?? latest;
    const round = (x: number | undefined) => x === undefined ? undefined : Math.round(x);
    const [end, lo, hi] = [round(cell('end_median').value), round(cell('end_lower').value),
      round(cell('end_upper').value)];
    const planEnd = plan?.predicted.studyEnd ?? null;
    const looks = new Set(String(df.getTag(LOOKS_TAG) ?? '').split(/[,;\s]+/).filter((x) => x).map(Number));
    const missing = Array.from({length: Math.max(0, now)}, (_, d) => d).filter((d) => planned(d) && !measured.has(d));
    const last = Math.max(latest, today ?? 0, hi ?? end ?? 0, planEnd ?? 0);

    const lines: HTMLElement[] = [];
    const todayMark = mark('Today', CAL.blue, true);
    if (today === null) {
      lines.push(marked(todayMark, `Day 0 has no date yet: today is taken as the last day measured, day ` +
        `${latest}. Set the date from any day's menu (right click) and it follows today.`));
    } else if (today < 0)
      lines.push(marked(todayMark, `Day 0 is in ${-today} days, on ${dateText(dateOf(m.start, 0)!, true)}.`));
    else {
      lines.push(marked(todayMark, `Day ${today}, ${dateText(dateOf(m.start, today)!, true)}: ${today} days in` +
        (end === undefined ? '.' : today <= end ? `, about ${end - today} to the expected end (day ${end}` +
          (lo !== undefined && hi !== undefined ? `, likely ${lo} to ${hi}).` : ').') :
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
          mark('Measure', CAL.slate, false));
      }
      if (looks.has(d))
        marks.push(mark('Review', CAL.dark, false));
      marks.push(...own);
      if (d === end)
        marks.push(mark('Expected end', CAL.dark, true));
      if (planEnd !== null && d === planEnd && d !== end)
        marks.push(mark('Planned end', CAL.dark, false));
      const detail = d === today ? 'today' : measured.has(d) ? 'measured' : planned(d) ? (d < now ?
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
      ...d === end ? [[d, 'Expected end of study'] as [number, string]] : []])};
    calendarBox.replaceChildren(weeksGrid({meta: m, startWeekday, last, day: look, onMenu: calendarDays.menu,
      onNote: calendarDays.addNote}));
  });

  // Decisions
  const decisionBody = div([], 'pmx-monitor-decisions');
  const logBox = div([], 'pmx-monitor-log');
  const record = (what: string) => {
    const df = studyTable();
    if (df) {
      addLog(df, what);
      logVersion.value++;
    }
  };
  const reveal = () => record(REVEALED);
  form.effect(() => {
    const r = result.value;
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
    fill(decisionBody, d, (a) => a.arm, [['Group', (a) => swatch(a.arm, colorsOf(r)[a.arm])],
      // only the low end in words: on the atlas's real arms a low chance held, a high one was too hopeful
      // (MONITOR_DESIGN.md §8), so a high figure is not shown
      ['Chance of a clear result', (a) => a.success < FUTILITY ? `Low (under ${pct(FUTILITY)})` : 'Still open'],
      ['TGI if it continues', (a) =>
        `${Math.round(a.tgi_median)}% (${Math.round(a.tgi_lower)}–${Math.round(a.tgi_upper)}%)`],
      ['Suggests', (a) => badge(a.success < FUTILITY ? 'Consider stopping' : 'Continue',
        {variant: a.success < FUTILITY ? 'warning' : 'success'})],
      ['You decide', (a) => divH([button('Continue', () => record(`${a.arm}: continue`)),
        button('Stop', () => record(`${a.arm}: stop`))])]]);
  });
  form.effect(() => {
    void logVersion.value;
    void result.value;
    const df = studyTable();
    const log = df ? logOf(df) : [];
    if (log.length === 0)
      logBox.replaceChildren(span('Nothing recorded yet.'));
    else {
      fill(logBox, log, (l) => l.when + l.what, [['When', (l) => l.when.slice(0, 16).replace('T', ' ')],
        ['Who', (l) => l.who], ['What', (l) => l.what]]);
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

  const saveButton = buttonWithIcon('Save study', 'save', () => {
    const df = studyTable();
    if (!df)
      return;
    remember();
    setMetaOn(df, withId(metaOf(df) ?? emptyMeta()));
    void saveStudy(studyOf(df, call.getParamValue('humane') ?? null, call.getParamValue('weightLimit') ?? null))
      .then((p) => grok.shell.info(`Study saved: ${p}`), (e) => error.value = String(e?.message ?? e));
  });
  saveButton.title = 'Save the plan, the measurements and the log to your files';
  const analyzeButton = buttonWithIcon('Analyze this study', 'chart-line', () => {
    const df = studyTable();
    if (!df)
      return;
    // Analyze shows the groups: that is revealing them
    if (!revealed.peek())
      reveal();
    // in place of an Analyze view already open
    [...grok.shell.views].find((v) => stageOf(v) === 'analyze')?.close();
    void analyzeView({study: df}).then((v) => showStage('analyze', v));
  });
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
  form.effect(() => {
    if (result.value)
      dataLine.refresh();
  });
  const panels = panelTabs([
    {id: 'today', title: 'Needs attention', hint: 'Who needs attention today; a person decides what to do',
      content: [limits, entered, attention, noteBox]},
    {id: 'growth', title: 'Growth', hint: 'Each group so far and what is expected up to the end',
      content: [growthHidden, growthAxis.row, sized(growth, 340)]},
    {id: 'calendar', title: 'Study calendar', hint: 'The study by its days: today, the days measured and to measure, ' +
      'doses, reviews and the expected end. Double-click a day to note what happened, right-click it for more',
    content: [calendarSummary, calendarBox]},
    {id: 'sheet', title: 'Measurement sheet', hint: 'The range the next value of each animal should fall in: check ' +
      'each value as you measure it', content: [sheet]},
    {id: 'decisions', title: 'Decisions', hint: 'Whether each treated group is still worth continuing, and what was ' +
      'decided', content: [decisionBody, block('Log', 'Who did what, and when', logBox), divH([saveButton,
      analyzeButton])]},
    {id: 'welfare', title: 'Welfare', hint: 'Whether the treatments are being tolerated, from the body weights so far',
      content: [welfareHidden, welfareBody]},
    {id: 'enrollment', title: 'Enrollment', hint: 'Random groups alike in their starting tumors',
      content: [divH([proposeButton, sheetButton]), allocationTable]},
  ], MONITOR_TABS, 'pmx-monitor-panels', 'monitor');
  exports.push(['Layout…', () => editLayout(panels.root), 'View']);

  const content = divV([dataLine.root, cards, panels.root], 'pmx-monitor-content');
  const needStudy = () => studyTable() ?? (() => {
    throw new Error('Open a study first');
  })();
  exports.push(
    ['Report (PDF)…', async () => {
      const df = needStudy();
      printSheet(await panels.all(() => reportHtml({title: df.name, stage: 'Monitor', meta: metaOf(df), content,
        log: logOf(df)})));
    }],
    ['Study file (.pmxstudy)', () => {
      const df = needStudy();
      remember();
      download(`${df.name}.pmxstudy`, JSON.stringify(studyOf(df, call.getParamValue('humane') ?? null,
        call.getParamValue('weightLimit') ?? null)), 'application/json');
    }],
    ['Measurements (CSV)', () => download(`${needStudy().name}.csv`, needStudy().toCsv(), 'text/csv')],
    ['Calendar (.ics)', () => {
      const df = needStudy();
      download(`${df.name}.ics`, icsOf(df.name, meta.peek(), monitorEvents.last, monitorEvents.events),
        'text/calendar');
    }]);
  Object.assign(content.style, {gap: 'var(--dg-space-l)'});
  form.effect(() => {
    content.hidden = !result.value;
    start.hidden = !!result.value;
  });
  // the first run under u2's loader; every later one under the thin bar, what was shown left in view
  const overlay = busyOverlay('Forecasting the study…');
  const bar = busyBar();
  const failure = div([], 'pmx-monitor-error');
  Object.assign(failure.style, {color: 'var(--dg-error-color, #c62828)'});
  form.effect(() => {
    overlay.hidden = !busy.value || result.value !== undefined;
    bar.hidden = !busy.value || result.value === undefined;
  });
  form.effect(() => {
    failure.hidden = error.value === null;
    failure.replaceChildren(...(error.value === null ? [] : [noteLine(error.value, 'error')]));
  });
  // the view follows the form by itself: what it shows is behind only while a field is missing or invalid
  const outdated = staleNote('A field is missing or invalid: the forecasts shown are for the last complete inputs');
  form.effect(() => {
    outdated.hidden = !stale.value || form.missingRequired.value.length + form.invalidFields.value.length === 0;
    content.style.opacity = stale.value ? '0.5' : '';
  });

  const panes = new Splitter([declarationPane(form.root),
    resultsFrame([failure, outdated, start, content], overlay, bar)],
  {direction: 'horizontal', sizes: [0.2, 0.8], minSize: DECLARATION_MIN});
  // a session's values pasted or typed into the table: the forecasts follow, after a pause (one value, then the next)
  self.own(onTableEdits(form, computed(() => [result.value, error.value]),
    studyTable, () => auto.changed(TABLE_EDITS)));
  // the table's own subscriptions are not built under the host: dropped with it
  self.own(stopWatching);
  const head = studyHeader(form, 'monitor', computed(() => {
    const r = result.value;
    if (!r)
      return undefined;
    const colors = colorsOf(r);
    const day = cell('day').value;
    const df = studyTable();
    return {name: df?.name ?? '', groups: revealed.value ?
      Object.entries(colors).map(([name, color]) => ({name, color})) : [],
    status: [r.status.get('line', 0), day === undefined ? '' : dayText(day),
      revealed.value ? '' : 'groups hidden'].filter((x) => x).join(' · '),
    meta: metaOf(df), setMeta: df ? (m) => {
      const was = !!metaOf(df)?.blinded;
      setMeta(m);
      if (m.blinded !== was)
        record(m.blinded ? 'Blinded entry on' : 'Blinded entry off: the groups are shown');
    } : undefined, groupMenu: (g, e) => studyMenu(studyNamesOf(df), 'group', g, {causedBy: e})};
  }), (s) => void openStage(s));
  return {
    root: stagePage(head.row, panes),
    ribbon: [[head.steps]],
    path: '/monitor',
    status: computed(() => busy.value ? 'Forecasting the study…' : result.value ? String(status.value?.get('line', 0)) :
      ''),
  };
}

/** Plan's "Follow this study": a new measurements table carrying the plan, opened in Monitor. */
export async function followPlan(plan: SavedPlan): Promise<void> {
  await showStage('monitor', await monitorView(newStudyTable(plan)));
}

/** A .pmxstudy file opened from Files or a drop. */
export async function openStudyFile(content: string): Promise<DG.DataFrame[]> {
  await showStage('monitor', await monitorView(tableOf(parseStudy(content))));
  return [];
}
