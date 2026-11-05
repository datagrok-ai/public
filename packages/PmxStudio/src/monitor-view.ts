/* The Monitor view, built as Analyze is: the form generated from pmxMonitorStudy's signature, the results from its
   outputs. The groups stay hidden until someone reveals them, and that goes in the study's log. A study's measurements
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
import {armColors} from './frames';
import {allocateGroups, MonitorFrames, monitorFrames, MonitorInputs, outsideBand} from './monitor';
import {PLAN_TAG, SavedPlan, parsePlan} from './plan-file';
import {planInput} from './plan-vs-observed';
import {NAMES, openStage, showStage} from './stages';
import {LIMITS_TAG, LOOKS_TAG, WEIGHT_TAG, addLog, logOf, newStudyTable, parseStudy, restoreFields, saveStudy,
  studyOf, tableOf} from './study-file';
import {panelGrid, Slot} from './panels';
import {studyHeader} from './study-header';
import {emptyStudy, studyDataLine, tagStudy, withoutExcluded} from './study-table';
import {emptyMeta, metaOf, setMetaOn, withId} from './study-meta';
import {studyInputsOf} from './explore';
import {tolerabilityTab} from './tolerability';
import {askedOf, autoRun, block, busyBar, busyOverlay, callInputs, cardGrid, dayText, declarationPane, DECLARATION_MIN,
  infoNote,
  onTableEdits, TABLE_EDITS,
  isBodyWeight, keepNamedColumns, noteLine, notes,
  notesBox, num, pct, RESULTS_MIN, resultsFrame, settleForm, sized, staleNote,
  signed, statCard, swatch, volumeAxis, wrapCells} from './ui';

/** Monitor's panels, in order and width (twelfths): today's work first, then the groups, then the start. */
const MONITOR_LAYOUT: Slot[] = [{id: 'today', span: 4}, {id: 'growth', span: 8}, {id: 'sheet', span: 12},
  {id: 'decisions', span: 6}, {id: 'welfare', span: 6}, {id: 'enrollment', span: 12}];

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
  return ViewHost.create({name: 'PMX Studio · Monitor', content: (self) => monitorContent(self, initial)}).view;
}

/** What the view shows, built under its host, which owns it all (the form, the charts, the panes). */
function monitorContent(self: ViewHost, initial?: DG.DataFrame): ViewContent {
  const call = DG.Func.byName('PmxStudio:pmxMonitorStudy').prepare();
  const result = signal<MonitorFrames | undefined>(undefined);
  const error = signal<string | null>(null);
  const stale = signal(false);
  const busy = signal(false);
  const revealed = signal(false);
  let latest = 0;
  // no Run: each change forecasts again once the changes stop
  let preparing = false;
  const auto = autoRun(() => form, () => preparing ? undefined : run());
  const form = funcForm(call, {layout: 'tall', onInputChanged: () => {
    stale.value = result.peek() !== undefined;
    auto.changed();
  }});
  settleForm(form, {control: 'Auto (the vehicle group)', evaluationDay: 'Auto (endpoint reached)'});
  keepNamedColumns(form, {bodyWeight: isBodyWeight});

  const studyTable = (): DG.DataFrame | null => {
    const v = call.getParamValue('study');
    return v instanceof DG.DataFrame ? v : null;
  };
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
    // ids that say nothing of the group, as a blinded study's
    const ids = [...new Set(df.col('animal')!.toList() as string[])].sort(() => Math.random() - 0.5);
    const coded = new Map(ids.map((a, i) => [a, `M${String(i + 1).padStart(3, '0')}`]));
    df.col('animal')!.init((i) => coded.get(df.get('animal', i))!);
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
    infoNote('<b>The groups stay hidden.</b> Monitor shows each animal without its group, so that whoever measures ' +
      'is not led by what they expect. Revealing them is a click, and it goes in the study\'s log.'),
  ], 'pmx-monitor-start');
  Object.assign(start.style, {gap: 'var(--dg-space-m)', maxWidth: '760px'});
  (start.children[1] as HTMLElement).style.gap = 'var(--dg-space-m)';
  (start.children[1] as HTMLElement).style.alignItems = 'center';

  const status = computed(() => result.value?.status);
  const cell = <T = number>(name: string) => computed(() => {
    const s = status.value;
    return (s && !s.col(name)?.isNone(0) ? s.get(name, 0) ?? undefined : undefined) as T | undefined;
  });
  const card = (label: string, value: ReadonlySignal<StatValue | undefined>, tip: () => string,
    format?: (x: StatValue) => string) => statCard(form, label, value, tip, format);
  const cards = cardGrid([
    card('Animals needing attention', cell('attention'), () => 'Animals listed under Today, each with its reason'),
    card('Next measurement', cell('next_day'), () => 'The next day at the study\'s usual interval', dayText),
    card('Expected end', cell('end_median'), () => {
      const [a, b] = [cell('end_lower').value, cell('end_upper').value];
      // not "8 in 10": halfway through, the atlas's MC38 studies fell in this range about half the time (§8)
      return `The day the median control reaches the endpoint, likely between ${dayText(a)} and ` +
        `${b === undefined ? 'later than the days forecast' : dayText(b)}`;
    }, dayText),
    card('End in the plan', cell('planned_end'), () => 'The last day the plan predicted (empty: no plan)', dayText),
  ], 'pmx-monitor-cards', 150);

  // Today
  const attention = div([], 'pmx-monitor-attention');
  const sheet = div([], 'pmx-monitor-sheet');
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

  // Growth: by group, once revealed
  const growthAxis = volumeAxis();
  const colorsOf = (r: MonitorFrames) => {
    const names = [...new Set(r.days.col('group')!.toList() as string[])];
    const colors = armColors(names);
    colors[r.status.get('control', 0)] = colors['Control'];
    return Object.fromEntries(names.map((n) => [n, colors[n]]));
  };
  const growthData = computed(() => {
    const r = result.value;
    if (!r || !revealed.value)
      return undefined;
    r.days.col('group')!.meta.colors.setCategorical(colorsOf(r));
    return r.days;
  });
  const growth = new BandChart(growthData, {x: 'day', lower: 'lower', upper: 'upper', line: 'median', split: 'group',
    xTitle: 'Study day', yTitle: 'Median volume (mm³)', log: growthAxis.log, markers: true});
  const growthHidden = infoNote('<b>The groups are hidden.</b> Reveal them to see each group\'s growth.');
  form.effect(() => {
    growth.root.hidden = !revealed.value;
    growthAxis.row.hidden = !revealed.value;
    growthHidden.hidden = revealed.value;
  });

  // Decisions
  const decisionBody = div([], 'pmx-monitor-decisions');
  const logBox = div([], 'pmx-monitor-log');
  const logVersion = signal(0);
  const record = (what: string) => {
    const df = studyTable();
    if (df) {
      addLog(df, what);
      logVersion.value++;
    }
  };
  const reveal = () => {
    revealed.value = true;
    record('Groups revealed');
  };
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
    fill(decisionBody, d, (a) => a.arm, [['Group', (a) => swatch(a.arm, armColors([a.arm])[a.arm])],
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
    if (!revealed.value)
      reveal();
    // in place of an Analyze view already open
    [...grok.shell.views].find((v) => v.name === NAMES.analyze)?.close();
    void analyzeView({study: df}).then((v) => showStage('analyze', v));
  });
  analyzeButton.title = 'Analyze this table as planned; this reveals the groups, and the log says so';

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
  const panels = panelGrid([
    {id: 'today', title: 'Needs attention', hint: 'Who needs attention today; a person decides what to do',
      content: [limits, entered, attention, noteBox]},
    {id: 'growth', title: 'Growth', hint: 'Each group so far and what is expected up to the end',
      content: [growthHidden, growthAxis.row, sized(growth, 340)]},
    {id: 'sheet', title: 'Measurement sheet', hint: 'The range the next value of each animal should fall in: check ' +
      'each value as you measure it', content: [sheet]},
    {id: 'decisions', title: 'Decisions', hint: 'Whether each treated group is still worth continuing, and what was ' +
      'decided', content: [decisionBody, block('Log', 'Who did what, and when', logBox), divH([saveButton,
      analyzeButton])]},
    {id: 'welfare', title: 'Welfare', hint: 'Whether the treatments are being tolerated, from the body weights so far',
      content: [welfareHidden, welfareBody]},
    {id: 'enrollment', title: 'Enrollment', hint: 'Random groups alike in their starting tumors',
      content: [divH([proposeButton]), allocationTable]},
  ], MONITOR_LAYOUT, 'pmx-monitor-panels');

  const content = divV([dataLine.root, cards, panels], 'pmx-monitor-content');
  Object.assign(content.style, {gap: 'var(--dg-space-l)', minWidth: `${RESULTS_MIN}px`});
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
  {direction: 'horizontal', sizes: [0.3, 0.7], minSize: DECLARATION_MIN});
  // a session's values pasted or typed into the table: the forecasts follow, after a pause (one value, then the next)
  self.own(onTableEdits(form, result, studyTable, () => auto.changed(TABLE_EDITS)));
  // the table's own subscriptions are not built under the host: dropped with it
  self.own(stopWatching);
  return {
    root: panes,
    ribbon: [[studyHeader(form, 'monitor', computed(() => {
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
      meta: metaOf(df), setMeta: df ? (m) => setMetaOn(df, m) : undefined};
    }), (s) => void openStage(s))]],
    path: '/monitor',
    status: computed(() => busy.value ? 'Forecasting the study…' : result.value ? String(status.value?.get('line', 0)) :
      ''),
  };
}

/** Plan's "Follow this study": a new measurements table carrying the plan, opened in Monitor. */
export async function followPlan(plan: SavedPlan): Promise<void> {
  const view = grok.shell.addView(await monitorView(newStudyTable(plan)));
  view.path = '/apps/PmxStudio/monitor';
}

/** A .pmxstudy file opened from Files or a drop. */
export async function openStudyFile(content: string): Promise<DG.DataFrame[]> {
  const view = grok.shell.addView(await monitorView(tableOf(parseStudy(content))));
  view.path = '/apps/PmxStudio/monitor';
  return [];
}
