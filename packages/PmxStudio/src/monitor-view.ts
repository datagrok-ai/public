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

import {exampleStudy} from './analyze';
import {BandChart} from './band-chart';
import {armColors} from './frames';
import {allocateGroups, MonitorFrames} from './monitor';
import {PLAN_TAG, SavedPlan, parsePlan} from './plan-file';
import {planInput} from './plan-vs-observed';
import {openStage} from './stages';
import {LOOKS_TAG, addLog, logOf, newStudyTable, parseStudy, saveStudy, studyOf, tableOf} from './study-file';
import {addRunButton, block, busyOverlay, cardGrid, dayText, declarationPane, DECLARATION_MIN, infoNote, noteLine, notes,
  notesBox, num, outputsOf, pct, RESULTS_MIN, resultsFrame, resultTabs, settleForm, sized, stageStrip, staleNote,
  statCard, swatch, volumeAxis, wrapCells} from './ui';

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
  const call = DG.Func.byName('PmxStudio:pmxMonitorStudy').prepare();
  const result = signal<MonitorFrames | undefined>(undefined);
  const error = signal<string | null>(null);
  const stale = signal(false);
  const busy = signal(false);
  const revealed = signal(false);
  let latest = 0;
  const form = funcForm(call, {layout: 'tall', onInputChanged: () => stale.value = result.peek() !== undefined});
  settleForm(form, {control: 'Auto (the vehicle group)', evaluationDay: 'Auto (endpoint reached)'});

  const studyTable = (): DG.DataFrame | null => {
    const v = call.getParamValue('study');
    return v instanceof DG.DataFrame ? v : null;
  };
  const run = async (prepare?: () => Promise<void>) => {
    busy.value = true;
    const mine = ++latest;
    try {
      await prepare?.();
      const c = await call.call();
      if (mine !== latest)
        return;
      error.value = null;
      stale.value = false;
      // the reviews stay with the study from its first run
      const df = studyTable();
      if (df && !df.getTag(LOOKS_TAG))
        df.setTag(LOOKS_TAG, String(call.getParamValue('looks') ?? ''));
      result.value = outputsOf<MonitorFrames>(c);
    } catch (e: any) {
      error.value = String(e?.message ?? e);
    } finally {
      busy.value = false;
    }
  };
  addRunButton(form, () => run(), busy);

  const setStudy = async (df: DG.DataFrame) => {
    grok.shell.addTableView(df);
    grok.shell.v = view;
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
    await setStudy(df);
  });
  void form.settled.then(() => {
    const study = form.getInput('study') as unknown as {actions?: {id: string, icon: string, tooltip: string,
      run: () => void}[]} | null;
    if (study?.actions) {
      study.actions = [...study.actions, {id: 'example', icon: 'flask', tooltip: 'Follow an example study: MC38 as ' +
        'Plan predicts it, halfway through', run: () => void openExample()}];
    }
    if (initial)
      void setStudy(initial);
  });
  const planField = planInput(form, call, {stale, error,
    result: computed(() => result.value && {analysis: result.value.status})});

  const start = divV([
    span('Open the study\'s measurements (one row per animal and day, day 0 the randomization), check the columns ' +
      'and run. Paste each session\'s values into the table and run again.'),
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
  const groupCols = (): [string, (r: AnimalRow) => string][] => revealed.value ? [['Group', (a) => a.group]] : [];
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
        ['Next expected', (a) => a.next_median === null ? '—' : `${mm3(a.next_lower)}–${mm3(a.next_upper)} mm³`]]);
    }
    fill(sheet, all, (a) => a.animal, [['Animal', (a) => a.animal], ...groupCols(),
      ['Last volume', (a) => mm3(a.last_volume)], ['Next expected', (a) => a.next_median === null ? '—' :
        `${mm3(a.next_lower)}–${mm3(a.next_upper)}`],
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
  const growthData = computed(() => {
    const r = result.value;
    if (!r || !revealed.value)
      return undefined;
    const names = [...new Set(r.days.col('group')!.toList() as string[])];
    const colors = armColors(names);
    colors[r.status.get('control', 0)] = colors['Control'];
    r.days.col('group')!.meta.colors.setCategorical(Object.fromEntries(names.map((n) => [n, colors[n]])));
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
      decisionBody.replaceChildren(infoNote('<b>The groups are hidden.</b> The chance of success tells how the groups ' +
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
  const saveButton = buttonWithIcon('Save study', 'save', () => {
    const df = studyTable();
    if (!df)
      return;
    void saveStudy(studyOf(df, call.getParamValue('humane') ?? null)).then((p) => grok.shell.info(`Study saved: ${p}`),
      (e) => error.value = String(e?.message ?? e));
  });
  saveButton.title = 'Save the plan, the measurements and the log to your files';
  const analyzeButton = buttonWithIcon('Analyze this study', 'chart-line', () => void openStage('analyze'));
  analyzeButton.title = 'Open Analyze, where this table can be analyzed as planned';

  const tabs = resultTabs('monitor', [
    {id: 'today', label: 'Today', intro: 'Who needs attention, and the range each animal\'s next value should fall in.',
      content: [block('Needs attention', 'A person decides what to do', attention), noteBox,
        block('Measurement sheet', 'Check each value as you measure it', sheet)]},
    {id: 'enrollment', label: 'Enrollment', intro: 'Random groups alike in their starting tumors.', content: [
      divH([proposeButton]), allocationTable]},
    {id: 'growth', label: 'Growth', intro: 'Each group so far and what is expected up to the end.', content: [
      growthHidden, growthAxis.row, sized(growth, 340)]},
    {id: 'decisions', label: 'Decisions', intro: 'Whether each treated group is still worth continuing, and what was ' +
      'decided.', content: [decisionBody, block('Log', 'Who did what, and when', logBox),
      divH([saveButton, analyzeButton])]},
  ]);

  const content = divV([cards, tabs.root], 'pmx-monitor-content');
  Object.assign(content.style, {gap: 'var(--dg-space-l)', minWidth: `${RESULTS_MIN}px`});
  form.effect(() => {
    content.hidden = !result.value;
    start.hidden = !!result.value;
  });
  const overlay = busyOverlay('Forecasting the study…');
  const failure = div([], 'pmx-monitor-error');
  Object.assign(failure.style, {color: 'var(--dg-error-color, #c62828)'});
  form.effect(() => overlay.hidden = !busy.value);
  form.effect(() => {
    failure.hidden = error.value === null;
    failure.replaceChildren(...(error.value === null ? [] : [noteLine(error.value, 'error')]));
  });
  const outdated = staleNote('The inputs changed: run to update the forecasts');
  form.effect(() => {
    outdated.hidden = !stale.value;
    content.style.opacity = stale.value ? '0.5' : '';
  });

  const panes = new Splitter([declarationPane(form.root), resultsFrame([failure, outdated, start, content], overlay)],
    {direction: 'horizontal', sizes: [0.3, 0.7], minSize: DECLARATION_MIN});
  const host = ViewHost.create({name: 'PMX Studio · Monitor', content: () => ({
    root: panes,
    ribbon: [[stageStrip('monitor', (s) => void openStage(s))]],
    path: '/monitor',
    status: computed(() => busy.value ? 'Forecasting the study…' : result.value ? String(status.value?.get('line', 0)) :
      ''),
  })});
  [form, planField, growthAxis.input, growth, tabs].forEach((c) => host.own(() => c.dispose()));
  const view: DG.ViewBase = host.view;
  return view;
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
