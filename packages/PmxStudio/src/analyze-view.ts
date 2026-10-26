/* The Analyze view, built as the Plan view is: the declaration form generated from pmxAnalyzeStudy's own signature
   (u2 funcForm: the study table, its columns and the analysis settings), the results from the function's outputs, in
   the order a scientist asks: did each treatment work (the figures, then the groups with their result), how did the
   tumors grow (band-chart.ts), which animals responded (each treated animal against its untreated twin, the
   animals listed by effect, the one selected drawn against its twin), then, folded, what serves a closer look: the
   statistics behind each result, the control check and the analysis protocol; the study against its plan has a tab
   of its own (plan-vs-observed.ts). With no table open, an example study predicted by Plan's model is opened and
   analyzed. */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {badge, BasicTable, button, computed, div, divH, divV, dot, signal, span, Splitter, tableFromMap,
  Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal, StatValue} from '@datagrok-libraries/u2';
import {appView, funcForm, viewers} from '@datagrok-libraries/u2/src/dg/index.js';

import {AnalysisFrames, exampleStudy, median} from './analyze';
import {armColors} from './frames';
import {CorrectionInput} from './arms-input';
import {axisAt, block, busyOverlay, cardGrid, copyRow, dayText, declarationPane, DECLARATION_MIN, formRunning, guide,
  infoNote, noteLine, notesBox, num, outputsOf, resultTabs, pct, pValue, RESULTS_MIN, resultsFrame, settleForm, signed,
  sized, staleNote, stageStrip, statCard, swatch, volumeAxis, wrapCells} from './ui';
import {openStage} from './stages';
import {noEffectNote, observedCell, planInput, planTab} from './plan-vs-observed';
import {nciReadings, NciReading, tolerabilityTab} from './tolerability';
import {BandChart} from './band-chart';

const TWIN = 'Untreated twin';
const MEASURED = 'Measured';

/** One treated animal's measurements and its twin's band, as the band chart draws two series. */
function animalFrame(twins: DG.DataFrame, animal: string, color: string): DG.DataFrame {
  const out: {day: number, series: string, lower: number, upper: number, median: number}[] = [];
  for (let i = 0; i < twins.rowCount; i++) {
    if (twins.get('animal', i) !== animal)
      continue;
    const day = twins.get('day', i);
    out.push({day, series: TWIN, lower: twins.get('twin_q05', i), upper: twins.get('twin_q95', i),
      median: twins.get('twin_median', i)});
    const v = twins.get('observed', i);
    if (v !== null && v !== undefined && !Number.isNaN(v))
      out.push({day, series: MEASURED, lower: v, upper: v, median: v});
  }
  const df = DG.DataFrame.fromObjects(out)!;
  // the twin in the control's slate, the animal in its group's colour
  df.col('series')!.meta.colors.setCategorical({[TWIN]: '#9aa5b1', [MEASURED]: color});
  df.col('day')!.meta.units = 'd';
  for (const q of ['lower', 'upper', 'median']) {
    df.col(q)!.meta.units = 'mm³';
    df.col(q)!.meta.format = '0.0';
  }
  return df;
}

/** A treated animal against its twin on the evaluation day, in words: smaller (or larger) by how much, and its 90%
 *  range. */
function effectWords(effect: number, lower: number, upper: number): string {
  const side = (x: number) => x >= 0 ? `${pct(x / 100)} smaller` : `${pct(-x / 100)} larger`;
  return `${side(effect)} than its untreated twin (90% range: ${side(lower)} to ${side(upper)})`;
}

/** An animal's effect and its 90% range as a bar on a shared axis from `lo` to 100% (the twin's size at 0): the dot
 *  filled for a responder, hollow otherwise. */
function effectBar(effect: number, lower: number, upper: number, lo: number, color: string): HTMLElement {
  const track = div([], 'pmx-effect-bar');
  Object.assign(track.style, {position: 'relative', height: '16px'});
  const range = div([]);
  Object.assign(range.style, {position: 'absolute', left: axisAt(lower, lo),
    width: `calc(${axisAt(upper, lo)} - ${axisAt(lower, lo)})`, top: '7px', height: '2px', background: color});
  const point = div([]);
  Object.assign(point.style, {position: 'absolute', left: `calc(${axisAt(effect, lo)} - 5px)`, top: '3px',
    width: '10px', height: '10px', borderRadius: '50%', boxSizing: 'border-box', border: `2px solid ${color}`,
    background: lower > 0 ? color : 'var(--dg-background, #fff)'});
  track.append(guide(0, lo), guide(50, lo), guide(100, lo), range, point);
  track.setAttribute('role', 'img');
  track.setAttribute('aria-label', effectWords(effect, lower, upper));
  return track;
}

/** The effect axis over the bars: a label over a tick at the left end (larger than the twin), at the twin (0), and at
 *  half its size and no tumor at all, in words; each tick over its guide in the bars. */
function effectAxis(lo: number): HTMLElement {
  const axis = div([], 'pmx-effect-axis');
  Object.assign(axis.style, {position: 'relative', height: '22px', fontSize: '11px', fontWeight: '400',
    color: 'var(--dg-text-color-muted, #68727f)'});
  for (const [x, text, shift] of [[lo, 'Larger', '0'], [0, 'Same as twin', '-50%'], [50, 'Half size', '-50%'],
    [100, 'No tumor', '-100%']] as [number, string, string][]) {
    const label = span(text);
    Object.assign(label.style, {position: 'absolute', top: '0', left: axisAt(x, lo), transform: `translateX(${shift})`,
      whiteSpace: 'nowrap', lineHeight: '14px'});
    const tick = div([]);
    Object.assign(tick.style, {position: 'absolute', bottom: '0', left: axisAt(x, lo), width: '1px', height: '5px',
      background: 'var(--dg-text-color-muted, #9aa3ad)'});
    axis.append(label, tick);
  }
  return axis;
}

const correctionWords = (corr: unknown, k: number) => k < 2 ? 'one comparison' : corr === 'dunnett' ?
  'Dunnett\'s correction' : corr === 'bonferroni' ? 'Bonferroni\'s correction' : 'no correction';

/** A group's result, with a red "!" after it when the group is over the NCI limits for tolerability: the criteria on
 *  hover, the Tolerability tab on click. */
function withNci(el: HTMLElement, x: NciReading | undefined, open: () => void): HTMLElement {
  if (x?.kind !== 'error')
    return el;
  const mark = span('!', 'pmx-nci-mark');
  const words = `Over the NCI limits for tolerability: ${x.criteria}.`;
  Object.assign(mark.style, {display: 'inline-flex', alignItems: 'center', justifyContent: 'center', width: '14px',
    height: '14px', borderRadius: '50%', flex: '0 0 auto', background: 'var(--dg-error, #EB6767)', color: '#fff',
    fontSize: '10px', fontWeight: '700', lineHeight: '1', cursor: 'pointer'});
  mark.setAttribute('role', 'button');
  mark.setAttribute('aria-label', `${words} Open Tolerability`);
  mark.tabIndex = 0;
  Tooltip.bind(mark, () => `${words} Click to see Tolerability.`);
  mark.addEventListener('click', open);
  mark.addEventListener('keydown', (ev) => {
    if (ev.key === 'Enter' || ev.key === ' ') {
      ev.preventDefault();
      open();
    }
  });
  const both = divH([el, mark]);
  Object.assign(both.style, {alignItems: 'center', gap: 'var(--dg-space-s)'});
  return both;
}

/** The results in words, a line each: the study, each treated group's effect and result, the test, the body weight. */
function summaryLines(r: AnalysisFrames, nci: Map<string, NciReading>): string[] {
  const v = (k: string) => r.analysis.get(k, 0);
  const c = r.comparisons;
  const control = v('control');
  const lines = [`${v('line')} study evaluated on day ${v('evaluation_day')}` + (v('from_endpoint') ?
    `, the first measurement day the median control reached ${v('endpoint')} mm³` : '') +
    `; control: ${control}, ${v('n_control')} animals.`];
  for (let i = 0; i < c.rowCount; i++) {
    const g = (k: string) => c.get(k, i);
    const p = pValue(g('p_value'));
    const tc = g('tc') === null ? '' : `, T/C ${pct(g('tc') / 100)}` + (g('tc_lower') === null ? '' :
      ` (95% CI ${Math.round(g('tc_lower'))}–${Math.round(g('tc_upper'))}%)`);
    lines.push(`${g('arm')} (${g('n')} animals): TGI ${pct(g('tgi') / 100)}${tc}; ${g('activity').toLowerCase()} ` +
      `against the control (p ${p.startsWith('<') ? p : `= ${p}`}); ${g('responders')} of ${g('twin_animals')} ` +
      'animals responded against their untreated twin.');
  }
  lines.push(`Tumor growth of each treated group against the control, two-sided α = ${num(v('alpha'), 2)} with ` +
    `${correctionWords(v('correction'), v('treated_arms'))}.`);
  if (nci.size) {
    const order = [control, ...c.col('arm')!.toList() as string[]].filter((g) => nci.has(g));
    const over = order.filter((g) => nci.get(g)!.kind === 'error');
    const nadirs = order.map((g) => {
      const x = nci.get(g)!.nadir;
      return `${g} ${x === null ? '—' : signed(x)}`;
    });
    lines.push(`Body weight, lowest mean change from day 0: ${nadirs.join(', ')}. ` + (over.length ?
      over.map((g) => `${g} over the NCI limits for tolerability: ${nci.get(g)!.criteria}.`).join(' ') :
      'No group over the NCI limits for tolerability.') + ' (NCI limits: Plowman et al. 1997.)');
  }
  return lines;
}

export async function analyzeView(): Promise<DG.ViewBase> {
  const call = DG.Func.byName('PmxStudio:pmxAnalyzeStudy').prepare();
  const read = (c: DG.FuncCall) => outputsOf<AnalysisFrames>(c);
  const result = signal<AnalysisFrames | undefined>(undefined);
  // the newest run: the form's Run and the example share the call, and a run that finishes after a newer one started
  // is dropped
  let latest = 0;
  const error = signal<string | null>(null);
  // a result that no longer matches the form is dimmed until the next run
  const stale = signal(false);
  // the multiple comparisons as Plan asks for them, over the treated groups of the last analysis (two until then)
  const treatedGroups = signal(2);
  const correction = new CorrectionInput({label: 'Multiple comparisons', arms: treatedGroups});
  const form = funcForm(call, {showRun: true, layout: 'tall', overrides: {correction: {input: correction}},
    onInputChanged: () => stale.value = result.peek() !== undefined,
    onRun: (c) => {
      latest++;
      error.value = null;
      stale.value = false;
      result.value = read(c as DG.FuncCall);
    }});
  settleForm(form, {control: 'Auto (the vehicle group)', evaluationDay: 'Auto (endpoint reached)',
    density: 'Auto (the cell line’s in the atlas)'});
  // u2 auto-picks a column for every column param, the optional ones too (the nearest name): an optional column keeps
  // the pick only when its name says what it is, so a study without body weights does not read its volumes as weights
  const optional: {[input: string]: (column: string) => boolean} = {
    bodyWeight: (c) => /weight|\bbw\b/i.test(c) && !/tum/i.test(c),
    tumorWeight: (c) => /tum/i.test(c) && /weight|mass/i.test(c),
    removalReason: (c) => /reason|removal|fate|cause|status/i.test(c),
  };
  form.effect(() => {
    void form.getInput('study')?.value.value;
    void form.settled.then(() => {
      for (const [name, named] of Object.entries(optional)) {
        const input = form.getInput(name);
        const column = input?.value.peek();
        if (input && typeof column === 'string' && !named(column))
          input.value.value = null;
      }
    });
  });
  // busy while the form runs (funcForm says so in aria-busy) and while the example study is prepared
  const running = formRunning(form);
  const preparing = signal(false);
  const busy = computed(() => running.value || preparing.value);

  // the example: Plan's study opened in the workspace like any table (the form's table picker then offers it), this
  // view brought back and the study analyzed at once
  const openExample = async () => {
    preparing.value = true;
    const mine = ++latest;
    try {
      const df = await exampleStudy();
      grok.shell.addTableView(df);
      grok.shell.v = view;
      const study = form.getInput('study');
      if (study)
        study.value.value = df.name;
      await form.settled;
      // the optional body and tumor weights: the example's own columns
      for (const [input, column] of [['bodyWeight', 'body_weight'], ['tumorWeight', 'tumor_weight']]) {
        const w = form.getInput(input);
        if (w)
          w.value.value = column;
      }
      await form.settled;
      const c = await call.call();
      if (mine !== latest)
        return;
      error.value = null;
      stale.value = false;
      result.value = read(c);
    } catch (e: any) {
      error.value = String(e?.message ?? e);
    } finally {
      preparing.value = false;
    }
  };
  // the example also among the Study field's own actions (beside open file, Files, database), within reach once a
  // study is open
  void form.settled.then(() => {
    const study = form.getInput('study') as unknown as {actions?: {id: string, icon: string, tooltip: string,
      run: () => void}[]} | null;
    if (study?.actions) {
      study.actions = [...study.actions, {id: 'example', icon: 'flask', tooltip: 'Analyze an example study: MC38 ' +
        'as Plan predicts it, a vehicle and two treated groups', run: () => void openExample()}];
    }
  });
  const planField = planInput(form, call, {stale, error, result});

  const start = divV([
    span('Open your study table (one row per animal and measurement day, day 0 the randomization), check its ' +
      'columns and run.'),
    divH([button('Analyze an example study', () => void openExample()),
      span('An MC38 study as Plan predicts it: a vehicle and two treated groups at TGI 40% and 70%, 8 animals ' +
        'each.')]),
    infoNote('<b>Analyzed as planned.</b> Analyze applies the analysis Plan pre-specifies: the study is evaluated ' +
      'when the median control reaches the endpoint, each treated group is compared with the control on tumor ' +
      'growth, and each treated animal with its untreated twin: what it would measure without treatment, predicted ' +
      'from its own tumor at randomization, this study\'s controls and the cell line\'s growth.'),
  ], 'pmx-analyze-start');
  Object.assign(start.style, {gap: 'var(--dg-space-m)', maxWidth: '760px'});
  (start.children[1] as HTMLElement).style.gap = 'var(--dg-space-m)';
  (start.children[1] as HTMLElement).style.alignItems = 'center';

  const analysis = computed(() => result.value?.analysis);
  const cell = <T = number>(name: string) =>
    computed(() => (analysis.value?.get(name, 0) ?? undefined) as T | undefined);
  const comparisons = computed(() => result.value?.comparisons);
  const significant = computed(() => {
    const c = comparisons.value;
    return c ? `${c.col('activity')!.toList().filter((x) => x === 'Slows growth').length} of ${c.rowCount}` :
      undefined;
  });
  const treatedAnalyzed = computed(() => {
    const c = comparisons.value;
    return c ? (c.col('n')!.toList() as number[]).reduce((t, x) => t + x, 0) : undefined;
  });
  const analyzed = computed(() => treatedAnalyzed.value === undefined ? undefined :
    treatedAnalyzed.value + (cell('n_control').value ?? 0));
  // the study's answers first, what each figure covers in its tooltip
  const card = (label: string, value: ReadonlySignal<StatValue | undefined>, tip: () => string,
    format?: (x: StatValue) => string) => statCard(form, label, value, tip, format);
  const cards = cardGrid([
    card('Treated groups that slow growth', significant, () => 'Treated groups whose tumors grew slower than the ' +
      `control's up to day ${cell('evaluation_day').value}`),
    card('Evaluation day', cell('evaluation_day'), () => cell('from_endpoint').value ? 'The first measurement day ' +
      `the median control reached ${cell('endpoint').value} mm³, as Plan defines the end of the study` :
      'Chosen, or the last measurement day if the median control never reached the endpoint', dayText),
    card('Animals analyzed', analyzed, () => `${cell('n_control').value} controls and ${treatedAnalyzed.value} ` +
      'treated animals with at least two measurements'),
    card('Controls on track', cell('control_coverage'), () => 'Control measurements inside the range the other ' +
      'controls predict for them (about 90% is as expected)', (x) => pct(x)),
  ], 'pmx-analyze-cards', 150);

  const groupsTable = div([]);
  const groupNotes = notesBox();
  const statsTable = div([]);
  const controlCheck = div([]);
  const protocol = div([]);
  let protocolText = '';
  const summary = divV([], 'pmx-results-summary');
  summary.style.gap = 'var(--dg-space-s)';
  let summaryText = '';
  // the Tolerability tab, at the groups against the NCI limits (the tabs are built below)
  const showNci = () => {
    tabs.activeTab.value = 'tolerability';
    requestAnimationFrame(() => tolerability.content.map((el) => el.querySelector('.pmx-nci')).find((el) => el)
      ?.scrollIntoView({behavior: 'smooth', block: 'center'}));
  };
  const colorsOf = (r: AnalysisFrames) => {
    const colors = armColors(r.comparisons.col('arm')!.toList());
    colors[r.analysis.get('control', 0)] = colors['Control'];
    return colors;
  };
  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    const v = (c: string) => r.analysis.get(c, 0);
    treatedGroups.value = v('treated_arms');
    const control = v('control');
    const colors = colorsOf(r);
    const allNotes = String(v('notes') ?? '').split(' | ').filter((n) => n);
    const odd = allNotes.filter((n) => n.startsWith('Control '));
    groupNotes.replaceChildren(...allNotes.filter((n) => !odd.includes(n)).map((n) => noteLine(n)));

    const nci = nciReadings(r);

    // the control first, then each treated group: its effect on the evaluation day, its animals against their
    // twins and its result (a table built in this run is disposed with it)
    const c = r.comparisons;
    type Row = {name: string, cells: (string | HTMLElement)[]};
    const rows: Row[] = [{name: control, cells: [swatch(control, colors[control]), String(v('n_control')), '—', '—',
      '—', '—', '—', withNci(span('Reference'), nci.get(control), showNci)]}];
    for (let i = 0; i < c.rowCount; i++) {
      const g = (k: string) => c.get(k, i);
      const tc = g('tc');
      const verdict = divH([dot(g('activity') === 'Slows growth' ? 'success' : 'warning'), span(g('activity'))]);
      Object.assign(verdict.style, {alignItems: 'center', gap: 'var(--dg-space-s)'});
      // no clear effect: what the study could have detected with its animals
      const why = noEffectNote((k) => r.versusPlan.get(k, i), v);
      if (why)
        Tooltip.bind(verdict, why);
      rows.push({name: g('arm'), cells: [swatch(g('arm'), colors[g('arm')]), String(g('n')), pct(g('tgi') / 100),
        tc === null ? '—' : observedCell(pct(tc / 100), g('tc_lower') === null ? '' :
          `(${Math.round(g('tc_lower'))}–${Math.round(g('tc_upper'))}%)`),
        g('paired_effect') === null ? '—' : pct(g('paired_effect') / 100),
        `${g('responders')} of ${g('twin_animals')}`, pValue(g('p_value')).replace(' ', ' '),
        withNci(verdict, nci.get(g('arm')), showNci)]});
    }
    const headers = ['Group', 'Analyzed', 'TGI', 'T/C (95% CI)', 'Smaller than its twin', 'Responders', 'p',
      'Result'];
    const table = new BasicTable<Row>({items: rows, keyOf: (x) => x.name, columns: headers.map((header, j) => ({header,
      render: (x: Row) => x.cells[j], align: j >= 1 && j <= 6 ? 'center' : 'left'}))});
    table.root.style.width = '100%';
    wrapCells(table.root);
    groupsTable.replaceChildren(table.root);

    // the statistics behind each result (the growth rates the test compares and their difference)
    const reached = (d: unknown) => d === null || d === undefined ? 'Not reached' : dayText(d);
    type Stat = {name: string, cells: string[]};
    const stats: Stat[] = [{name: control, cells: [control, String(num(c.get('slope_control_per_day', 0))), '—', '—',
      reached(v('control_median_tte'))]}];
    for (let i = 0; i < c.rowCount; i++) {
      const g = (k: string) => c.get(k, i);
      stats.push({name: g('arm'), cells: [g('arm'), String(num(g('slope_per_day'))),
        `${num(g('slope_difference'), 2)} (${num(g('slope_difference_lower'), 2)} to ` +
        `${num(g('slope_difference_upper'), 2)})`, pValue(g('p_value')), reached(g('median_time_to_endpoint'))]});
    }
    const statHeaders = ['Group', 'Growth rate (ln V per day)', 'Difference from the control (95% CI)', 'p',
      `Median day to ${v('endpoint')} mm³`];
    const st = new BasicTable<Stat>({items: stats, keyOf: (x) => x.name, columns: statHeaders.map((header, j) =>
      ({header, render: (x: Stat) => x.cells[j], align: j >= 1 ? 'center' : 'left'}))});
    st.root.style.width = '100%';
    wrapCells(st.root);
    statsTable.replaceChildren(st.root);

    // the controls against each other: how many measurements fall in their range, and the ones that grow apart
    const coverage = v('control_coverage');
    controlCheck.replaceChildren(divV([
      span(coverage === null ? 'Too few controls to check them against each other.' :
        `${pct(coverage)} of the control measurements fall in the range the other controls predict for them ` +
        `(about 90% when the controls grow as the ${v('line')} model expects).`),
      ...(odd.length ? odd.map((n) => noteLine(n)) : [noteLine('No control grows apart from the others.', 'success')]),
    ], 'pmx-control-check'));
    (controlCheck.firstElementChild as HTMLElement).style.gap = 'var(--dg-space-s)';

    // the analysis protocol, as Plan pre-specifies it
    const k = v('treated_arms');
    const corr = v('correction');
    const rowsOf: {[row: string]: string} = {
      'Primary comparison': 'Each treated group against the control on tumor growth: the per-animal slope of log ' +
        'volume from randomization to the evaluation day, by a two-sample t test',
      'Multiple comparisons': k < 2 ? 'One treated group: none needed' : corr === 'dunnett' ?
        `Dunnett's many-to-one test over ${k} comparisons` : corr === 'bonferroni' ?
          `Bonferroni over ${k} comparisons` : `None: each of the ${k} comparisons at α = ${num(v('alpha'), 2)} ` +
          '(exploratory)',
      'Significance': `α = ${num(v('alpha'), 2)}, two-sided; the intervals are adjusted as the p-values`,
      'Evaluation day': v('from_endpoint') ? `Day ${v('evaluation_day')}, the first measurement day the median ` +
        `control reached ${v('endpoint')} mm³, as Plan defines the end of the study` : `Day ${v('evaluation_day')}`,
      'Effects': 'TGI = 1 − ΔT/ΔC, each animal\'s change from its own volume at randomization, and T/C of the ' +
        'mean volumes, on the evaluation day (NCI definitions, as declared in Plan), each with its 95% interval',
      'Untreated twin': `What each treated animal would measure untreated, predicted from its own tumor at ` +
        `randomization, this study's controls and the ${v('line')} growth model. Its effect: 1 − measured / twin on ` +
        'the evaluation day; a responder\'s 90% range lies above 0. Supporting the primary comparison, ' +
        'not replacing it',
      'Small tumors': `Volumes below ${num(v('v_lod'), 2)} mm³, the smallest a caliper records, are counted at it`,
    };
    protocol.replaceChildren(tableFromMap(rowsOf));
    protocolText = [`Analysis protocol (${r.analysis.name})`, ...Object.entries(rowsOf).map(([key, t]) =>
      `${key}: ${t}`)].join('\n');

    const lines = summaryLines(r, nci);
    summary.replaceChildren(...lines.map((l) => span(l)));
    summaryText = lines.join('\n');
  });

  const versus = planTab(form, result, colorsOf);
  const tolerability = tolerabilityTab(form, result, colorsOf);

  const growthAxis = volumeAxis();
  const growth = new BandChart(computed(() => result.value?.volumes), {x: 'day', lower: 'lower', upper: 'upper',
    line: 'median', split: 'group', xTitle: 'Day after randomization', yTitle: 'Tumor volume', log: growthAxis.log,
    markers: true});

  // each treated animal against its twin: the animals listed by group, the strongest first, the one selected drawn
  type AnimalRow = {animal: string, group: string, effect: number, lower: number, upper: number};
  const picked = signal<string | undefined>(undefined);
  // what the badges mean, on each badge and under the list
  const day = () => cell('evaluation_day').value;
  const responderText = () => `Responder: on day ${day()} its tumor is below the range its untreated twin would ` +
    'reach';
  const largerText = () => `On day ${day()} its tumor is above the range its untreated twin would reach`;
  const explained = (el: HTMLElement, tip: () => string) => {
    Tooltip.bind(el, tip);
    return el;
  };
  // the legend, each line after the mark it explains, as the Plan calendar's summary
  const legendLine = (mark: HTMLElement, text: HTMLElement) => {
    const line = divH([mark, text]);
    Object.assign(line.style, {alignItems: 'baseline', gap: 'var(--dg-space-s)'});
    return line;
  };
  const responderWords = span('');
  form.effect(() => responderWords.textContent = result.value ? `On day ${day()} its tumor is below the range its ` +
    'untreated twin would reach: its 90% range lies right of the twin.' : '');
  const legendDot = (filled: boolean) => {
    const d = span('');
    Object.assign(d.style, {display: 'inline-block', width: '10px', height: '10px', borderRadius: '50%',
      boxSizing: 'border-box', border: '2px solid var(--dg-text-color-muted, #68727f)', flex: '0 0 auto',
      background: filled ? 'var(--dg-text-color-muted, #68727f)' : 'transparent'});
    return d;
  };
  const dots = divH([legendDot(true), legendDot(false)]);
  dots.style.gap = '3px';
  const responderNote = divV([
    legendLine(form.runInScope(() => explained(badge('Responder', {variant: 'success'}), responderText)),
      responderWords),
    legendLine(dots, span(`Each bar: how much smaller the treated tumor is than its twin's (the dot, filled for a ` +
      `responder) and its 90% range; the line at "Same as twin" is the tumor's untreated size.`)),
  ], 'pmx-responder-legend');
  Object.assign(responderNote.style, {gap: 'var(--dg-space-s)', margin: 'var(--dg-space-m) 0 0', fontSize: '12px',
    color: 'var(--dg-text-color-muted, #68727f)', lineHeight: '1.45'});
  // the rows' columns: the animal, its effect in words, the bar, its badge; one grid the rows share (subgrid), so
  // every row's bar sits under the axis and the badge column fits its widest badge, beside the chart too
  const animalList = div([], 'pmx-animal-list');
  Object.assign(animalList.style, {display: 'grid', columnGap: 'var(--dg-space-m)',
    gridTemplateColumns: 'auto auto minmax(140px, 1fr) max-content'});
  const gridRow = (cells: (HTMLElement | string)[], cls: string) => {
    const row = div(cells.map((c) => typeof c === 'string' ? span(c) : c), cls);
    (row.lastElementChild as HTMLElement).style.justifySelf = 'start';
    Object.assign(row.style, {display: 'grid', gridTemplateColumns: 'subgrid', gridColumn: '1 / -1',
      alignItems: 'center', padding: '3px var(--dg-space-s)'});
    return row;
  };
  const sideOf = (x: number) => x >= 0 ? `${pct(x / 100)} smaller` : `${pct(-x / 100)} larger`;
  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    const e = r.effects;
    const colors = colorsOf(r);
    const all: AnimalRow[] = Array.from({length: e.rowCount}, (_, i) => ({animal: e.get('animal', i),
      group: e.get('group', i), effect: e.get('effect', i), lower: e.get('effect_lower', i),
      upper: e.get('effect_upper', i)}));
    // the axis from 50% larger: the room goes to 0-100% where the effects are; a longer range is cut at the edge
    const lo = -50;
    const head = gridRow(['Animal', 'Against its twin', 'Tumor size against its twin, 90% range', ''],
      'pmx-forest-head');
    head.style.fontWeight = '500';
    const scale = gridRow(['', '', effectAxis(lo), ''], 'pmx-forest-axis');
    Object.assign(scale.style, {paddingTop: '0', borderBottom: '1px solid var(--dg-border-color, #e1e6ec)'});
    const rows: HTMLElement[] = [head, scale];
    const buttons = new Map<string, HTMLElement>();
    const mark = (animal: string | undefined) => buttons.forEach((b, a) => {
      b.setAttribute('aria-pressed', String(a === animal));
      b.style.background = a === animal ? 'color-mix(in srgb, var(--dg-accent, #2083d5) 12%, transparent)' : '';
    });
    for (const g of r.comparisons.col('arm')!.toList() as string[]) {
      const own = all.filter((a) => a.group === g).sort((a, b) => b.effect - a.effect);
      if (own.length === 0)
        continue;
      // the group's line: how many respond and its median animal against its twin
      const responders = own.filter((a) => a.lower > 0).length;
      const mid = median(own.map((a) => a.effect));
      const summary = gridRow([swatch(g, colors[g]), sideOf(mid), span(`${responders} of ${own.length} respond`), ''],
        'pmx-forest-group');
      Object.assign(summary.style, {fontWeight: '600', marginTop: 'var(--dg-space-s)',
        background: 'var(--dg-surface-variant, #f4f6f8)'});
      rows.push(summary);
      for (const a of own) {
        const badgeEl = a.lower > 0 ? explained(badge('Responder', {variant: 'success'}), responderText) :
          a.upper < 0 ? explained(badge('Larger than its twin', {variant: 'warning'}), largerText) : span('');
        const row = gridRow([span(a.animal), sideOf(a.effect), effectBar(a.effect, a.lower, a.upper, lo,
          colors[a.group]), badgeEl], 'pmx-forest-row');
        // a button: picked by click or keyboard, the one shown marked
        row.setAttribute('role', 'button');
        row.tabIndex = 0;
        Object.assign(row.style, {cursor: 'pointer', borderRadius: '3px'});
        const pick = () => {
          picked.value = a.animal;
          mark(a.animal);
        };
        row.addEventListener('click', pick);
        row.addEventListener('keydown', (ev) => {
          if (ev.key === 'Enter' || ev.key === ' ') {
            ev.preventDefault();
            pick();
          }
        });
        buttons.set(a.animal, row);
        rows.push(row);
      }
    }
    animalList.replaceChildren(...rows);
    // the strongest animal of the first treated group, until one is picked
    const keep = picked.peek();
    picked.value = keep !== undefined && buttons.has(keep) ? keep : buttons.keys().next().value;
    mark(picked.peek());
  });
  const pickedRow = computed(() => {
    const e = result.value?.effects;
    const a = picked.value;
    const i = e && a ? e.col('animal')!.toList().indexOf(a) : -1;
    return i < 0 ? undefined : {animal: a!, group: e!.get('group', i) as string, effect: e!.get('effect', i),
      lower: e!.get('effect_lower', i), upper: e!.get('effect_upper', i)};
  });
  const twinTable = computed(() => {
    const r = result.value;
    const p = pickedRow.value;
    return r && p ? animalFrame(r.twins, p.animal, colorsOf(r)[p.group]) : undefined;
  });
  const twinAxis = volumeAxis();
  const animalTwin = new BandChart(twinTable, {x: 'day', lower: 'lower', upper: 'upper', line: 'median',
    split: 'series', xTitle: 'Day after randomization', yTitle: 'Tumor volume', log: twinAxis.log, markers: true});
  const effectLine = span('');
  form.effect(() => {
    const p = pickedRow.value;
    effectLine.textContent = !p ? '' : `${p.animal}, day ${cell('evaluation_day').value}: ` +
      effectWords(p.effect, p.lower, p.upper) + '.';
  });
  const twinHead = divV([effectLine, twinAxis.row]);
  Object.assign(twinHead.style, {gap: 'var(--dg-space-s)', margin: '0 0 var(--dg-space-s)', fontWeight: '500'});
  const twinSide = divV([twinHead, sized(animalTwin, 320)], 'pmx-twin-side');
  Object.assign(twinSide.style, {position: 'sticky', top: '0', alignSelf: 'start'});
  const animals = div([divV([animalList, responderNote]), twinSide], 'pmx-animals');
  Object.assign(animals.style, {display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(460px, 1fr))',
    gap: 'var(--dg-space-l)', alignItems: 'start'});

  // the growth rate of each animal, the quantity the test compares
  const slopes = viewers.boxPlot(computed(() => {
    const r = result.value;
    if (!r)
      return undefined;
    const a = r.animals;
    const df = DG.DataFrame.fromObjects(Array.from({length: a.rowCount}, (_, i) => ({'Group': a.get('group', i),
      'Growth rate (ln V per day)': a.get('slope_per_day', i)})).filter((x) => x['Growth rate (ln V per day)'] !==
      null))!;
    df.col('Group')!.meta.colors.setCategorical(colorsOf(r));
    // the control first, then the treated groups in their order
    df.col('Group')!.setCategoryOrder([r.analysis.get('control', 0), ...r.comparisons.col('arm')!.toList()]);
    return df;
  }), {valueColumnName: 'Growth rate (ln V per day)', categoryColumnNames: ['Group'], showPValue: false,
    showStatistics: false});

  // the figures answer the study's question; each tab, one closer question
  const tabs = resultTabs('analyze', [
    {id: 'efficacy', label: 'Efficacy', intro: 'Whether each treatment slowed tumor growth against the control, and ' +
      'by how much.', content: [
      block('Experimental groups', 'Each treated group against the control on the evaluation day', groupsTable,
        groupNotes),
      block('Observed tumor growth', 'Median and range (smallest to largest tumor) of each group, to the evaluation ' +
        'day', growthAxis.row, sized(growth, 380))]},
    {id: 'animals', label: 'Each animal', intro: 'How each treated tumor compares with what it would have measured ' +
      'untreated, and which animals responded.', content: [
      block('Each animal against its untreated twin', 'Select an animal to see its growth', animals)]},
    {id: 'tolerability', label: 'Tolerability', intro: 'Whether the treatments were well tolerated.',
      content: tolerability.content},
    {id: 'plan', label: 'Plan vs observed', intro: 'Whether the study went as it was planned.',
      content: versus.content},
    {id: 'statistics', label: 'Statistics', intro: 'The numbers behind each result, for a closer look or a report.',
      content: [
        block('Results summary', 'The study’s results in words, ready to paste',
          copyRow(() => summaryText, 'Results summary'), summary),
        block('Growth rates', 'The per-animal slope of log volume each test compares', statsTable, sized(slopes, 260)),
        block('Control check', 'Each control against the range the other controls predict', controlCheck),
        block('Analysis protocol', 'As Plan pre-specifies it, ready to paste',
          copyRow(() => protocolText, 'Analysis protocol'), protocol)]},
  ]);
  const content = divV([cards, tabs.root], 'pmx-analyze-content');
  Object.assign(content.style, {gap: 'var(--dg-space-l)', minWidth: `${RESULTS_MIN}px`});
  form.effect(() => {
    content.hidden = !result.value;
    start.hidden = !!result.value;
  });

  const overlay = busyOverlay('Analyzing the study…');
  const failure = div([], 'pmx-analyze-error');
  Object.assign(failure.style, {color: 'var(--dg-error-color, #c62828)'});
  form.effect(() => overlay.hidden = !busy.value);
  form.effect(() => {
    failure.hidden = error.value === null;
    failure.replaceChildren(...(error.value === null ? [] : [noteLine(error.value, 'error')]));
  });
  const outdated = staleNote('The inputs changed: run to update the analysis');
  form.effect(() => {
    outdated.hidden = !stale.value;
    content.style.opacity = stale.value ? '0.5' : '';
  });

  const panes = new Splitter([declarationPane(form.root),
    resultsFrame([failure, outdated, start, content], overlay)],
  {direction: 'horizontal', sizes: [0.3, 0.7], minSize: DECLARATION_MIN});
  const view: DG.ViewBase = appView({
    name: 'PMX Studio · Analyze',
    content: panes,
    ribbon: [[stageStrip('analyze', (s) => void openStage(s))]],
    path: '/analyze',
    own: [form, correction, planField, growthAxis.input, twinAxis.input, growth, animalTwin, ...versus.own,
      ...tolerability.own, tabs],
    status: computed(() => busy.value ? 'Analyzing the study…' : result.value ? result.value.analysis.name : ''),
  });
  return view;
}
