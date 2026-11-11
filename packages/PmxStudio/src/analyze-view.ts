/* The Analyze view, built as the Plan view is: the declaration form generated from pmxAnalyzeStudy's own signature
   (u2 funcForm: the study table, its columns and the analysis settings), the results from the function's outputs, in
   the order a scientist asks: did each treatment work (the figures, then the groups with their result), how did the
   tumors grow (band-chart.ts), which animals responded (each treated animal against its untreated twin, the
   animals listed by effect, the one selected drawn against its twin), then, folded, what serves a closer look: the
   statistics behind each result, the control check and the analysis protocol; the study against its plan has a tab
   of its own (plan-vs-observed.ts). */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {badge, BasicTable, computed, div, divH, divV, dot, loader, signal, span, tableFromMap,
  Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal, StatValue} from '@datagrok-libraries/u2';
import {funcForm, ViewHost, viewers} from '@datagrok-libraries/u2/dg';
import type {ViewContent} from '@datagrok-libraries/u2/dg';

import {AnalysisFrames, AnalysisInputs, PREDICTED, analysisDeclarationOf, analysisFrames, warmTwin}
  from './analyze';
import {median} from './engine';
import {groupColors} from './frames';
import {CorrectionInput} from './arms-input';
import {askedOf, autoRun, axisAt, callInputs, cardGrid, dayText, keepSame, onTableEdits, only, TABLE_EDITS, guide,
  exampleLine, infoNote, noteLine, notesBox, num, pct, pValue, settleForm, signed, isBodyWeight, keepNamedColumns,
  sized, markStage,
  statCard, swatch, volumeAxis, wrapCells, stagePanes} from './ui';
import {logOf, restoreFields} from './study-file';
import {openStage} from './stages';
import {noEffectNote, observedCell, planInput, planTab} from './plan-vs-observed';
import {nciReadings, NciReading, tolerabilityTab} from './tolerability';
import {BandChart} from './band-chart';
import {explore, exploreContent, StudyInputs, studyInputsOf} from './explore';
import {editLayout, panelTabs, Tab} from './panels';
import {stagePage, studyHeader} from './study-header';
import {AnimalSelection, exclusionsOf, mirrorSelection, studyDataLine, studyMenu, studyOf, tagStudy,
  withoutExcluded} from './study-table';
import {metaOf, setMetaOn} from './study-meta';
import {Exports, download, exportMenu, reportHtml} from './report';
import {printSheet} from './sheet';
import {SpiderChart, WaterfallChart, spiderLines} from './spider-chart';
import {LineInput} from './line-input';
import {Competing, Survival, competingOf, survivalOf} from './survival';
import {KmChart} from './km-chart';
import {synergyDialog} from './synergy';
import {arriveDialog, arriveItems} from './arrive';

/** Analyze's tabs and their panels, in order and width (twelfths): the answer and how the groups grew; each animal;
 *  how the animals bore it; whether the study went as planned. */
const ANALYZE_TABS: Tab[] = [
  {name: 'Results', layout: [{id: 'groups', span: 12}, {id: 'growth', span: 8}, {id: 'rates', span: 4},
    {id: 'km', span: 8}, {id: 'explore', span: 4}]},
  {name: 'Each animal', layout: [{id: 'waterfall', span: 6}, {id: 'spider', span: 6}, {id: 'animals', span: 12}]},
  {name: 'Tolerability', layout: [{id: 'tolerability', span: 12}]},
  {name: 'Against the plan', layout: [{id: 'plan', span: 12}, {id: 'controls', span: 12}]}];

/** Under a growth chart: how many animals each group had on each measurement day, so a thinning group shows. */
function animalsPerDay(volumes: DG.DataFrame, groups: string[], colors: {[g: string]: string}): HTMLElement {
  const days = [...new Set(volumes.col('day')!.toList() as number[])].sort((a, b) => a - b);
  const n = new Map<string, number>();
  for (let i = 0; i < volumes.rowCount; i++)
    n.set(`${volumes.get('group', i)}|${volumes.get('day', i)}`, volumes.get('n', i));
  const table = document.createElement('table');
  table.className = 'pmx-n-table';
  const head = table.createTHead().insertRow();
  const th = (text: string) => Object.assign(document.createElement('th'), {textContent: text});
  head.append(th('Animals by day'), ...days.map((d) => th(String(d))));
  const body = table.createTBody();
  for (const g of groups) {
    const row = body.insertRow();
    row.insertCell().append(swatch(g, colors[g]));
    for (const d of days)
      row.insertCell().textContent = String(n.get(`${g}|${d}`) ?? '');
  }
  return div([table], 'pmx-n-scroll');
}

const TWIN = 'Untreated twin';

/** The inputs the twin reads (the study's volumes to the evaluation day, the line's model, the control): a change of
 *  any other (the correction, α, the weights, the plan) finds the twin in the engine's cache. */
const TWIN_INPUTS = ['study', 'animal', 'group', 'day', 'volume', 'line', 'control', 'endpoint', 'evaluationDay'];
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

const correctionWords = (corr: unknown, k: number) => k < 2 ? 'one comparison' : corr === 'bonferroni' ?
  'Bonferroni\'s correction' : corr === 'none' ? 'no correction' : 'Holm\'s correction';

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
  const twinned = r.analysis.get('twinned', 0);
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
      `against the control (p ${p.startsWith('<') ? p : `= ${p}`})` + (twinned ? `; ${g('responders')} of ` +
      `${g('twin_animals')} animals responded against their untreated twin.` : '.'));
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

/** What Analyze opens with: a study table (with its plan in the `pmx.plan` tag, if any). */
export interface AnalyzeStart {study?: DG.DataFrame}

/** The Analyze view; given a study, it analyzes it at once. */
export async function analyzeView(opening?: AnalyzeStart): Promise<DG.ViewBase> {
  const exports: Exports = [];
  const host = ViewHost.create({name: 'PMX Studio', content: (self) => analyzeContent(self, opening,
    exports)});
  exportMenu(host.view, exports);
  return markStage(host.view, 'analyze');
}

/** What the view shows, built under its host, which owns it all. */
function analyzeContent(self: ViewHost, opening: AnalyzeStart | undefined, exports: Exports): ViewContent {
  // the twin's shaders compiled in Crux's worker while the study is picked, once a page (a failure only leaves the
  // first twin slower); here and not when the app opens, which starts no engine
  void warmTwin().catch(() => null);
  const call = DG.Func.byName('PmxStudio:pmxAnalyzeStudy').prepare();
  // the engine itself on the form's inputs, not the registered function's FuncCall (see callInputs)
  const run = (extra: object) =>
    analysisFrames(analysisDeclarationOf(withoutExcluded(callInputs<AnalysisInputs>(call, extra))));
  const result = signal<AnalysisFrames | undefined>(undefined);
  // the newest run: the form's changes and an opened study share the call, and a run that finishes after a newer one
  // started is dropped
  let latest = 0;
  const error = signal<string | null>(null);
  // a result that no longer matches the form is dimmed until the next run
  const stale = signal(false);
  // the multiple comparisons as Plan asks for them, over the treated groups of the last analysis (two until then)
  const treatedGroups = signal(2);
  const correction = new CorrectionInput({label: 'Multiple comparisons', arms: treatedGroups});
  // no Run: each change analyzes again once the changes stop
  let preparing = false;
  const auto = autoRun(() => form, () => preparing ? undefined : analyze());
  const overrides = {correction: {input: correction}, line: {input: new LineInput({label: 'Tumor model'})}};
  const form = funcForm(call, {layout: 'tall', overrides, onInputChanged: () => {
    stale.value = result.peek() !== undefined;
    auto.changed();
  }});
  settleForm(form, {control: 'Auto (the vehicle group)', evaluationDay: 'Auto (endpoint reached)',
    density: 'Auto (the cell line’s in the atlas)'});
  keepNamedColumns(form, {
    bodyWeight: isBodyWeight,
    tumorWeight: (c) => /tum/i.test(c) && /weight|mass/i.test(c),
    removalReason: (c) => /reason|removal|fate|cause|status/i.test(c),
  });
  // busy while the study is analyzed, and before that while an opened study is prepared
  const busy = signal(false);
  const planState = {stale, error, result, reading: Promise.resolve() as Promise<unknown>,
    changed: () => auto.changed()};
  // the study and columns the shown result read (Your data opens them)
  let resultStudy: StudyInputs | null = null;
  // the same without its excluded rows: what the analysis read, and so what KM, Fine-Gray, synergy and Your data read
  let keptStudy: StudyInputs | null = null;
  const shown = (reads: StudyInputs | null) => {
    resultStudy = reads;
    keptStudy = reads && withoutExcluded(reads);
    selection.attach(reads);
  };
  // what the last analysis read of the study (its table as text, the columns' roles, the removal column): the panels
  // that draw the study's own rows (each animal's line, the days to the endpoint) follow it, since a value edited may
  // leave every table of the result as it was (the medians by day) and still move them
  const dataRead = signal('');
  // the animals selected, through the study table's own selection: what one chart picks, every chart marks
  const selection = new AnimalSelection();
  self.own(() => selection.dispose());
  // the twin inputs of the last analyses that got their twin, as many as the engine keeps (twinOf's 8)
  const twinsCached: string[] = [];
  const analyze = async (prepare?: () => Promise<void>) => {
    busy.value = true;
    const mine = ++latest;
    try {
      // a study the view fills in changes the form as it goes: none of those changes runs by itself
      preparing = prepare !== undefined;
      await prepare?.();
      preparing = false;
      auto.cancel();
      await planState.reading;
      const sent = askedOf(call);
      const reads = studyInputsOf(call);
      tagStudy(reads);
      const twinSent = askedOf(call, TWIN_INPUTS);
      const removal = call.getParamValue('removalReason');
      const data = `${twinSent}\u0000${removal instanceof DG.Column ? removal.name : removal ?? ''}`;
      // the twin unchanged (cached): one call, with it
      if (twinsCached.includes(twinSent)) {
        const frames = await run({withTwin: true});
        if (mine !== latest)
          return;
        error.value = null;
        stale.value = askedOf(call) !== sent;
        shown(reads);
        result.value = keepSame(result.peek(), frames);
        dataRead.value = data;
        return;
      }
      // else the analysis first; then the same inputs with the twins (seconds), whatever the form holds by then
      const withTwin = withoutExcluded(callInputs<AnalysisInputs>(call, {withTwin: true}));
      const frames = await run({withTwin: false});
      if (mine !== latest)
        return;
      error.value = null;
      stale.value = askedOf(call) !== sent;
      shown(reads);
      result.value = keepSame(result.peek(), frames);
      dataRead.value = data;
      void analysisFrames(analysisDeclarationOf(withTwin)).then((done) => {
        if (mine !== latest)
          return;
        if (!twinsCached.includes(twinSent) && twinsCached.push(twinSent) > 8)
          twinsCached.shift();
        result.value = keepSame(result.peek(), done);
      }, (e) => {
        if (mine === latest)
          error.value = String(e?.message ?? e);
      });
    } catch (e: any) {
      // only the last run speaks: an older one failing late would write over a newer result's place
      if (mine === latest) {
        error.value = String(e?.message ?? e);
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

  // a study brought from elsewhere (the hub, Monitor's "Analyze this study"): opened in the workspace like any table
  // (the form's table picker then offers it), this view brought back, picked in the form and analyzed at once; its
  // columns as the form picks them by name
  const openStudy = (getStudy: () => Promise<DG.DataFrame>) => analyze(async () => {
    const df = await getStudy();
    if (!grok.shell.tables.includes(df))
      grok.shell.addTableView(df);
    // a view not yet shown is left to showStage: making it current would add it, and showStage add it a second time
    if ([...grok.shell.views].some((v) => v.dart === self.view.dart))
      grok.shell.v = self.view;
    const study = form.getInput('study');
    if (study)
      study.value = df.name;
    await form.settled;
    // a table from Monitor brings its body weight column and limits
    restoreFields(form, df);
    await form.settled;
  });
  void form.settled.then(() => {
    if (opening?.study)
      void openStudy(async () => opening.study!);
  });
  planInput(form, call, planState);

  const start = divV([
    span('Open your study table (one row per animal and measurement day, day 0 the randomization), check its ' +
      'columns: the analysis follows them, and any change you make to the table.'),
    exampleLine(),
    infoNote('<b>Analyzed as planned.</b> Analyze applies the analysis Plan pre-specifies: the study is evaluated ' +
      'when the median control reaches the endpoint, each treated group is compared with the control on tumor ' +
      'growth, and each treated animal with its untreated twin: what it would measure without treatment, predicted ' +
      'from its own tumor at randomization, this study\'s controls and the cell line\'s growth.'),
  ], 'pmx-analyze-start');
  Object.assign(start.style, {gap: 'var(--dg-space-m)', maxWidth: '760px'});

  const analysis = computed(() => result.value?.analysis);
  const cell = <T = number>(name: string) =>
    computed(() => (analysis.value?.get(name, 0) ?? undefined) as T | undefined);
  const comparisons = computed(() => result.value?.comparisons);
  // what depends on the twins waits for them: a loader, or a dash if they failed
  const twinned = computed(() => !!analysis.value?.get('twinned', 0));
  const waitingTwin = (words?: string) => error.peek() !== null ? span('—') : loader(words);
  const coverage = cell('control_coverage');
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
    card('Controls on track', computed(() => twinned.value ? coverage.value : undefined), () => 'Control ' +
      'measurements inside the range the other controls predict for them (about 90% is as expected)', (x) => pct(x)),
  ], 'pmx-analyze-cards', 150);

  const groupsTable = div([]);
  const groupNotes = notesBox();
  const statsTable = div([]);
  const controlCheck = div([]);
  const protocol = div([]);
  let protocolText = '';
  let statsCsv = '';
  const summary = divV([], 'pmx-results-summary');
  summary.style.gap = 'var(--dg-space-s)';
  let summaryText = '';
  // the Tolerability panel, at the groups against the NCI limits (the panels are built below)
  const showNci = () => {
    panels.show('tolerability');
    requestAnimationFrame(() => tolerability.content.map((el) => el.querySelector('.pmx-nci')).find((el) => el)
      ?.scrollIntoView({behavior: 'smooth', block: 'center'}));
  };
  // what a panel reads besides its own tables: the control and the treated groups (their colours follow from them)
  const groupsOf = (r: AnalysisFrames) => [r.analysis.get('control', 0),
    r.comparisons.col('arm')!.toList().join('\u0000')];
  const colorsOf = (r: AnalysisFrames) => {
    return groupColors(r.comparisons.col('arm')!.toList(), r.analysis.get('control', 0));
  };
  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    const v = (c: string) => r.analysis.get(c, 0);
    const done = v('twinned');
    void error.value;
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
        !done ? waitingTwin() : g('paired_effect') === null ? '—' : pct(g('paired_effect') / 100),
        !done ? waitingTwin() : `${g('responders')} of ${g('twin_animals')}`, pValue(g('p_value')).replace(' ', ' '),
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
    statsCsv = DG.DataFrame.fromObjects(stats.map((s) => Object.fromEntries(statHeaders.map((h, j) =>
      [h, s.cells[j]]))))!.toCsv();

    // the controls against each other: how many measurements fall in their range, and the ones that grow apart
    const inBand = v('control_coverage');
    controlCheck.replaceChildren(!done ? waitingTwin('Checking the controls against each other…') : divV([
      span(inBand === null ? 'Too few controls to check them against each other.' :
        `${pct(inBand)} of the control measurements fall in the range the other controls predict for them ` +
        `(about 90% when the controls grow as the ${v('line')} model expects).`),
      ...(odd.length ? odd.map((n) => noteLine(n)) : [noteLine('No control grows apart from the others.', 'success')]),
    ], 'pmx-control-check'));
    if (done)
      (controlCheck.firstElementChild as HTMLElement).style.gap = 'var(--dg-space-s)';

    // the analysis protocol, as Plan pre-specifies it
    const k = v('treated_arms');
    const corr = v('correction');
    const rowsOf: {[row: string]: string} = {
      'Primary comparison': 'Each treated group against the control on tumor growth: the per-animal slope of log ' +
        'volume from randomization to the evaluation day, by Welch’s t test (each group its own variance)',
      'Multiple comparisons': k < 2 ? 'One treated group: none needed' : corr === 'none' ?
        `None: each of the ${k} comparisons at α = ${num(v('alpha'), 2)} (exploratory)` :
        `${corr === 'bonferroni' ? 'Bonferroni' : 'Holm'} over ${k} comparisons`,
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
    const excluded = exclusionsOf(studyInputsOf(call));
    if (excluded)
      rowsOf['Exclusions'] = excluded;
    protocol.replaceChildren(tableFromMap(rowsOf));
    protocolText = [`Analysis protocol (${r.analysis.name})`, ...Object.entries(rowsOf).map(([key, t]) =>
      `${key}: ${t}`)].join('\n');

    const lines = summaryLines(r, nci);
    summary.replaceChildren(...lines.map((l) => span(l)));
    summaryText = lines.join('\n');
  });

  const versus = planTab(form, result, colorsOf);
  versus.own.forEach((o) => self.own(() => o.dispose()));
  const tolerability = tolerabilityTab(form, only(result, ['tolerability', 'weights', 'volumes'], groupsOf),
    colorsOf);

  const growthAxis = volumeAxis();
  // the plan over what was observed: the median control a study run as planned shows (the plan's, else the
  // line's atlas), a grey band under the groups
  const showPredicted = signal(true);
  const predictedInput = ui.input.bool('Predicted controls', {value: true,
    tooltipText: 'The range of a study\'s median control run as planned (95 of 100 labs): the plan\'s, or the ' +
      'line\'s atlas without one', onValueChanged: (v) => showPredicted.value = v});
  growthAxis.row.prepend(predictedInput.root);
  const growthData = computed(() => {
    const r = result.value;
    if (!r)
      return undefined;
    const p = r.predicted;
    const ahead = Array.from({length: p.rowCount}, (_, i) => i).filter((i) => p.get('series', i) === PREDICTED);
    if (!showPredicted.value || ahead.length === 0)
      return r.volumes;
    const v = r.volumes;
    const row = (df: DG.DataFrame, i: number, group: string) => ({group, day: df.get('day', i),
      median: df.get('median', i), lower: df.get('lower', i), upper: df.get('upper', i)});
    const df = DG.DataFrame.fromObjects([...ahead.map((i) => row(p, i, PREDICTED)),
      ...Array.from({length: v.rowCount}, (_, i) => row(v, i, v.get('group', i)))])!;
    for (const c of ['day', 'median', 'lower', 'upper']) {
      df.col(c)!.meta.units = v.col(c)?.meta.units ?? '';
      df.col(c)!.meta.format = v.col(c)?.meta.format ?? '';
    }
    df.col('group')!.meta.colors.setCategorical({[PREDICTED]: '#9aa5b1', ...colorsOf(r)});
    df.meta.formulaLines.addAll(v.meta.formulaLines.items);
    return df;
  });
  const growth = new BandChart(growthData, {x: 'day', lower: 'lower', upper: 'upper',
    line: 'median', split: 'group', xTitle: 'Day after randomization', yTitle: 'Tumor volume', log: growthAxis.log,
    markers: true});

  // each treated animal against its twin: the animals listed by group, the strongest first, the one selected drawn
  type AnimalRow = {animal: string, group: string, effect: number, lower: number, upper: number};
  const picked = signal<string | undefined>(undefined);
  // a control animal selected elsewhere (the waterfall, a line): it has no twin, which the twin's side says
  const untwinned = signal<string | undefined>(undefined);
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
    legendLine(form.build(() => explained(badge('Responder', {variant: 'success'}), responderText)),
      responderWords),
    legendLine(dots, span(`Each bar: how much smaller the treated tumor is than its twin's (the dot, filled for a ` +
      `responder) and its 90% range; the line at "Same as twin" is the tumor's untreated size.`)),
  ], 'pmx-responder-legend');
  Object.assign(responderNote.style, {gap: 'var(--dg-space-s)', margin: 'var(--dg-space-m) 0 0', fontSize: '12px',
    color: 'var(--dg-text-color-muted, #68727f)', lineHeight: '1.45'});
  // the rows' columns: the animal, its effect in words, the bar, its badge; one grid the rows share (subgrid), so
  // every row's bar sits under the axis and the badge column fits its widest badge, beside the chart too
  const twinSide = divV([], 'pmx-twin-side');
  Object.assign(twinSide.style, {position: 'sticky', top: '0', alignSelf: 'start'});
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
    // (a pick stays until the twins are in)
    const done = r.analysis.get('twinned', 0);
    void error.value;
    responderNote.hidden = !done;
    twinSide.hidden = !done;
    if (!done) {
      const waiting = waitingTwin('Comparing each animal with its untreated twin…');
      waiting.style.gridColumn = '1 / -1';
      animalList.replaceChildren(waiting);
      return;
    }
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
        row.dataset.animal = a.animal;
        // a plain click selects it (never clears it: picked from the waterfall, it stays); Ctrl, Cmd or Shift toggle
        const choose = (ev: MouseEvent | KeyboardEvent) => {
          pick();
          if (ev.ctrlKey || ev.metaKey || ev.shiftKey)
            selection.pick(a.animal, ev);
          else
            selection.select(new Set([a.animal]));
        };
        row.addEventListener('click', choose);
        row.addEventListener('contextmenu', (ev) => {
          ev.preventDefault();
          studyMenu(studyOf(call.getParamValue('study') as DG.DataFrame | null), 'animal', a.animal,
            {causedBy: ev}, [['Show against its twin', pick]]);
        });
        row.addEventListener('keydown', (ev) => {
          if (ev.key === 'Enter' || ev.key === ' ') {
            ev.preventDefault();
            choose(ev);
          }
        });
        buttons.set(a.animal, row);
        rows.push(row);
      }
    }
    animalList.replaceChildren(...rows);
    markSelected(selection.animals.peek());
    // the strongest animal of the first treated group, until one is picked
    const keep = picked.peek();
    picked.value = keep !== undefined && buttons.has(keep) ? keep : buttons.keys().next().value;
    mark(picked.peek());
  });
  // the animals selected anywhere, outlined in the list
  const markSelected = (on: ReadonlySet<string>) => animalList.querySelectorAll<HTMLElement>('.pmx-forest-row')
    .forEach((row) => row.classList.toggle('is-selected', on.has(row.dataset.animal!)));
  // one animal selected anywhere (a bar, a line): drawn against its twin too
  form.effect(() => {
    const on = selection.animals.value;
    markSelected(on);
    const [one] = on;
    const listed = on.size === 1 && !!animalList.querySelector(`[data-animal="${CSS.escape(one)}"]`);
    untwinned.value = on.size === 1 && !listed && animalList.querySelector('.pmx-forest-row') ? one : undefined;
    if (untwinned.peek() !== undefined)
      picked.value = undefined;
    if ((listed && one !== picked.peek()) || untwinned.peek() !== undefined) {
      if (listed)
        picked.value = one;
      animalList.querySelectorAll<HTMLElement>('.pmx-forest-row').forEach((b) => {
        const at = listed && b.dataset.animal === one;
        b.setAttribute('aria-pressed', String(at));
        b.style.background = at ? 'color-mix(in srgb, var(--dg-accent, #2083d5) 12%, transparent)' : '';
      });
    }
  });
  const twinsOnly = only(result, ['twins', 'effects'], groupsOf);
  const pickedRow = computed(() => {
    const e = twinsOnly.value?.effects;
    const a = picked.value;
    const i = e && a ? e.col('animal')?.toList().indexOf(a) ?? -1 : -1;
    return i < 0 ? undefined : {animal: a!, group: e!.get('group', i) as string, effect: e!.get('effect', i),
      lower: e!.get('effect_lower', i), upper: e!.get('effect_upper', i)};
  });
  const twinTable = computed(() => {
    const r = twinsOnly.value;
    const p = pickedRow.value;
    return r && p ? animalFrame(r.twins, p.animal, colorsOf(r)[p.group]) : undefined;
  });
  const twinAxis = volumeAxis();
  const animalTwin = new BandChart(twinTable, {x: 'day', lower: 'lower', upper: 'upper', line: 'median',
    split: 'series', xTitle: 'Day after randomization', yTitle: 'Tumor volume', log: twinAxis.log, markers: true});
  const sizedTwin = sized(animalTwin, 320);
  const effectLine = span('');
  form.effect(() => {
    const p = pickedRow.value;
    const alone = untwinned.value;
    effectLine.textContent = p ? `${p.animal}, day ${cell('evaluation_day').value}: ` +
      effectWords(p.effect, p.lower, p.upper) + '.' : alone ? `${alone} has no twin: only a treated animal is ` +
      'compared with an untreated twin (a control is the untreated growth itself). Pick a treated one, in the list ' +
      'or the waterfall.' : '';
    sizedTwin.style.display = twinAxis.row.style.display = alone && !p ? 'none' : '';
  });
  const twinHead = divV([effectLine, twinAxis.row]);
  Object.assign(twinHead.style, {gap: 'var(--dg-space-s)', margin: '0 0 var(--dg-space-s)', fontWeight: '500'});
  twinSide.append(twinHead, sizedTwin);
  const animals = div([divV([animalList, responderNote]), twinSide], 'pmx-animals');
  Object.assign(animals.style, {display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(460px, 1fr))',
    gap: 'var(--dg-space-l)', alignItems: 'start'});

  // the growth rate of each animal, the quantity the test compares
  // under the growth chart, the animals behind each point
  const perDay = div([], 'pmx-n-per-day');
  const volumesOnly = only(result, ['volumes'], groupsOf);
  form.effect(() => {
    const r = volumesOnly.value;
    perDay.replaceChildren(...(r ? [animalsPerDay(r.volumes, [r.analysis.get('control', 0),
      ...r.comparisons.col('arm')!.toList()], colorsOf(r))] : []));
  });

  // each animal as a bar in its group's colour, on what the scientist compares them by: its growth rate (the slope of
  // log volume the test compares), its best response (mRECIST's, % from randomization) or how much smaller it is than
  // its untreated twin; sorted so the better animals end on the right
  const WATERFALL = {
    'Growth rate': {col: 'slope_per_day', title: 'Growth rate (log volume per day)', order: 'desc', control: true,
      none: '<b>No growth rates.</b> Each animal needs two measurements or more.'},
    'Best response': {col: 'best_response', title: 'Best response (%)', order: 'desc', control: true,
      none: '<b>No best response.</b> Each animal needs measurements from day 10 on for its mRECIST class.'},
    'Against its twin': {col: 'effect', title: 'Smaller than its twin (%)', order: 'asc', control: false,
      none: '<b>No twins.</b> The twin needs the cell line’s growth and each animal’s volume at ' +
        'randomization.'},
  } as const;
  type WaterfallBy = keyof typeof WATERFALL;
  const waterfallBy = signal<WaterfallBy>('Growth rate');
  const byInput = ui.input.choice('Compare by', {value: 'Growth rate', items: Object.keys(WATERFALL),
    onValueChanged: (v) => waterfallBy.value = (v ?? 'Growth rate') as WaterfallBy,
    tooltipText: 'Growth rate: the slope of log volume each test compares. Best response: mRECIST’s smallest ' +
      'volume from day 10 against randomization. Against its twin: how much smaller than its untreated twin.'});
  const animalsOnly = only(result, ['animals', 'effects'], groupsOf);
  const waterfallData = computed(() => {
    const r = animalsOnly.value;
    const by = WATERFALL[waterfallBy.value];
    const t = by.col === 'effect' ? r?.effects : r?.animals;
    if (!r || !t?.col(by.col))
      return undefined;
    const control = r.analysis.get('control', 0);
    const bars = Array.from({length: t.rowCount}, (_, i) => ({animal: String(t.get('animal', i)),
      group: String(t.get('group', i)), value: t.get(by.col, i) as number}))
      // a missing value may come back as null, undefined or NaN: only numbers make a bar
      .filter((x) => (by.control || x.group !== control) && Number.isFinite(x.value));
    const groups = [control, ...r.comparisons.col('arm')!.toList()].filter((g) => bars.some((b) => b.group === g));
    return bars.length ? {bars, colors: colorsOf(r), title: by.title, order: by.order, groups} : undefined;
  });
  const waterfall = new WaterfallChart(waterfallData, (animal, e) =>
    studyMenu(studyOf(call.getParamValue('study') as DG.DataFrame | null), 'animal', animal, {causedBy: e}),
  selection);
  const waterfallMissing = infoNote('');
  form.effect(() => {
    const has = !!waterfallData.value;
    waterfall.root.hidden = !has;
    waterfallMissing.hidden = has || !result.value;
    waterfallMissing.querySelector('span')!.innerHTML = WATERFALL[waterfallBy.value].none;
  });

  // each animal's growth as a % change from randomization, in its group's colour (the study's own rows, as analyzed)
  const spider = new SpiderChart(computed(() => {
    const r = volumesOnly.value;
    void dataRead.value;
    const s = resultStudy;
    return r && s ? {lines: spiderLines(s), colors: colorsOf(r)} : undefined;
  }), selection);

  // the strict Kaplan-Meier: days to the endpoint as measured, an animal that never reached it censored
  // projected, by its slider: at the last day measured, what was measured; further right, each animal still under
  // the endpoint carried on at the rate of its last three measurements to that day
  const survival = signal<Survival | undefined>(undefined);
  const kmOnly = only(result, ['volumes'], (r) => [...groupsOf(r), r.analysis.get('endpoint', 0), dataRead.value]);
  const horizon = signal<number | null>(null);
  // one computation at a time, of the latest day asked for: a slider dragged fast asks for many, and answers
  // arriving out of order made the curves jump
  let computing = false;
  const compute = async () => {
    if (computing)
      return;
    computing = true;
    try {
      for (;;) {
        const r = kmOnly.peek();
        const s = keptStudy;
        const h = horizon.peek();
        if (!r || !s) {
          survival.value = undefined;
          return;
        }
        const k = await survivalOf(s, r.analysis.get('control', 0), r.analysis.get('endpoint', 0), colorsOf(r),
          h === null ? undefined : {horizon: h}).catch(() => undefined);
        if (kmOnly.peek() === r && horizon.peek() === h) {
          survival.value = k;
          return;
        }
      }
    } finally {
      computing = false;
    }
  };
  form.effect(() => {
    void kmOnly.value;
    void horizon.value;
    void compute();
  });
  // the animals removed for another reason: competing with the endpoint (Fine-Gray), when the study has any
  const competing = signal<Competing[] | null>(null);
  form.effect(() => {
    const r = kmOnly.value;
    const s = keptStudy;
    const removal = call.getParamValue('removalReason');
    if (!r || !s) {
      competing.value = null;
      return;
    }
    void competingOf(s, removal instanceof DG.Column ? removal.name : null, r.analysis.get('control', 0),
      r.analysis.get('endpoint', 0)).then((c) => {
      if (kmOnly.peek() === r)
        competing.value = c;
    }, () => competing.value = null);
  });
  const km = new KmChart(survival, horizon);
  const kmWords = div([], 'pmx-km-words');
  form.effect(() => {
    const k = survival.value;
    const p = new Map((k?.logrank ?? []).map((t) => [t.group, t.p_value]));
    // projected: no test (the tests are the measured study's), every figure said to be projected
    const on = k?.horizon !== undefined;
    kmWords.replaceChildren(...(on ? [noteLine(`Projected beyond day ${k!.lastDay} to day ${k!.horizon}: each ` +
      'animal still under the endpoint carried on at the rate of its last three measurements. Not a result: the ' +
      'tests and the report\'s conclusions are the measured study\'s.', 'warning')] : []),
    ...(k?.medians ?? []).map((m) => noteLine(`${m.group}: ${on ? 'projected ' : ''}median ${m.median === null ?
      'not reached' : dayText(m.median)}, ${m.events} of ${m.n} ${on ? 'would reach' : 'reached'} the endpoint` +
      (!on && p.has(m.group) ? `; log-rank against the control ${pValue(p.get(m.group))}` : ''))),
    ...(!on && competing.value?.length ? [noteLine('With the animals removed for another reason as a competing ' +
      'event (Fine-Gray): ' + competing.value.map((c) => `${c.group}, subdistribution hazard ratio ` +
      `${c.shr.toFixed(2)} (${c.lower.toFixed(2)}–${c.upper.toFixed(2)}), p ${pValue(c.p)}, ${c.competing} removed`)
      .join('; '), 'success')] : []));
  });

  const slopesData = computed(() => {
    const r = animalsOnly.value;
    if (!r)
      return undefined;
    const a = r.animals;
    const df = DG.DataFrame.fromObjects(Array.from({length: a.rowCount}, (_, i) => ({'Group': a.get('group', i),
      'Animal': String(a.get('animal', i)), 'Growth rate (ln V per day)': a.get('slope_per_day', i)}))
      .filter((x) => x['Growth rate (ln V per day)'] !== null))!;
    df.col('Group')!.meta.colors.setCategorical(colorsOf(r));
    // the control first, then the treated groups in their order
    df.col('Group')!.setCategoryOrder([r.analysis.get('control', 0), ...r.comparisons.col('arm')!.toList()]);
    return df;
  });
  // its points selected with the study's animals, and a lasso on it selects them in the study
  let unmirror = () => {};
  form.effect(() => {
    const df = slopesData.value;
    void selection.animals.value;
    unmirror();
    unmirror = df ? mirrorSelection(df, 'Animal', selection) : () => {};
  });
  self.own(() => unmirror());
  const slopes = viewers.boxPlot(slopesData, {valueColumnName: 'Growth rate (ln V per day)',
    categoryColumnNames: ['Group'], showPValue: false,
    showStatistics: false});

  // the study's data in Datagrok's own table view; here, one line on it and the way there
  const dataLine = studyDataLine(() => studyInputsOf(call));
  self.own(dataLine.dispose);
  form.effect(() => {
    if (result.value)
      dataLine.refresh();
  });

  // the figures answer the study's question; each panel, one closer question
  const panels = panelTabs([
    {id: 'groups', title: 'Experimental groups', hint: 'Each treated group against the control on the evaluation day',
      content: [groupsTable, groupNotes]},
    {id: 'growth', title: 'Observed tumor growth', hint: 'Median and range (smallest to largest tumor) ' +
      'of each group, to the evaluation day; in grey, the median control predicted by the plan',
    content: [growthAxis.row, sized(growth, 380), perDay]},
    {id: 'rates', title: 'Growth rates', hint: 'The per-animal slope of log volume each test compares',
      content: [sized(slopes, 380)]},
    {id: 'waterfall', title: 'Each animal compared', hint: 'Every animal as a bar in its group\u2019s colour, on its ' +
      'growth rate, its best response or its effect against its twin; the better animals on the right',
    content: [byInput.root, waterfallMissing, sized(waterfall, 300)]},
    {id: 'spider', title: 'Each tumor over time', hint: 'Every animal’s volume as a % change from ' +
      'randomization, in its group’s colour; click a group in the legend to hide or show its animals',
    content: [sized(spider, 300)]},
    {id: 'km', title: 'Time to endpoint', hint: 'Kaplan-Meier: the animals of each group still below the ' +
      'endpoint volume, day by day, as measured; an animal that never reached it counts until its last measurement',
    content: [sized(km, 320), kmWords]},
    {id: 'animals', title: 'Each animal against its untreated twin', hint: 'How each treated tumor compares with ' +
      'what it would have measured untreated. Select an animal to see its growth', content: [animals]},
    {id: 'tolerability', title: 'Tolerability', hint: 'Whether the treatments were well tolerated',
      content: tolerability.content},
    {id: 'plan', title: 'Plan vs observed', hint: 'Whether the study went as it was planned',
      content: versus.content},
    {id: 'controls', title: 'Control check', hint: 'Each control against the range the other controls predict',
      content: [controlCheck]},
    {id: 'explore', title: 'Your data', hint: 'The study’s animals and measurements as Datagrok tables, to chart, ' +
      'filter, group and compare as you like', content: exploreContent((which) => {
      const r = result.peek();
      if (r && keptStudy)
        explore(which, r, keptStudy, colorsOf(r), self.view);
    })},
  ], ANALYZE_TABS, 'pmx-analyze-panels', 'analyze');
  exports.push(['Layout…', () => editLayout(panels.root), 'View']);
  const content = divV([dataLine.root, cards, panels.root], 'pmx-analyze-content');
  const copied = (text: string, what: string) => {
    studyNow();
    return navigator.clipboard.writeText(text).then(() => grok.shell.info(`${what} copied`));
  };
  const studyNow = () => {
    const df = call.getParamValue('study') as DG.DataFrame | null;
    const r = result.value;
    if (!df || !r)
      throw new Error('Analyze a study first');
    // an export is what is on screen: not while that is behind the table or the form
    if (busy.peek() || stale.peek())
      throw new Error('The analysis is not up to date yet: export once it is');
    return {df, r};
  };
  exports.push(
    ['Report (PDF)…', async () => {
      const {df, r} = studyNow();
      printSheet(await panels.all(() => reportHtml({title: r.analysis.name, stage: 'Analyze', meta: metaOf(df),
        content, log: logOf(df), annex: [['Results', summary], ['Statistics', statsTable],
          ['Analysis protocol', protocol]]})));
    }],
    ['Copy the results summary', () => copied(summaryText, 'Results summary')],
    ['ARRIVE checklist…', () => {
      const {df, r} = studyNow();
      const c = r.comparisons;
      const plan = r.analysis.get('plan', 0);
      arriveDialog(arriveItems({control: r.analysis.get('control', 0), arms: c.col('arm')!.toList() as string[],
        plan: plan ? String(plan) : null, exclusions: exclusionsOf(studyOf(df)), meta: metaOf(df),
        // from the log alone: a plan says how many animals, not that the groups were formed at random
        randomized: logOf(df).some((e) => /randomi[sz]|allocat/i.test(e.what)),
        endpoint: r.analysis.get('endpoint', 0), evaluationDay: r.analysis.get('evaluation_day', 0),
        protocol: protocolText, analyzed: [{group: r.analysis.get('control', 0), n: r.analysis.get('n_control', 0)},
          ...Array.from({length: c.rowCount}, (_, i) => ({group: c.get('arm', i), n: c.get('n', i)}))]}));
    }, 'Analyze'],
    ['Combination synergy…', () => {
      const {r} = studyNow();
      if (!keptStudy)
        throw new Error('Analyze a study first');
      synergyDialog(keptStudy, r.analysis.get('control', 0), r.comparisons.col('arm')!.toList() as string[],
        r.analysis.get('evaluation_day', 0));
    }, 'Analyze'],
    ['Copy the analysis protocol', () => copied(protocolText, 'Analysis protocol')],
    ['Statistics (CSV)', () => download(`${studyNow().df.name} - statistics.csv`, statsCsv, 'text/csv')],
    ['Measurements (CSV)', () => download(`${studyNow().df.name}.csv`, studyNow().df.toCsv(), 'text/csv')],
    ['Per-animal results (CSV)', () => download(`${studyNow().df.name} - animals.csv`, studyNow().r.animals.toCsv(),
      'text/csv')]);
  Object.assign(content.style, {gap: 'var(--dg-space-l)'});
  form.effect(() => {
    content.hidden = !result.value;
    start.hidden = !!result.value;
  });

  const panes = stagePanes(form, content, {busy, result, error, stale, before: [start],
    busyText: 'Analyzing the study…', staleText: 'the analysis shown is for the last complete inputs'});
  // an edit of the study table (a value corrected, rows added): analyzed again, after a pause
  self.own(onTableEdits(form, computed(() => [result.value, error.value]),
    () => call.getParamValue('study'), () => auto.changed(TABLE_EDITS)));
  const head = studyHeader(form, 'analyze', computed(() => {
    const r = result.value;
    if (!r)
      return undefined;
    const colors = colorsOf(r);
    const control = r.analysis.get('control', 0) as string;
    const arms = r.comparisons.col('arm')!.toList() as string[];
    const df = call.getParamValue('study') as DG.DataFrame | null;
    // the study's own name, as Monitor's header (the analysis table's name says only what was compared)
    return {name: df?.name || r.analysis.name, groups: [control, ...arms].map((name) => ({name, color: colors[name]})),
      status: `${r.analysis.get('line', 0)} · evaluated on ${dayText(r.analysis.get('evaluation_day', 0))}`,
      meta: () => metaOf(df), setMeta: df ? (m) => setMetaOn(df, m) : undefined,
      groupMenu: (g, e) => studyMenu(studyOf(df), 'group', g, {causedBy: e})};
  }), (s) => void openStage(s));
  return {
    root: stagePage(head.row, panes),
    ribbon: [[head.steps]],
    path: '/analyze',
    status: computed(() => busy.value ? 'Analyzing the study…' : result.value ?
      (call.getParamValue('study') as DG.DataFrame | null)?.name || result.value.analysis.name : ''),
  };
}
