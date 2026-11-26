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
import {badge, BasicTable, buttonWithIcon, computed, div, divH, divV, dot, effect, loader, Section, signal, span,
  tableFromMap, Tooltip} from '@datagrok-libraries/u2';
import {offerCopy, saveGate, tableSavers} from './native';
import {domainStorageOn, isReadOnly, markTableFromFile, newRowNote, open, rowOfTable, save, saveCopy, saveError,
  setRowOn, tableFileOwner, withCode} from './store';
import {FINISHED_TAG, LINE_TAG, STUDY_TAG, setHistoryOn, studyOf as savedStudyOf, stamped} from './study-file';
import {roleOf} from './study-reader';
import type {ReadonlySignal, StatValue} from '@datagrok-libraries/u2';
import {funcForm, TableInput, ViewHost, viewers} from '@datagrok-libraries/u2/dg';
import type {ViewContent} from '@datagrok-libraries/u2/dg';

import {AnalysisFrames, AnalysisInputs, PREDICTED, analysisDeclarationOf, analysisFrames, doublingOf, intervalNote,
  warmTwin} from './analyze';
import {median} from './engine';
import {EmptyState, NoControl} from './crux';
import {StudyLink, follow, groupColors, hasOwnColors, linkStudy, paint, resetColors} from './colors';
import {CorrectionInput} from './arms-input';
import {askedOf, autoRun, axisAt, callInputs, cardGrid, dayText, keepSame, onTableEdits, only, TABLE_EDITS, guide,
  exampleLine, infoNote, noteLine, notesBox, num, pct, plainError, pValue, settleForm, signed, isBodyWeight,
  keepNamedColumns,
  sized, markStage,
  statCard, swatch, volumeAxis, wrapCells, stagePanes, whenSettled, hintLink, PanelFact, PanelSection,
  showDepth, textLines, factList} from './ui';
import {groupsShown, lineOf, logOf, planAtlas, restoreFields} from './study-file';
import {TWIN_NOTE, savedHistory} from './history-models';
import {markLine} from './ui';
import {startMismatch} from './calendar';
import {PLAN_TAG, planKeyOf} from './plan-file';
import {openStage} from './stages';
import {noEffectNote, observedCell, planInput, planTab} from './plan-vs-observed';
import {nciReadings, NciReading, tolerabilityTab} from './tolerability';
import {BandChart} from './band-chart';
import {explore, StudyInputs, studyInputsOf} from './explore';
import {editLayout, panelTabs, Tab} from './panels';
import {stagePage, studyHeader} from './study-header';
import {AnimalSelection, exclusionsOf, mirrorSelection, studyDataLine, studyMenu, studyOf, tagStudy,
  keepStudyField, withoutExcluded, animalLine, animalMarks, animalPanel as panelOfAnimal} from './study-table';
import {armParts, emptyMeta, metaOf, parseLine, setMetaOn} from './study-meta';
import {ExposureSection, GroupDetails, detailsElement, exposureSection, groupDetails, sectionElement}
  from './translation-report';
import {Exports, download, exportMenu, reportHtml} from './report';
import {esc, printSheet} from './sheet';
import {SpiderChart, WaterfallChart, spiderLines} from './spider-chart';
import {LineInput} from './line-input';
import {Competing, Survival, competingOf, survivalOf} from './survival';
import {KmChart} from './km-chart';
import {CombinationStep, combinationStep, notCombined, stepWhy, stepWords, synergyAlpha, synergyDialog}
  from './synergy';
import {findCombinations, readProblem} from './arm-name';

/** Each combination's step on the ladder, in the words shown under it, for the regulatory table. */
const stepsInWords = (steps: Map<string, {step: CombinationStep}>) =>
  new Map([...steps].map(([ab, x]) => [ab, stepWords(x.step)]));
import {ARRIVE_CREDIT, arriveDialog, arriveItems, checklistBlocks, checklistDialog} from './arrive';
import {CTD_TITLE, M4S_URL, ctdBlocks, ctdItems, ctdStudyOf} from './ctd';
import {DOCX_TYPE, docx} from './docx';
import {CURVE_BAND, DrugExposure, ExposureAxis, ExposureCurve, HumanExposure, NcaByArm, ReferenceRow,
  SavedExposure,
  aucUsed, concentrationMapping, cvOk, drugMultiples,
  exposureAxis, exposureCurve, exposureHint, exposureJson, groupColumnOf, humanMultiple, medianDay, medianTgi,
  ncaByArm,
  ncaInstalled, parseExposure, readout, referenceRow, runNca}
  from './exposure';
import {EXPOSURE_TAG, exposureOfTag, roleColumns} from './study-file';
import {humanReference, linePrior, modelName} from './priors';
import {humanExposureDialog, humanExposureSummary} from './human-exposure';
import {PkSource, modellingDatasetDialog} from './modelling-dataset';

/** Analyze's tabs and their panels, in order and width (twelfths): the answer and how the groups grew; each animal;
 *  how the animals bore it; whether the study went as planned. */
const ANALYZE_TABS: Tab[] = [
  {name: 'Results', layout: [{id: 'groups', span: 12}, {id: 'growth', span: 12}, {id: 'km', span: 12}]},
  {name: 'Each animal', layout: [{id: 'waterfall', span: 6}, {id: 'spider', span: 6},
    {id: 'animals', span: 12}]},
  {name: 'Tolerability', layout: [{id: 'tolerability', span: 12}]},
  {name: 'Against the plan', layout: [{id: 'plan', span: 12}, {id: 'controls', span: 12}]}];

/** Under a growth chart: how many animals each group had on each measurement day, said only when a group thins (with
 *  no dropouts it is a grid of one number). */
function animalsPerDay(volumes: DG.DataFrame, groups: string[], colors: {[g: string]: string}): HTMLElement | null {
  const days = [...new Set(volumes.col('day')!.toList() as number[])].sort((a, b) => a - b);
  const n = new Map<string, number>();
  for (let i = 0; i < volumes.rowCount; i++)
    n.set(`${volumes.get('group', i)}|${volumes.get('day', i)}`, volumes.get('n', i));
  const lost = groups.map((g) => {
    const counts = days.map((d) => n.get(`${g}|${d}`)).filter((x): x is number => typeof x === 'number');
    // left by the end, not the max minus the min: a group measured from a later day has not lost any
    return {g, n: counts.length ? Math.max(...counts) - counts[counts.length - 1] : 0};
  }).filter((x) => x.n > 0);
  if (!lost.length)
    return null;
  const table = document.createElement('table');
  table.className = 'pmx-n-table';
  const head = table.createTHead().insertRow();
  const th = (text: string) => Object.assign(document.createElement('th'), {textContent: text});
  head.append(th(''), ...days.map((d) => th(String(d))));
  const body = table.createTBody();
  for (const g of groups) {
    const row = body.insertRow();
    row.insertCell().append(swatch(g, colors[g]));
    for (const d of days)
      row.insertCell().textContent = String(n.get(`${g}|${d}`) ?? '');
  }
  const lead = div([], 'pmx-n-lead');
  lead.textContent = `Animals by day: ${lost.map((x) => `${x.g} lost ${x.n}`).join(', ')}`;
  return div([lead, div([table], 'pmx-n-scroll')]);
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
  return `${side(effect)} than it would be untreated (90% range: ${side(lower)} to ${side(upper)})`;
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
 *  half its size and no tumor at all, in words; each tick over its guide in the bars. The labels take two heights in
 *  turn, so neighbours never overlap in a narrow column. */
function effectAxis(lo: number): HTMLElement {
  const axis = div([], 'pmx-effect-axis');
  Object.assign(axis.style, {position: 'relative', height: '34px', fontSize: '11px', fontWeight: '400',
    color: 'var(--dg-text-color-muted, #68727f)'});
  for (const [i, [x, text, shift]] of ([[lo, 'Larger', '0'], [0, 'Same as twin', '-50%'], [50, 'Half size', '-50%'],
    [100, 'No tumor', '-100%']] as [number, string, string][]).entries()) {
    const label = span(text);
    Object.assign(label.style, {position: 'absolute', top: i % 2 ? '13px' : '0', left: axisAt(x, lo),
      transform: `translateX(${shift})`,
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

/** A group's result, with a red "!" after it when the group is over the welfare limits (NCI's): the criteria on
 *  hover, the Tolerability tab on click. */
function withNci(el: HTMLElement, x: NciReading | undefined, open: () => void): HTMLElement {
  if (x?.kind !== 'error')
    return el;
  const mark = span('!', 'pmx-nci-mark');
  const words = `Over welfare limits: ${x.criteria}.`;
  Object.assign(mark.style, {display: 'inline-flex', alignItems: 'center', justifyContent: 'center', width: '14px',
    height: '14px', borderRadius: '50%', flex: '0 0 auto', background: 'var(--dg-error, #EB6767)', color: '#fff',
    fontSize: '10px', fontWeight: '700', lineHeight: '1', cursor: 'pointer'});
  mark.setAttribute('role', 'button');
  mark.setAttribute('aria-label', `${words} Open Tolerability`);
  mark.tabIndex = 0;
  Tooltip.bind(mark, () => `${words} Click to see Tolerability.`);
  mark.addEventListener('click', (ev) => {
    ev.stopPropagation(); // not the row's click too (a treated group's details)
    open();
  });
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

/** The note on a result that its T/C interval seems to contradict (analyze.ts intervalNote), or null. */
function groupsNote(r: AnalysisFrames): string | null {
  const c = r.comparisons;
  return intervalNote(Array.from({length: c.rowCount}, (_, i) => ({activity: c.get('activity', i) as string,
    tcLower: c.get('tc_lower', i) as number | null, tcUpper: c.get('tc_upper', i) as number | null})),
  r.analysis.get('evaluation_day', 0));
}

/** The results in words, a line each: the study, each treated group's effect and result, the test, the body weight. */
function summaryLines(r: AnalysisFrames, nci: Map<string, NciReading>): string[] {
  const twinned = r.analysis.get('twinned', 0);
  const v = (k: string) => r.analysis.get(k, 0);
  const c = r.comparisons;
  const control = v('control');
  const lines = [`${modelName(v('line'))} study evaluated on day ${v('evaluation_day')}` + (v('from_endpoint') ?
    `, the first measurement day the median control reached ${v('endpoint')} mm³` : '') +
    `; control: ${control}, ${v('n_control')} animals.`];
  for (let i = 0; i < c.rowCount; i++) {
    const g = (k: string) => c.get(k, i);
    const p = pValue(g('p_value'));
    const tc = g('tc') === null ? '' : `, T/C ${pct(g('tc') / 100)}` + (g('tc_lower') === null ? '' :
      ` (95% CI ${Math.round(g('tc_lower'))}–${Math.round(g('tc_upper'))}%)`);
    lines.push(`${g('arm')} (${g('n')} animals): TGI ${pct(g('tgi') / 100)}${tc}; ${g('activity').toLowerCase()} ` +
      `against the control (p ${p.startsWith('<') ? p : `= ${p}`})` + (twinned ? `; ${g('responders')} of ` +
      `${g('twin_animals')} animals responded (below the range they would reach untreated).` : '.'));
  }
  const note = groupsNote(r);
  if (note)
    lines.push(note);
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
      over.map((g) => `${g} over welfare limits: ${nci.get(g)!.criteria}.`).join(' ') :
      'No group over welfare limits.') + ' (Welfare limits: the NCI limits, Plowman et al. 1997.)');
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
  // a failure as the user reads it: what is missing as it is said, else what failed and what to try
  const AGAIN = 'Try again; if it fails again, reload the page.';
  const failure = (e: unknown, fallback: string) =>
    e instanceof EmptyState || e instanceof NoControl ? e.message : plainError(e, fallback);
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
  // the study and columns the shown result read (the tables opened from Analyze's menu read them)
  let resultStudy: StudyInputs | null = null;
  // the same without its excluded rows: what the analysis read, and so what KM, Fine-Gray, synergy and the tables read
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
      const blinded = call.getParamValue('study');
      if (blinded instanceof DG.DataFrame && !groupsShown(blinded)) {
        error.value = 'Reveal the groups in Monitor first.';
        return;
      }
      const sent = askedOf(call);
      const reads = studyInputsOf(call);
      tagStudy(reads);
      // the correction's count follows the table now picked, before (and whether or not) it can be analyzed: a
      // one-group table failing to analyze kept the last study's "over 2 comparisons"
      if (reads) {
        const groups = new Set(reads.study.col(reads.group)?.toList().filter((g) => g !== null && g !== '') ?? []);
        treatedGroups.value = Math.max(0, groups.size - 1);
      }
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
          error.value = failure(e, `Each animal could not be compared with its twin. ${AGAIN}`);
      });
    } catch (e: any) {
      // only the last run speaks: an older one failing late would write over a newer result's place
      if (mine === latest) {
        error.value = failure(e, `The analysis could not be computed. ${AGAIN}`);
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
  whenSettled(form, () => {
    const studyInput = form.getInput('study') as unknown as TableInput | null;
    keepStudyField(studyInput, () => call.getParamValue('study') ?? null);
    // the picker offers the open tables that look like a study (tagged as one by a stage, or with a volume or caliper
    // column by its name), and the one picked: not every table open (a concentrations table opened for NCA is none)
    if (studyInput) {
      const offer = studyInput.setItems.bind(studyInput);
      const looksLikeStudy = (name: string) => {
        const t = grok.shell.tables.find((x) => x.name === name);
        return !t || name === studyInput.value || !!t.getTag(STUDY_TAG) ||
          t.columns.names().some((c) => ['volume', 'length'].includes(roleOf(c)?.role ?? ''));
      };
      studyInput.setItems = (items: string[]) => offer(items.filter(looksLikeStudy));
      studyInput.setItems(grok.shell.tables.map((t) => t.name));
    }
    if (opening?.study)
      void openStudy(async () => opening.study!);
  });
  planInput(form, call, planState);

  // exposure from the NCA package, when installed: a concentration table picked here, or opened while the field is
  // empty, gives each group NCA's AUC (AUCinf, else AUClast) and Cmax beside its effect (a description, kept with
  // the study: never analyzed)
  const ncaRun = signal<Awaited<ReturnType<typeof runNca>> | null>(null);
  const ncaSays = signal<string | null>(null);
  // the concentrations NCA read (its table and mapping), for the modelling dataset
  let ncaSource: PkSource | null = null;
  const exposureField = ncaInstalled() ? form.build(() => new TableInput({label: 'Exposure (NCA)', nullable: true,
    tooltipText: 'Plasma concentrations over time, one row per animal and sample: NCA gives each group its AUC ' +
      '(AUCinf where NCA gives it, else AUClast) and Cmax, shown beside its effect and kept with the study (not the ' +
      'concentrations). The analysis does not change.'})) : null;
  // the patients' exposure, by drug, typed under Advanced (empty until typed, or kept with the study): each group's
  // mouse exposure is then also shown as a multiple of it
  const humans = signal(new Map<string, HumanExposure>());
  // what the study kept of NCA's exposure, shown until a concentration table is read again
  const savedNca = signal<NcaByArm | null>(null);
  // NCA's exposure shown now (read again or kept), to keep with the patients' when they change
  let shownNca: NcaByArm | null = null;
  const referenceRows = signal<ReferenceRow[]>([]);
  void humanReference().then((t) => referenceRows.value = t.drugs).catch(() => {});
  const humanField = exposureField ? ui.input.string('Patients\' exposure', {value: '',
    placeholder: 'None: click to enter',
    tooltipText: 'The patients\' exposure, measured or projected, for each drug of the groups: the groups table then ' +
      'shows the mouse exposure as a multiple of it. It is kept with the study, and the analysis does not change.'}) :
    null;
  // the exposure on the study's table (EXPOSURE_TAG), so that saving the study keeps it
  let sectionNow: () => ExposureSection | null = () => null;
  // the group details shown in the context panel, cleared with the view
  let lastDetails: HTMLElement | null = null;
  let lastGroup = '';
  // a group's animals for its details (set with the animal list)
  let animalsOf: (group: string, color: string) => PanelSection[] = () => [];
  // each group's statistics, as the report's Statistics table has them
  let statsOf = new Map<string, PanelFact[]>();
  const groupSections = (group: string, color: string): PanelSection[] => [...animalsOf(group, color),
    ...statsOf.has(group) ? [{title: 'Statistics', content: () => factList(statsOf.get(group)!)}] : []];
  self.own(() => {
    if (lastDetails && grok.shell.o === lastDetails)
      grok.shell.o = null;
  });
  let detailsNow: (arm: string) => GroupDetails = (arm) => groupDetails({group: arm, read: '', tgi: null, day: null,
    x: undefined, table: null, live: false, multiple: null, human: undefined, curve: null});
  // the analysis's effect, kept with the exposure for Compare (only while it is this study's: its groups in it)
  let shownEffect: SavedExposure['effect'] | undefined;
  // the patients' exposure on the study's tag as Analyze last read or wrote it: one set there since (Compare's "Add
  // patients' exposure…") is taken in, not written over
  let humansSeen = '';
  const humansOnTag = (df: DG.DataFrame | null) => JSON.stringify(exposureOfTag(df)?.humans ?? null);
  const keepExposure = (nca: NcaByArm | null, typedHere = false) => {
    const df = call.getParamValue('study') as DG.DataFrame | null;
    if (!typedHere && df && humansOnTag(df) !== humansSeen) {
      humans.value = parseExposure(exposureOfTag(df)).humans;
      if (humanField)
        humanField.value = humanExposureSummary(humans.peek());
    }
    let groups = new Set<string>();
    try {
      groups = new Set(df ? roleColumns(df).group.categories : []);
    } catch {/* no group column: no effect kept */}
    const effect = shownEffect && Object.keys(shownEffect.tgi).every((a) => groups.has(a)) ? shownEffect : undefined;
    const exposure = exposureJson(humans.peek(), nca);
    const kept = effect ? {...exposure ?? {}, effect} : exposure;
    if (!df || (!kept && !df.getTag(EXPOSURE_TAG)))
      return;
    const text = JSON.stringify(kept ?? {});
    if (df.getTag(EXPOSURE_TAG) !== text)
      df.setTag(EXPOSURE_TAG, text);
    humansSeen = humansOnTag(df);
  };
  const openHumans = () => {
    const r = result.peek();
    if (!r) {
      grok.shell.info('Analyze a study first: the drugs are read from its groups');
      return;
    }
    void humanExposureDialog(r.comparisons.col('arm')!.toList() as string[], humans.peek(), (m) => {
      humans.value = m;
      if (humanField)
        humanField.value = humanExposureSummary(m);
      keepExposure(shownNca, true);
    }, metaOf(call.getParamValue('study') as DG.DataFrame | null))
      .catch((e) => grok.shell.error(plainError(e, 'The patients\' exposure could not be kept. Try again.')));
  };
  if (humanField) {
    humanField.readOnly = true;
    humanField.input.style.cursor = 'pointer';
    humanField.input.addEventListener('click', openHumans);
  }
  // a study opened with exposure it kept: shown at once, without asking NCA again
  let exposureStudy: DG.DataFrame | null | undefined;
  // a concentration table was read in the field (so emptying it is the user's clear)
  let fieldLoaded = false;
  form.effect(() => {
    if (!result.value)
      return;
    const df = call.getParamValue('study') as DG.DataFrame | null;
    if (df === exposureStudy)
      return;
    // another study's NCA run is not this one's: the field emptied (not a clear: what this study kept stays); the
    // first study keeps a table picked before it was analyzed
    if (exposureStudy !== undefined) {
      fieldLoaded = false;
      if (exposureField?.table)
        exposureField.table = null;
      ncaRun.value = null;
      ncaSource = null;
    }
    exposureStudy = df;
    const kept = parseExposure(exposureOfTag(df));
    // only what this study kept: a drug's exposure typed for another study is not carried into this one
    humans.value = kept.humans;
    humansSeen = humansOnTag(df);
    savedNca.value = kept.nca;
    if (humanField)
      humanField.value = humanExposureSummary(humans.peek());
  });
  if (exposureField) {
    whenSettled(form, () => {
      form.getInput('density')?.root.after(exposureField.root);
      if (humanField)
        exposureField.root.after(humanField.root);
    });
    let reading = 0;
    form.effect(() => {
      const t = exposureField.table;
      const mine = ++reading;
      ncaRun.value = null;
      ncaSays.value = null;
      ncaSource = null;
      if (!t) {
        // the field cleared: what the study kept goes too
        if (fieldLoaded) {
          fieldLoaded = false;
          savedNca.value = null;
          keepExposure(null);
        }
        return;
      }
      fieldLoaded = true;
      const missed: string[] = [];
      void concentrationMapping(t, missed).then(async (m) => {
        if (mine !== reading)
          return;
        if (!m) {
          ncaSays.value = `NCA does not read ${t.name} as concentrations over time` +
            `${missed.length ? ` (${missed.join('; ')})` : ''}: it reads the units from the column names, as in ` +
            'TIME_H, CONC_NG_ML and DOSE_MG_KG.';
          return;
        }
        const run = await runNca(t, m);
        if (mine === reading) {
          ncaSource = {table: t, mapping: m, groupCol: run.groupCol};
          ncaRun.value = run;
        }
      }).catch((e) => {
        if (mine === reading)
          ncaSays.value = `NCA could not compute exposure: ${e?.message ?? e}`;
      });
    });
    self.subscribe(grok.events.onTableAdded, (e: DG.EventData<DG.DataFrameArgs>) => {
      const t = e.args.dataFrame;
      if (exposureField.value || t.name === form.getInput('study')?.value)
        return;
      void concentrationMapping(t).then((m) => {
        if (m && !exposureField.value)
          exposureField.table = t;
      });
    });
  }

  const start = divV([
    span('Open your study table (one row per animal and measurement day, day 0 the randomization), check its ' +
      'columns: the analysis follows them, and any change you make to the table.'),
    exampleLine(),
    infoNote('<b>Analyzed as planned.</b> Analyze applies the analysis Plan pre-specifies: the study is evaluated ' +
      'when the median control reaches the endpoint, each treated group is compared with the control on tumor ' +
      'growth, and each treated animal with what it would measure without treatment, predicted ' +
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
  const openSynergy = () => {
    const r = result.peek();
    const df = call.getParamValue('study') as DG.DataFrame | null;
    if (!r || !keptStudy)
      return;
    // the plan's α if a plan is attached (the one picked in the form, else the study's), else Analyze's
    const planText = call.getParamValue('plan')?.trim() || df?.getTag(PLAN_TAG);
    synergyDialog(keptStudy, r.analysis.get('control', 0), r.comparisons.col('arm')!.toList() as string[],
      r.analysis.get('evaluation_day', 0), synergyAlpha(planText, r.analysis.get('alpha', 0)), notCombined(planText));
  };
  // each combination found in the groups' names (unless the plan unticked it): where it stands against its two drugs,
  // under its result; computed after the analysis, shown when it arrives
  const ladder = signal(new Map<string, {step: CombinationStep, why: string}>());
  // what the shown ladder was computed on: the study read (a new one each analysis run, the same when only the twin
  // arrives) and the analysis' choices, the arms the plan unticked among them
  let laddered: {s: StudyInputs | null, key: string} = {s: null, key: ''};
  // a combination left without its step against its drugs, and why
  const ladderNotes = notesBox();
  form.effect(() => {
    const r = result.value;
    const s = keptStudy;
    const treated = (r?.comparisons.col('arm')!.toList() ?? []) as string[];
    const planText = call.getParamValue('plan')?.trim() ||
      (call.getParamValue('study') as DG.DataFrame | null)?.getTag(PLAN_TAG);
    const off = notCombined(planText);
    const key = r && s ? [...treated, '|', ...off, r.analysis.get('control', 0), r.analysis.get('line', 0),
      r.analysis.get('evaluation_day', 0)].join('\u0000') : '';
    if (s === laddered.s && key === laddered.key)
      return;
    const mine = laddered = {s, key};
    if (ladder.peek().size)
      ladder.value = new Map();
    ladderNotes.replaceChildren();
    const say = (lines: string[]) => {
      if (laddered === mine)
        ladderNotes.replaceChildren(...lines.map((l) => noteLine(l)));
    };
    if (!r || !s)
      return;
    const combos = findCombinations(treated).filter((c) => !off.includes(treated[c.ab]));
    if (!combos.length)
      return;
    const control: string = r.analysis.get('control', 0);
    const lastDay: number = r.analysis.get('evaluation_day', 0) ?? Infinity;
    const line: string = r.analysis.get('line', 0);
    void linePrior(line).then(async (lp) => {
      const logBeta = lp.study_prior?.mean_log_beta;
      if (logBeta === undefined) {
        say(combos.map((c) => `${treated[c.ab]}: not compared with its two drugs (${line} has no reference ` +
          'growth rate).'));
        return;
      }
      const found = new Map<string, {step: CombinationStep, why: string}>();
      const missing: string[] = [];
      for (const c of combos) {
        const [a, b, ab] = [treated[c.a], treated[c.b], treated[c.ab]];
        const x = await combinationStep(s, {control, a, b, ab}, Math.exp(logBeta), lastDay).then((v) => v ?? 'few',
          (e) => plainError(e, `${ab}: the comparison with its two drugs could not be computed. ${AGAIN}`));
        if (x === 'few')
          missing.push(`${ab}: too few animals to compare it with its two drugs.`);
        else if (typeof x === 'string')
          missing.push(x);
        else
          found.set(ab, {step: x.step, why: stepWhy(a, b, x.floored)});
      }
      if (laddered === mine && found.size)
        ladder.value = found;
      say(missing);
    }).catch((e) => say([plainError(e, `The combinations could not be compared with their drugs just now. ${AGAIN}`)]));
  });
  // above the groups, one sentence on exposure and effect; under them, a muted line on what would add to it
  // (the report gives the sentence once, in its section on exposure and patients)
  const readoutLine = div([], 'pmx-noprint');
  const hintLine = div([]);
  // a result its T/C interval seems to contradict: one quiet line on why, under the groups
  const intervalLine = div([]);
  Object.assign(hintLine.style, {display: 'flex', alignItems: 'center', gap: 'var(--dg-space-s)'});
  const muted = (text: string) => {
    const e = span(text);
    e.style.color = 'var(--grey-4)';
    return e;
  };
  // the concentrations' picker, as the Exposure (NCA) field's own button opens it
  const addConcentrations = () => {
    const open = exposureField?.actions.find((a) => a.id === 'open-file');
    if (open) {
      void Promise.resolve(open.run(exposureField!)).catch((e) =>
        grok.shell.error(plainError(e, 'The concentrations could not be opened. Try again.')));
    } else {
      exposureField?.root.scrollIntoView({behavior: 'smooth', block: 'center'});
      exposureField?.root.querySelector<HTMLElement>('input, select, [tabindex]')?.focus();
    }
  };
  const groupNotes = notesBox();
  // NCA's exposure: why a group has none, and exposure against effect when the groups span 3 doses or more
  const exposureNotes = notesBox();
  const exposurePoints = signal<{arm: string, tgi: number, auc: number, aucUnit: string, dose: string, drug: string,
    route: string, color: string}[] | null>(null);
  // one drug by one route at 3 doses or more: exposures of different drugs share no axis
  const exposureShown = computed(() => {
    const p = exposurePoints.value;
    return p && new Set(p.map((x) => `${x.drug}|${x.route}`)).size === 1 && new Set(p.map((x) => x.dose)).size >= 3 ?
      p : null;
  });
  // the groups on one axis with their 1 − T/C of medians (exposure.ts exposureAxis, medianTgi: the measure the curve
  // and its band were validated on, not the groups table's TGI), read at once; then the effect expected between the
  // doses tested (exposureCurve), or why there is none
  const exposureGroups = signal<{axis: ExposureAxis, unit: string, day: number | null,
    points: {arm: string, x: number, tgi: number, color: string}[]} | null>(null);
  const exposureFit = signal<ExposureCurve | {why: string} | null>(null);
  let fitting = '';
  form.effect(() => {
    const p = exposureShown.value;
    const r = result.peek();
    const s = keptStudy;
    const control: string = r?.analysis.get('control', 0);
    const lastDay: number | null = r?.analysis.get('evaluation_day', 0) ?? null;
    const key = p && r && s ? [JSON.stringify(p), dataRead.peek(), control, lastDay].join('\u0000') : '';
    if (key === fitting)
      return;
    fitting = key;
    exposureFit.value = null;
    exposureGroups.value = null;
    if (!p || !r || !s)
      return;
    const [animal, group, day, volume] = [s.animal, s.group, s.day, s.volume].map((c) => s.study.col(c)!);
    const rows = Array.from({length: s.study.rowCount}, (_, i) => i)
      .filter((i) => !animal.isNone(i) && !group.isNone(i) && !day.isNone(i) && !volume.isNone(i))
      .map((i) => ({animal: String(animal.get(i)), group: String(group.get(i)), day: day.get(i),
        volume: volume.get(i)}));
    const groups = p.map((x) => ({arm: x.arm, auc: x.auc}));
    // each group's dose from its details first, then its name
    const studyMeta = metaOf(call.getParamValue('study') as DG.DataFrame | null);
    const axis = exposureAxis(groups, studyMeta);
    const points = p.flatMap((x, i) => {
      const tgi = medianTgi(rows, x.arm, control, lastDay);
      return tgi === null ? [] : [{arm: x.arm, x: axis.of(groups[i]), tgi, color: x.color}];
    });
    // the day the medians are read on, of the groups plotted: named on the axis only when they share it
    const days = new Set(points.map((x) => medianDay(rows, x.arm, control, lastDay)));
    exposureGroups.value = {axis, unit: axis.axis === 'dose' ? axis.unit : p[0].aucUnit,
      day: days.size === 1 ? [...days][0] : null, points};
    void exposureCurve(groups, rows, control, lastDay, studyMeta).then((c) => {
      if (fitting === key)
        exposureFit.value = c;
    }, (e) => {
      if (fitting === key)
        exposureFit.value = {why: plainError(e, 'it could not be computed just now; try again').replace(/\.$/, '')};
    });
  });
  // the groups as points (a band of no width each), the curve dashed with its band; on AUClast, or on the dose when
  // exposure does not rise in step with it (a chart for each axis, the one with data shown)
  const CURVE = 'Expected between the doses tested';
  const exposureData = computed(() => {
    const g = exposureGroups.value;
    if (!exposureShown.value || !g)
      return undefined;
    const f = exposureFit.value;
    const xName = g.axis.axis === 'dose' ? 'Dose' : 'AUClast';
    const df = DG.DataFrame.fromObjects([
      ...g.points.map((x) => ({series: x.arm, [xName]: x.x, lower: x.tgi, upper: x.tgi, tgi: x.tgi})),
      ...(f && 'curve' in f ? f.curve : []).map((c) => ({series: CURVE, [xName]: c.x, lower: c.lo, upper: c.hi,
        tgi: c.tgi}))])!;
    if (!df.rowCount)
      return undefined;
    df.col(xName)!.meta.units = g.unit;
    df.col('tgi')!.meta.units = '%';
    exposureOptions.forEach((o) => o.yTitle = `TGI (medians${g.day !== null ? `, day ${g.day}` : ''})`);
    df.col('series')!.meta.colors.setCategorical({[CURVE]: '#9aa5b1',
      ...Object.fromEntries(g.points.map((x) => [x.arm, x.color]))});
    return df;
  });
  // the Y axis and the tooltip name the measure, with its day when the groups share one (set as the data is built)
  const exposureOptions = ['AUClast', 'Dose'].map((x) => ({x, lower: 'lower', upper: 'upper', line: 'tgi',
    split: 'series', yTitle: 'TGI (medians)', yInTooltip: true, markers: true, dashed: [CURVE]}));
  const exposureCharts = exposureOptions.map((o) => new BandChart(computed(() =>
    exposureData.value?.col(o.x) ? exposureData.value : undefined), o));
  const exposureText = span('');
  form.effect(() => {
    const f = exposureFit.value;
    const day = exposureGroups.value?.day;
    const tgi = `Each point is one group, at its TGI (medians${day != null ? `, day ${day}` : ''}): 1 − T/C of the ` +
      'median tumor volumes on the last day up to the evaluation day with 3 animals in the group and the control. ' +
      'The curve and its band were validated on this measure; the groups table gives the analysis\'s TGI.';
    exposureText.textContent = `${tgi} ${f && 'curve' in f ? 'The dashed line is the effect expected between the ' +
      `doses tested, never beyond them.${f.note ? ` ${f.note}` : ''} The band is ±${CURVE_BAND} points: in the ` +
      'validation, 4 doses left out in 5 were within it (12 doses, 9 studies); some missed by 20 to 30.' :
      f ? `No curve: ${f.why}.` : ''}`;
    exposureCharts.forEach((c) => c.root.style.display = exposureData.value?.col(c === exposureCharts[0] ?
      'AUClast' : 'Dose') ? '' : 'none');
  });
  const exposureBox = divV([exposureText, ...exposureCharts.map((c) => sized(c, 260))]);
  form.effect(() => exposureBox.style.display = exposureData.value ? '' : 'none');
  const statsTable = div([]);
  const controlCheck = div([]);
  const protocol = div([]);
  let protocolText = '';
  let statsCsv = '';
  // each group's measured Kaplan-Meier median (set by the Time to endpoint panel, never its projection): the one
  // "median day to the endpoint" the statistics and their file give, as that panel does
  const kmMedians = signal<Survival['medians'] | null>(null);
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
    const m = metaOf(call.getParamValue('study') as DG.DataFrame | null);
    return groupColors(r.comparisons.col('arm')!.toList(), r.analysis.get('control', 0), m?.colors, m?.slots);
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
    const steps = ladder.value;
    const allNotes = String(v('notes') ?? '').split(' | ').filter((n) => n);
    const odd = allNotes.filter((n) => n.startsWith('Control '));
    groupNotes.replaceChildren(...allNotes.filter((n) => !odd.includes(n)).map((n) => noteLine(n)));

    const nci = nciReadings(r);
    // each treated group's exposure, matched to it by NCA's group name or its drug, dose and route
    // (NCA read again replaces what the study kept, and is kept in its place)
    const run = ncaRun.value;
    const meta = metaOf(call.getParamValue('study') as DG.DataFrame | null);
    const nca = run ? ncaByArm(r.comparisons.col('arm')!.toList(), run, ncaSource?.table.name ?? '', meta) :
      savedNca.value;
    shownNca = nca;
    shownEffect = {day: v('evaluation_day') ?? null, tgi: Object.fromEntries((r.comparisons.col('arm')!.toList() as
      string[]).map((a, i) => [a, r.comparisons.get('tgi', i)]).filter(([, t]) => Number.isFinite(t))),
    ctd: ctdStudyOf(r, meta, r.analysis.name, nca, stepsInWords(steps))};
    keepExposure(nca);
    const exposure = nca?.arms ?? null;
    // a combination's drugs, each with its own exposure (none: no NCA, or nothing for any of its drugs)
    const drugX = nca?.drugs ?? new Map<string, DrugExposure[]>();
    // each treated group's exposure as a multiple of the human's, once a human exposure is typed for its drug
    const human = humans.value;
    const multiples = new Map<string, ReturnType<typeof humanMultiple>>();
    if (nca && [...human.values()].some((h) => h.auc !== null)) {
      for (const arm of nca.arms.keys())
        multiples.set(arm, humanMultiple(arm, nca.arms.get(arm)?.exposure, nca.pk.get(arm) ?? null, human, meta));
    }
    // a combination's drugs, each against its own patients' exposure
    const drugMs = new Map<string, ReturnType<typeof drugMultiples>>();
    if ([...human.values()].some((h) => h.auc !== null)) {
      for (const [arm, ds] of drugX)
        drugMs.set(arm, drugMultiples(arm, ds, human, meta));
    }
    const drugWhys = [...drugMs.values()].flat().map((d) => d.m?.why).filter((w): w is string => !!w);
    const withHuman = [...multiples.values()].some((m) => m !== null) ||
      [...drugMs.values()].some((ds) => ds.some((d) => d.m !== null));
    // exposure kept under a group name the table no longer has (renamed before Monitor saved the study): said
    const unmatched = !run && exposure ? [...exposure].filter(([a, x]) => x.exposure && a !== control &&
      !(r.comparisons.col('arm')!.toList() as string[]).includes(a)).map(([a]) => a) : [];
    exposureNotes.replaceChildren(...ncaSays.value ? [noteLine(ncaSays.value)] : [], ...!exposure ? [] :
      [...!run && nca ? [muted(`Exposure from NCA, saved with the study (table ${nca.table || 'not named'})`)] : [],
        ...unmatched.length ? [noteLine(`Exposure kept for ${unmatched.join(', ')}: no group has that name now, so ` +
          'it is not shown')] : [],
        ...[...exposure.values(), ...[...drugX.values()].flat().map((d) => d.x)].filter((x) => x.why)
          .map((x) => noteLine(x.why!)),
        ...[...multiples.values()].filter((m) => m?.why).map((m) => noteLine(m!.why!)),
        ...drugWhys.map((w) => noteLine(w))]);
    // one sentence on exposure and effect, and the one thing that would add to it
    const arms = r.comparisons.col('arm')!.toList() as string[];
    const said = readout(arms.map((arm, i) => {
      const x = exposure?.get(arm)?.exposure;
      const m = multiples.get(arm);
      const own = armParts(arm, meta);
      const h = human.get(own.length === 1 ? own[0].key : '');
      const a = x ? aucUsed(x) : null;
      return {arm, tgi: r.comparisons.get('tgi', i), auc: a?.auc ?? null, aucUnit: a?.unit ?? '',
        multiple: m?.ratio !== undefined ? {ratio: m.ratio, range: m.range, cvKnown: cvOk(h?.cv), weekly: !!m.weekly} :
          null};
    }), meta);
    readoutLine.replaceChildren(...said ? [span(said)] : []);
    // the report's section on exposure and patients, as the groups are now (the curve as drawn when it is asked)
    sectionNow = () => exposureSection({sentence: said, curve: !!exposureFit.peek() && 'curve' in exposureFit.peek()!,
      table: nca ? nca.table : null, rows: arms.flatMap((arm, i) => {
        const parts = armParts(arm, meta);
        const tgi = r.comparisons.get('tgi', i);
        // a combination: a row per drug, "A + B (A)"
        const ds = drugX.get(arm);
        if (ds) {
          return parts.map((p) => ({group: `${arm} (${p.drug})`, tgi,
            x: ds.find((d) => d.key === p.key)?.x, human: human.get(p.key),
            multiple: drugMs.get(arm)?.find((d) => d.key === p.key)?.m ?? null}));
        }
        return [{group: arm, tgi, x: exposure?.get(arm), multiple: multiples.get(arm) ?? null,
          human: parts.length === 1 ? human.get(parts[0].key) : undefined}];
      })});
    const drugs = arms.flatMap((arm) => armParts(arm, meta));
    const hint = exposureHint({ncaInstalled: !!exposureField, mousePk: !!exposure && [...exposure.values(),
      ...[...drugX.values()].flat().map((d) => d.x)].some((x) => x.exposure),
    patients: drugs.some((p) => human.get(p.key)?.auc != null), placed: withHuman,
    why: [...multiples.values()].find((m) => m?.why)?.why ?? drugWhys[0], inTable: [...new Set(drugs
      .filter((p) => referenceRow(p.key, referenceRows.value)).map((p) => p.drug))]});
    // a treated group's details for the context panel, as computed here (the curve as drawn when it is asked)
    detailsNow = (arm) => {
      const i = arms.indexOf(arm);
      const parts = armParts(arm, meta);
      const f = exposureFit.peek();
      const pts = f && 'curve' in f ? [...exposureGroups.peek()?.points ?? []].sort((a, b) => a.x - b.x) : null;
      const at = pts ? pts.findIndex((p) => p.arm === arm) : -1;
      const g = (k: string) => i < 0 ? null : r.comparisons.get(k, i);
      const step = steps.get(arm);
      return groupDetails({group: arm, read: parseLine(arm, meta), tgi: g('tgi'), color: colors[arm],
        n: arm === control ? v('n_control') : g('n'),
        tc: g('tc') === null ? null : {value: g('tc'), lower: g('tc_lower'), upper: g('tc_upper')},
        result: g('activity'), step: step ? {words: stepWords(step.step), why: step.why} : null,
        responders: done && g('twin_animals') ? `${g('responders')} of ${g('twin_animals')}` : null,
        tgiCi: {lower: g('tgi_lower'), upper: g('tgi_upper')}, nadir: nci.get(arm)?.nadir ?? null,
        day: v('evaluation_day') ?? null, x: exposure?.get(arm), table: nca?.table ?? null, live: !!run,
        multiple: multiples.get(arm) ?? null, human: parts.length === 1 ? human.get(parts[0].key) : undefined,
        curve: !pts || arm === control ? null : at < 0 ? 'off' : at === 0 ? 'lowest' :
          at === pts.length - 1 ? 'highest' : 'between',
        // a combination: each drug's exposure apart
        ...drugX.has(arm) ? {drugs: parts.map((p) => ({drug: p.drug, x: drugX.get(arm)!.find((d) => d.key === p.key)?.x,
          human: human.get(p.key), multiple: drugMs.get(arm)?.find((d) => d.key === p.key)?.m ?? null}))} : {}});
    };
    // a group's details open in the context panel: built again from the analysis and exposure as they are now
    if (lastDetails && grok.shell.o === lastDetails && (arms.includes(lastGroup) || lastGroup === control)) {
      grok.shell.o = lastDetails = detailsElement(detailsNow(lastGroup),
        groupSections(lastGroup, colors[lastGroup]));
    }
    hintLine.replaceChildren(...hint ? [muted(hint.text), ...hint.add ? [ui.button('Add…', () =>
      hint.add === 'patients' ? openHumans() : addConcentrations())] : []] : []);
    exposurePoints.value = !exposure ? null : Array.from({length: r.comparisons.rowCount}, (_, i) => i)
      .filter((i) => exposure.get(r.comparisons.get('arm', i))?.exposure)
      .map((i) => ({arm: r.comparisons.get('arm', i) as string, tgi: r.comparisons.get('tgi', i) as number,
        color: colors[r.comparisons.get('arm', i)], ...exposure.get(r.comparisons.get('arm', i))!.exposure!}));

    // the control first, then each treated group: its effect on the evaluation day, its animals against their
    // twins and its result (a table built in this run is disposed with it)
    const c = r.comparisons;
    type Row = {name: string, cells: (string | HTMLElement)[]};
    const rows: Row[] = [{name: control, cells: [swatch(control, colors[control]),
      String(v('n_control')), '—', '—', '—', '—', '—', withNci(span('Reference'), nci.get(control), showNci)]}];
    for (let i = 0; i < c.rowCount; i++) {
      const g = (k: string) => c.get(k, i);
      const tc = g('tc');
      const verdict = divH([dot(g('activity') === 'Slows growth' ? 'success' : 'warning'), span(g('activity'))]);
      Object.assign(verdict.style, {alignItems: 'center', gap: 'var(--dg-space-s)'});
      // no clear effect: what the study could have detected with its animals
      const why = noEffectNote((k) => r.versusPlan.get(k, i), v);
      if (why)
        Tooltip.bind(verdict, why);
      // a combination: under its result, where it stands against its two drugs
      const step = steps.get(g('arm'));
      const stepLine = step ? muted(stepWords(step.step)) : null;
      if (step) {
        Tooltip.bind(stepLine!, step.why);
        // the same detail dialog the menu item used to open
        Object.assign(stepLine!.style, {cursor: 'pointer', textDecoration: 'underline'});
        // only the dialog: not the row's click too, which would open the context panel beside it
        stepLine!.onclick = (e) => {
          e.stopPropagation();
          openSynergy();
        };
      }
      const resultCell = stepLine ? divV([verdict, stepLine]) : verdict;
      // what its details and name say it received, on hover; under the name only what could not be read
      const name = swatch(g('arm'), colors[g('arm')]);
      const read = parseLine(g('arm'), meta);
      if (read)
        Tooltip.bind(name, read);
      const problem = readProblem(armParts(g('arm'), meta));
      const nameCell = problem ? divV([name, muted(`⚠ ${problem}`)]) : name;
      rows.push({name: g('arm'), cells: [nameCell, String(g('n')),
        pct(g('tgi') / 100),
        tc === null ? '—' : observedCell(pct(tc / 100), g('tc_lower') === null ? '' :
          `(${Math.round(g('tc_lower'))}–${Math.round(g('tc_upper'))}%)`),
        !done ? waitingTwin() : g('paired_effect') === null ? '—' : pct(g('paired_effect') / 100),
        !done ? waitingTwin() : `${g('responders')} of ${g('twin_animals')}`, pValue(g('p_value')).replace(' ', ' '),
        withNci(resultCell, nci.get(g('arm')), showNci)]});
    }
    const headers = ['Group', 'Analyzed', 'TGI', 'T/C (95% CI)', 'Below its predicted untreated size', 'Responders',
      'p', 'Result'];
    const table = new BasicTable<Row>({items: rows, keyOf: (x) => x.name, columns: headers.map((header, j) => ({header,
      render: (x: Row) => x.cells[j], align: j >= 1 && j < headers.length - 1 ? 'center' : 'left'})),
    // a group's details in the context panel
    onRowClick: (x) => {
      lastDetails = showDepth(detailsElement(detailsNow(lastGroup = x.name), groupSections(x.name, colors[x.name])));
    }});
    table.root.style.width = '100%';
    wrapCells(table.root);
    groupsTable.replaceChildren(table.root);

    // the statistics behind each result (the growth rates the test compares and their difference)
    // the median day to the endpoint is the Kaplan-Meier's, as Time to endpoint shows it (not tgi's T−C, a median
    // over the animals that reached it, which disagreed with the panel); '…' until the curves are computed
    const km = new Map((kmMedians.value ?? []).map((m) => [m.group, m.median]));
    const reached = (g: string) => !kmMedians.value ? '…' : !km.has(g) ? '—' :
      km.get(g) === null ? 'Not reached' : dayText(km.get(g));
    const later = (g: string) => {
      const [t, k] = [km.get(g), km.get(control)];
      return typeof t === 'number' && typeof k === 'number' ? (t - k).toFixed(1) : '—';
    };
    type Stat = {name: string, cells: string[]};
    const doubling = (s: number) => {
      const d = doublingOf(s);
      return d === null ? 'Does not grow' : d.toFixed(1);
    };
    const stats: Stat[] = [{name: control, cells: [control, String(num(c.get('slope_control_per_day', 0))),
      doubling(c.get('slope_control_per_day', 0)), '—', '—', reached(control), '—']}];
    for (let i = 0; i < c.rowCount; i++) {
      const g = (k: string) => c.get(k, i);
      stats.push({name: g('arm'), cells: [g('arm'), String(num(g('slope_per_day'))), doubling(g('slope_per_day')),
        `${num(g('slope_difference'), 2)} (${num(g('slope_difference_lower'), 2)} to ` +
        `${num(g('slope_difference_upper'), 2)})`, pValue(g('p_value')), reached(g('arm')), later(g('arm'))]});
    }
    const statHeaders = ['Group', 'Growth rate (ln V per day)', 'Doubling time (days)',
      'Difference from the control (95% CI)', 'p', `Median day to ${v('endpoint')} mm³`, 'Days later than the control'];
    const st = new BasicTable<Stat>({items: stats, keyOf: (x) => x.name, columns: statHeaders.map((header, j) =>
      ({header, render: (x: Stat) => x.cells[j], align: j >= 1 ? 'center' : 'left'}))});
    st.root.style.width = '100%';
    wrapCells(st.root);
    statsOf = new Map(stats.map((x) => [x.name,
      statHeaders.slice(1).map((h, j) => ({label: h, value: x.cells[j + 1]}))]));
    statsTable.replaceChildren(st.root);
    // the file leads with the headline numbers the groups table shows, TGI and T/C with their 95% CI (none for the
    // control), as numbers
    const headline = ['tgi', 'tgi_lower', 'tgi_upper', 'tc', 'tc_lower', 'tc_upper'];
    const headlineHeaders = ['TGI (%)', 'TGI 95% CI lower', 'TGI 95% CI upper', 'T/C (%)', 'T/C 95% CI lower',
      'T/C 95% CI upper'];
    const armRow = new Map((c.col('arm')!.toList() as string[]).map((a, i) => [a, i]));
    statsCsv = DG.DataFrame.fromObjects(stats.map((s) => {
      const i = armRow.get(s.name);
      return {[statHeaders[0]]: s.cells[0], ...Object.fromEntries(headline.map((k, j) => {
        const x = i === undefined ? null : c.get(k, i);
        return [headlineHeaders[j], typeof x === 'number' && Number.isFinite(x) ? num(x) : ''];
      })), ...Object.fromEntries(statHeaders.slice(1).map((h, j) => [h, s.cells[j + 1]]))};
    }))!.toCsv();

    // the controls against each other: how many measurements fall in their range, and the ones that grow apart
    const inBand = v('control_coverage');
    controlCheck.replaceChildren(!done ? waitingTwin('Checking the controls against each other…') : divV([
      span(inBand === null ? 'Too few controls to check them against each other.' :
        `${pct(inBand)} of the control measurements fall in the range the other controls predict for them ` +
        `(about 90% when the controls grow as the ${modelName(v('line'))} model expects).`),
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
        `randomization, this study's controls and the ${modelName(v('line'))} growth model. Its effect: ` +
        '1 − measured / twin on the evaluation day; a responder\'s 90% range lies above 0. Supporting the primary ' +
        'comparison, not replacing it',
      'Small tumors': `Volumes below ${num(v('v_lod'), 2)} mm³, the smallest a caliper records, are counted at it`,
    };
    const excluded = exclusionsOf(studyInputsOf(call));
    if (excluded)
      rowsOf['Exclusions'] = excluded;
    protocol.replaceChildren(tableFromMap(rowsOf));
    protocolText = [`Analysis protocol (${r.analysis.name})`, ...Object.entries(rowsOf).map(([key, t]) =>
      `${key}: ${t}`)].join('\n');

    const note = groupsNote(r);
    intervalLine.replaceChildren(...note ? [muted(note)] : []);
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
    paint(df, 'group', {[PREDICTED]: '#9aa5b1', ...colorsOf(r)});
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
  const responderText = () => `Responder: on day ${day()} its tumor is below the range it would reach untreated`;
  const largerText = () => `On day ${day()} its tumor is above the range it would reach untreated`;
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
  form.effect(() => responderWords.textContent = result.value ? `On day ${day()} its tumor is below the range it ` +
    'would reach untreated: its 90% range lies right of "Same as twin".' : '');
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
  // an animal's details in the context panel: how it ends against its twin, the reading folded
  // each treated animal against its twin, as the list below last read them
  const twins = new Map<string, AnimalRow>();
  const verdict = (a: AnimalRow) => a.lower > 0 ? 'Responder' : a.upper < 0 ? 'Larger than its twin' :
    'Within its twin\'s range';
  // an animal in the context panel: its curve, readings and marks, then how it ends against its twin
  const showAnimal = (animal: string, color: string) => {
    const n = studyOf(call.getParamValue('study') as DG.DataFrame | null);
    if (!n)
      return;
    const a = twins.get(animal);
    const own = result.peek()?.animals;
    const at = own ? (own.col('animal')!.toList() as unknown[]).findIndex((x) => String(x) === animal) : -1;
    const o = (k: string) => at < 0 || !own!.col(k) || own!.col(k)!.isNone(at) ? null : own!.get(k, at);
    const rate = o('slope_per_day');
    const mine: PanelFact[] = [
      ...rate !== null ? [{label: 'Growth rate', value: `${num(rate)} ln V/day` +
        (o('doubling_time') !== null ? `, doubles in ${num(o('doubling_time'))} d` : '')}] : [],
      ...o('best_response') !== null ? [{label: 'Best response', value: `${signed(o('best_response'))}` +
        (o('mrecist') ? ` (${o('mrecist')})` : '')}] : [],
      ...o('time_to_endpoint') !== null ? [{label: 'Endpoint', value: o('endpoint_reached') ?
        `day ${num(o('time_to_endpoint'))}` : `not reached by day ${num(o('time_to_endpoint'))}`}] : []];
    lastDetails = showDepth(panelOfAnimal(n, animal, {color, showGroup: true, extra: [...mine, ...!a ? [] : [
      {label: 'Against its twin', value: sideOf(a.effect), tip: `90% range: ${sideOf(a.lower)} to ${sideOf(a.upper)}`},
      {label: 'Result', value: verdict(a)}]]}));
    lastGroup = '';
  };
  // a group's animals for its context panel: each one's curve, last volume and status; a click shows the animal
  const MARK_WORDS: {[m: string]: string} = {excluded: 'Excluded', partly: 'Some excluded', flagged: 'Flagged'};
  animalsOf = (group, color) => {
    const n = studyOf(call.getParamValue('study') as DG.DataFrame | null);
    if (!n)
      return [];
    const [g, an] = [n.study.col(n.group)!, n.study.col(n.animal)!];
    const names = [...new Set(Array.from({length: n.study.rowCount}, (_, i) => i)
      .filter((i) => String(g.get(i)) === group).map((i) => String(an.get(i))))]
      .sort((x, y) => x.localeCompare(y, undefined, {numeric: true}));
    if (!names.length)
      return [];
    const marks = animalMarks(n);
    // an animal removed early, and why, as Tolerability lists it ("T3-04 (day 14, weight)"): without body weight
    // there is no such list
    const tol = result.value?.tolerability;
    const removedAs = new Map<string, string>();
    for (let i = 0; tol?.col('removed_animals') && i < tol.rowCount; i++) {
      for (const m of String(tol.get('removed_animals', i) ?? '').matchAll(/(?:^|, )(.+?) \(day (\d+), ([^)]+)\)/g))
        removedAs.set(m[1], `Removed on day ${m[2]} (${m[3]})`);
    }
    const status = (x: string) => [MARK_WORDS[marks.get(x) ?? ''], removedAs.get(x),
      twins.has(x) ? verdict(twins.get(x)!) : '']
      .filter((w) => w && !w.startsWith('Within')).join(' · ');
    return [{title: 'Animals', open: true, content: () => textLines(names.map((x) =>
      animalLine(n, x, color, status(x), () => showAnimal(x, color))))}];
  };
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
      const waiting = waitingTwin('Comparing each animal with what it would measure untreated…');
      waiting.style.gridColumn = '1 / -1';
      animalList.replaceChildren(waiting);
      return;
    }
    const e = r.effects;
    const colors = colorsOf(r);
    const all: AnimalRow[] = Array.from({length: e.rowCount}, (_, i) => ({animal: e.get('animal', i),
      group: e.get('group', i), effect: e.get('effect', i), lower: e.get('effect_lower', i),
      upper: e.get('effect_upper', i)}));
    twins.clear();
    all.forEach((a) => twins.set(a.animal, a));
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
          showAnimal(a.animal, colors[a.group]);
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
  // two lines kept for it: an animal whose sentence takes one line, picked after one whose takes two, does not pull
  // the chart up
  const effectLine = span('');
  Object.assign(effectLine.style, {lineHeight: '1.45', minHeight: '2.9em'});
  form.effect(() => {
    const p = pickedRow.value;
    const alone = untwinned.value;
    effectLine.textContent = p ? `${p.animal}, day ${cell('evaluation_day').value}: ` +
      effectWords(p.effect, p.lower, p.upper) + '.' : alone ? `${alone} has no twin: only a treated animal is ` +
      'compared with what it would measure untreated (a control is the untreated growth itself). Pick a treated ' +
      'one, in the list or the waterfall.' : '';
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
  // the study's animals or measurements as Datagrok tables (to chart, filter, cluster there); the animals picked in
  // Analyze come marked in them
  const openTable = (which: 'animals' | 'measurements') => {
    const r = result.peek();
    if (r && keptStudy)
      explore(which, r, keptStudy, colorsOf(r), self.view, selection.animals.peek());
  };
  const openAnimalsLink = div([hintLink('Open as table', () => openTable('animals'))], 'pmx-noprint');
  const openMeasurementsLink = div([hintLink('Open measurements', () => openTable('measurements'))],
    'pmx-noprint');
  const volumesOnly = only(result, ['volumes'], groupsOf);
  form.effect(() => {
    const r = volumesOnly.value;
    const t = r && animalsPerDay(r.volumes, [r.analysis.get('control', 0), ...r.comparisons.col('arm')!.toList()],
      colorsOf(r));
    perDay.replaceChildren(...t ? [t] : []);
  });

  // each animal as a bar in its group's colour, on what the scientist compares them by: its growth rate (the slope of
  // log volume the test compares), its best response (mRECIST's, % from randomization) or how much smaller it is than
  // its untreated twin; sorted so the better animals end on the right
  const WATERFALL = {
    'Growth rate': {col: 'slope_per_day', title: 'Growth rate (log volume per day)', order: 'desc', control: true,
      none: '<b>No growth rates.</b> Each animal needs two measurements or more.'},
    'Best response': {col: 'best_response', title: 'Best response (%)', order: 'desc', control: true,
      none: '<b>No best response.</b> Each animal needs measurements from day 10 on for its mRECIST class.'},
    'Against its twin': {col: 'effect', title: 'Below its predicted untreated size (%)', order: 'asc', control: false,
      none: '<b>No twins.</b> The twin needs the cell line’s growth and each animal’s volume at ' +
        'randomization.'},
  } as const;
  type WaterfallBy = keyof typeof WATERFALL;
  const waterfallBy = signal<WaterfallBy>('Growth rate');
  const byInput = ui.input.choice('Compare by', {value: 'Growth rate', items: Object.keys(WATERFALL),
    onValueChanged: (v) => waterfallBy.value = (v ?? 'Growth rate') as WaterfallBy,
    tooltipText: 'Growth rate: the slope of log volume each test compares. Best response: mRECIST’s smallest ' +
      'volume from day 10 against randomization. Against its twin: how much smaller than it would be untreated.'});
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
  // why the curves, or the animals removed for other reasons, are not shown (null: they are, or need not be)
  const kmFailed = signal<string | null>(null);
  const competingFailed = signal<string | null>(null);
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
        let failed: string | null = null;
        const k = await survivalOf(s, r.analysis.get('control', 0), r.analysis.get('endpoint', 0), colorsOf(r),
          h === null ? undefined : {horizon: h}).catch((e) => {
          failed = plainError(e, `The survival curves could not be drawn. ${AGAIN}`);
          return undefined;
        });
        if (kmOnly.peek() === r && horizon.peek() === h) {
          kmFailed.value = failed;
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
  // the measured medians for the statistics, whatever day the slider projects to (one more 2 ms call)
  form.effect(() => {
    const r = kmOnly.value;
    const s = keptStudy;
    kmMedians.value = null;
    if (r && s) {
      void survivalOf(s, r.analysis.get('control', 0), r.analysis.get('endpoint', 0), colorsOf(r)).then((k) => {
        if (kmOnly.peek() === r)
          kmMedians.value = k.medians;
      }, () => {
        if (kmOnly.peek() === r)
          kmMedians.value = []; // the panel says why; the statistics give '—'
      });
    }
  });
  // the animals removed for another reason: competing with the endpoint (Fine-Gray), when the study has any
  const competing = signal<Competing[] | null>(null);
  form.effect(() => {
    const r = kmOnly.value;
    const s = keptStudy;
    const removal = call.getParamValue('removalReason');
    competingFailed.value = null;
    if (!r || !s) {
      competing.value = null;
      return;
    }
    void competingOf(s, removal instanceof DG.Column ? removal.name : null, r.analysis.get('control', 0),
      r.analysis.get('endpoint', 0)).then((c) => {
      if (kmOnly.peek() === r)
        competing.value = c;
    }, (e) => {
      competing.value = null;
      if (kmOnly.peek() === r)
        competingFailed.value = plainError(e, `Animals removed for other reasons could not be counted apart. ${AGAIN}`);
    });
  });
  // the curves in the study's colours as they are now (the survival keeps those it was computed with)
  const km = new KmChart(computed(() => {
    const k = survival.value;
    const r = result.value;
    return k && r ? {...k, colors: colorsOf(r)} : k;
  }), horizon);
  const kmWords = div([], 'pmx-km-words');
  form.effect(() => {
    const k = survival.value;
    const p = new Map((k?.logrank ?? []).map((t) => [t.group, t.p_value]));
    // projected: no test (the tests are the measured study's), every figure said to be projected
    const on = k?.horizon !== undefined;
    kmWords.replaceChildren(...[kmFailed.value, competingFailed.value].filter((x): x is string => !!x)
      .map((x) => noteLine(x)), ...(on ? [noteLine(`Projected beyond day ${k!.lastDay} to day ${k!.horizon}: each ` +
      'animal still under the endpoint carried on at the rate of its last three measurements. Not a result: the ' +
      'tests and the report\'s conclusions are the measured study\'s.', 'warning')] : []),
    // no animal at the endpoint: one line says so, not a "median not reached, 0 of n" per group
    ...k?.medians.length && k.medians.every((m) => m.events === 0) ? [noteLine(`No animal ${on ? 'would reach' :
      'reached'} the endpoint${on ? ` by day ${k.horizon}` : ''}: no group has a median time to it.`)] :
      (k?.medians ?? []).map((m) => noteLine(`${m.group}: ${on ? 'projected ' : ''}median ${m.median === null ?
      'not reached' : dayText(m.median)}, ${m.events} of ${m.n} ${on ? 'would reach' : 'reached'} the endpoint` +
      (!on && p.has(m.group) ? `; log-rank against the control ${pValue(p.get(m.group))}` : ''))),
    ...(!on && competing.value?.length ? [noteLine('Counting the animals removed for other reasons apart: ' +
      competing.value.map((c) => c.none ? `${c.group}: not estimable, no ${c.none === 'group' ? 'animal of the ' +
        'group' : 'control'} reached the endpoint, ${c.competing} removed` :
        `${c.group} has ${c.shr.toFixed(2)} times the control's rate of reaching the endpoint (95% range ` +
      `${c.lower.toFixed(2)}–${c.upper.toFixed(2)}), p ${pValue(c.p)}, ${c.competing} removed`)
        .join('; '), 'success')] : []));
  });

  // short enough for the value axis' title; the unit, ln V, in the section's hint
  const RATE = 'Growth rate per day';
  const slopesData = computed(() => {
    const r = animalsOnly.value;
    if (!r)
      return undefined;
    const a = r.animals;
    const df = DG.DataFrame.fromObjects(Array.from({length: a.rowCount}, (_, i) => ({'Group': a.get('group', i),
      'Animal': String(a.get('animal', i)), [RATE]: a.get('slope_per_day', i)}))
      .filter((x) => x[RATE] !== null))!;
    paint(df, 'Group', colorsOf(r));
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
  // a figure, not a tool: each box named under it, level (the platform shortens a long name to its box; the box's
  // tooltip names it whole), its legend for the platform's colour picker (right click), no column or colour selectors
  const slopes = viewers.boxPlot(slopesData, {valueColumnName: RATE,
    categoryColumnNames: ['Group'], showPValue: false, showStatistics: false, showCategoryAxis: true,
    labelOrientation: 'Horz', markerColorColumnName: 'Group', legendVisibility: 'Always', showCategorySelector: false,
    showValueSelector: false, showColorSelector: false, showSizeSelector: false});
  // the growth rates the tests compare, folded under the results table until asked; laid out hidden, the chart
  // follows its width once it opens
  const rates = new Section({title: 'Growth rates', hint: 'Each animal’s slope of log volume (ln V per day)',
    expanded: false});
  rates.body.append(sized(slopes, 380));
  self.own(() => rates.dispose());
  rates.effect(() => {
    if (rates.expanded.value)
      window.dispatchEvent(new Event('resize'));
  });

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
      content: [readoutLine, groupsTable, intervalLine, hintLine, groupNotes, ladderNotes, rates.root, exposureNotes,
        exposureBox]},
    {id: 'growth', title: 'Observed tumor growth', hint: 'Median and range (smallest to largest tumor) ' +
      'of each group, to the evaluation day; in grey, the median control predicted (by the plan, else by the ' +
      'atlas for this model)',
    content: [growthAxis.row, openMeasurementsLink, sized(growth, 380), perDay]},
    {id: 'waterfall', title: 'Each animal compared', hint: 'Every animal as a bar in its group\u2019s colour, on its ' +
      'growth rate, its best response or its effect against its twin; the better animals on the right',
    content: [byInput.root, waterfallMissing, sized(waterfall, 300), openAnimalsLink]},
    {id: 'spider', title: 'Each tumor over time', hint: 'Every animal’s volume as a % change from ' +
      'randomization, in its group’s colour; click a group in the legend to hide or show its animals',
    content: [sized(spider, 300)]},
    {id: 'km', title: 'Time to endpoint', hint: 'Kaplan-Meier: the animals of each group still below the ' +
      'endpoint volume, day by day, as measured; an animal that never reached it counts until its last measurement',
    content: [sized(km, 320), kmWords]},
    {id: 'animals', title: 'Each animal against its twin', hint: 'Its twin: what each treated tumor would have ' +
      'measured untreated, predicted from its size at randomization, the study\'s controls and the cell line\'s ' +
      'growth. Select an animal to see its growth', content: [animals]},
    {id: 'tolerability', title: 'Tolerability', hint: 'Whether the treatments were well tolerated',
      content: tolerability.content},
    {id: 'plan', title: 'Plan vs observed', hint: 'Whether the study went as it was planned',
      content: versus.content},
    {id: 'controls', title: 'Control check', hint: 'Each control against the range the other controls predict',
      content: [controlCheck]},
  ], ANALYZE_TABS, 'pmx-analyze-panels', 'analyze');
  exports.push(['Layout…', () => editLayout(panels.root), 'View']);
  // the study's colours: a colour changed in a chart is the study's, in every stage open on it
  const studyDf = () => call.getParamValue('study') as DG.DataFrame | null;
  const link: StudyLink = {meta: () => metaOf(studyDf()), setMeta: (m) => {
    const df = studyDf();
    if (df)
      setMetaOn(df, m);
  }, colors: () => {
    const r = result.peek();
    return r ? colorsOf(r) : {};
  }, table: studyDf, plan: () => planKeyOf(studyDf()),
  control: () => result.peek()?.analysis.get('control', 0) ?? null};
  self.own(linkStudy(link));
  form.effect(() => {
    const r = result.value;
    if (r)
      follow(link, Object.values(r).filter((x): x is DG.DataFrame => x instanceof DG.DataFrame));
  });
  exports.push(['Reset colors', () => resetColors(link), 'View',
    () => hasOwnColors(link) ? null : 'The colors are the automatic ones']);
  // a model planned from its history: its twin is not validated for it (D694)
  const twinLine = markLine(TWIN_NOTE);
  form.effect(() => {
    void result.value;
    const df = studyDf();
    twinLine.hidden = !savedHistory(String(call.getParamValue('line') || (df ? lineOf(df) : '') || ''));
  });
  // day 0's date against the table, as Monitor's calendar checks it: a date after the measurements is not printed
  const lastMeasured = (): number | null => {
    const df = studyDf();
    try {
      const max = df ? roleColumns(df).day.stats.max : NaN;
      return Number.isFinite(max) ? max : null;
    } catch {
      return null; // no day column picked yet
    }
  };
  const dateLine = markLine('');
  form.effect(() => {
    void result.value;
    const wrong = startMismatch(metaOf(studyDf())?.start ?? null, lastMeasured());
    dateLine.hidden = !wrong;
    dateLine.lastElementChild!.textContent = wrong ? `${wrong} Set day 0's date in the study's details (click the ` +
      'study\'s name).' : '';
  });
  const content = divV([dataLine.root, twinLine, dateLine, cards, panels.root], 'pmx-analyze-content');
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
      const section = sectionNow();
      printSheet(await panels.all(() => reportHtml({title: r.analysis.name, stage: 'Analyze', meta: metaOf(df),
        groups: [r.analysis.get('control', 0), ...r.comparisons.col('arm')!.toList() as string[]],
        // the plan the analysis used: the one picked in the form, else the study's (analyze.ts)
        content, log: logOf(df), planAtlas: planAtlas(call.getParamValue('plan')?.trim() || df?.getTag(PLAN_TAG)),
        annex: [['Results', summary], ...section ? [['Exposure and patients', sectionElement(section)] as
          [string, HTMLElement]] : [], ['Statistics', statsTable], ['Analysis protocol', protocol]],
        lastDay: lastMeasured()})));
    }],
    // the Essential 10, or ARRIVE 2.0 whole (the Recommended Set added); Download (.docx) is the author checklist
    ['ARRIVE checklist…', function arrive(full = false) {
      const {df, r} = studyNow();
      const c = r.comparisons;
      const plan = r.analysis.get('plan', 0);
      const items = arriveItems({control: r.analysis.get('control', 0), arms: c.col('arm')!.toList() as string[],
        plan: plan ? String(plan) : null, exclusions: exclusionsOf(studyOf(df)), meta: metaOf(df),
        // randomised: the app allocated the groups (the log says so), or the study was followed from its plan, which
        // balances them on day 0's volume (as Plan's checklist says)
        randomized: !!df?.getTag(PLAN_TAG) || logOf(df).some((e) => /randomi[sz]|allocat/i.test(e.what)),
        endpoint: r.analysis.get('endpoint', 0), evaluationDay: r.analysis.get('evaluation_day', 0),
        protocol: protocolText, alpha: r.analysis.get('alpha', 0), correction: r.analysis.get('correction', 0),
        analyzed: [{group: r.analysis.get('control', 0), n: r.analysis.get('n_control', 0)},
          ...Array.from({length: c.rowCount}, (_, i) => ({group: c.get('arm', i), n: c.get('n', i)}))]}, full);
      const which = ui.input.choice('Checklist', {value: full ? 'ARRIVE 2.0 (21 items)' : 'Essential 10',
        items: ['Essential 10', 'ARRIVE 2.0 (21 items)'], nullable: false, onValueChanged: (x) => {
          d.close();
          arrive(x !== 'Essential 10');
        }});
      const title = full ? 'ARRIVE 2.0 checklist' : 'ARRIVE Essential 10 checklist';
      const sub = [metaOf(df)?.id, `From the analysis of ${df.name}`,
        `prepared ${new Date().toISOString().slice(0, 10)}`].filter((x) => x).join(' · ');
      const d = arriveDialog(items, {full, details: () => head.details(() => arrive(full)), options: which.root,
        download: async () => download(`${df.name} - ${title}.docx`, await docx(checklistBlocks({title, sub, items,
          credit: ARRIVE_CREDIT, last: ['Section, line', 'where']}), {title, landscape: true}), DOCX_TYPE)});
    }],
    // the study as a row of the nonclinical dossier (CTD), each group a sub-row; its checklist is the preview
    ['Regulatory (CTD) table…', function ctd() {
      const {df, r} = studyNow();
      const s = ctdStudyOf(r, metaOf(df), r.analysis.name, shownNca, stepsInWords(ladder.peek()));
      const sub = [s.id, `From the analysis of ${df.name}`, `prepared ${new Date().toISOString().slice(0, 10)}`]
        .filter((x) => x).join(' · ');
      checklistDialog({title: 'Regulatory (CTD) table', guide: 'ICH M4S(R2)', url: M4S_URL, asker: 'The table',
        download: async () => download(`${df.name} - CTD table.docx`,
          await docx(ctdBlocks([s], {sub}), {title: CTD_TITLE, landscape: true}), DOCX_TYPE)},
      ctdItems(s), {details: () => head.details(ctd)});
    }],
    ['Animals', () => openTable('animals'), 'Data > Open as table'],
    ['Measurements', () => openTable('measurements'), 'Data > Open as table'],
    ['Statistics', () => download(`${studyNow().df.name} - statistics.csv`, statsCsv, 'text/csv'), 'Export > CSV'],
    ['Measurements', () => download(`${studyNow().df.name}.csv`, studyNow().df.toCsv(), 'text/csv'),
      'Export > CSV'],
    ['Per-animal results', () => download(`${studyNow().df.name} - animals.csv`, studyNow().r.animals.toCsv(),
      'text/csv'), 'Export > CSV'],
    ['Modelling dataset (NONMEM)…', async () => {
      const {df} = studyNow();
      const s = studyInputsOf(call);
      if (!s)
        throw new Error('Analyze a study first');
      // no table picked since the study was opened: the one NCA read when it was kept, if it is open
      const recorded = savedNca.peek()?.table || null;
      let pk = ncaSource;
      const t = !pk && recorded ? grok.shell.tables.find((x) => x.name === recorded) : undefined;
      const m = t ? await concentrationMapping(t) : null;
      if (t && m)
        pk = {table: t, mapping: m, groupCol: groupColumnOf(t, m)};
      modellingDatasetDialog(s, metaOf(df), pk, recorded);
    }],
    ['Copy the results summary', () => copied(summaryText, 'Results summary')],
    ['Copy the analysis protocol', () => copied(protocolText, 'Analysis protocol')]);
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
      status: `${modelName(r.analysis.get('line', 0))} · evaluated on ${dayText(r.analysis.get('evaluation_day', 0))}`,
      meta: () => metaOf(df), setMeta: df ? (m) => setMetaOn(df, m) : undefined,
      groupMenu: (g, e) => studyMenu(studyOf(df), 'group', g, {causedBy: e})};
    // Plan opens with the plan the analysis uses: the one picked, else the study's own
  }), (s) => void openStage(s, call.getParamValue('plan')?.trim() ||
    (call.getParamValue('study') as DG.DataFrame | null)?.getTag(PLAN_TAG)));
  // Save (and Ctrl+S): what Analyze changes on the study's table (the patients' exposure, the effect kept, exclusions).
  // A table followed in Monitor is saved by that view; any other (from a file, a drop), with rows on, is saved here as
  // its row: the one it came from (its limits kept), else a new one, the first time it is saved
  const savedAs = new WeakMap<DG.DataFrame, string>();
  const snapshot = (df: DG.DataFrame) =>
    [df.toCsv(), df.getTag(EXPOSURE_TAG) ?? '', df.getTag('pmx.meta') ?? ''].join('|');
  const saveHere = async (df: DG.DataFrame): Promise<void> => {
    const at = rowOfTable(df);
    const was = at ? (await open('study', at.id)).doc : null;
    const m = await withCode(metaOf(df) ?? emptyMeta());
    setMetaOn(df, m);
    const s = await stamped(savedStudyOf(df, was?.humane ?? null, was?.weightLimit ?? null));
    try {
      const saved = await save('study', s, at, {fileOwner: tableFileOwner(df)});
      setRowOn(df, saved);
      markTableFromFile(df, '');
      if (saved.meta) {
        setMetaOn(df, saved.meta);
        grok.shell.info(newRowNote('study', saved, m.id));
      }
      // (a first save gives the study its ID, from withCode above: the header shows it)
      head.refresh();
    } catch (e) {
      if (!at || !isReadOnly(e))
        throw e;
      return offerCopy(async () => {
        const {saved, doc} = await saveCopy('study', s, at);
        setRowOn(df, saved);
        df.name = doc.name;
        setMetaOn(df, doc.meta!);
        head.refresh();
        savedAs.set(df, snapshot(df));
        grok.shell.info(`Saved as your copy: ${esc(doc.name)}`); // (a toast reads HTML: the name is never markup)
      });
    }
    savedAs.set(df, snapshot(df));
    grok.shell.info('Study saved');
  };
  const saver = (): (() => Promise<void>) | null => {
    const df = studyDf();
    if (!df)
      return null;
    // (Monitor's save may give the study a new ID: the header here shows it too)
    const viaMonitor = tableSavers.get(df);
    const monitor = viaMonitor && (() => {
      finished(df);
      return viaMonitor().then(() => head.refresh());
    });
    return monitor ?? (domainStorageOn() ? () => {
      finished(df);
      return saveHere(df).catch((e) => {
        grok.shell.error(saveError(e, 'study'));
        throw e;
      });
    } : null);
  };
  // saved from Analyze, the study is finished (the Overview lists it so, whichever view saves it next), with the
  // model it was analyzed with: a table opened straight here has none of its own, and its row would list no model
  const finished = (df: DG.DataFrame) => {
    df.setTag(FINISHED_TAG, 'true');
    const line = String(call.getParamValue('line') || '');
    if (line && !lineOf(df)) {
      df.setTag(LINE_TAG, line);
      setHistoryOn(df, line);
    }
  };
  const saveButton = buttonWithIcon('Save study', 'save', () => void saver()?.().catch(() => {}));
  saveButton.title = 'Save the study with what Analyze added (the patients’ exposure, exclusions)';
  effect(() => {
    void result.value;
    const df = studyDf();
    saveButton.hidden = !saver();
    if (df && !savedAs.has(df))
      savedAs.set(df, snapshot(df));
  });
  // (closing asks only for a table saved here: one followed in Monitor is asked about there)
  queueMicrotask(() => saveGate(self, 'study', () => saver()?.() ?? Promise.resolve(), () => {
    const df = studyDf();
    return !!df && domainStorageOn() && !tableSavers.has(df) && savedAs.has(df) && savedAs.get(df) !== snapshot(df);
  }, () => {
    const df = studyDf();
    return !!df && !rowOfTable(df);
  }));
  return {
    root: stagePage(head.row, panes),
    ribbon: [[head.steps], [saveButton]],
    path: '/analyze',
    status: computed(() => busy.value ? 'Analyzing the study…' : result.value ?
      (call.getParamValue('study') as DG.DataFrame | null)?.name || result.value.analysis.name : ''),
  };
}
