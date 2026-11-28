/* The Analyze view, built as the Plan view is: the declaration form generated from pmxAnalyzeStudy's own signature
   (u2 funcForm: the study table, its columns and the analysis settings), the results from the function's outputs, in
   the order a scientist asks: did each treatment work (the figures, then the groups with their result), how did the
   tumors grow (band-chart.ts), which animals responded (each treated animal against its untreated twin, the
   animals listed by effect, the one selected drawn against its twin), then, folded, what serves a closer look: the
   statistics behind each result, the control check and the analysis protocol; the study against its plan has a tab
   of its own (plan-vs-observed.ts). */
import {maxOf} from './format';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {badge, BasicTable, buttonWithIcon, computed, div, divH, divV, dot, effect, loader, signal, span,
  tableFromMap, Tooltip} from '@datagrok-libraries/u2';
import {AnimalPanels, AnimalView, growthFrame, ordinal, placeWords} from './animal-panel';
import {holdTable, leaveAsks, offerCopy, saveGate, saveStudyAs, savedTables, studyPickList, tableSavers}
  from './native';
import {applyKept, keptOfCall, recallForm, rememberForm} from './form-memory';
import {domainStorageOn, isReadOnly, markTableFromFile, newRowNote, open, rowOfTable, save, saveCopy, saveError,
  setRowOn, tableFileOwner, withCode} from './store';
import {FINISHED_TAG, LIMITS_TAG, LINE_TAG, LOOKS_TAG, STUDY_TAG, setAnalysisOn, setHistoryOn, studyOf as savedStudyOf,
  stamped}
  from './study-file';
import {roleOf} from './study-reader';
import {funcForm, TableInput, ViewHost, viewers} from '@datagrok-libraries/u2/dg';
import type {ViewContent} from '@datagrok-libraries/u2/dg';

import {AnalysisFrames, AnalysisInputs, PREDICTED, analysisDeclarationOf, analysisFrames, doublingOf, inMm3,
  GROWTH_RATE_TIP, TOO_FEW, comparisonsMade, tooFewWhy, growthRateWords, growthVsControl, intervalWords, oneDecimal,
  pctPerDay, warmTwin}
  from './analyze';
import {CORRECTION_NAME, correctionOf, correctionOver, isFiniteNumber, median} from './engine';
import {EmptyState, NoControl} from './crux';
import {PREDICTED_COLOR, StudyLink, follow, groupColors, hasOwnColors, linkStudy, paint, resetColors} from './colors';
import {CorrectionInput} from './arms-input';
import {askedOf, autoRun, axisAt, callInputs, cardGrid, dayText, keepSame, onTableEdits, only, TABLE_EDITS, guide,
  exampleLine, infoNote, noteLine, notesBox, num, pct, plainError, pValue, settleForm, signed, isBodyWeight,
  keepNamedColumns,
  sized, markStage,
  statCard, swatch, volumeAxis, wrapCells, stagePanes, whenSettled, hintLink, PanelFact, PanelSection, PanelVerdict,
  showDepth, textLines, factList, docSubtitle, tgiText, unadjustedP} from './ui';
import {analysisFieldsOf, groupsShown, lineOf, logOf, planAtlas, restoreFields, studyNameOf} from './study-file';
import {TWIN_NOTE, savedHistory} from './history-models';
import {markLine} from './ui';
import {icsOf, startMismatch} from './calendar';
import {exportName} from './study-file';
import {PLAN_TAG, planKeyOf, planOfStudy} from './plan-file';
import {canReturn, keepOnLeaving, linkRow, openStage} from './stages';
import {noEffectNote, planInput, planTab} from './plan-vs-observed';
import {nciReadings, NciReading, tolerabilityTab} from './tolerability';
import {BandChart} from './band-chart';
import {explore, StudyInputs, studyInputsOf} from './explore';
import {customizePage, panelTabs, Tab} from './panels';
import {stagePage, studyHeader} from './study-header';
import {compareWith, openStudyDialog, stageLinks, takesFiles} from './hub-view';
import {AnimalSelection, exclusionsOf, mirrorSelection, studyDataLine, studyOf, tagStudy,
  keepStudyField, withoutExcluded, animalLine, animalMarks, animalReadings, gain} from './study-table';
import {META_TAG, armParts, emptyMeta, metaOf, parseLine, setMetaOn} from './study-meta';
import type {StudyMeta} from './study-meta';
import {ExposureSection, GroupDetails, detailsElement, exposureSection, groupDetails, sectionElement}
  from './translation-report';
import {Exports, download, htmlTable, reasonOf, reportHtml, stageMenu, stageToolbox} from './report';
import {esc, printSheet} from './sheet';
import {SpiderChart, WaterfallChart, spiderLines} from './spider-chart';
import {LineInput} from './line-input';
import {Competing, Survival, competingOf, survivalOf} from './survival';
import {KmChart} from './km-chart';
import {CombinationStep, combinationSteps, notCombined, stepWords, synergyAlpha, synergyDialog}
  from './synergy';
import {findCombinations, readProblem} from './arm-name';
import {PmxGroup, bindObject, popupFor} from './handlers';

/** The exposure panel's hint, for the axis it shows. */
const axisHint = (x: 'AUClast' | 'Dose') => `Each group’s TGI against its ${x === 'Dose' ? 'dose' :
  'mouse exposure (NCA’s AUClast)'}, and the effect expected between the doses tested; shown when one drug was ` +
  'given at 3 doses or more';

/** The statistics file's rows (Export > CSV > Statistics): the control, then each treated group with its test on the
 *  slopes (the one interval, at `slope_difference_confidence`) and its TGI and T/C as the groups table gives them;
 *  `kmDay`: a group's measured median day to the endpoint (null: none). */
export function statisticsRows(r: AnalysisFrames, kmDay: (group: string) => number | null): {[k: string]: unknown}[] {
  const [c, v] = [r.comparisons, (k: string) => r.analysis.get(k, 0)];
  const [control, done] = [v('control') as string, v('twinned')];
  const fin = (x: unknown) => isFiniteNumber(x) ? x : null;
  return [{group: control, role: 'control', evaluation_day: v('evaluation_day'),
    n_test: v('n_control'), n_evaluation_day: v('n_control_day'),
    growth_rate_ln_per_day: fin(c.get('slope_control_per_day', 0)),
    doubling_time_days: doublingOf(c.get('slope_control_per_day', 0)),
    median_days_to_endpoint: kmDay(control)},
  ...Array.from({length: c.rowCount}, (_, i) => {
    const g = (k: string) => fin(c.get(k, i));
    const arm = c.get('arm', i) as string;
    const [t, k] = [kmDay(arm), kmDay(control)];
    return {group: arm, role: 'treated', evaluation_day: v('evaluation_day'), n_test: g('n'),
      n_evaluation_day: g('n_day'), growth_rate_ln_per_day: g('slope_per_day'),
      doubling_time_days: doublingOf(c.get('slope_per_day', i)), slope_difference: g('slope_difference'),
      slope_difference_lower: g('slope_difference_lower'), slope_difference_upper: g('slope_difference_upper'),
      slope_difference_confidence: g('slope_difference_confidence'),
      // the same difference as the groups table shows it, in %/day (round 3, analyst W9)
      ...pctOf(c, i), p_value: c.get('activity', i) === TOO_FEW ? null : g('p_value'),
      p_raw: c.get('activity', i) === TOO_FEW ? null : g('p_raw'),
      p_adjustment: v('correction'), tgi_pct: g('tgi'), tc_pct: g('tc'),
      median_days_to_endpoint: t, days_later_than_control: t !== null && k !== null ? t - k : null,
      result: c.get('activity', i), twin_animals: done ? g('twin_animals') : null,
      responders: done ? g('responders') : null};
  })];
}

/** The statistics file's growth against the control in %/day (analyze.ts pctPerDay), as the groups table shows it;
 *  nulls: no test. */
function pctOf(c: DG.DataFrame, i: number): {[k: string]: number | null} {
  const tested = c.get('activity', i) !== TOO_FEW;
  const at = (k: string) => {
    const x = c.get(k, i);
    return tested && isFiniteNumber(x) ? Math.round(10 * pctPerDay(x)) / 10 : null;
  };
  return {growth_vs_control_pct_per_day: at('slope_difference'),
    growth_vs_control_pct_per_day_lower: at('slope_difference_lower'),
    growth_vs_control_pct_per_day_upper: at('slope_difference_upper')};
}

/** A group on the exposure chart: its TGI as the groups table gives it, its exposure and colour. */
export interface ExposurePoint {arm: string, tgi: number,
  auc: number, aucUnit: string, dose: string, drug: string, route: string, color: string}

/** The treated groups with NCA's exposure, at the groups table's TGI (none without a TGI). */
export function exposurePointsOf(r: AnalysisFrames, exposure: Map<string, ArmExposure>,
  colors: {[group: string]: string}): ExposurePoint[] {
  const c = r.comparisons;
  const ci = (k: string, i: number) => Number.isFinite(c.get(k, i)) ? c.get(k, i) as number : null;
  return (c.col('arm')!.toList() as string[]).flatMap((arm, i) => {
    const x = exposure.get(arm)?.exposure;
    const tgi = ci('tgi', i);
    return x && tgi !== null ? [{arm, tgi, color: colors[arm], ...x}] : [];
  });
}

/** The exposure chart's series of the curve expected between the doses tested (dashed, with its band). */
export const CURVE = 'Expected between the doses tested';

/** The exposure chart's rows: each group at its TGI (no interval: one day's description), on AUClast or on the dose
 *  (exposureAxis), and the curve expected between the doses tested when there is one (exposureCurve, on the same
 *  TGI); undefined: no group. */
export function exposureFrame(p: ExposurePoint[], meta: StudyMeta | null,
  curve: ExposureCurve['curve'] = []): DG.DataFrame | undefined {
  if (!p.length)
    return undefined;
  const groups = p.map((x) => ({arm: x.arm, auc: x.auc}));
  const axis = exposureAxis(groups, meta);
  const xName = axis.axis === 'dose' ? 'Dose' : 'AUClast';
  const df = DG.DataFrame.fromObjects([...p.map((x, i) => ({series: x.arm, [xName]: axis.of(groups[i]),
    lower: x.tgi, upper: x.tgi, tgi: x.tgi})),
  ...curve.map((c) => ({series: CURVE, [xName]: c.x, lower: c.lo, upper: c.hi, tgi: c.tgi}))])!;
  df.col(xName)!.meta.units = axis.axis === 'dose' ? axis.unit : p[0].aucUnit;
  df.col('tgi')!.meta.units = '%';
  df.col('series')!.meta.colors.setCategorical({[CURVE]: '#9aa5b1',
    ...Object.fromEntries(p.map((x) => [x.arm, x.color]))});
  return df;
}

/** Each combination's step on the ladder, in the words shown under it, for the regulatory table. */
const stepsInWords = (steps: Map<string, {step: CombinationStep}>) =>
  new Map([...steps].map(([ab, x]) => [ab, stepWords(x.step)]));

/** What a study keeps of its analysis `r` (Compare reads it): each treated group's TGI on the evaluation day, the
 *  regulatory table's study, and the control's growth rate with its predicted range. The demo studies Compare saves
 *  keep it the same way. */
export function keptEffect(r: AnalysisFrames, meta: StudyMeta | null, nca: NcaByArm | null,
  steps = new Map<string, {step: CombinationStep}>()): NonNullable<SavedExposure['effect']> {
  const v = (c: string) => r.analysis.get(c, 0);
  const arms = r.comparisons.col('arm')!.toList() as string[];
  return {day: v('evaluation_day') ?? null, tgi: Object.fromEntries(arms
    .map((a, i) => [a, r.comparisons.get('tgi', i)]).filter(([, t]) => Number.isFinite(t))),
  ctd: ctdStudyOf(r, meta, r.analysis.name, nca, stepsInWords(steps)),
  control: {slope: v('control_slope'), q025: v('control_slope_q025'), median: v('control_slope_median'),
    q975: v('control_slope_q975'), endpoint: v('endpoint')}};
}
import {ARRIVE_CREDIT, arriveDialog, arriveItems, checklistBlocks, checklistDialog} from './arrive';
import {CTD_TITLE, M4S_URL, ctdBlocks, ctdItems, ctdStudyOf} from './ctd';
import {DOCX_TYPE, docx} from './docx';
import {ArmExposure, CURVE_BAND, DrugExposure, ExposureCurve, HumanExposure, NcaByArm, ReferenceRow,
  SavedExposure,
  aucUsed, concentrationMapping, cvOk, drugMultiples,
  exposureAxis, exposureCurve, exposureHint, exposureJson, groupColumnOf, humanMultiple,
  ncaByArm,
  ncaInstalled, parseExposure, readout, referenceRow, runNca, shownMultiple}
  from './exposure';
import {EXPOSURE_TAG, exposureOfTag, roleColumns} from './study-file';
import {humanReference, linePrior, modelName} from './priors';
import {humanExposureDialog, humanExposureSummary} from './human-exposure';
import {PkSource, modellingDatasetDialog} from './modelling-dataset';
import {mm3Text} from './monitor';

/** Rows as CSV: numbers as numbers to `digits` significant figures (float noise off), empty for none, text quoted
 *  where it must be. */
export function csvText(rows: {[k: string]: unknown}[], digits = 6): string {
  const keys = [...new Set(rows.flatMap((r) => Object.keys(r)))];
  const cell = (x: unknown) => {
    if (x === null || x === undefined || (typeof x === 'number' && !Number.isFinite(x)))
      return '';
    const t = typeof x === 'number' ? String(Number(x.toPrecision(digits))) : String(x);
    return /[",\n\r]/.test(t) ? `"${t.replace(/"/g, '""')}"` : t;
  };
  return [keys.map(cell).join(','), ...rows.map((r) => keys.map((k) => cell(r[k])).join(','))].join('\n') + '\n';
}

/** Whether `df` reads as a study: tagged as one by a stage, or with a volume or caliper column by its name (a
 *  concentrations table opened for NCA is none). */
export const readsAsStudy = (df: DG.DataFrame) => !!df.getTag(STUDY_TAG) ||
  [...Array(df.columns.length).keys()].some((j) =>
    ['volume', 'length'].includes(roleOf(df.columns.byIndex(j).name)?.role ?? ''));
export const NOT_A_STUDY = 'Not a study table: no animal, day or volume column. Pick the study’s measurements, or ' +
  'open them with File > Open (Ctrl+O).';

/** A table as CSV with its values as stored, not as its columns display them (the platform's toCsv writes the
 *  column's format: 19.26 as 19.3, 83 as 83.0): a 32-bit column to the 7 figures it holds, else to `digits`. */
export function tableCsv(df: DG.DataFrame, digits = 15): string {
  const cols = df.columns.toList();
  const fig = cols.map((c) => c.getRawData() instanceof Float32Array ? Math.min(7, digits) : digits);
  const rows = Array.from({length: df.rowCount}, (_, i) => Object.fromEntries(cols.map((c, j) => [c.name,
    c.isNone(i) ? null : c.type === DG.TYPE.DATE_TIME ? c.getString(i) : typeof c.get(i) === 'number' ?
      Number(c.get(i).toPrecision(fig[j])) : c.get(i)])));
  return rows.length ? csvText(rows, 15) : cols.map((c) => c.name).join(',') + '\n';
}

/** Analyze's tabs and their panels, in order and width (twelfths): the answer and how the groups grew; each animal;
 *  how the animals bore it; whether the study went as planned. */
export const ANALYZE_TABS: Tab[] = [
  {name: 'Results', layout: [{id: 'groups', span: 12}, {id: 'growth', span: 12},
    {id: 'km', span: 12, folded: true}, {id: 'exposure', span: 12, folded: true}, {id: 'rates', span: 12, folded: true},
    {id: 'notes', span: 12, folded: true}]},
  {name: 'Each animal', layout: [{id: 'waterfall', span: 12}, {id: 'spider', span: 12, folded: true},
    {id: 'animals', span: 12, folded: true}]},
  {name: 'Tolerability', layout: [{id: 'tolerability', span: 12}]},
  {name: 'Against the plan', layout: [{id: 'plan', span: 12}, {id: 'controls', span: 12}]}];

/** How many animals each group had on each measurement day (the volumes' `n`): the days, and each group's counts. */
export interface AnimalsByDay {days: number[], groups: string[], n: Map<string, Map<number, number>>}

export function animalsByDay(volumes: DG.DataFrame, groups: string[]): AnimalsByDay {
  const days = [...new Set(volumes.col('day')!.toList() as number[])].sort((a, b) => a - b);
  const n = new Map(groups.map((g) => [g, new Map<number, number>()]));
  for (let i = 0; i < volumes.rowCount; i++)
    n.get(volumes.get('group', i))?.set(volumes.get('day', i), volumes.get('n', i));
  return {days, groups, n};
}

/** The one line under the growth chart when a group thins ("Animals by day: Vehicle lost 3, …"); null: none lost
 *  (the counts by day are in each group's details, the report and Open as table). */
export function lostLine(a: AnimalsByDay): string | null {
  const lost = a.groups.map((g) => {
    const counts = a.days.map((d) => a.n.get(g)?.get(d)).filter((x): x is number => typeof x === 'number');
    // left by the end, not the max minus the min: a group measured from a later day has not lost any
    return {g, n: counts.length ? maxOf(counts) - counts[counts.length - 1] : 0};
  }).filter((x) => x.n > 0);
  return lost.length ? `Animals by day: ${lost.map((x) => `${x.g} lost ${x.n}`).join(', ')}` : null;
}

/** The counts as a table: a row per group, a column per day (the report's, and Open as table's). */
export const animalsByDayFrame = (a: AnimalsByDay) => DG.DataFrame.fromObjects(a.groups.map((g) => ({Group: g,
  ...Object.fromEntries(a.days.map((d) => [`Day ${d}`, a.n.get(g)?.get(d) ?? null]))})))!;

/** The same for the report. */
function byDayTable(a: AnimalsByDay): HTMLElement {
  const el = div([]);
  el.innerHTML = htmlTable(['Group', ...a.days.map((d) => `Day ${d}`)],
    a.groups.map((g) => [g, ...a.days.map((d) => a.n.get(g)?.get(d) ?? '')]));
  return el;
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
  inMm3(df, ['lower', 'upper', 'median']);
  return df;
}

/** Each animal removed early, by name: its day and why, as Tolerability lists them ("T3-04 (day 14, weight)");
 *  without body weight there is no such list. */
export function removedAnimals(tol: DG.DataFrame | undefined): Map<string, {day: number, why: string}> {
  const out = new Map<string, {day: number, why: string}>();
  for (let i = 0; tol?.col('removed_animals') && i < tol.rowCount; i++) {
    for (const m of String(tol.get('removed_animals', i) ?? '').matchAll(/(?:^|, )(.+?) \(day (\d+), ([^)]+)\)/g))
      out.set(m[1], {day: Number(m[2]), why: m[3]});
  }
  return out;
}

/** A group's animals removed early, in a few words ("2 removed on day 7 (weight)"; on several days or for
 *  several reasons: "3 removed early"), from Tolerability's list; null: none, or no body weight to list them. */
export function removedOf(tol: DG.DataFrame | undefined, group: string): string | null {
  const i = tol?.col('removed_animals') && tol.col('group') ?
    (tol.col('group')!.toList() as unknown[]).findIndex((x) => String(x) === group) : -1;
  if (i < 0)
    return null;
  const all = [...String(tol!.get('removed_animals', i) ?? '').matchAll(/(?:^|, )(.+?) \(day (\d+), ([^)]+)\)/g)];
  if (!all.length)
    return null;
  const same = all.every((m) => m[2] === all[0][2] && m[3] === all[0][3]);
  return `${all.length} removed ${same ? `on day ${all[0][2]} (${all[0][3]})` : 'early'}`;
}

/** What Analyze says where the twins and the atlas's checks would be, with no tumor model picked. */
export const NO_MODEL = 'Pick the tumor model to see each animal against its twin';
const noModelNote = (cell = false) => {
  // (a table's cell: a dash, the words on hover)
  const el = span(cell ? '—' : `${NO_MODEL}.`);
  if (cell)
    Tooltip.bind(el, NO_MODEL);
  else
    el.style.color = 'var(--dg-text-color-muted, #68727f)';
  return el;
};

/** An effect against the twin (%), in words: how much smaller, or larger. */
/** A percent whole, with a decimal under 1 (a responder's lower end of 0.3% read "0%": round 3, analyst W6). */
export const effectPct = (x: number) => {
  const v = Math.abs(x);
  return v > 0 && v < 1 ? v.toFixed(1) : String(Math.round(v));
};
const sideOf = (x: number) => `${effectPct(x)}% ${x >= 0 ? 'smaller' : 'larger'}`;

/** A treated animal against its twin on the evaluation day, in words: smaller (or larger) by how much, and its 90%
 *  range, low to high ("60–85% smaller"; across the twin "10% larger to 30% smaller"). */
export function effectWords(effect: number, lower: number, upper: number): string {
  const range = lower >= 0 ? `${effectPct(lower)}–${effectPct(upper)}% smaller` : upper <= 0 ?
    `${effectPct(upper)}–${effectPct(lower)}% larger` : `${sideOf(lower)} to ${sideOf(upper)}`;
  return `${sideOf(effect)} than it would be untreated (90% range: ${range})`;
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

/** The optional columns Analyze takes by their names (the form's, and the demo studies' saved analysis). */
export const OPTIONAL_COLUMNS: {[input: string]: (name: string) => boolean} = {
  bodyWeight: isBodyWeight,
  tumorWeight: (c) => /tum/i.test(c) && /weight|mass/i.test(c),
  removalReason: (c) => /reason|removal|fate|cause|status/i.test(c),
};

/** The tested group with the largest TGI (a group too few to test is not counted: round 4, analyst A4); null: none. */
export function bestTgi(c: DG.DataFrame): {arm: string, tgi: number} | null {
  let best: {arm: string, tgi: number} | null = null;
  for (let i = 0; i < c.rowCount; i++) {
    const x = c.get('tgi', i);
    if (c.get('activity', i) !== TOO_FEW && c.get('n_day', i) > 0 && isFiniteNumber(x) && (!best || x > best.tgi))
      best = {arm: c.get('arm', i), tgi: x};
  }
  return best;
}

/** A cell's text with its reason on hover ("—" and why: round 4, analyst W7). */
function tipped(text: string, tip: string): HTMLElement {
  const el = span(text);
  Tooltip.bind(el, tip);
  return el;
}

/** The Analyzed cell: the group's animals measured on the evaluation day, "7 of 10" when the growth test has more,
 *  why on hover (round 4, analyst U1/W8: in place of a sentence under the table). */
export function analyzedText(group: string, nDay: number, n: number, day: number): {text: string, tip: string} {
  return nDay < n ? {text: `${nDay} of ${n}`, tip: `${group}: ${nDay} of ${n} measured on day ${day} (TGI, T/C); ` +
    `growth test uses all ${n}.`} : {text: String(nDay), tip: `${group}: ${nDay} measured on day ${day}, all in the ` +
    'growth test.'};
}
const analyzedCell = (group: string, nDay: number, n: number, day: number) => {
  const a = analyzedText(group, nDay, n, day);
  return tipped(a.text, a.tip);
};

/** A p-value cell, never broken over two lines ("< 0.001": round 4, analyst B1). */
/** The adjusted p, the unadjusted one on hover (unadjustedP), also where the correction left it unchanged. */
function pCell(p: unknown, raw?: unknown): HTMLElement {
  const el = span(pValue(p));
  el.style.whiteSpace = 'nowrap';
  const hover = unadjustedP(raw);
  if (hover)
    Tooltip.bind(el, hover);
  return el;
}

const correctionWords = (corr: unknown, k: number) => k < 2 ? 'one comparison' : corr === 'none' ? 'no correction' :
  `${CORRECTION_NAME[correctionOf(corr)]}'s correction`;

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

/** A treated group's test on the slopes in words (analyze.ts growthRateWords), row `i` of the comparisons; null:
 *  none. */
const growthOf = (c: DG.DataFrame, i: number) => c.get('activity', i) === TOO_FEW ? null :
  growthRateWords(c.get('slope_difference', i), c.get('slope_difference_lower', i),
    c.get('slope_difference_upper', i), c.get('p_value', i), c.get('slope_difference_confidence', i));

/** The groups table's growth cell (analyze.ts growthVsControl), "−7.1 (−9.6 to −4.5)" %/day; null: no test (a group
 *  with too few animals: its stored NaNs may read back as 0). */
export const growthCellOf = (c: DG.DataFrame, i: number) => c.get('activity', i) === TOO_FEW ? null :
  growthVsControl(c.get('slope_difference', i), c.get('slope_difference_lower', i), c.get('slope_difference_upper', i));

/** The results in words, a line each: the study, each treated group's effect and result, the test, the body weight. */
function summaryLines(r: AnalysisFrames, nci: Map<string, NciReading>): string[] {
  const twinned = r.analysis.get('twinned', 0);
  const v = (k: string) => r.analysis.get(k, 0);
  const c = r.comparisons;
  const control = v('control');
  const lines = [`${v('line') ? `${modelName(v('line'))} study` : 'Study'} evaluated on day ${v('evaluation_day')}` + (v('from_endpoint') ?
    `, the first measurement day the median control reached ${v('endpoint')} mm³` : '') +
    `; control: ${control}, ${v('n_control_day')} animals on that day.`];
  for (let i = 0; i < c.rowCount; i++) {
    const g = (k: string) => c.get(k, i);
    const growth = growthOf(c, i);
    const tc = g('tc') === null ? '' : ` and T/C ${pct(g('tc') / 100)}`;
    lines.push(`${g('arm')} (${g('n_day')} animals): ${g('activity').toLowerCase()} against the control` +
      (growth ? `; ${growth}` : '') + `; TGI ${tgiText(g('tgi'))}${tc} on day ` +
      `${v('evaluation_day')}` + (twinned ? `; ${g('responders')} of ` +
      `${g('twin_animals')} animals responded (below the range they would reach untreated).` : '.'));
  }
  lines.push(`Tumor growth of each treated group against the control, two-sided α = ${num(v('alpha'), 2)} with ` +
    `${correctionWords(v('correction'), v('tested_arms'))}.`);
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

/** What Analyze opens with: a study table (with its plan in the `pmx.plan` tag, if any) and, with the NCA package,
 *  a concentration table for its exposure. */
/** `evaluationDay`: the day to analyze on, in place of Auto (Monitor's last complete day, when the latest is under
 *  way). */
export interface AnalyzeStart {study?: DG.DataFrame, exposure?: DG.DataFrame, evaluationDay?: number}

/** The Analyze view; given a study, it analyzes it at once. */
export async function analyzeView(opening?: AnalyzeStart): Promise<DG.ViewBase> {
  const exports: Exports = [];
  const host = ViewHost.create({name: 'PMX Studio', content: (self) => analyzeContent(self, opening,
    exports)});
  return markStage(host.view, 'analyze');
}

/** What the view shows, built under its host, which owns it all. */
function analyzeContent(self: ViewHost, opening: AnalyzeStart | undefined, exports: Exports): ViewContent {
  // the twin's shaders compiled in Crux's worker while the study is picked, once a page (a failure only leaves the
  // first twin slower); here and not when the app opens, which starts no engine
  void warmTwin().catch(() => null);
  const call = DG.Func.byName('PmxStudio:pmxAnalyzeStudy').prepare();
  const studyDf = () => call.getParamValue('study') as DG.DataFrame | null;
  // a group's object acts on its study only while that study is the one on screen (QA 0.1.3 round 3, analyst 1)
  let closed = false;
  self.own(() => closed = true);
  const liveNow = () => {
    const df = studyDf();
    return () => !closed && !!df && studyDf()?.dart === df.dart;
  };
  // the plan the analysis uses, as text: the one picked in the form, else the study's own
  const planText = () => call.getParamValue('plan')?.trim() || studyDf()?.getTag(PLAN_TAG);
  // the control as its field shows it: emptied, it is automatic again (the call kept the last name typed, and a
  // cleared field showing "Auto" still analyzed with it)
  const controlNow = (): string | null => {
    const box = form.getInput('control')?.root.querySelector('input');
    const v = box ? box.value : call.getParamValue('control');
    return typeof v === 'string' && v.trim() ? v.trim() : null;
  };
  const inputsNow = (extra: object) => withoutExcluded(callInputs<AnalysisInputs>(call,
    {control: controlNow(), ...extra}));
  // the engine itself on the form's inputs, not the registered function's FuncCall (see callInputs)
  const run = (extra: object) => analysisFrames(analysisDeclarationOf(inputsNow(extra)));
  const result = signal<AnalysisFrames | undefined>(undefined);
  // the newest run: the form's changes and an opened study share the call, and a run that finishes after a newer one
  // started is dropped
  let latest = 0;
  const error = signal<string | null>(null);
  // a failure as the user reads it: what is missing as it is said, else what failed and what to try
  const AGAIN = 'Try again; if it fails again, reload the page.';
  const failure = (e: unknown, fallback: string) => {
    // a control all of whose rows are excluded is in the table, just not analyzed: said so
    const s = studyInputsOf(call);
    const inTable = e instanceof NoControl && !!e.control && !!s?.study.col(s.group)?.toList().includes(e.control);
    if (e instanceof NoControl && inTable) {
      return `The control group "${e.control}" is excluded from the analysis: include it again, or pick another ` +
        'one in Control group.';
    }
    return e instanceof EmptyState || e instanceof NoControl ? e.message : plainError(e, fallback);
  };
  // a result that no longer matches the form is dimmed until the next run
  const stale = signal(false);
  // the multiple comparisons as Plan asks for them, over the treated groups of the last analysis (two until then)
  const treatedGroups = signal(2);
  const correction = new CorrectionInput({label: 'Multiple comparisons', arms: treatedGroups});
  // no Run: each change analyzes again once the changes stop
  let preparing = false;
  const auto = autoRun(() => form, () => preparing ? undefined : analyze());
  const overrides = {correction: {input: correction}, line: {input: new LineInput({label: 'Tumor model'})}};
  // the study the page took in (openStudy), and whether the form has settled (a study field changing from then on is
  // the user's pick)
  let shownStudy: DG.DataFrame | null = null;
  let ready = false;
  const form = funcForm(call, {layout: 'tall', overrides, onInputChanged: (name, value) => {
    // a model picked for a table that had none is the study's from now on (no longer "assumed")
    const df = studyDf();
    // another study picked in the field: taken in as an opened one is, its fields as last set here
    if (name === 'study' && df && ready && !preparing && df !== shownStudy)
      return void switchTo(df);
    // (a value the user set for this study: kept for it, given back when the page shows it again; QA 0.1.3 round 3
    // chaos 2)
    if (name !== 'study' && !preparing && df && df === shownStudy)
      rememberForm(df, 'analyze', keptOfCall(call));
    // a study picked: its own model (its plan's), else none, never the last study's nor a default (QA 0.1.3 A2)
    if (name === 'study' && df)
      queueMicrotask(() => form.getInput('line') && (form.getInput('line')!.value = lineOf(df)));
    if (name === 'line' && typeof value === 'string' && value && df && !lineOf(df))
      df.setTag(LINE_TAG, value);
    stale.value = result.peek() !== undefined;
    auto.changed();
  }});
  settleForm(form, {control: 'Auto (the vehicle group)', evaluationDay: 'Auto (endpoint reached)',
    density: 'Auto (the cell line’s in the atlas)'});
  keepNamedColumns(form, OPTIONAL_COLUMNS);
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
  // a failed run leaves no result behind it: what the form or the table now say is not what was shown (another
  // study's figures, its name in the header, a change "from the previous result" across two studies)
  const clearResult = () => {
    result.value = undefined;
    stale.value = false;
    shown(null);
  };
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
        clearResult();
        return;
      }
      // a table picked that is no study (a PK table): said so, not "No control group found"
      if (blinded instanceof DG.DataFrame && !readsAsStudy(blinded)) {
        error.value = NOT_A_STUDY;
        clearResult();
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
      // the twin unchanged (cached): one call, with it; else the analysis first, then the same inputs with the twins
      // (seconds), whatever the form holds by then
      const cached = twinsCached.includes(twinSent);
      const withTwin = cached ? null : inputsNow({withTwin: true});
      const frames = await run({withTwin: cached});
      if (mine !== latest)
        return;
      error.value = null;
      stale.value = askedOf(call) !== sent;
      shown(reads);
      result.value = keepSame(result.peek(), frames);
      dataRead.value = data;
      if (!withTwin)
        return;
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
        clearResult();
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
    shownStudy = df;
    // (named "X (2)" while the page it replaces held its own "X": its name back as that one goes, picked again)
    holdTable(self, df, () => {
      const field = form.getInput('study');
      if (field)
        field.value = df.name;
    });
    // a view not yet shown is left to showStage: making it current would add it, and showStage add it a second time
    if ([...grok.shell.views].some((v) => v.dart === self.view.dart)) {
      grok.shell.v = self.view;
      // (the URL names this study at once: Ctrl+O of a new file no longer keeps the last one's ?study=, QA 0.1.3
      // round 4 chaos N6)
      linkRow(self.view, 'study', rowOfTable(df)?.id);
    }
    const study = form.getInput('study');
    if (study)
      study.value = df.name;
    await form.settled;
    // a table from Monitor brings its body weight column and limits, and the day to read when its latest is under way
    restoreFields(form, df);
    // and what it was analyzed with when saved here: its method, endpoint and evaluation day (QA 0.1.3 r4: a Holm
    // study reopened as Bonferroni)
    for (const [name, value] of Object.entries(analysisFieldsOf(df))) {
      const input = form.getInput(name);
      if (input)
        input.value = value;
    }
    // (what the user set for it on this page before: Back, Forward, the field picking it again)
    applyKept(form, recallForm(df, 'analyze'));
    const day = form.getInput('evaluationDay');
    if (df === opening?.study && opening.evaluationDay !== undefined && day)
      day.value = opening.evaluationDay;
    await form.settled;
  });
  whenSettled(form, () => {
    const studyInput = form.getInput('study') as unknown as TableInput | null;
    keepStudyField(studyInput, () => call.getParamValue('study') ?? null);
    // the picker offers the open tables that look like a study (tagged as one by a stage, or with a volume or caliper
    // column by its name), and the one picked: not every table open (a concentrations table opened for NCA is none)
    if (studyInput) {
      // one way in for new data: File > Open (Ctrl+O) and a drop, through the importer's checks (columns, units,
      // calipers, wide tables, dates); the field only picks among the studies already open. The platform's own
      // open-file, Files and database icons read a table raw (calipers taken as volumes, a wide table as nonsense)
      studyInput.actions = [];
      const offer = studyInput.setItems.bind(studyInput);
      const looksLikeStudy = (name: string) => {
        const t = grok.shell.tables.find((x) => x.name === name);
        return !t || name === studyInput.value || readsAsStudy(t);
      };
      studyInput.setItems = (items: string[]) => offer(studyPickList(items.filter(looksLikeStudy), studyInput.value));
      studyInput.setItems(grok.shell.tables.map((t) => t.name));
    }
    if (opening?.study)
      void openStudy(async () => opening.study!);
    ready = true;
  });
  // another study picked in the Study field while this one has changes not saved: the field goes back to it and asks
  // first (Save / Discard / Cancel), before the page takes the new one; Cancel keeps it, its fields as they were (QA
  // 0.1.3 round 3 analyst 8). ponytail: the form gives no hook before a change, so the field is set back at once
  const switchTo = async (df: DG.DataFrame) => {
    const old = shownStudy;
    // (the view read only with a study already shown: the first one comes while the page is built, before its view)
    const ask = old ? leaveAsks.get(self.view.dart) : undefined;
    if (old && ask && domainStorageOn() && grok.shell.tables.includes(old) && !tableSavers.has(old) &&
      asSaved.changed(old)) {
      (form.getInput('study') as unknown as TableInput | null)!.table = old;
      if (!await ask(true))
        return void openStudy(async () => old);
    }
    return void openStudy(async () => df);
  };
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
  // a study table is no concentration table: never offered as one (round 4, analyst U6). After the input's own
  // following of the open tables, so this list is the last word
  if (exposureField) {
    const pkTables = () => grok.shell.tables.filter((t) => t.name && !t.getTag(STUDY_TAG) && !t.col('volume'))
      .map((t) => t.name);
    exposureField.setItems([...new Set(pkTables())]);
    for (const e of [grok.events.onTableAdded, grok.events.onTableRemoved])
      self.subscribe(e, () => exposureField.setItems([...new Set(pkTables())]));
  }
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
  let lastDetails: object | null = null;
  let lastGroup = '';
  // a group's animals for its details (set with the animal list)
  let animalsOf: (group: string, color: string) => PanelSection[] = () => [];
  // each group's statistics, as the report's Statistics table has them
  let statsOf = new Map<string, PanelFact[]>();
  // the animals each group had on each measurement day (set with the growth chart)
  let byDay: AnimalsByDay | null = null;
  const groupSections = (group: string, color: string): PanelSection[] => [...animalsOf(group, color),
    ...statsOf.has(group) ? [{title: 'Statistics', content: () => factList(statsOf.get(group)!)}] : [],
    ...byDay?.n.get(group)?.size ? [{title: 'Animals by day', content: () => factList(byDay!.days
      .filter((d) => byDay!.n.get(group)!.has(d)).map((d) => ({label: `Day ${d}`,
        value: String(byDay!.n.get(group)!.get(d))})))}] : []];
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
    const df = studyDf();
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
    }, metaOf(studyDf()))
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
    const df = studyDf();
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
    // opened with its concentrations (the demo's mouse PK): read as a table picked in the field
    // (opened in the workspace first, as a table the user picks is: the field offers open tables only)
    if (opening?.exposure) {
      if (!grok.shell.tables.includes(opening.exposure))
        grok.shell.addTable(opening.exposure);
      exposureField.table = opening.exposure;
    }
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
  // no tumor model picked: no twin to wait for (QA 0.1.3 A2), said where it would be, a dash in a table's cell
  const waitingTwin = (words?: string) => error.peek() !== null ? span('—') :
    result.peek()?.analysis.get('line', 0) === '' ? noModelNote(!words) :
      loader(words);
  const coverage = cell('control_coverage');
  const significant = computed(() => {
    const c = comparisons.value;
    // (over the groups tested: one too few to test is no group that fails to slow growth, round 4 A4)
    const act = c?.col('activity')!.toList() ?? [];
    return c ? `${act.filter((x) => x === 'Slows growth').length} of ${act.filter((x) => x !== TOO_FEW).length}` :
      undefined;
  });
  const treatedAnalyzed = computed(() => {
    const c = comparisons.value;
    return c ? (c.col('n_day')!.toList() as number[]).reduce((t, x) => t + x, 0) : undefined;
  });
  const analyzed = computed(() => treatedAnalyzed.value === undefined ? undefined :
    treatedAnalyzed.value + (cell('n_control_day').value ?? 0));
  // the study's answers first, what each figure covers in its tooltip
  const card = statCard.bind(null, form);
  // two figures, the study's answer; the evaluation day, the animals and the controls are in the status line
  const largest = computed(() => {
    const c = comparisons.value;
    if (!c || !c.col('tgi'))
      return undefined;
    return bestTgi(c)?.tgi;
  });
  // (a change is shown against the previous result of the same study only)
  const cards = cardGrid([
    card('Treated groups that slow growth', significant, () => 'Treated groups whose tumors grew slower than the ' +
      `control's up to day ${cell('evaluation_day').value}, of those tested`, undefined, true, studyDf),
    // (the table's rounding, its group and day: round 4, analyst W10)
    card('Best TGI', largest, () => 'The largest TGI among the tested groups on the evaluation day', (x) => {
      const b = comparisons.peek() ? bestTgi(comparisons.peek()!) : null;
      return `${tgiText(x)}${b ? ` (${b.arm}, day ${cell('evaluation_day').peek()})` : ''}`;
    }, true, studyDf),
  ], 'pmx-analyze-cards', 150);
  const statusTip = () => (cell('from_endpoint').peek() ? 'Evaluated on the first measurement day the median ' +
    `control reached ${mm3Text(cell('endpoint').peek() as number)} mm³, as Plan defines the end of the study.` : 'Evaluated on the day ' +
    'chosen, or the last measurement day if the median control never reached the endpoint.') +
    ` ${cell('n_control_day').peek()} controls and ${treatedAnalyzed.peek()} treated animals measured on the ` +
    'evaluation day (the growth test takes every animal measured twice or more: see Analyzed). Controls on ' +
    'track: control measurements inside the range the other controls predict for them (about 90% is as expected).';

  const groupsTable = div([]);
  const openSynergy = () => {
    const r = result.peek();
    if (!r || !keptStudy)
      return;
    // the plan's α if a plan is attached, else Analyze's
    const plan = planText();
    synergyDialog(keptStudy, r.analysis.get('control', 0), r.comparisons.col('arm')!.toList() as string[],
      r.analysis.get('evaluation_day', 0), synergyAlpha(plan, r.analysis.get('alpha', 0)), notCombined(plan));
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
    const off = notCombined(planText());
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
    if (!findCombinations(treated).some((c) => !off.includes(treated[c.ab])))
      return;
    const control: string = r.analysis.get('control', 0);
    const lastDay: number = r.analysis.get('evaluation_day', 0) ?? Infinity;
    const line: string = r.analysis.get('line', 0);
    void (line ? linePrior(line) : Promise.resolve(null)).then(async (lp) => {
      const logBeta = lp?.study_prior?.mean_log_beta;
      const {found, missing} = await combinationSteps(s, treated, off, control, line, lastDay,
        logBeta === undefined ? undefined : Math.exp(logBeta),
        (e, ab) => plainError(e, `${ab}: the comparison with its two drugs could not be computed. ${AGAIN}`));
      if (laddered === mine && found.size)
        ladder.value = found;
      say(missing);
    }).catch((e) => say([plainError(e, `The combinations could not be compared with their drugs just now. ${AGAIN}`)]));
  });
  // above the groups, one sentence on exposure and effect; under them, a muted line on what would add to it
  // (the report gives the sentence once, in its section on exposure and patients)
  const readoutLine = div([], 'pmx-noprint');
  const hintLine = div([]);
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
  const exposurePoints = signal<ExposurePoint[] | null>(null);
  // one drug by one route at 3 doses or more: exposures of different drugs share no axis
  const exposureShown = computed(() => {
    const p = exposurePoints.value;
    return p && new Set(p.map((x) => `${x.drug}|${x.route}`)).size === 1 && new Set(p.map((x) => x.dose)).size >= 3 ?
      p : null;
  });
  // the effect expected between the doses tested (exposureCurve, on the same TGI), or why there is none
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
    if (!p || !r || !s)
      return;
    const [animal, group, day, volume] = [s.animal, s.group, s.day, s.volume].map((c) => s.study.col(c)!);
    const rows = Array.from({length: s.study.rowCount}, (_, i) => i)
      .filter((i) => !animal.isNone(i) && !group.isNone(i) && !day.isNone(i) && !volume.isNone(i))
      .map((i) => ({animal: String(animal.get(i)), group: String(group.get(i)), day: day.get(i),
        volume: volume.get(i)}));
    void exposureCurve(p.map((x) => ({arm: x.arm, auc: x.auc})), rows, control, lastDay, metaOf(studyDf()))
      .then((c) => {
        if (fitting === key)
          exposureFit.value = c;
      }, (e) => {
        if (fitting === key)
          exposureFit.value = {why: plainError(e, 'it could not be computed just now; try again').replace(/\.$/, '')};
      });
  });
  // the groups on one axis at the groups table's TGI (the one TGI, 1 − ΔT/ΔC on the evaluation day),
  // and the curve between the doses tested, dashed with its band
  const exposureData = computed(() => {
    const f = exposureFit.value;
    // each group's dose from its details first, then its name
    const df = exposureFrame(exposureShown.value ?? [], metaOf(studyDf()), f && 'curve' in f ? f.curve : []);
    const day = result.peek()?.analysis.get('evaluation_day', 0);
    exposureOptions.forEach((o) => o.yTitle = `TGI${day != null ? `, day ${day}` : ''}`);
    return df;
  });
  // the Y axis and the tooltip name the measure, with its day (set as the data is built)
  const exposureOptions = ['AUClast', 'Dose'].map((x) => ({x, lower: 'lower', upper: 'upper', line: 'tgi',
    split: 'series', yTitle: 'TGI', yInTooltip: true, markers: true, dashed: [CURVE], xPad: true}));
  const exposureCharts = exposureOptions.map((o) => new BandChart(computed(() =>
    exposureData.value?.col(o.x) ? exposureData.value : undefined), o));
  const exposureText = span('');
  form.effect(() => {
    const f = exposureFit.value;
    exposureText.textContent = 'Each point is one group at its TGI on the evaluation day, as the groups table gives ' +
      `it. ${f && 'curve' in f ? 'The dashed line is the effect expected between the doses tested, never beyond ' +
      `them, on the same TGI.${f.note ? ` ${f.note}` : ''} The band is ±${CURVE_BAND} points.` :
      f ? `No curve: ${f.why}.` : ''}`.trim();
  });
  form.effect(() => {
    exposureCharts.forEach((c) => c.root.style.display = exposureData.value?.col(c === exposureCharts[0] ?
      'AUClast' : 'Dose') ? '' : 'none');
  });
  // its own panel, after growth and time to endpoint: the groups panel keeps the headline (each group's multiple of
  // the patients' exposure, the sentence), the charts would push the growth down
  const exposureBox = divV([exposureText, ...exposureCharts.map((c) => sized(c, 260))]);
  const statsTable = div([]);
  const controlCheck = div([]);
  const protocol = div([]);
  let protocolText = '';
  let statsCsv = '';
  let statsShownCsv = '';
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
    const m = metaOf(studyDf());
    return groupColors(r.comparisons.col('arm')!.toList(), r.analysis.get('control', 0), m?.colors, m?.slots);
  };
  form.effect(() => {
    const r = result.value;
    if (!r)
      return;
    const v = (c: string) => r.analysis.get(c, 0);
    const done = v('twinned');
    void error.value;
    // the selector's hint counts the comparisons made, not the treated groups
    treatedGroups.value = comparisonsMade(r.analysis);
    const control = v('control');
    const arms = r.comparisons.col('arm')!.toList() as string[];
    const colors = colorsOf(r);
    const steps = ladder.value;
    const allNotes = String(v('notes') ?? '').split(' | ').filter((n) => n);
    const odd = allNotes.filter((n) => n.startsWith('Control '));
    groupNotes.replaceChildren(...allNotes.filter((n) => !odd.includes(n)).map((n) => noteLine(n)));

    const nci = nciReadings(r);
    // each treated group's exposure, matched to it by NCA's group name or its drug, dose and route
    // (NCA read again replaces what the study kept, and is kept in its place)
    const run = ncaRun.value;
    const meta = metaOf(studyDf());
    const nca = run ? ncaByArm(arms, run, ncaSource?.table.name ?? '', meta) : savedNca.value;
    shownNca = nca;
    shownEffect = keptEffect(r, meta, nca, steps);
    keepExposure(nca);
    const exposure = nca?.arms ?? null;
    // a combination's drugs, each with its own exposure (none: no NCA, or nothing for any of its drugs)
    const drugX = nca?.drugs ?? new Map<string, DrugExposure[]>();
    // each treated group's exposure as a multiple of the human's, once a human exposure is typed for its drug
    const human = humans.value;
    const typed = [...human.values()].some((h) => h.auc !== null);
    const multiples = new Map<string, ReturnType<typeof humanMultiple>>();
    if (nca && typed) {
      for (const arm of nca.arms.keys())
        multiples.set(arm, humanMultiple(arm, nca.arms.get(arm)?.exposure, nca.pk.get(arm) ?? null, human, meta));
    }
    // a combination's drugs, each against its own patients' exposure
    const drugMs = new Map<string, ReturnType<typeof drugMultiples>>();
    if (typed) {
      for (const [arm, ds] of drugX)
        drugMs.set(arm, drugMultiples(arm, ds, human, meta));
    }
    const drugWhys = [...drugMs.values()].flat().map((d) => d.m?.why).filter((w): w is string => !!w);
    const withHuman = [...multiples.values()].some((m) => m !== null) ||
      [...drugMs.values()].some((ds) => ds.some((d) => d.m !== null));
    // exposure kept under a group name the table no longer has (renamed before Monitor saved the study): said
    const unmatched = !run && exposure ? [...exposure].filter(([a, x]) => x.exposure && a !== control &&
      !arms.includes(a)).map(([a]) => a) : [];
    exposureNotes.replaceChildren(...ncaSays.value ? [noteLine(ncaSays.value)] : [], ...!exposure ? [] :
      [...!run && nca ? [muted(`Exposure from NCA, saved with the study (table ${nca.table || 'not named'})`)] : [],
        ...unmatched.length ? [noteLine(`Exposure kept for ${unmatched.join(', ')}: no group has that name now, so ` +
          'it is not shown')] : [],
        ...[...exposure.values(), ...[...drugX.values()].flat().map((d) => d.x)].filter((x) => x.why)
          .map((x) => noteLine(x.why!)),
        ...[...multiples.values()].filter((m) => m?.why).map((m) => noteLine(m!.why!)),
        ...drugWhys.map((w) => noteLine(w))]);
    // one sentence on exposure and effect, and the one thing that would add to it
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
    // the curve as drawn when it is asked: whether there is one, and where a group sits among the doses drawn
    const drawn = () => {
      const f = exposureFit.peek();
      return !!f && 'curve' in f;
    };
    const curvePlace = (arm: string): 'lowest' | 'highest' | 'between' | 'off' | null => {
      const p = exposureShown.peek();
      if (!drawn() || !p || arm === control)
        return null;
      const axis = exposureAxis(p, meta);
      const at = [...p].sort((a, b) => axis.of(a) - axis.of(b)).findIndex((x) => x.arm === arm);
      return at < 0 ? 'off' : at === 0 ? 'lowest' : at === p.length - 1 ? 'highest' : 'between';
    };
    // the report's section on exposure and patients, as the groups are now
    sectionNow = () => exposureSection({sentence: said, curve: drawn(),
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
    // a treated group's details for the context panel, as computed here
    detailsNow = (arm) => {
      const i = arms.indexOf(arm);
      const parts = armParts(arm, meta);
      const g = (k: string) => i < 0 ? null : r.comparisons.get(k, i);
      const step = steps.get(arm);
      return groupDetails({group: arm, read: parseLine(arm, meta), tgi: g('tgi'), color: colors[arm],
        n: arm === control ? v('n_control_day') : g('n_day'),
        tc: g('tc'), growth: i < 0 ? null : growthOf(r.comparisons, i),
        result: g('activity'), step: step ? {words: stepWords(step.step), why: step.why} : null,
        responders: done && g('twin_animals') ? `${g('responders')} of ${g('twin_animals')}` : null,
        nadir: nci.get(arm)?.nadir ?? null,
        p: g('activity') === TOO_FEW ? null : g('p_value'), removed: removedOf(r.tolerability, arm),
        day: v('evaluation_day') ?? null, x: exposure?.get(arm), table: nca?.table ?? null, live: !!run,
        multiple: multiples.get(arm) ?? null, human: parts.length === 1 ? human.get(parts[0].key) : undefined,
        curve: curvePlace(arm),
        // a combination: each drug's exposure apart
        ...drugX.has(arm) ? {drugs: parts.map((p) => ({drug: p.drug, x: drugX.get(arm)!.find((d) => d.key === p.key)?.x,
          human: human.get(p.key), multiple: drugMs.get(arm)?.find((d) => d.key === p.key)?.m ?? null}))} : {}});
    };
    // a group's details open in the context panel: built again from the analysis and exposure as they are now
    if (lastDetails && grok.shell.o === lastDetails && (arms.includes(lastGroup) || lastGroup === control)) {
      const group = lastGroup;
      grok.shell.o = lastDetails = new PmxGroup(studyOf(studyDf()), group, colors[group], {details: () =>
        detailsElement(detailsNow(group), groupSections(group, colors[group])), live: liveNow()});
    }
    hintLine.replaceChildren(...hint ? [muted(hint.text), ...hint.add ? [ui.button('Add…', () =>
      hint.add === 'patients' ? openHumans() : addConcentrations())] : []] : []);
    exposurePoints.value = !exposure ? null : exposurePointsOf(r, exposure, colors);

    // the control first, then each treated group: its effect on the evaluation day, its animals against their
    // twins and its result (a table built in this run is disposed with it)
    const c = r.comparisons;
    type Row = {name: string, cells: (string | HTMLElement)[]};
    const day = v('evaluation_day');
    const rows: Row[] = [{name: control, cells: [swatch(control, colors[control]),
      analyzedCell(control, v('n_control_day'), v('n_control'), day), '—', '—', '—', '—', '—', '—',
      withNci(span('Reference'), nci.get(control), showNci)]}];
    for (let i = 0; i < c.rowCount; i++) {
      const g = (k: string) => c.get(k, i);
      const tc = g('tc');
      // (none measured on the day: no TGI nor T/C; too few animals: no p nor verdict, chaos 4)
      const tooFew = g('activity') === TOO_FEW;
      const none = !(g('n_day') > 0);
      const verdict = tooFew ? tipped('—', `Not tested: ${tooFewWhy(g('n'))}`) :
        divH([dot(g('activity') === 'Slows growth' ? 'success' : 'warning'), span(g('activity'))]);
      Object.assign(verdict.style, {alignItems: 'center', gap: 'var(--dg-space-s)'});
      // no clear effect: what the study could have detected with its animals
      // why, on hover: what the study could have detected
      const why = noEffectNote((k) => r.versusPlan.get(k, i), v);
      if (why)
        Tooltip.bind(verdict, why);
      // the test the result is read on, its one interval: numbers in the cell, the sentence on hover (round 3 W1)
      const growth = growthOf(c, i);
      const growthCell = span(growthCellOf(c, i) ?? '—');
      Tooltip.bind(growthCell, growth ? `${growth}. ${GROWTH_RATE_TIP}` : tooFew ? `Not tested: ${tooFewWhy(g('n'))}` :
        GROWTH_RATE_TIP);
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
      // the group's name stands for it: its right click offers the group's actions
      const name = bindObject(new PmxGroup(studyOf(studyDf()), g('arm'), colors[g('arm')], {details: () =>
        detailsElement(detailsNow(g('arm')), groupSections(g('arm'), colors[g('arm')])), live: liveNow()}),
      swatch(g('arm'), colors[g('arm')]));
      const read = parseLine(g('arm'), meta);
      if (read)
        Tooltip.bind(name, read);
      const problem = readProblem(armParts(g('arm'), meta));
      // (grey, no warning sign: what could not be read, why on hover: round 4, analyst W9)
      const problemLine = problem ? muted(problem) : null;
      if (problemLine) {
        Tooltip.bind(problemLine, `Read from the group's name: ${problem}. Exposure and dose-response need a dose in ` +
          'mg/kg: add it to the name or in Study details.');
      }
      const nameCell = problemLine ? divV([name, problemLine]) : name;
      rows.push({name: g('arm'), cells: [nameCell, analyzedCell(g('arm'), g('n_day'), g('n'), day),
        none || !isFiniteNumber(g('tgi')) ? tipped('—', `No animal of ${g('arm')} measured on day ${day}`) :
          tgiText(g('tgi')),
        none || !isFiniteNumber(tc) ? '—' : pct(tc / 100),
        !done ? waitingTwin() : !isFiniteNumber(g('paired_effect')) || !g('twin_animals') ? '—' :
          pct(g('paired_effect') / 100),
        !done ? waitingTwin() : !g('twin_animals') ? '—' : `${g('responders')} of ${g('twin_animals')}`, growthCell,
        tooFew ? tipped('—', `Not tested: ${tooFewWhy(g('n'))}`) : pCell(g('p_value'), g('p_raw')),
        withNci(resultCell, nci.get(g('arm')), showNci)]});
    }
    // each treated group's exposure as a multiple of the patients', once one is placed (a combination: each drug's)
    if (withHuman) {
      const tipped = (text: string, tip?: string) => {
        const el = span(text);
        if (tip)
          Tooltip.bind(el, tip);
        return el;
      };
      const multipleCell = (arm: string): string | HTMLElement => {
        const ds = drugMs.get(arm);
        if (ds?.some((d) => d.m?.ratio !== undefined)) {
          return divV(ds.map((d) => d.m?.ratio !== undefined ? tipped(`${d.drug} ${shownMultiple(d.m)}`, d.m.tip) :
            tipped(`${d.drug} —`, d.m?.why)));
        }
        const m = multiples.get(arm);
        return m?.ratio !== undefined ? tipped(shownMultiple(m), m.tip) : m?.why ? tipped('—', m.why) : '—';
      };
      for (const row of rows)
        row.cells.splice(row.cells.length - 1, 0, row.name === control ? '—' : multipleCell(row.name));
    }
    // the interval's level and why on the growth header's hover (W2); p: as the correction makes it
    // (the level and k of the comparisons made, the cells' own: round 4, analyst A3/A4)
    const ci = intervalWords(v('confidence'), v('correction'), v('tested_arms'));
    const GROWTH = `Growth vs control, %/day (${ci.header})`;
    const headers = ['Group', 'Analyzed', 'TGI', 'T/C', 'Smaller than its twin (median)', 'Below its twin', GROWTH, ci.p,
      ...withHuman ? ['× patients’ exposure'] : [], 'Result'];
    const tips: {[h: string]: string} = {
      // (the sign said first, QA 0.1.3 round 4 technician wording)
      [GROWTH]: `Negative = grows slower than control. How much more (+) or less (−) the tumors grow each day ` +
        `than the control's, with its ${ci.tip} ${GROWTH_RATE_TIP}`,
      [ci.p]: `Welch's t test on the growth rates, ${correctionWords(v('correction'), v('tested_arms'))}.`,
      'Smaller than its twin (median)': 'The median, over the group’s animals, of how much smaller each tumor is ' +
        'than its digital twin predicts it would be untreated, on the evaluation day.',
      'Below its twin': 'Animals whose tumor is below the 90% range its digital twin predicts untreated, on the ' +
        'evaluation day (not mRECIST responses: the regulatory table counts those as CR and PR).'};
    const table = new BasicTable<Row>({items: rows, keyOf: (x) => x.name, columns: headers.map((header, j) => ({header,
      tooltip: tips[header], render: (x: Row) => x.cells[j],
      align: j >= 1 && j < headers.length - 1 ? 'center' : 'left'})),
    // a group's details in the context panel
    onRowClick: (x) => {
      lastDetails = showDepth(new PmxGroup(studyOf(studyDf()), lastGroup = x.name, colors[x.name], {details: () =>
        detailsElement(detailsNow(x.name), groupSections(x.name, colors[x.name])), live: liveNow()}));
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
    // in %/day as the groups table (the slopes' ln V per day stay in the statistics file: round 3, analyst 11)
    const daily = (s: unknown) => isFiniteNumber(s) ? oneDecimal(pctPerDay(s)) : '—';
    const stats: Stat[] = [{name: control, cells: [control, daily(c.get('slope_control_per_day', 0)),
      doubling(c.get('slope_control_per_day', 0)), '—', '—', reached(control), '—']}];
    for (let i = 0; i < c.rowCount; i++) {
      const g = (k: string) => c.get(k, i);
      stats.push({name: g('arm'), cells: [g('arm'), daily(g('slope_per_day')), doubling(g('slope_per_day')),
        // (never broken over two lines: B1)
        growthCellOf(c, i) ?? '—', g('activity') === TOO_FEW ? '—' : pValue(g('p_value')).replace(' ', ' '),
        reached(g('arm')), later(g('arm'))]});
    }
    // the interval the result is read on: 95%, or wider with a correction for the comparisons
    const statHeaders = ['Group', 'Tumor growth per day (%)', 'Doubling time (days)',
      `Growth vs control, %/day (${ci.header})`, ci.p, `Median day to ${v('endpoint')} mm³`,
      'Days later than the control'];
    const st = new BasicTable<Stat>({items: stats, keyOf: (x) => x.name, columns: statHeaders.map((header, j) =>
      ({header, render: (x: Stat) => x.cells[j], align: j >= 1 ? 'center' : 'left'}))});
    st.root.style.width = '100%';
    wrapCells(st.root);
    statsOf = new Map(stats.map((x) => [x.name,
      statHeaders.slice(1).map((h, j) => ({label: h, value: x.cells[j + 1]}))]));
    statsTable.replaceChildren(st.root);
    // the file leads with the headline numbers the groups table shows, TGI and T/C on the evaluation day (none for
    // the control), as numbers
    const headline = ['tgi', 'tc'];
    const headlineHeaders = ['TGI (%)', 'T/C (%)'];
    const armRow = new Map((c.col('arm')!.toList() as string[]).map((a, i) => [a, i]));
    // (the text as shown, never a typed table: the platform's toCsv wrote a shown 129 as "129.0")
    statsShownCsv = csvText(stats.map((s) => {
      const i = armRow.get(s.name);
      return {[statHeaders[0]]: s.cells[0], ...Object.fromEntries(headline.map((k, j) => {
        // (as the groups table shows them: whole percent, round 3 analyst 11)
        const x = i === undefined || !(c.get('n_day', i) > 0) ? null : c.get(k, i);
        return [headlineHeaders[j], isFiniteNumber(x) ? Math.round(x) : ''];
      })), ...Object.fromEntries(statHeaders.slice(1).map((h, j) => [h, s.cells[j + 1]]))};
    }));
    // the same to analyze further: a column per number (each bound of an interval its own), numbers as numbers, the
    // animals of each group in the test and on the evaluation day
    statsCsv = csvText(statisticsRows(r, (g) => kmMedians.value && km.has(g) ? km.get(g) ?? null : null));

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
    const k = v('tested_arms');
    const corr = v('correction');
    const rowsOf: {[row: string]: string} = {
      'Primary comparison': 'Each treated group against the control on tumor growth: the per-animal slope of log ' +
        'volume from randomization to the evaluation day, by Welch’s t test (each group its own variance)',
      'Multiple comparisons': k < 2 ? 'One treated group: none needed' : corr === 'none' ?
        `None: each of the ${k} comparisons at α = ${num(v('alpha'), 2)} (exploratory)` :
        correctionOver(corr, k),
      'Significance': `α = ${num(v('alpha'), 2)}, two-sided; ${ci.tip}`,
      'Evaluation day': v('from_endpoint') ? `Day ${v('evaluation_day')}, the first measurement day the median ` +
        `control reached ${v('endpoint')} mm³, as Plan defines the end of the study` : `Day ${v('evaluation_day')}`,
      'Effects': 'TGI = 1 − ΔT/ΔC, each animal\'s change from its own volume at randomization, and T/C of the ' +
        'mean volumes, on the evaluation day (NCI definitions, as declared in Plan): descriptions of that day, ' +
        'without an interval; the result and its interval come from the primary comparison alone',
      // (no model: none computed, never described as if it were: round 4, analyst W14)
      'Untreated twin': !v('line') ? 'Not computed (no tumor model chosen).' : `What each treated animal would ` +
        `measure untreated, predicted from its own tumor at randomization, this study's controls and the ` +
        `${modelName(v('line'))} growth model. Its effect: ` +
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
  // nothing predicted (no tumor model, no plan): no switch for it (round 4, analyst U5)
  effect(() => {
    const p = result.value?.predicted;
    predictedInput.root.style.display = p && Array.from({length: p.rowCount}, (_, i) => i)
      .some((i) => p.get('series', i) === PREDICTED) ? '' : 'none';
  });
  const growthData = computed(() => {
    const r = result.value;
    if (!r)
      return undefined;
    const p = r.predicted;
    const ahead = Array.from({length: p.rowCount}, (_, i) => i).filter((i) => p.get('series', i) === PREDICTED);
    if (!showPredicted.value || ahead.length === 0)
      return r.volumes;
    const v = r.volumes;
    // (the animals behind each point, for the tooltip and the thin stretches; none for the prediction)
    const row = (df: DG.DataFrame, i: number, group: string) => ({group, day: df.get('day', i),
      median: df.get('median', i), lower: df.get('lower', i), upper: df.get('upper', i),
      n: df.col('n') && !df.col('n')!.isNone(i) ? df.get('n', i) : NaN});
    const df = DG.DataFrame.fromObjects([...ahead.map((i) => row(p, i, PREDICTED)),
      ...Array.from({length: v.rowCount}, (_, i) => row(v, i, v.get('group', i)))])!;
    for (const c of ['day', 'median', 'lower', 'upper']) {
      df.col(c)!.meta.units = v.col(c)?.meta.units ?? '';
      df.col(c)!.meta.format = v.col(c)?.meta.format ?? '';
    }
    paint(df, 'group', {[PREDICTED]: PREDICTED_COLOR, ...colorsOf(r)});
    df.meta.formulaLines.addAll(v.meta.formulaLines.items);
    return df;
  });
  const growth = new BandChart(growthData, {x: 'day', lower: 'lower', upper: 'upper',
    line: 'median', split: 'group', xTitle: 'Day after randomization', yTitle: 'Tumor volume', log: growthAxis.log,
    markers: true, count: 'n'});

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
  // the row of the animal drawn against its twin, marked (the others not)
  const press = (row: HTMLElement, on: boolean) => {
    row.setAttribute('aria-pressed', String(on));
    row.style.background = on ? 'color-mix(in srgb, var(--dg-accent, #2083d5) 12%, transparent)' : '';
  };
  // an animal's details in the context panel: how it ends against its twin, the reading folded
  // each treated animal against its twin, as the list below last read them
  const twins = new Map<string, AnimalRow>();
  const verdict = (a: AnimalRow) => a.lower > 0 ? 'Responder' : a.upper < 0 ? 'Larger than its twin' :
    'Within its twin\'s range';
  // an animal in the context panel, the same from the list, the waterfall, the curves, a group's animals and the
  // animals table: its answer and effect bar, its curve (against its twin), its place in its group, its weight, its
  // details; the charts there made once and given each animal's data, let go with the view
  const groupOfAnimal = (animal: string): string | undefined => {
    const own = result.peek()?.animals;
    const at = own ? (own.col('animal')!.toList() as unknown[]).findIndex((x) => String(x) === animal) : -1;
    return at < 0 ? undefined : String(own!.get('group', at));
  };
  const colorOfAnimal = (animal: string): string | undefined => {
    const r = result.peek();
    const g = groupOfAnimal(animal);
    return r && g !== undefined ? colorsOf(r)[g] : undefined;
  };
  const twinAxis = volumeAxis();
  const animalPanels = new AnimalPanels({names: () => studyOf(studyDf()), color: colorOfAnimal,
    view: (animal) => animalView(animal), log: twinAxis.log,
    // marked on the page: selected alone (the list, the waterfall, the curves and the twin chart follow)
    onShow: (animal) => selection.select(new Set([animal])),
    // another study opened, a value corrected, a group renamed: the panel shown follows (A1, T11, D18); another
    // study picked in the form clears it at once, not once its analysis is back (`stale`: round 3, analyst 1)
    changed: () => [result.value, dataRead.value, error.value, stale.value]});
  self.own(() => animalPanels.dispose());
  // every entry point's
  const showAnimal = (animal: string) => {
    if (!studyOf(studyDf()))
      return;
    lastDetails = animalPanels.show(animal);
    lastGroup = '';
  };
  // a click on an animal anywhere: with Ctrl, Cmd or Shift only the selection changes
  const clickAnimal = (animal: string, e?: MouseEvent | KeyboardEvent) => {
    if (e && (e.ctrlKey || e.metaKey || e.shiftKey))
      selection.pick(animal, e);
    else
      showAnimal(animal);
  };
  const animalView = (animal: string): AnimalView => {
    const n = studyOf(studyDf());
    const r = result.peek();
    if (!n || !r)
      return {};
    const color = colorOfAnimal(animal);
    const a = twins.get(animal);
    const own = r.animals;
    const at = (own.col('animal')!.toList() as unknown[]).findIndex((x) => String(x) === animal);
    const o = (k: string) => at < 0 || !own.col(k) || own.col(k)!.isNone(at) ? null : own.get(k, at);
    const rate = o('slope_per_day');
    const best = o('best_response');
    const an = result.peek()?.analysis;
    const [day, endpoint] = [an?.get('evaluation_day', 0) ?? null, an?.get('endpoint', 0) ?? null];
    // the growth in plain words, folded: the rate in ln V per day and mRECIST's class on hover
    const growth: PanelFact[] = [
      ...o('doubling_time') !== null ? [{label: 'Doubles', value: `every ${num(o('doubling_time'))} days`,
        ...rate !== null ? {tip: `Tumor growth: ${oneDecimal(pctPerDay(rate))}% per day`} : {}}] : [],
      ...best !== null ? [{label: 'Best response', value: `${pct(Math.abs(best) / 100)} ${best <= 0 ? 'smaller' :
        'larger'} than on day 0`, ...o('mrecist') ? {tip: `mRECIST: ${o('mrecist')}`} : {}}] : []];
    // (the endpoint without a day of its own: "not reached by day 16" read against a volume of day 21)
    const mine: PanelFact[] = o('time_to_endpoint') === null ? [] : [{label: 'Endpoint',
      value: o('endpoint_reached') ? `Reached on day ${num(o('time_to_endpoint'))}` :
        `Not reached${endpoint !== null ? ` (${num(endpoint)} mm³)` : ''}`}];
    // the answer, as the list says it, and under it the list's bar (its range and day on hover)
    const answer: PanelVerdict | null = !a ? null : {text: a.lower > 0 ? `Responder: ${sideOf(a.effect)} than its ` +
      'twin' : verdict(a), tone: a.lower > 0 ? 'ok' : a.upper < 0 ? 'warn' : 'neutral'};
    const bar = a ? div([effectAxis(-50), effectBar(a.effect, a.lower, a.upper, -50, color ?? '#9aa5b1')],
      'pmx-animal-effect') : null;
    if (bar)
      ui.tooltip.bind(bar, `${day !== null ? `Day ${day}: ` : ''}${effectWords(a!.effect, a!.lower, a!.upper)}`);
    // its curve: against its twin; a control (or a twin not in yet), its own growth
    const readings = animalReadings(n, animal);
    const itself = growthFrame(readings.days, readings.volumes, color);
    const curve = a && r.twins?.rowCount ?
      {title: 'Against its twin', frame: animalFrame(r.twins, animal, color ?? '')} :
      itself ? {title: 'Growth', frame: itself} : null;
    // its place in its group by best response, the page's waterfall order (the best on the right)
    // a treated animal with its twin: by its effect against its twin, the list's order (the most smaller first: the
    // panel's place and the list's agree, QA 0.1.3 A11); else (a control) by best response, the smallest growth or
    // the largest named as such
    const g = groupOfAnimal(animal);
    const byTwin = a ? [...twins.values()].filter((t) => t.group === g && Number.isFinite(t.effect))
      .map((t) => ({animal: t.animal, group: t.group, value: t.effect})) : [];
    const bars = byTwin.length ? byTwin : Array.from({length: own.rowCount}, (_, i) =>
      ({animal: String(own.get('animal', i)), group: String(own.get('group', i)),
        value: own.col('best_response') ? own.get('best_response', i) : NaN}))
      .filter((b) => b.group === g && Number.isFinite(b.value));
    const place = byTwin.length ? byTwin.filter((t) => t.value > a!.effect).length + 1 : 0;
    const words = byTwin.length ? `${ordinal(place)} of ${bars.length} against its twin, ${sideOf(a!.effect)}` :
      best !== null && placeWords(bars, animal) ? `${placeWords(bars, animal)}, ${gain(best)} best` : '';
    const group = words ? {bars, order: 'desc' as const, summary: words} : null;
    return {verdict: answer, bar, curve, group, extra: [...growth, ...mine],
      removedDay: removedAnimals(r.tolerability).get(animal)?.day ?? null,
      weightLimit: call.getParamValue('weightLimit') ?? 20};
  };
  // a group's animals for its context panel: each one's curve, last volume and status; a click shows the animal
  const MARK_WORDS: {[m: string]: string} = {excluded: 'Excluded', partly: 'Some excluded', flagged: 'Flagged'};
  animalsOf = (group, color) => {
    const n = studyOf(studyDf());
    if (!n)
      return [];
    const [g, an] = [n.study.col(n.group)!, n.study.col(n.animal)!];
    const names = [...new Set(Array.from({length: n.study.rowCount}, (_, i) => i)
      .filter((i) => String(g.get(i)) === group).map((i) => String(an.get(i))))]
      .sort((x, y) => x.localeCompare(y, undefined, {numeric: true}));
    if (!names.length)
      return [];
    const marks = animalMarks(n);
    // an animal removed early, and why, as Tolerability lists it
    const removed = removedAnimals(result.value?.tolerability);
    const removedAs = (x: string) => removed.has(x) ?
      `Removed on day ${removed.get(x)!.day} (${removed.get(x)!.why})` : '';
    const status = (x: string) => [MARK_WORDS[marks.get(x) ?? ''], removedAs(x),
      twins.has(x) ? verdict(twins.get(x)!) : '']
      .filter((w) => w && !w.startsWith('Within')).join(' · ');
    // (folded, its count in its header: open, its lines pushed the other panes off the panel)
    return [{title: `Animals ${names.length}`, content: () => textLines(names.map((x) =>
      animalLine(n, x, color, status(x), () => showAnimal(x))))}];
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
    const mark = (animal: string | undefined) => buttons.forEach((b, a) => press(b, a === animal));
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
        row.dataset.animal = a.animal;
        // a plain click shows it and selects it (never clears it: picked from the waterfall, it stays); Ctrl, Cmd or
        // Shift toggle it in the selection
        const choose = (ev: MouseEvent | KeyboardEvent) => clickAnimal(a.animal, ev);
        row.addEventListener('click', choose);
        // its right click (and the context panel's Actions): the animal's actions; from the menu it becomes the
        // current object too
        bindObject(animalPanels.objectOf(a.animal), row);
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
      animalList.querySelectorAll<HTMLElement>('.pmx-forest-row')
        .forEach((b) => press(b, listed && b.dataset.animal === one));
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
  // under the growth chart, one line when a group loses animals (each point's count is in its tooltip)
  const perDay = div([], 'pmx-n-per-day');
  // the study's animals or measurements as Datagrok tables (to chart, filter, cluster there); the animals picked in
  // Analyze come marked in them
  const openTable = (which: 'animals' | 'measurements') => {
    const r = result.peek();
    if (r && keptStudy)
      explore(which, r, keptStudy, colorsOf(r), self.view, selection.animals.peek(), showAnimal);
  };
  const openAnimalsLink = div([hintLink('Open as table', () => openTable('animals'))], 'pmx-noprint');
  const openMeasurementsLink = div([hintLink('Open measurements', () => openTable('measurements'))],
    'pmx-noprint');
  const volumesOnly = only(result, ['volumes'], groupsOf);
  form.effect(() => {
    const r = volumesOnly.value;
    byDay = r ? animalsByDay(r.volumes, [r.analysis.get('control', 0), ...r.comparisons.col('arm')!.toList()]) : null;
    const line = byDay && lostLine(byDay);
    perDay.replaceChildren(...line ? [noteLine(line)] : []);
  });

  // each animal as a bar in its group's colour, on what the scientist compares them by: its growth rate (the slope of
  // log volume the test compares), its best response (mRECIST's, % from randomization) or how much smaller it is than
  // its untreated twin; sorted so the better animals end on the right
  const WATERFALL = {
    'Growth rate': {col: 'slope_per_day', title: 'Tumor growth per day (%)', order: 'desc', control: true,
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
      'volume from day 10 against randomization. Against its twin: how much smaller than it would be untreated.'});
  const animalsOnly = only(result, ['animals', 'effects'], groupsOf);
  const waterfallData = computed(() => {
    const r = animalsOnly.value;
    const by = WATERFALL[waterfallBy.value];
    const t = by.col === 'effect' ? r?.effects : r?.animals;
    if (!r || !t?.col(by.col))
      return undefined;
    const control = r.analysis.get('control', 0);
    // (a growth rate in %/day, as the groups table: designer round 3 wording 3)
    const value = (i: number) => {
      const x = t.get(by.col, i) as number;
      return by.col === 'slope_per_day' && Number.isFinite(x) ? pctPerDay(x) : x;
    };
    const bars = Array.from({length: t.rowCount}, (_, i) => ({animal: String(t.get('animal', i)),
      group: String(t.get('group', i)), value: value(i)}))
      // a missing value may come back as null, undefined or NaN: only numbers make a bar
      .filter((x) => (by.control || x.group !== control) && Number.isFinite(x.value));
    const groups = [control, ...r.comparisons.col('arm')!.toList()].filter((g) => bars.some((b) => b.group === g));
    return bars.length ? {bars, colors: colorsOf(r), title: by.title, order: by.order, groups} : undefined;
  });
  // a bar's right click: the animal's actions (a canvas cannot be bound: a menu of the same actions)
  // a bar's click: the animal selected everywhere and its card in the context panel (not the platform's "n selected
  // rows"); with Ctrl, Cmd or Shift only the selection changes
  const waterfall = new WaterfallChart(waterfallData,
    (animal, e) => popupFor(animalPanels.objectOf(animal), {causedBy: e}),
    {animals: selection.animals, pick: clickAnimal});
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
  }), {animals: selection.animals, pick: clickAnimal});

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

  // short enough for the value axis' title: in %/day, as the groups table
  const RATE = 'Tumor growth per day (%)';
  const slopesData = computed(() => {
    const r = animalsOnly.value;
    if (!r)
      return undefined;
    const a = r.animals;
    const df = DG.DataFrame.fromObjects(Array.from({length: a.rowCount}, (_, i) => ({'Group': a.get('group', i),
      'Animal': String(a.get('animal', i)), [RATE]: isFiniteNumber(a.get('slope_per_day', i)) ?
        pctPerDay(a.get('slope_per_day', i)) : null}))
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

  // the study's data in Datagrok's own table view; here, one line on it and the way there
  const dataLine = studyDataLine(() => studyInputsOf(call));
  self.own(dataLine.dispose);
  // after a failed run too: the control's rows excluded fail the analysis, and the line still counted them (#19)
  form.effect(() => {
    if (result.value || error.value)
      dataLine.refresh();
  });

  // the figures answer the study's question; each panel, one closer question
  // (the exposure panel's hint says the axis it shows, set below)
  const axisNote = span(axisHint('AUClast'));
  const panels = panelTabs([
    {id: 'groups', title: 'Experimental groups', hint: 'Each treated group against the control on the evaluation day',
      content: [readoutLine, groupsTable, groupNotes]},
    {id: 'notes', title: 'Notes on the results', hint: 'The combinations and the exposure, in words',
      content: [hintLine, ladderNotes, exposureNotes], summary: computed(() => result.value ?
        'the combinations against their drugs, the exposure' : '')},
    {id: 'rates', title: 'Growth rates', hint: 'Each animal’s tumor growth per day (%), from the slope of its log ' +
      'volume: what the test compares', content: [sized(slopes, 380)], summary: computed(() => result.value ?
      'each animal’s growth per day, by group' : '')},
    {id: 'growth', title: 'Observed tumor growth', hint: 'Median and range (smallest to largest tumor) ' +
      'of each group, to the evaluation day; in grey, the median control predicted (by the plan, else by the ' +
      'atlas for this model)',
    content: [divH([predictedInput.root, openMeasurementsLink], 'pmx-growth-tools'), growthAxis.row,
      sized(growth, 380), perDay]},
    {id: 'waterfall', title: 'Each animal compared', hint: 'Every animal as a bar in its group\u2019s color, on its ' +
      'growth rate, its best response or its effect against its twin; the better animals on the right',
    content: [byInput.root, waterfallMissing, sized(waterfall, 300), openAnimalsLink]},
    {id: 'spider', title: 'Each tumor over time', hint: 'Every animal’s volume as a % change from ' +
      'randomization, in its group’s color; click a group in the legend to hide or show its animals',
    content: [sized(spider, 300)], summary: computed(() => result.value ? 'every animal’s volume as a % change' : '')},
    {id: 'km', title: 'Time to endpoint', hint: 'Kaplan-Meier: the animals of each group still below the ' +
      'endpoint volume, day by day, as measured; an animal that never reached it counts until its last measurement',
    content: [sized(km, 320), kmWords], summary: computed(() => result.value ? 'the animals of each group still ' +
      'below the endpoint, day by day' : '')},
    {id: 'exposure', title: 'Exposure and effect', hint: axisNote, content: [exposureBox],
      summary: computed(() => exposureData.value ? `each group’s effect against its ${exposureData.value
        .col('Dose') ? 'dose' : 'exposure'}` : '')},
    {id: 'animals', title: 'Each animal against its twin', hint: 'Its twin: what each treated tumor would have ' +
      'measured untreated, predicted from its size at randomization, the study\'s controls and the cell line\'s ' +
      'growth. Click an animal to see its curve against its twin', content: [animals],
    summary: computed(() => result.value ? 'click an animal for its curve against its twin' : '')},
    {id: 'tolerability', title: 'Tolerability', hint: 'Whether the treatments were well tolerated',
      content: tolerability.content},
    {id: 'plan', title: 'Plan vs observed', hint: 'Whether the study went as it was planned',
      content: versus.content},
    {id: 'controls', title: 'Control check', hint: 'Each control against the range the other controls predict',
      content: [controlCheck]},
  ], ANALYZE_TABS, 'pmx-analyze-panels', 'analyze', form);
  exports.push(['Customize this page', () => customizePage(panels.root), 'View']);
  // the exposure panel only when there is a chart to show (the layout editor still lists it)
  form.effect(() => {
    const cell = exposureBox.closest<HTMLElement>('.pmx-panel');
    if (cell)
      cell.style.display = exposureData.value ? '' : 'none';
    // the hint says the axis shown: the mouse exposure, or the dose when NCA read none
    if (exposureData.value)
      axisNote.textContent = axisHint(exposureData.value.col('Dose') ? 'Dose' : 'AUClast');
  });
  // the study's colours: a colour changed in a chart is the study's, in every stage open on it
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
  // (a report says the mismatch in its header: the line, an instruction, is for the screen)
  const dateLine = markLine('');
  dateLine.classList.add('pmx-noprint');
  form.effect(() => {
    void result.value;
    const wrong = startMismatch(metaOf(studyDf())?.start ?? null, lastMeasured());
    dateLine.hidden = !wrong;
    dateLine.lastElementChild!.textContent = wrong ? `${wrong} Set day 0's date in the study's details (click the ` +
      'study\'s name).' : '';
  });
  // a table with no tumor model of its own (no plan, never given one): the one in the form is assumed, and said
  const assumedLine = markLine('');
  form.effect(() => {
    void result.value;
    const df = studyDf();
    const line = String(call.getParamValue('line') || '');
    const assumed = !!result.value && !!df && !lineOf(df) && !planText();
    assumedLine.hidden = !assumed;
    // (a table read from a file names no model: none is preset, so the field says "Not set" until one is picked)
    assumedLine.lastElementChild!.textContent = !assumed ? '' : line ? `Tumor model assumed: ${modelName(line)}. ` +
      'Pick the study\'s own in Tumor model.' : 'Tumor model: not set. Pick the study\'s in Tumor model for the twins ' +
      'and the atlas checks (none until then: no model is assumed).';
  });
  const content = divV([dataLine.root, assumedLine, twinLine, dateLine, cards, panels.root], 'pmx-analyze-content');
  const copied = (text: string, what: string) => {
    studyNow();
    return navigator.clipboard.writeText(text).then(() => grok.shell.info(`${what} copied`));
  };
  const studyNow = () => {
    const df = studyDf();
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
      // titled by the study, its ID after it, as its header (the analysis' table name read "Analysis · 6 groups":
      // round 3, W11)
      const id = metaOf(df)?.id;
      printSheet(await panels.all(() => reportHtml({title: `${df.name}${id ? ` ${id}` : ''}`, stage: 'Analyze',
        meta: metaOf(df),
        groups: [r.analysis.get('control', 0), ...r.comparisons.col('arm')!.toList() as string[]],
        // the plan the analysis used (analyze.ts)
        content, log: logOf(df), planAtlas: planAtlas(planText()),
        annex: [['Results', summary], ...section ? [['Exposure and patients', sectionElement(section)] as
          [string, HTMLElement]] : [], ['Statistics', statsTable], ...byDay && lostLine(byDay) ?
          [['Animals by day', byDayTable(byDay)] as [string, HTMLElement]] : [], ['Analysis protocol', protocol]],
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
      const sub = docSubtitle(metaOf(df)?.id, `From the analysis of ${df.name}`);
      const d = arriveDialog(items, {full, details: () => head.details(() => arrive(full)), options: which.root,
        download: async () => download(`${df.name} - ${title}.docx`, await docx(checklistBlocks({title, sub, items,
          credit: ARRIVE_CREDIT, last: ['Section, line', 'where']}), {title, landscape: true}), DOCX_TYPE)});
    }, 'Documents'],
    // the study as a row of the nonclinical dossier (CTD), each group a sub-row; its checklist is the preview
    ['Regulatory (CTD) table…', function ctd() {
      const {df, r} = studyNow();
      const s = ctdStudyOf(r, metaOf(df), r.analysis.name, shownNca, stepsInWords(ladder.peek()));
      const sub = docSubtitle(s.id, `From the analysis of ${df.name}`);
      checklistDialog({title: `Regulatory (CTD) table${s.id ? `: ${s.id}` : ''}`, guide: 'ICH M4S(R2)', url: M4S_URL, asker: 'The table',
        download: async () => download(`${df.name} - CTD table.docx`,
          await docx(ctdBlocks([s], {sub}), {title: CTD_TITLE, landscape: true}), DOCX_TYPE)},
      ctdItems(s), {details: () => head.details(ctd)});
    }, 'Documents'],
    // the study as Monitor exports it: its file (with what Analyze added: the exposure, exclusions) and its calendar
    ['Study file (.pmxstudy)', async () => {
      const {df} = studyNow();
      const limits = JSON.parse(df.getTag(LIMITS_TAG) || '{}');
      const s = await stamped(savedStudyOf(df, limits.humane ?? null, limits.weightLimit ?? null));
      download(`${exportName(df)}.pmxstudy`, JSON.stringify(s), 'application/json');
    }],
    ['Calendar (.ics)', () => {
      const {df, r} = studyNow();
      const days = [...new Set(roleColumns(df).day.toList().filter(isFiniteNumber))].sort((a, b) => a - b);
      const evaluation = r.analysis.get('evaluation_day', 0);
      const events = [...days.map((d): [number, string] => [d, 'Tumors measured']),
        ...isFiniteNumber(evaluation) ? [[evaluation, 'Evaluation day'] as [number, string]] : []];
      download(`${exportName(df)}.ics`, icsOf(exportName(df), metaOf(df) ?? emptyMeta(), lastMeasured() ?? 0, events),
        'text/calendar');
    }],
    ['Animals', () => openTable('animals'), 'Data > Open as table'],
    ['Animals by day', () => {
      if (!byDay)
        throw new Error('Analyze a study first');
      const t = animalsByDayFrame(byDay);
      t.name = `${studyNow().df.name} - animals by day`;
      grok.shell.addTableView(t);
    }, 'Data > Open as table'],
    ['Measurements', () => openTable('measurements'), 'Data > Open as table'],
    ['Statistics', () => download(`${studyNow().df.name} - statistics.csv`, statsCsv, 'text/csv'), 'Export > CSV'],
    ['Statistics as shown', () => download(`${studyNow().df.name} - statistics as shown.csv`, statsShownCsv,
      'text/csv'), 'Export > CSV'],
    // the values as stored, not as the columns display them (19.26, not 19.3)
    ['Measurements', () => download(`${studyNow().df.name}.csv`, tableCsv(studyNow().df), 'text/csv'),
      'Export > CSV'],
    ['Per-animal results', () => download(`${studyNow().df.name} - animals.csv`, tableCsv(studyNow().r.animals, 6),
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
    const df = studyDf();
    // the study's own name, as Monitor's header (the analysis table's name says only what was compared)
    return {name: (df && studyNameOf(df)) || r.analysis.name, groups: [control, ...arms].map((name) => ({name, color: colors[name]})),
      status: [modelName(r.analysis.get('line', 0)),
        `evaluated on ${dayText(r.analysis.get('evaluation_day', 0)).toLowerCase()}`,
        analyzed.value === undefined ? '' : `${analyzed.value} animals analyzed`,
        twinned.value && isFiniteNumber(coverage.value) ? `controls on track ${pct(coverage.value)}` : '']
        .filter((x) => x).join(' · '), statusTip,
      meta: () => metaOf(df), setMeta: df ? (m) => setMetaOn(df, m) : undefined,
      groupObject: (g, color) => new PmxGroup(studyOf(df), g, color, {details: () =>
        detailsElement(detailsNow(g), groupSections(g, color)), live: liveNow()})};
    // Plan opens with the plan the analysis uses: the one picked, else the study's own
  }), (s) => void openStage(s, planOfStudy(planText(), metaOf(studyDf())?.id)), computed(() => ({monitor:
    !result.value && !canReturn('monitor') ? 'Needs a study: open one, or follow a plan' : null})));
  // Save (and Ctrl+S): what Analyze changes on the study's table (the patients' exposure, the effect kept, exclusions).
  // A table followed in Monitor is saved by that view; any other (from a file, a drop), with rows on, is saved here as
  // its row: the one it came from (its limits kept), else a new one, the first time it is saved
  // (the effect Analyze keeps on the exposure tag is its own reading, not a change of the user's: opening a saved
  // study left it unsaved, so it stayed open beside the next page; the table's text taken only once it changed)
  const asSaved = savedTables(self, (df) => df.toCsv(), (df) =>
    [JSON.stringify({...exposureOfTag(df) ?? {}, effect: undefined}), ...[META_TAG, LOOKS_TAG, PLAN_TAG]
      .map((t) => df.getTag(t) ?? '')].join('|'));
  const saveHere = async (df: DG.DataFrame): Promise<void> => {
    const at = rowOfTable(df);
    // (a row gone to the trash since: its limits unread, the save below offers "Save as new" as every page's does;
    // read first, it said "could not be saved" three times, QA 0.1.3 round 4 chaos N2, N11)
    const was = at ? await open('study', at.id).then((r) => r.doc, () => null) : null;
    const m = await withCode(metaOf(df) ?? emptyMeta());
    setMetaOn(df, m);
    const s = await stamped(savedStudyOf(df, was?.humane ?? null, was?.weightLimit ?? null));
    try {
      const saved = await save('study', s, at, {fileOwner: tableFileOwner(df)});
      setRowOn(df, saved);
      linkRow(self.view, 'study', saved.id);
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
        linkRow(self.view, 'study', saved.id);
        df.name = doc.name;
        setMetaOn(df, doc.meta!);
        head.refresh();
        asSaved.set(df);
        grok.shell.info(`Saved as your copy: ${esc(doc.name)}`); // (a toast reads HTML: the name is never markup)
      }, e);
    }
    asSaved.set(df);
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
        // (kept until dismissed: a save refused is read and acted on, QA 0.1.3 technician F12)
        grok.shell.error(saveError(e, 'study'), {autoHide: false});
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
    // and what it was read with, set again when it is reopened (restoreFields)
    const p = (name: string) => call.getParamValue(name);
    setAnalysisOn(df, {correction: correctionOf(p('correction')), endpoint: p('endpoint'),
      evaluationDay: p('evaluationDay')});
  };
  const saveButton = buttonWithIcon('Save study', 'save', () => void saver()?.().catch(() => {}));
  saveButton.title = 'Save the study with what Analyze added (the patients’ exposure, exclusions)';
  // Save as…, as Monitor's: the study under a name of its own, a new study (QA 0.1.3 round 3 designer U3)
  if (domainStorageOn()) {
    exports.push(['Save as…', () => {
      const df = studyDf();
      if (df)
        saveStudyAs(df, (m) => setMetaOn(df, m), () => (saver()?.() ?? Promise.resolve()).then(() => head.refresh()));
    }, 'File', () => studyDf() ? null : 'Open a study first']);
  }
  effect(() => {
    void result.value;
    const df = studyDf();
    saveButton.hidden = !saver();
    const id = rowOfTable(df)?.id;
    // (the view is made once the content is: the first run waits for it)
    queueMicrotask(() => linkRow(self.view, 'study', id));
    if (df && !asSaved.has(df))
      asSaved.set(df);
  });
  // (closing asks only for a table saved here: one followed in Monitor is asked about there)
  queueMicrotask(() => saveGate(self, 'study', () => saver()?.() ?? Promise.resolve(), () => {
    const df = studyDf();
    return !!df && domainStorageOn() && !tableSavers.has(df) && asSaved.changed(df);
  }, () => {
    const df = studyDf();
    return !!df && !rowOfTable(df);
  }, asSaved.changes));
  const links = stageLinks({stage: 'analyze'}, () => head.details());
  // (Compare is no stage: opened from here on this study, or from the Overview)
  links.study.push(['Compare with other studies…', () => void compareWith(rowOfTable(studyDf())?.id),
    'This study beside your other saved studies with measurements']);
  // another stage shown in this view's place: this study, to come back to
  keepOnLeaving(self, async () => {
    const df = studyDf();
    return df ? () => analyzeView({study: df}) : null;
  });
  queueMicrotask(() => stageMenu(self.view, exports, {open: () => void openStudyDialog(undefined, {stage: 'analyze'}),
    save: {text: () => saveButton.textContent?.trim() || 'Save study', run: () => saveButton.click(),
      why: () => saveButton.hidden ? 'Open a study first' : reasonOf(saveButton)}}, links, (d) => self.own(d)));
  const toolbox = stageToolbox(exports, links.study, links.help, panels.root);
  // its "On this page" only while the page shows its sections, as Monitor's (QA 0.1.3 round 4 chaos N12)
  effect(() => {
    const pane = toolbox.querySelector('.pmx-index')?.closest<HTMLElement>('.d4-accordion-pane');
    if (pane)
      pane.style.display = result.value ? '' : 'none';
  });
  return {
    root: takesFiles(stagePage(head.row, panes), {stage: 'analyze'}),
    ribbon: [[head.steps]],
    toolbox,
    path: '/analyze',
    status: computed(() => busy.value ? 'Analyzing the study…' : result.value ?
      studyDf()?.name || result.value.analysis.name : ''),
  };
}
