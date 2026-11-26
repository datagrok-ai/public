/* Compare: studies with measurements side by side, each treated group's effect against the patients' exposure.
   What each study kept is read as saved (its groups, its NCA and patients' exposure, and the analysis's TGI on its
   evaluation day, pmx.exposure); a study analyzed before that was kept gives 1 − T/C of the median volumes against its
   own control (exposure.ts medianTgi), said so. It describes what was measured: no direction, no cause, nothing
   beyond the groups. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BasicTable, RadioInput, Signal, computed, signal} from '@datagrok-libraries/u2';
import type {StatValue} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

import {HumanExposure, HumanMultiple, cvOk, drugMultiples, fuOk, humanMultiple, humansInto, medianTgi, multipleText,
  parseExposure, patientsRange, savedEffect, shownMultiple, sourceText} from './exposure';
import type {ArmExposure, CurveRow} from './exposure';
import {PreviousStudy, controlOf, readDayOf} from './study-effect';
import {armParts, detailsName, metaOf} from './study-meta';
import {describeArm, plainArm} from './arm-name';
import {previousOfSaved} from './previous-studies';
import {entryText, open, save, saveError} from './store';
import {EXPOSURE_TAG, exposureOfTag, parseStudy} from './study-file';
import {humanExposureDialog} from './human-exposure';
import type {LibraryEntry} from './library';
import {boundText, exposureSection, sectionElement} from './translation-report';
import {Exports, exportMenu, reportHtml, stageToolbox} from './report';
import {entered, esc, printSheet, sig} from './sheet';
import {initQuietly} from './band-chart';
import {PREDICTED_COLOR, groupColors} from './colors';
import {STANDS_ALONE, openStage} from './stages';
import {compareHelp} from './hub-view';
import {steps} from './study-header';
import {cardGrid, helpPage, sized, statCard} from './ui';
import type {Stage} from './ui';
import {BandChart} from './band-chart';
import {editLayout, panelGrid} from './panels';
import {modelName} from './priors';
import {atlasControls, controlsTo, doublingOf, slopeSide} from './analyze';
import {ENDPOINT, median} from './engine';
import type {Range} from './plan-checks';

const NO_MOUSE = 'the study kept no mouse exposure (NCA)';
/** Under a combination's multiples, and in the report: whose AUC they use. */
export const COMBO_AUC = 'each drug at its AUC in this combination group\'s mice';

/** One treated group of one study, as Compare shows it. */
export interface CompareRow {
  study: string, model: string, group: string,
  /** the drug (a combination's drugs joined by ' + ') and its key ('' for a combination: never placed) */
  drug: string, key: string, treatment: string, dose: string | null,
  /** a drug alone: its dose as a number and its unit (the dose chart); null: no dose read, or a combination */
  doseValue: number | null, doseUnit: string,
  /** the analysis's result for the group (Slows growth, No clear effect…), as the study kept it; null: not kept */
  result: string | null,
  /** TGI (%), on `day`; `analyzed`: the analysis's (1 − ΔT/ΔC), else 1 − T/C of the medians */
  tgi: number | null, day: number | null, analyzed: boolean,
  multiple: HumanMultiple | null,
  /** where the patients' exposure came from; null: not given */
  source: string | null, cv: number | null, unbound: boolean,
  /** the group's NCA exposure and its drug's patients' exposure, as kept (for the report) */
  x?: ArmExposure, human?: HumanExposure,
  /** a combination with NCA's exposure kept for its drugs: each drug against its patients' exposure (not charted) */
  drugs?: {drug: string, m: HumanMultiple | null, x?: ArmExposure, human?: HumanExposure}[],
}

/** A study's groups, in the order of its rows. */
const groupsIn = (s: PreviousStudy) => [...new Set(s.rows.map((r) => r.group).filter((g): g is string => !!g))];

/** A study's treated groups: each with the analysis's TGI when the study kept it (else 1 − T/C of the medians on
 *  the last day with 3 animals in it and its control), and its exposure as a multiple of the patients'. */
export function compareRows(s: PreviousStudy, model: string): CompareRow[] {
  const groups = groupsIn(s);
  const control = controlOf(groups, s.meta);
  const rows = s.rows.filter((r): r is CurveRow => r.group !== null && r.day !== null && r.volume !== null);
  const {humans, nca} = parseExposure(s.exposure);
  const effect = savedEffect(s.exposure);
  const ctd = s.exposure?.effect?.ctd;
  return groups.filter((g) => g !== control).map((g) => {
    const parts = armParts(g, s.meta);
    const one = parts.length === 1 ? parts[0] : null;
    const h = one ? humans.get(one.key) : undefined;
    const x = nca?.arms.get(g)?.exposure;
    const multiple = one && h && h.auc !== null ? nca && (x || nca.pk.get(g)) ?
      humanMultiple(g, x, nca.pk.get(g) ?? null, humans, s.meta) : {why: NO_MOUSE} : null;
    const analyzed = effect?.tgi.get(g);
    // a combination's source: each of its drugs that has one, named (its multiples are each drug's)
    const given = parts.flatMap((p) => humans.get(p.key) ? [`${p.drug}: ${sourceText(humans.get(p.key)!)}`] : []);
    return {study: s.name, model, group: g, drug: parts.map((p) => p.drug).join(' + ') || g, key: one?.key ?? '',
      treatment: describeArm(parts) || g, dose: one?.dose ? `${entered(one.dose.value)} ${one.dose.unit}` : null,
      doseValue: one?.dose?.value ?? null, doseUnit: one?.dose?.unit ?? '',
      result: ctd?.groups.find((x) => x.name === g)?.activity ?? null,
      ...analyzed !== undefined ? {tgi: analyzed, day: effect!.day, analyzed: true} :
        {tgi: control ? medianTgi(rows, g, control, null) : null, day: control ? readDayOf(s, g, control) : null,
          analyzed: false},
      multiple, source: one ? h ? sourceText(h) : null : given.length ? given.join('; ') : null,
      cv: h && cvOk(h.cv) ? h.cv : null,
      unbound: !!h && fuOk(h.fuMouse) && fuOk(h.fuHuman), x: nca?.arms.get(g), human: h,
      ...parts.length > 1 && nca?.drugs?.get(g) ?
        {drugs: drugMultiples(g, nca.drugs.get(g)!, humans, s.meta).map(({drug, key, m}) => ({drug, m,
          x: nca.drugs!.get(g)!.find((d) => d.key === key)?.x, human: humans.get(key)}))} : {}};
  });
}

/** The rows the chart shows: those with a TGI and a multiple of the patients' exposure. */
export const chartRows = (rows: CompareRow[]) =>
  rows.filter((r) => r.tgi !== null && Number.isFinite(r.tgi) && r.multiple?.ratio !== undefined);

/** The one line when the studies are of more than one model; null: one model. */
export const modelWarning = (models: string[]) => new Set(models.filter((m) => m)).size > 1 ?
  'Different models: tumour growth inhibition is not directly comparable across models.' : null;

/** One hint per drug no study gives the patients' exposure of (its first study, where it would be kept). */
export function hintRows(rows: CompareRow[]): {drug: string, key: string, study: string}[] {
  const given = new Set(rows.filter((r) => r.source).map((r) => r.key));
  const out = new Map<string, {drug: string, key: string, study: string}>();
  for (const r of rows) {
    if (r.key && !given.has(r.key) && !out.has(r.key))
      out.set(r.key, {drug: r.drug, key: r.key, study: r.study});
  }
  return [...out.values()];
}

/** What a TGI is, when it is not the analysis's (shown with an asterisk). */
export const tgiNote = (r: CompareRow) => r.analyzed ? null :
  `1 − T/C of medians${r.day !== null ? `, day ${r.day}` : ''} (not from the analysis)`;

/** The footnote of the TGIs marked with an asterisk; null: every TGI shown is the analysis's. */
export const tgiFootnote = (rows: CompareRow[]): string | null => rows.some((r) => !r.analyzed && r.tgi !== null) ?
  '* The study kept no analysis: 1 − T/C of the median volumes, on the last day with 3 animals in the group and ' +
  'its control.' : null;

const pct = (x: number) => `${Math.round(x)}%`;
const rangeWord = (m: HumanMultiple) => m.range === 'above' ? ', above the patients\' range' :
  m.range === 'below' ? ', below the patients\' range' : '';

/** One clause per drug on what its groups gave, from the numbers alone (the study named when the drug is in more
 *  than one); then the drugs not placed, and why. null: nothing to say. */
export function compareSentence(rows: CompareRow[]): string | null {
  const placed = chartRows(rows);
  const drugs = [...new Set(placed.map((r) => r.drug))];
  const said = drugs.map((d) => {
    const gs = placed.filter((r) => r.drug === d).sort((a, b) => a.multiple!.ratio! - b.multiple!.ratio!);
    const studies = new Set(rows.filter((r) => r.drug === d).map((r) => r.study)).size > 1;
    const each = gs.map((r) => `${pct(r.tgi!)}${r.dose ? ` at ${r.dose}` : ''}${studies ? ` in ${r.study}` : ''} ` +
      `(${multipleText(r.multiple!.ratio!)}${rangeWord(r.multiple!)}` +
      `${r.multiple!.weekly ? ', averaged over the week' : ''})`);
    const list = each.length > 1 ? `${each.slice(0, -1).join(', ')} and ${each[each.length - 1]}` : each[0];
    return `${d} gave a TGI of ${list}`;
  });
  const placedKeys = new Set(placed.map((r) => r.key));
  const noMouse = [...new Set(rows.filter((r) => r.key && !placedKeys.has(r.key) && r.multiple?.why === NO_MOUSE)
    .map((r) => r.drug))].map((d) => `${d} has no mouse exposure`);
  const missing = hintRows(rows).map((h) => `${h.drug} has no patients' exposure yet`);
  const out = [...said, ...noMouse, ...missing];
  return out.length ? `${out.join('; ')}.`.replace(/^./, (c) => c.toUpperCase()) : null;
}

/** What the exposures are and where they come from. */
export function compareFooter(rows: CompareRow[]): string {
  // (a combination's source, its drugs' each named, stays on its row: here it would repeat the single drugs')
  const sources = [...new Set(rows.filter((r) => r.source && r.key).map((r) => r.source!))];
  return `${boundText(rows.map((r) => r.human))} Mouse AUC from each study's NCA` +
    `${sources.length ? `; patients' exposure from ${sources.join('; ')}` : ''}.`;
}

/** A study file's text with the patients' exposure set to `humans`, all else as the file has it now (read again
 *  just before writing: the study may have been saved since Compare read it). */
export function withHumans(text: string, humans: Map<string, HumanExposure>): string {
  const raw = JSON.parse(text);
  raw.exposure = humansInto(raw.exposure, humans);
  return JSON.stringify(raw);
}

/** The study's tables open now (its ID and its name: a copy of the table is left alone) given the same patients'
 *  exposure as its file: a Save from Monitor then keeps it, and Analyze takes it in (analyze-view keepExposure). */
export function humansToOpen(id: string, name: string, humans: Map<string, HumanExposure>,
  tables: DG.DataFrame[] = grok.shell.tables): void {
  for (const df of id ? tables.filter((t) => metaOf(t)?.id === id && t.name === name) : [])
    df.setTag(EXPOSURE_TAG, JSON.stringify(humansInto(exposureOfTag(df), humans)));
}

/** The text of each entry: the library lists studies without reading them, so each is read now from its file; one
 *  that cannot be read (deleted, damaged) is skipped, by name. */
export async function compareTexts(entries: LibraryEntry[],
  read: (f: DG.FileInfo) => Promise<string> = (f) => grok.dapi.files.readAsText(f)):
Promise<{texts: Map<LibraryEntry, string>, skipped: string[]}> {
  const read1 = async (e: LibraryEntry): Promise<string | null> => {
    try {
      const text = e.text || (e.file ? await read(e.file) : e.row ? await entryText(e) : '');
      parseStudy(text);
      return text;
    } catch {
      return null;
    }
  };
  const texts = await Promise.all(entries.map(read1));
  return {texts: new Map(entries.flatMap((e, i): [LibraryEntry, string][] =>
    texts[i] === null ? [] : [[e, texts[i]!]])),
  skipped: entries.filter((_, i) => texts[i] === null).map((e) => e.name || e.id || 'A study')};
}

/** Each study's name as Compare shows it: with its ID when another study shown has the same name, and "running, day
 *  n" when it is still followed (its results are not final); unique, so the chart keys its series by it. */
export function studyLabels(entries: Pick<LibraryEntry, 'name' | 'id' | 'stage' | 'day'>[]): string[] {
  const count = new Map<string, number>();
  for (const e of entries)
    count.set(e.name, (count.get(e.name) ?? 0) + 1);
  const out = entries.map((e) => {
    const more = [...(count.get(e.name) ?? 0) > 1 && e.id ? [e.id] : [],
      ...e.stage === 2 ? [`running, day ${e.day ?? 0}`] : []];
    return `${e.name || e.id || 'A study'}${more.length ? ` (${more.join(', ')})` : ''}`;
  });
  // (the same name and no ID to tell them apart: numbered)
  return out.map((l, i) => out.indexOf(l) === i ? l : `${l} #${out.slice(0, i + 1).filter((x) => x === l).length}`);
}

/** The report's section for the compared studies: the sentence, the groups, the models' warning, what it does
 *  not say and the sources (translation-report.ts). */
export const compareSection = (rows: CompareRow[], models: string[]) => exposureSection({sentence:
  compareSentence(rows), warning: modelWarning(models), curve: false, rows: rows.flatMap((r) => {
  const common = {tgi: r.tgi, study: r.study, model: r.model, ...r.analyzed ? {} : {medians: {day: r.day}}};
  // a combination: a row per drug, "A + B (A)", with the multiple the table shows (as Analyze's report does)
  return r.drugs?.some((d) => d.m) ? r.drugs.map((d) => ({...common, group: `${r.treatment} (${d.drug})`, x: d.x,
    multiple: d.m, human: d.human})) : [{...common, group: r.treatment, x: r.x, multiple: r.multiple,
    human: r.human}];
})});

/** The chart: each group at its TGI against its multiple of the patients' exposure, on a log axis ticked at each
 *  doubling, one colour per study (its name in the legend, above the plot), each drug's patients' range shaded
 *  behind. echarts, as the band chart (the platform's scatter plot gave no log ticks at these few points). */
function compareChart(shown: CompareRow[], y: string): {root: HTMLElement, dispose: () => void} {
  const studies = [...new Set(shown.map((r) => r.study))];
  const ranges = [...new Set(shown.map((r) => r.cv))].map((cv) => ({cv, ...patientsRange(cv),
    drugs: [...new Set(shown.filter((r) => r.cv === cv).map((r) => r.drug))].join(', ')}));
  const xs = [...shown.map((r) => r.multiple!.ratio!), ...ranges.flatMap((b) => [b.lo, b.hi])];
  const fmt = (v: number) => `${sig(v, 2)}×`;
  const [lo, hi] = [Math.floor(Math.log2(Math.min(...xs))), Math.ceil(Math.log2(Math.max(...xs)))];
  // a point in the right third of the axis is labelled on its left: on its right the label ran off the chart
  const leftOf = (ratio: number) => (Math.log2(ratio) - lo) / Math.max(hi - lo, 1) > 2 / 3;
  return drawnChart(320, {
    animation: false,
    grid: {left: 56, right: 32, top: 40, bottom: 44},
    legend: {top: 0, data: studies},
    // (echarts reads the tooltip as HTML: the study's and the group's names are escaped)
    tooltip: {trigger: 'item', formatter: (p: any) => `${esc(p.seriesName)}<br>${esc(p.name)}: ` +
      `TGI ${Math.round(p.value[1])}% at ${fmt(p.value[0])} the patients' exposure`},
    xAxis: {type: 'log', logBase: 2, min: 2 ** lo, max: 2 ** hi, name: 'Multiple of the patients\' exposure',
      nameLocation: 'middle', nameGap: 28, axisLabel: {formatter: fmt}},
    yAxis: {type: 'value', name: y, nameLocation: 'middle', nameGap: 40},
    series: studies.map((s, i) => ({type: 'scatter', name: s, symbolSize: 10,
      color: DG.Color.toHtml(DG.Color.getCategoricalColor(i)),
      data: shown.filter((r) => r.study === s).map((r) => ({name: `${r.drug} ${r.dose ?? ''}`.trim(),
        value: [r.multiple!.ratio!, r.tgi!], ...leftOf(r.multiple!.ratio!) ? {label: {position: 'left'}} : {}})),
      label: {show: true, position: 'right', formatter: '{b}', fontSize: 10},
      // two groups at nearly one point (one dose in two studies): their labels stacked, not drawn over each other
      labelLayout: {moveOverlap: 'shiftY'},
      // the patients' ranges, behind the points (on the first series: drawn once)
      ...i === 0 ? {markArea: {silent: true, itemStyle: {color: 'rgba(154, 165, 177, 0.18)'},
        label: {position: 'insideTop', color: '#717581', fontSize: 10},
        data: ranges.map((b) => [{xAxis: b.lo, name: `${b.drugs}: patients' range` +
          (cvOk(b.cv) ? ` (CV ${sig(b.cv!)}%)` : ' (3-fold convention)')}, {xAxis: b.hi}])}} : {}})),
  });
}

/** An echarts chart of `height` drawn with `option`, following its panel's width; dispose lets it go. */
function drawnChart(height: number, option: {[k: string]: any}): {root: HTMLElement, dispose: () => void} {
  const root = ui.div([], {style: {height: `${height}px`, width: '100%'}});
  const chart = initQuietly(root);
  chart.setOption(option);
  const resized = ui.onSizeChanged(root).subscribe(() => chart.resize());
  return {root, dispose: () => {
    resized.unsubscribe();
    chart.dispose();
  }};
}

// ---- the growth of every group, the doses, and each control against the atlas ----

/** Every group's median volume and range (smallest to largest tumor) by day, a series per group (with its study when
 *  there are several), in its study's own colour; `relative`: each animal's volume over its first. */
export function growthFrame(studies: {s: PreviousStudy, label: string}[], relative: boolean):
  {df: DG.DataFrame, colors: {[series: string]: string}, dashed: string[]} | undefined {
  const out: {series: string, day: number, median: number, lower: number, upper: number, n: number}[] = [];
  const treatmentOf = new Map<string, string>();
  // each group in its study's own colour, as Analyze draws it (colors.ts groupColors: the study's own, else its slots)
  const colors: {[series: string]: string} = {};
  const dashed: string[] = [];
  studies.forEach(({s, label}, k) => {
    const groups = groupsIn(s);
    const control = controlOf(groups, s.meta);
    const own = groupColors(groups, control ?? '', s.meta?.colors, s.meta?.slots);
    const rows = s.rows.filter((r) => r.group !== null && r.day !== null && r.volume !== null && r.volume > 0)
      .sort((a, b) => a.day! - b.day!);
    const first = new Map<string, number>();
    for (const r of rows) {
      if (!first.has(r.animal))
        first.set(r.animal, r.volume!);
    }
    for (const g of groups) {
      const treatment = g === control ? 'Control' : describeArm(armParts(g, s.meta)) || g;
      const series = studies.length > 1 ? `${treatment} · ${label}` : treatment;
      treatmentOf.set(series, treatment);
      colors[series] = own[g];
      // each study after the first: its lines dashed (the colour is the treatment's)
      if (k > 0)
        dashed.push(series);
      const byDay = new Map<number, number[]>();
      for (const r of rows.filter((x) => x.group === g)) {
        const v = relative ? r.volume! / first.get(r.animal)! : r.volume!;
        byDay.set(r.day!, [...byDay.get(r.day!) ?? [], v]);
      }
      for (const [day, xs] of byDay)
        out.push({series, day, median: median(xs), lower: Math.min(...xs), upper: Math.max(...xs), n: xs.length});
    }
  });
  if (!out.length)
    return undefined;
  const df = DG.DataFrame.fromObjects(out)!;
  df.col('day')!.meta.units = 'd';
  for (const c of ['median', 'lower', 'upper']) {
    df.col(c)!.meta.units = relative ? '× day 0' : 'mm³';
    df.col(c)!.meta.format = relative ? '0.00' : '0.0';
  }
  df.col('series')!.meta.colors.setCategorical(colors);
  return {df, colors, dashed};
}

/** The treated groups of one drug alone, by dose (a drug in two units: two curves), the lowest dose first. */
export function doseGroups(rows: CompareRow[]): {name: string, rows: CompareRow[]}[] {
  const out = new Map<string, CompareRow[]>();
  for (const r of rows) {
    if (r.key && r.doseValue !== null && r.doseValue > 0 && r.tgi !== null && Number.isFinite(r.tgi)) {
      const name = `${r.drug} (${r.doseUnit})`;
      out.set(name, [...out.get(name) ?? [], r]);
    }
  }
  return [...out].map(([name, rs]) => ({name, rows: rs.sort((a, b) => a.doseValue! - b.doseValue!)}));
}

/** The chart of the doses: each drug's TGI at each dose, a colour per drug, each point named by its study. */
function doseChart(groups: {name: string, rows: CompareRow[]}[], several: boolean):
  {root: HTMLElement, dispose: () => void} {
  const doses = groups.flatMap((g) => g.rows.map((r) => r.doseValue!));
  const log = Math.max(...doses) / Math.min(...doses) >= 4;
  // a point in the right half of the axis is labelled on its left: on its right the label ran off the chart
  const at = (v: number) => log ? Math.log(v) : v;
  const [lo, hi] = [Math.min(...doses.map(at)), Math.max(...doses.map(at))];
  const leftOf = (v: number) => hi > lo && (at(v) - lo) / (hi - lo) > 1 / 2;
  // two studies at one dose with close effects (CT26 and MC38 at 3 mg/kg): one label above its point, one below,
  // not side by side on the same line
  const all = groups.flatMap((g) => g.rows);
  const near = (r: CompareRow) => all.filter((o) => o !== r && o.doseValue === r.doseValue &&
    Math.abs(o.tgi! - r.tgi!) < 12);
  // (an equal effect: the first drawn above)
  const above = (r: CompareRow, o: CompareRow) => o.tgi! < r.tgi! ||
    o.tgi === r.tgi && all.indexOf(o) > all.indexOf(r);
  const placed = (r: CompareRow) => near(r).length ? {label: {position: near(r).every((o) => above(r, o)) ? 'top' :
    'bottom'}} : leftOf(r.doseValue!) ? {label: {position: 'left'}} : {};
  return drawnChart(300, {
    animation: false,
    grid: {left: 56, right: 24, top: 40, bottom: 44},
    legend: {top: 0, data: groups.map((g) => g.name)},
    tooltip: {trigger: 'item', formatter: (p: any) => `${esc(p.seriesName)}<br>${esc(p.name)}: ` +
      `TGI ${Math.round(p.value[1])}% at ${sig(p.value[0])}`},
    // a log axis where the doses spread (a dose ladder), else a plain one
    xAxis: {type: log ? 'log' : 'value', name: 'Dose',
      nameLocation: 'middle', nameGap: 28, axisLabel: {formatter: (v: number) => sig(v)}},
    yAxis: {type: 'value', name: 'TGI (%)', nameLocation: 'middle', nameGap: 40},
    series: groups.map((g, i) => ({type: 'line', name: g.name, symbolSize: 9,
      color: DG.Color.toHtml(DG.Color.getCategoricalColor(i)),
      // (a line within one study's ladder only: across studies the points stand alone)
      lineStyle: {opacity: new Set(g.rows.map((r) => r.study)).size === 1 ? 1 : 0},
      data: g.rows.map((r) => ({name: r.study, value: [r.doseValue!, r.tgi!],
        ...placed(r)})),
      label: {show: several, position: 'right', formatter: '{b}', fontSize: 10},
      // (labels that still collide: the one drawn later hidden, its study in the tooltip)
      labelLayout: {moveOverlap: 'shiftY', hideOverlap: true}})),
  });
}

/** A study's control, read as Analyze reads it: its animals' mean growth rate of log volume (per day, least squares,
 *  each animal measured on 2 days or more) from day 0 to the evaluation day (the analysis's when it kept one, else the
 *  first day the median control reached `endpoint`, else the last day); its last day measured and its median volume
 *  on day 0 (null: none). */
export function controlGrowth(s: PreviousStudy, endpoint = ENDPOINT, evaluationDay: number | null = null):
  {slope: number | null, lastDay: number | null, toDay: number | null, vRand: number | null, n: number} {
  const control = controlOf(groupsIn(s), s.meta);
  const rows = s.rows.filter((r) => r.day !== null && r.volume !== null && r.volume > 0);
  const day0 = rows.filter((r) => r.day === 0).map((r) => r.volume!);
  const mine = rows.filter((r) => control !== null && r.group === control);
  const byDay = new Map<number, number[]>();
  for (const r of mine)
    byDay.set(r.day!, [...byDay.get(r.day!) ?? [], r.volume!]);
  const reached = [...byDay].sort((a, b) => a[0] - b[0]).find(([, xs]) => median(xs) >= endpoint)?.[0] ?? null;
  const lastDay = mine.length ? Math.max(...mine.map((r) => r.day!)) : null;
  const toDay = evaluationDay ?? reached ?? lastDay;
  // (before day 0 the tumor is not yet randomized; past the evaluation day a control slows near its limit)
  const used = mine.filter((r) => r.day! >= 0 && (toDay === null || r.day! <= toDay));
  const animals = [...new Set(mine.map((r) => r.animal))];
  const slopes = animals.flatMap((a) => {
    const pts = used.filter((r) => r.animal === a).map((r) => [r.day!, Math.log(r.volume!)]);
    if (new Set(pts.map((p) => p[0])).size < 2)
      return [];
    const [mx, my] = [pts.reduce((t, p) => t + p[0], 0) / pts.length, pts.reduce((t, p) => t + p[1], 0) / pts.length];
    return [pts.reduce((t, p) => t + (p[0] - mx) * (p[1] - my), 0) / pts.reduce((t, p) => t + (p[0] - mx) ** 2, 0)];
  });
  return {slope: slopes.length ? slopes.reduce((t, x) => t + x, 0) / slopes.length : null, lastDay, toDay,
    vRand: day0.length ? median(day0) : null, n: animals.length};
}

/** A control's doubling time against the atlas's range for its model (its growth rate's 95% range, as Analyze's grey
 *  band): within, or how many times faster or slower than the atlas's median; a control that did not grow is slower
 *  (its doubling and times null). */
export interface ControlCheck {doubling: number | null, lo: number | null, hi: number | null,
  side: 'within' | 'faster' | 'slower' | null, times: number | null}

const SIDES = {'As predicted': 'within', 'Faster than predicted': 'faster', 'Slower than predicted': 'slower'} as const;

export function controlVsAtlas(slope: number | null, r: Range | null): ControlCheck {
  if (!r || slope === null || !Number.isFinite(slope))
    return {doubling: doublingOf(slope), lo: null, hi: null, side: null, times: null};
  // the rate against the range first, as Analyze reads it: a control that did not grow is the furthest below
  const side = SIDES[slopeSide(slope, r)];
  const [doubling, mid] = [doublingOf(slope), doublingOf(r.median)];
  return {doubling, lo: doublingOf(r.q975), hi: doublingOf(r.q025), side,
    times: mid === null || doubling === null ? null : side === 'slower' ? doubling / mid : mid / doubling};
}

const days = (x: number | null) => x === null ? '—' : `${sig(x, 2)} d`;
/** A doubling time's atlas range in words ("2.1–3.4 d"; a range open at an end: "2.1 d or more", "up to 3.4 d"). */
export const rangeText = (c: ControlCheck) => c.lo === null && c.hi === null ? '—' : c.hi === null ?
  `${days(c.lo)} or more` : c.lo === null ? `up to ${days(c.hi)}` : `${sig(c.lo, 2)}–${days(c.hi)}`;
/** Where a control stands: "Within the atlas range", "1.6× faster", "2× slower" or "Did not grow". */
export const sideText = (c: ControlCheck) => c.side === null ? 'No atlas prediction' : c.side === 'within' ?
  'Within the atlas range' : c.times === null ? 'Did not grow' : `${sig(c.times, 2)}× ${c.side}`;

/** One plain sentence on the controls that grew outside the atlas's range; null: none did. */
export function driftSentence(rows: {study: string, model: string, check: ControlCheck}[]): string | null {
  const drift = rows.filter((r) => r.check.side === 'faster' || r.check.side === 'slower');
  if (!drift.length)
    return null;
  return drift.map((r) => `${r.study}: its control ${r.check.doubling === null ? 'did not grow, where' :
    `doubled every ${days(r.check.doubling)}, ${r.check.side} than`} the atlas predicts for ${modelName(r.model)} ` +
    `(${rangeText(r.check)})`).join('; ') +
    '. Its treated groups are read against a control unlike the atlas\'s: compare their TGI with care.';
}

/** A study's control against the atlas: the analysis's own reading when the study kept it (its rate and range, as
 *  Analyze's "Against the plan" shows them), else its animals' rate to the evaluation day against `atlas`. */
export function controlCheck(s: PreviousStudy, atlas: Range | null): {check: ControlCheck, fromAnalysis: boolean} {
  const kept = s.exposure?.effect?.control;
  if (kept && [kept.slope, kept.q025, kept.median].every((x) => typeof x === 'number' && Number.isFinite(x))) {
    return {check: controlVsAtlas(kept.slope, {q025: kept.q025!, median: kept.median!, q975: kept.q975 ?? null}),
      fromAnalysis: true};
  }
  const g = controlGrowth(s, kept?.endpoint ?? ENDPOINT, s.exposure?.effect?.day ?? null);
  return {check: controlVsAtlas(g.slope, atlas), fromAnalysis: false};
}

/** The Compare view of the selected studies with measurements (library entries, each read from its file: the
 *  library lists them without their text): panels in a grid the user can lay out, as each stage's. */
export async function compareView(entries: LibraryEntry[],
  read?: (f: DG.FileInfo) => Promise<string>): Promise<DG.View> {
  const {texts, skipped} = await compareTexts(entries, read);
  const shownEntries = entries.filter((e) => texts.has(e));
  const view = DG.View.create();
  view.name = 'PMX Studio · Compare';
  // a page of the app, quiet as the Overview is (stages.ts quietPanels)
  view.root.dataset[STANDS_ALONE] = 'true';
  view.helpUrl = helpPage();
  const host = ViewHost.of(view);
  // as every page of the app: the stages to go to (none of them current: Compare is no stage), the report in the
  // view's menu and the toolbox
  view.setRibbonPanels([[ui.div([steps('compare' as Stage, (s) => void openStage(s), (d) => host.own(d))],
    'pmx-head')]]);
  let report = () => {};
  const exports: Exports = [['Report (PDF)…', () => report()]];
  const muted = (t: string) => {
    const el = ui.divText(t);
    el.style.color = 'var(--grey-4)';
    return el;
  };
  const labels = studyLabels(shownEntries);
  const studiesNow = () => shownEntries.map((e, i) => ({e, s: previousOfSaved(parseStudy(texts.get(e)!)),
    label: labels[i]}));
  const models = shownEntries.map((e) => e.line);
  const warning = modelWarning(models);

  // (i) the summary: what is compared, filled as each part is read
  const nStudies = signal<StatValue | undefined>(shownEntries.length);
  const nTreated = signal<StatValue | undefined>(undefined);
  const nPlaced = signal<StatValue | undefined>(undefined);
  const nWithin = signal<StatValue | undefined>(undefined);
  // figures no result follows: no change shown, nothing to pin
  const card = (label: string, value: Signal<StatValue | undefined>, tip: () => string) =>
    statCard(host, label, value, tip, undefined, false);
  const cards = cardGrid([
    card('Studies', nStudies, () => `${shownEntries.length} studies, of ` +
      `${[...new Set(models.filter((m) => m).map(modelName))].join(', ') || 'no model given'}`),
    card('Treated groups', nTreated, () => 'Every treated group of the studies compared'),
    card('Against the patients’ exposure', nPlaced, () => 'Treated groups with both a TGI and a ' +
      'multiple of the patients’ exposure: the groups on the exposure chart'),
    card('Controls within the atlas', nWithin, () => 'The studies whose control grew within the atlas’s ' +
      'range for its model'),
  ], 'pmx-compare-cards', 140, false);

  // (ii) the growth of every group: Linear, Log or each animal over its first volume
  const axis = host.build(() => new RadioInput({items: ['Linear', 'Log', 'Relative V/V0'], buttons: true,
    value: 'Linear',
    itemTooltips: {'Linear': 'Linear volume axis', 'Log': 'Logarithmic volume axis',
      'Relative V/V0': 'Each animal’s volume over its first: studies started at different sizes, side by side'}}));
  const growth = growthFrame(studiesNow(), false);
  const relative = growthFrame(studiesNow(), true);
  const growthData = computed(() => (axis.value === 'Relative V/V0' ? relative : growth)?.df);
  const growthChart = growth ? host.build(() => new BandChart(growthData, {x: 'day', lower: 'lower', upper: 'upper',
    line: 'median', split: 'series', xTitle: 'Day after randomization', yTitle: 'Tumor volume', markers: true,
    log: computed(() => axis.value === 'Log'), dashed: growth.dashed, count: 'n'})) : null;
  const axisRow = ui.divH([ui.span(['Volume axis']), axis.root], 'pmx-noprint');
  Object.assign(axisRow.style, {justifyContent: 'flex-end', alignItems: 'center', gap: 'var(--dg-space-m)'});

  // (v) the controls against the atlas: filled once the atlas has predicted each study's
  const controlsBox = ui.divV([muted('Reading the atlas for each study’s model…')]);
  const exposureBox = ui.divV([]);
  const doseBox = ui.divV([]);
  const groupsBox = ui.divV([]);

  const panels = host.build(() => panelGrid([
    {id: 'summary', title: 'The studies', hint: 'What is compared', content: [...warning ? [muted(warning)] : [],
      ...skipped.map((n) => muted(`${n} could not be read (moved, deleted or damaged): left out.`)), cards]},
    {id: 'growth', title: 'Tumor growth', hint: 'Median and range (smallest to largest tumor) of each group, in its ' +
      'treatment’s colour; each study after the first dashed', content: growth ? [axisRow, sized(growthChart!, 360)] :
      [muted('No volumes to draw.')]},
    {id: 'dose', title: 'Dose and effect', hint: 'Each drug’s TGI at each dose it was given, across the studies',
      content: [doseBox]},
    {id: 'exposure', title: 'Exposure and effect', hint: 'Each treated group’s TGI against its ' +
      'multiple of the patients’ exposure; the patients’ range shaded', content: [exposureBox]},
    {id: 'controls', title: 'Controls against the atlas', hint: 'Each study’s median control over the atlas’s ' +
      'prediction for its model (its 95% range, as Analyze’s grey band)', content: [controlsBox]},
    {id: 'groups', title: 'Treated groups', hint: 'Every treated group: its TGI, its exposure and its result',
      content: [groupsBox]},
  ], [{id: 'summary', span: 12}, {id: 'growth', span: 6}, {id: 'dose', span: 6}, {id: 'exposure', span: 12},
    {id: 'controls', span: 12}, {id: 'groups', span: 12}], 'pmx-compare', 'compare'));
  exports.push(['Layout…', () => editLayout(panels), 'View']);
  exportMenu(view, exports);
  view.toolbox = stageToolbox(exports, [['Back to the Overview', () => void openStage('home'),
    'Your studies, to pick others to compare']], compareHelp());
  const page = ui.divV([panels], {style: {padding: '16px', overflow: 'auto'}});
  view.root.replaceChildren(page);

  // the charts drawn last, let go when drawn again or when the view closes
  let disposeCharts = () => {};
  host.own(() => disposeCharts());
  // (vi) and (iii), (iv): drawn again when a patients' exposure is added
  const draw = () => {
    const studies = studiesNow();
    // each study by its label (unique: two studies of one name told apart by their IDs), in the charts' series, the
    // table, the sentence and the report
    const rows = studies.flatMap(({e, s, label}) => compareRows(s, e.line).map((r) => ({...r, study: label})));
    const sentence = compareSentence(rows);
    const shown = chartRows(rows);
    nTreated.value = rows.length;
    nPlaced.value = `${shown.length} of ${rows.length}`;
    disposeCharts();
    const drawn: {root: HTMLElement, dispose: () => void}[] = [];
    if (shown.length) {
      const y = shown.every((r) => r.analyzed) ? 'TGI from the analysis (%)' : shown.some((r) => r.analyzed) ?
        'TGI from the analysis, * 1 − T/C of medians (%)' : 'TGI, * 1 − T/C of medians (%)';
      drawn.push(compareChart(shown, y));
    }
    const doses = doseGroups(rows);
    if (doses.length)
      drawn.push(doseChart(doses, studies.length > 1));
    disposeCharts = () => drawn.forEach((d) => d.dispose());
    const exposureChart = shown.length ? drawn[0].root : null;
    doseBox.replaceChildren(doses.length ? drawn[drawn.length - 1].root :
      muted('No drug alone with its dose: a dose is read from the group’s details, or its name.'));

    // built under the view's host: draw() also runs after an await (a patients' exposure added), when nothing builds
    // the page's width shared by the columns, long text wrapped in its cell: a long source would widen the table past
    // the window
    const wrap = (t: string) => ui.div([t], {style: {whiteSpace: 'normal', overflowWrap: 'anywhere'}});
    const table = host.build(() => new BasicTable<CompareRow>({items: rows, fixedLayout: true,
      keyOf: (r) => `${r.study}|${r.group}`, columns: [
        {header: 'Study', render: (r) => wrap(r.study)},
        {header: 'Model', render: (r) => r.model ? modelName(r.model) : '—'},
        // in plain words, the codes on hover (as the group's panel says it)
        {header: 'Group', render: (r) => {
          const el = wrap(plainArm(r.treatment));
          if (el.textContent !== r.treatment)
            ui.tooltip.bind(el, r.treatment);
          return el;
        }},
        {header: 'TGI', render: (r) => {
          if (r.tgi === null)
            return '—';
          const el = ui.span([`${pct(r.tgi)}${r.analyzed ? '' : ' *'}`]);
          const note = tgiNote(r);
          if (note)
            ui.tooltip.bind(el, note);
          return el;
        }},
        {header: 'Day', render: (r) => r.day === null ? '—' : String(r.day)},
        {header: '× patients’ exposure', render: (r) => {
          // a combination: each drug's multiple, its name first
          if (r.drugs?.some((d) => d.m)) {
            // each drug at the AUC NCA read in this group's own mice: not the same as the drug's group alone (1.9×
            // here, 2.2× there is two groups' exposure, not two computations)
            return ui.divV([...r.drugs.map((d) => {
              const el = ui.divText(`${d.drug} ${d.m?.ratio !== undefined ? shownMultiple(d.m) : '—'}`);
              if (d.m)
                ui.tooltip.bind(el, d.m.ratio !== undefined ? d.m.tip : d.m.why);
              return el;
            }), muted(COMBO_AUC)]);
          }
          const m = r.multiple;
          if (m?.ratio !== undefined) {
            const el = ui.span([shownMultiple(m)]);
            ui.tooltip.bind(el, m.tip);
            return el;
          }
          return ui.divV([ui.span(['—']), muted(m?.why ?? (r.key ? 'no patients\' exposure' : 'a combination'))]);
        }},
        {header: 'Result', render: (r) => r.result ?? '—'},
        {header: 'Patients\' exposure from', render: (r) => wrap(r.source ? `${r.source}${r.cv !== null ?
          ` · CV ${sig(r.cv)}%` : ''}` : 'Not given')},
      ]}));
    // one hint per drug without the patients' exposure: typed here, kept with its study
    const hints = hintRows(rows).map((hint) => {
      const add = ui.button('Add patients\' exposure…', () => {
        const {e, s} = studies.find((x) => x.label === hint.study)!;
        const {file, row} = e;
        if (!file && !row)
          return grok.shell.warning('This study is not saved in your files: open it in Analyze to keep it');
        const names = groupsIn(s).map((g) => detailsName(g, s.meta));
        // the study as it is now, changed and written back: its file, or its row at the version just read
        const update = async (change: (now: string) => string): Promise<string> => {
          if (file) {
            const text = change(await grok.dapi.files.readAsText(file));
            await grok.dapi.files.writeAsText(file, text);
            return text;
          }
          const {doc, saved} = await open('study', row!.id);
          const text = change(JSON.stringify(doc));
          await save('study', parseStudy(text), saved);
          return text;
        };
        void humanExposureDialog(names, parseExposure(s.exposure).humans, (m) => {
          void update((now) => withHumans(now, m)).then((text) => {
            texts.set(e, text);
            humansToOpen(e.id, e.name, m);
            draw();
          }).catch((err) => grok.shell.error(row ? saveError(err, 'study') :
            `The patients' exposure could not be saved with ${esc(hint.study)}: ${esc(err)}`));
        }).catch((err) => grok.shell.error(String(err)));
      });
      return ui.divH([muted(`${hint.drug}: add the patients' exposure (a projection is fine) to place it on the ` +
        'chart.'), add], {style: {alignItems: 'center', gap: '8px'}});
    });
    const footnote = tgiFootnote(rows);
    exposureBox.replaceChildren(...sentence ? [ui.divText(sentence)] : [],
      ...exposureChart ? [exposureChart] : [muted('No treated group has both a TGI and a multiple of the ' +
        'patients’ exposure yet.')], ...hints, muted(compareFooter(rows)));
    groupsBox.replaceChildren(table.root, ...footnote ? [muted(footnote)] : []);
    // the report: the page as drawn, then the section
    report = () => {
      const section = compareSection(rows, models);
      void reportHtml({title: 'Compare', stage: 'Compare', meta: null, content: page, annex: section ?
        [['Exposure and patients', sectionElement(section)]] : []}).then(printSheet,
        (err) => grok.shell.error(`The report could not be built: ${err}`));
    };
  };
  draw();
  void drawControls(studiesNow(), controlsBox, host, nWithin, muted);
  return view;
}

/** Each study's control against the atlas's prediction for its model: the chart, the table and the sentence. */
async function drawControls(studies: {e: LibraryEntry, s: PreviousStudy, label: string}[], box: HTMLElement,
  host: ViewHost, nWithin: Signal<StatValue | undefined>, muted: (t: string) => HTMLElement): Promise<void> {
  const read = await Promise.all(studies.map(async ({e, s, label}) => {
    // the study's own endpoint when its analysis kept it (the atlas's evaluation day is where its median reaches it)
    const endpoint = s.exposure?.effect?.control?.endpoint ?? ENDPOINT;
    const g = controlGrowth(s, endpoint, s.exposure?.effect?.day ?? null);
    // a study with no model or no control measured has no prediction: said on its row
    const atlas = e.line && g.lastDay !== null ? await atlasControls(e.line, g.vRand ?? 100, g.lastDay,
      Math.max(2, g.n), endpoint).catch(() => null) : null;
    return {label, model: e.line, start: s.meta?.start ?? null, g, atlas,
      ...controlCheck(s, atlas?.controlSlope ?? null), s};
  }));
  // (the view closed while the atlas predicted: nothing to draw)
  if (host.scope.isDisposed)
    return;
  nWithin.value = `${read.filter((r) => r.check.side === 'within').length} of ${read.length}`;
  const rows: object[] = [];
  read.forEach((r, i) => {
    const color = DG.Color.toHtml(DG.Color.getCategoricalColor(i));
    if (r.atlas && r.g.lastDay !== null) {
      for (const b of controlsTo(r.atlas.controls, r.g.lastDay)) {
        // the atlas in the predicted grey, as Analyze draws its predicted control (dashed: told apart by study)
        rows.push({series: `${r.label}: atlas`, color: PREDICTED_COLOR, day: b.day, lower: b.q025, upper: b.q975,
          median: b.median});
      }
    }
    const control = controlOf(groupsIn(r.s), r.s.meta);
    const byDay = new Map<number, number[]>();
    for (const x of r.s.rows.filter((x) => control !== null && x.group === control && x.day !== null &&
      x.volume !== null))
      byDay.set(x.day!, [...byDay.get(x.day!) ?? [], x.volume!]);
    for (const [day, xs] of [...byDay].sort((a, b) => a[0] - b[0])) {
      const m = median(xs);
      rows.push({series: `${r.label}: control`, color, day, lower: m, upper: m, median: m});
    }
  });
  const sentence = driftSentence(read.map((r) => ({study: r.label, model: r.model, check: r.check})));
  const noDay0 = read.filter((r) => r.atlas && !r.fromAnalysis && r.g.vRand === null).map((r) => r.label);
  const table = host.build(() => new BasicTable<typeof read[number]>({items: read, fixedLayout: true,
    keyOf: (r) => r.label, columns: [
      {header: 'Study', render: (r) => r.label},
      {header: 'Day 0', render: (r) => r.start ?? '—'},
      {header: 'Model', render: (r) => r.model ? modelName(r.model) : '—'},
      {header: 'Doubling time', render: (r) => days(r.check.doubling)},
      {header: 'Atlas range', render: (r) => rangeText(r.check)},
      {header: 'Against the atlas', render: (r) => sideText(r.check)},
    ]}));
  let chart: HTMLElement[] = [];
  if (rows.length) {
    const df = DG.DataFrame.fromObjects(rows.map(({color: _, ...x}: any) => x))!;
    df.col('day')!.meta.units = 'd';
    for (const c of ['median', 'lower', 'upper'])
      df.col(c)!.meta.units = 'mm³';
    df.col('series')!.meta.colors.setCategorical(Object.fromEntries(rows.map((x: any) => [x.series, x.color])));
    const atlasSeries = [...new Set(rows.map((x: any) => x.series as string))].filter((x) => x.endsWith(': atlas'));
    const band = host.build(() => new BandChart(signal(df), {x: 'day', lower: 'lower', upper: 'upper', line: 'median',
      split: 'series', xTitle: 'Day after randomization', yTitle: 'Median control volume', markers: true,
      dashed: atlasSeries}));
    chart = [sized(band, 320)];
  }
  box.replaceChildren(...sentence ? [ui.divText(sentence)] : [], ...chart, table.root,
    muted('Doubling time: of the control’s mean growth rate of log volume from day 0 to the evaluation day; the ' +
      'analysis’s own rate and range when the study kept them (as Analyze’s Against the plan shows them), else ' +
      'the control animals’ against the atlas: 95 of 100 studies of the model, randomized at the study’s day-0 ' +
      'median' + (noDay0.length ? ` (${noDay0.join(', ')} measured no day 0: its prediction starts at 100 mm³)` :
      '') + '.'));
}
