/* Compare: studies with measurements side by side, each treated group's effect against the patients' exposure.
   What each study kept is read as saved (its groups, its NCA and patients' exposure, and the analysis's TGI on its
   evaluation day, pmx.exposure); a study analyzed before that was kept gives the same TGI from its rows
   (study-effect.ts growthTgi). It describes what was measured: no direction, no cause, nothing beyond the groups. */
import {maxOf, minOf} from './format';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {BasicTable, Signal, computed, dot, signal} from '@datagrok-libraries/u2';
import type {StatValue} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

import {HumanExposure, HumanMultiple, cvOk, drugMultiples, fuOk, humanMultiple, humansInto, multipleText,
  parseExposure, patientsRange, savedEffect, shownMultiple, sourceText} from './exposure';
import type {ArmExposure} from './exposure';
import {PreviousStudy, controlOf, growthTgi, readDayOf} from './study-effect';
import {armParts, detailsName, metaOf} from './study-meta';
import {describeArm, plainArm} from './arm-name';
import {previousOfSaved} from './previous-studies';
import {entryText, open, save, saveError} from './store';
import {EXPOSURE_TAG, exposureOfTag, parseStudy} from './study-file';
import {humanExposureDialog} from './human-exposure';
import type {LibraryEntry} from './library';
import {boundText, exposureSection, sectionElement} from './translation-report';
import {Exports, reportHtml, stageMenu, stageToolbox} from './report';
import type {ToolLink} from './report';
import {entered, esc, printSheet, sig} from './sheet';
import {initQuietly} from './band-chart';
import {PREDICTED_COLOR, groupColors} from './colors';
import {STANDS_ALONE, keepOnLeaving, openStage} from './stages';
import {compareHelp, openStudyDialog} from './hub-view';
import {crumb, headBar} from './study-header';
import {axisSwitch, helpPage, noteLine, sized, swatch} from './ui';
import {BandChart} from './band-chart';
import {customizePage, panelGrid} from './panels';
import {modelName} from './priors';
import {atlasControls, controlsTo, doublingOf, slopeSide} from './analyze';
import {ENDPOINT, median} from './engine';
import {likelyTypos, typoKey} from './monitor';
import {findOutliers, rulesInForce} from './outlier-rules';
import {groupsByAnimal} from './animal-group';
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
  /** TGI (%, 1 − ΔT/ΔC), on `day`; `analyzed`: the analysis's, else computed from the study's rows (growthTgi) */
  tgi: number | null, day: number | null, analyzed: boolean,
  /** a TGI read from the rows (a running study's): the likely typos it read as they are, in the group and its
   *  control (Monitor's rule and its warning: QA 0.1.3 A7) */
  toCheck?: number,
  multiple: HumanMultiple | null,
  /** where the patients' exposure came from; null: not given */
  source: string | null, cv: number | null, unbound: boolean,
  /** the group's NCA exposure and its drug's patients' exposure, as kept (for the report) */
  x?: ArmExposure, human?: HumanExposure,
  /** a combination with NCA's exposure kept for its drugs: each drug against its patients' exposure (not charted) */
  drugs?: {drug: string, m: HumanMultiple | null, x?: ArmExposure, human?: HumanExposure}[],
}

/** The measurements Compare's growth and controls read: a volume over 0 on a day, neither flagged nor a likely typo
 *  (monitor.ts likelyTypos, Monitor's rule: a ×10 slip spiked its group's range). */
export function trustedRows(s: PreviousStudy): PreviousStudy['rows'] {
  const rows = s.rows.filter((r) => r.day !== null && r.volume !== null && r.volume > 0 && !r.flagged);
  const typos = likelyTypos(rows.map((r) => ({animal: r.animal, day: r.day!, volume: r.volume!})));
  return rows.filter((r) => !typos.has(typoKey(r.animal, r.day!)));
}

/** Each group's volumes to check: Monitor's likely typos (monitor.ts likelyTypos) and what the study's own "values to
 *  check" find (the outlier rules in force, as Analyze and the study table count them: a running study's TGI read on
 *  values its page asked to check went unmarked, QA 0.1.3 round 3 analyst 5); a value once each. */
export function typosByGroup(s: PreviousStudy): Map<string, number> {
  const rows = s.rows.filter((r) => r.day !== null && r.volume !== null && r.volume > 0 && !r.flagged);
  const typos = likelyTypos(rows.map((r) => ({animal: r.animal, day: r.day!, volume: r.volume!})));
  const animal = rows.map((r) => r.animal);
  const found = new Set(findOutliers({animal, group: groupsByAnimal(animal, rows.map((r) => r.group)),
    day: rows.map((r) => r.day), volume: rows.map((r) => r.volume), weight: null}, rulesInForce().rules)
    .filter((f) => f.on === 'volume').map((f) => f.row));
  const out = new Map<string, number>();
  rows.forEach((r, i) => {
    if (r.group && (typos.has(typoKey(r.animal, r.day!)) || found.has(i)))
      out.set(r.group, (out.get(r.group) ?? 0) + 1);
  });
  return out;
}

/** Whether a study has a tumor measured (a new study's empty sheet has none). */
export const measured = (s: PreviousStudy) => s.rows.some((r) => r.volume !== null && Number.isFinite(r.volume));

/** A study's groups, in the order of its rows. */
const groupsIn = (s: PreviousStudy) => [...new Set(s.rows.map((r) => r.group).filter((g): g is string => !!g))];

/** A study's treated groups: each with the analysis's TGI when the study kept it (else the same TGI from its rows on
 *  the last day with 3 animals in it and its control), and its exposure as a multiple of the patients'. */
export function compareRows(s: PreviousStudy, model: string): CompareRow[] {
  const groups = groupsIn(s);
  const control = controlOf(groups, s.meta);
  const day = (g: string) => control ? readDayOf(s, g, control) : null;
  const {humans, nca} = parseExposure(s.exposure);
  const effect = savedEffect(s.exposure);
  const ctd = s.exposure?.effect?.ctd;
  const typos = typosByGroup(s);
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
        {tgi: control ? growthTgi(s.rows, g, control, day(g)) : null, day: day(g), analyzed: false,
          toCheck: (typos.get(g) ?? 0) + (control ? typos.get(control) ?? 0 : 0)},
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
  'Different models: tumor growth inhibition is not directly comparable across models.' : null;

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

/** Where a TGI comes from, on its cell: the study, the day and the definition (said once). */
export const tgiSource = (r: CompareRow) => r.analyzed ? `${r.study}${r.day !== null ? `, day ${r.day}` : ''}: ` +
  'TGI = 1 − ΔT/ΔC, each animal\'s change from day 0, from the study\'s analysis' :
  // (round 4, analyst W16)
  `${r.study}: ${r.day !== null ? `day ${r.day} ` : ''}TGI, computed from the raw measurements (the study has no ` +
  'saved analysis)';

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
  const common = {tgi: r.tgi, study: r.study, model: r.model};
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
  const [lo, hi] = [Math.floor(Math.log2(minOf(xs))), Math.ceil(Math.log2(maxOf(xs)))];
  // a point in the right third of the axis is labelled on its left: on its right the label ran off the chart
  const leftOf = (ratio: number) => (Math.log2(ratio) - lo) / Math.max(hi - lo, 1) > 2 / 3;
  return drawnChart(320, {
    animation: false,
    grid: {left: 56, right: 32, top: 40, bottom: 44},
    // (many studies: one scrolled line, never over the plot: chaos round 3 #22)
    legend: {top: 0, type: 'scroll', data: studies},
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

/** Each study's groups' colours, one per study and group. One study: as Analyze draws them (colors.ts groupColors:
 *  the study's own, else its slots). Several: a hue per study (45° apart, the second round of 8 between them) and a
 *  shade per group, the control darkest; the palette's pairs of near colours (a cyan, two light greens) no longer
 *  meet across studies. */
export function compareColors(studies: {s: PreviousStudy}[]): {[group: string]: string}[] {
  return studies.map(({s}, k) => {
    const groups = groupsIn(s);
    const control = controlOf(groups, s.meta) ?? '';
    if (studies.length === 1)
      return groupColors(groups, control, s.meta?.colors, s.meta?.slots);
    const ordered = [...groups.filter((g) => g === control), ...groups.filter((g) => g !== control)];
    const hue = (210 + 45 * (k % 8) + 22.5 * Math.floor(k / 8)) % 360;
    return Object.fromEntries(ordered.map((g, j) =>
      [g, hslHex(hue, 65, ordered.length < 2 ? 45 : 28 + 47 * j / (ordered.length - 1))]));
  });
}

/** An HSL colour (degrees, percents) as #rrggbb. */
function hslHex(h: number, s: number, l: number): string {
  const [S, L] = [s / 100, l / 100];
  const f = (n: number) => {
    const k = (n + h / 30) % 12;
    const c = L - S * Math.min(L, 1 - L) * Math.max(-1, Math.min(k - 3, 9 - k, 1));
    return Math.round(255 * c).toString(16).padStart(2, '0');
  };
  return `#${f(0)}${f(8)}${f(4)}`;
}

/** Every group's median volume and range (smallest to largest tumor) by day, a series per group (with its study when
 *  there are several), each study and group a colour of its own (compareColors); `relative`: each animal's volume
 *  over its first. A flagged measurement is left out (a typo ×10 spiked its group's range). */
export function growthFrame(studies: {s: PreviousStudy, label: string}[], relative: boolean):
  {df: DG.DataFrame, colors: {[series: string]: string}, symbols: {[series: string]: string}} | undefined {
  const out: {series: string, day: number, median: number, lower: number, upper: number, n: number}[] = [];
  const treatmentOf = new Map<string, string>();
  const colors: {[series: string]: string} = {};
  const symbols: {[series: string]: string} = {};
  const studyColors = compareColors(studies);
  studies.forEach(({s, label}, k) => {
    const groups = groupsIn(s);
    const control = controlOf(groups, s.meta);
    const own = studyColors[k];
    const rows = trustedRows(s).filter((r) => r.group !== null).sort((a, b) => a.day! - b.day!);
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
      // each study its mark, as the dose chart's (every study after the first "dashed" told none apart: round 3,
      // analyst 13); a line dashed stays a stretch of few animals
      if (studies.length > 1)
        symbols[series] = STUDY_MARKS[k % STUDY_MARKS.length][0];
      const byDay = new Map<number, number[]>();
      for (const r of rows.filter((x) => x.group === g)) {
        const v = relative ? r.volume! / first.get(r.animal)! : r.volume!;
        byDay.set(r.day!, [...byDay.get(r.day!) ?? [], v]);
      }
      for (const [day, xs] of byDay)
        out.push({series, day, median: median(xs), lower: minOf(xs), upper: maxOf(xs), n: xs.length});
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
  return {df, colors, symbols};
}

/** The treated groups of one drug alone, by dose (a drug in two units: two curves), the lowest dose first; the
 *  studies in `leaveOut` (the running ones, whose effect is not final) left out. */
export function doseGroups(rows: CompareRow[], leaveOut: ReadonlySet<string> = new Set()):
  {name: string, rows: CompareRow[]}[] {
  const out = new Map<string, CompareRow[]>();
  for (const r of rows) {
    if (!leaveOut.has(r.study) && r.key && r.doseValue !== null && r.doseValue > 0 && r.tgi !== null &&
      Number.isFinite(r.tgi)) {
      const name = `${r.drug} (${r.doseUnit})`;
      out.set(name, [...out.get(name) ?? [], r]);
    }
  }
  return [...out].map(([name, rs]) => ({name, rows: rs.sort((a, b) => a.doseValue! - b.doseValue!)}));
}

/** The marks that tell the studies apart in the dose chart, the first study's first (8: unique to 8 studies). */
export const STUDY_MARKS: [string, string][] = [['circle', '●'], ['diamond', '◆'], ['triangle', '▲'], ['rect', '■'],
  ['emptyCircle', '○'], ['emptyDiamond', '◇'], ['emptyTriangle', '△'], ['emptyRect', '□']];

/** The chart of the doses: each drug's TGI at each dose, a colour per drug (each study's ladder a line), a mark per
 *  study; the key above the plot (never over it), every point's study, group and TGI on hover. */
function doseChart(groups: {name: string, rows: CompareRow[]}[], several: boolean):
  {root: HTMLElement, dispose: () => void} {
  const doses = groups.flatMap((g) => g.rows.map((r) => r.doseValue!));
  const log = maxOf(doses) / minOf(doses) >= 4;
  const studies = [...new Set(groups.flatMap((g) => g.rows.map((r) => r.study)))];
  const mark = (study: string) => STUDY_MARKS[studies.indexOf(study) % STUDY_MARKS.length];
  const colour = (i: number) => DG.Color.toHtml(DG.Color.getCategoricalColor(i));
  const units = [...new Set(groups.map((g) => g.rows[0].doseUnit))];
  const chart = drawnChart(300, {
    animation: false,
    grid: {left: 56, right: 24, top: 12, bottom: 48},
    tooltip: {trigger: 'item', formatter: (p: any) => `${esc(p.data.study)}<br>${esc(p.seriesName)} at ` +
      `${sig(p.value[0])}: TGI ${Math.round(p.value[1])}%`},
    // a log axis where the doses spread (a dose ladder), else a plain one
    xAxis: {type: log ? 'log' : 'value', name: units.length === 1 ? `Dose (${units[0]})` : 'Dose',
      nameLocation: 'middle', nameGap: 30, axisLabel: {formatter: (v: number) => sig(v)}},
    yAxis: {type: 'value', name: 'TGI (%)', nameLocation: 'middle', nameGap: 40},
    series: groups.flatMap((g, i) => studies.map((study) => ({study, rows: g.rows.filter((r) => r.study === study)}))
      .filter((x) => x.rows.length).map(({study, rows}) => ({type: 'line', name: g.name, color: colour(i),
        symbol: mark(study)[0], symbolSize: 10, lineStyle: {width: 2},
        data: rows.map((r) => ({study, value: [r.doseValue!, r.tgi!]}))}))),
  });
  // the key: each drug's colour, then each study's mark (only when there are several)
  const key = ui.divH([...groups.map((g, i) => swatch(g.name, colour(i))),
    ...several ? studies.map((st) => {
      const m = ui.span([mark(st)[1]], {style: {color: 'var(--grey-5)', marginRight: '6px'}});
      return ui.span([m, st], {style: {whiteSpace: 'nowrap'}});
    }) : []], {style: {flexWrap: 'wrap', gap: '4px 16px', fontSize: '12px', marginBottom: '4px', maxHeight: '72px',
    overflowY: 'auto'}});
  return {root: ui.divV([key, chart.root]), dispose: chart.dispose};
}

/** A study's control, read as Analyze reads it: its animals' mean growth rate of log volume (per day, least squares,
 *  each animal measured on 2 days or more) from day 0 to the evaluation day (the analysis's when it kept one, else the
 *  first day the median control reached `endpoint`, else the last day); its last day measured and its median volume
 *  on day 0 (null: none). */
export function controlGrowth(s: PreviousStudy, endpoint = ENDPOINT, evaluationDay: number | null = null):
  {slope: number | null, lastDay: number | null, toDay: number | null, vRand: number | null, n: number} {
  const control = controlOf(groupsIn(s), s.meta);
  const rows = trustedRows(s);
  const day0 = rows.filter((r) => r.day === 0).map((r) => r.volume!);
  const mine = rows.filter((r) => control !== null && r.group === control);
  const byDay = new Map<number, number[]>();
  for (const r of mine)
    byDay.set(r.day!, [...byDay.get(r.day!) ?? [], r.volume!]);
  const reached = [...byDay].sort((a, b) => a[0] - b[0]).find(([, xs]) => median(xs) >= endpoint)?.[0] ?? null;
  const lastDay = mine.length ? maxOf(mine.map((r) => r.day!)) : null;
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
  // a study with no volume yet (a new one, its sheet still empty) has nothing to compare: left out, said, not counted
  const empty = entries.filter((e) => texts.has(e) && !measured(previousOfSaved(parseStudy(texts.get(e)!))));
  const shownEntries = entries.filter((e) => texts.has(e) && !empty.includes(e));
  const view = DG.View.create();
  view.name = 'PMX Studio';
  // a page of the app, quiet as the Overview is (stages.ts quietPanels)
  view.root.dataset[STANDS_ALONE] = 'true';
  view.helpUrl = helpPage();
  const host = ViewHost.of(view);
  // left for a stage, Compare comes back on Back or Forward with the same studies
  keepOnLeaving(host, async () => () => compareView(entries, read));
  // as every page of the app: PMX Studio to the Overview (Compare is no stage: no steps), the report in the view's
  // menu and the toolbox
  let report = () => {};
  const exports: Exports = [['Report (PDF)…', () => report()]];
  view.setRibbonPanels([[headBar(ui.div([crumb(() => void openStage('home'))], 'pmx-head-center'),
    ui.span([`Compare · ${shownEntries.length} studies`], 'pmx-head-chip'))]]);
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
  // what is compared, in one line filled as each part is read (each figure's meaning on hover)
  const summary = ui.divText('', 'pmx-compare-summary');
  ui.tooltip.bind(summary, () => `${shownEntries.length} studies, of ` +
    `${[...new Set(models.filter((m) => m).map(modelName))].join(', ') || 'no model given'}; with the patients’ ` +
    'exposure: the treated groups with both a TGI and a multiple of it (the exposure chart); controls within the ' +
    'atlas: the studies whose control grew within the atlas’s range for its model');
  host.effect(() => {
    const said = (v: StatValue | undefined) => v === undefined ? '…' : String(v);
    summary.textContent = [`${said(nStudies.value)} studies`, `${said(nTreated.value)} treated groups`,
      `${said(nPlaced.value)} with the patients’ exposure`, `controls within the atlas: ${said(nWithin.value)}`]
      .join(' · ');
  });

  // (ii) the growth of every group: Linear, Log or each animal over its first volume, from the toggle in the chart
  const axisMode = signal('Linear');
  const axisRow = axisSwitch([['Linear', 'Volume on a linear axis'], ['Log', 'Volume on a logarithmic axis'],
    ['V/V0', 'Each animal’s volume over its first']], axisMode);
  const growth = growthFrame(studiesNow(), false);
  const relative = growthFrame(studiesNow(), true);
  const growthData = computed(() => (axisMode.value === 'V/V0' ? relative : growth)?.df);
  const several = shownEntries.length > 1;
  const growthChart = growth ? host.build(() => new BandChart(growthData, {x: 'day', lower: 'lower', upper: 'upper',
    line: 'median', split: 'series', xTitle: 'Day after randomization', yTitle: 'Tumor volume', markers: true,
    log: computed(() => axisMode.value === 'Log'), symbols: growth.symbols, count: 'n', legend: !several})) : null;
  // several studies: the key a line per study, its groups after it (the chart's own legend paged 36 entries)
  const keyColors = compareColors(studiesNow());
  const growthKey = several ? ui.divV(studiesNow().map(({s, label}, k) => {
    const groups = groupsIn(s);
    const control = controlOf(groups, s.meta);
    const own = keyColors[k];
    return ui.divH([ui.span([`${STUDY_MARKS[k % STUDY_MARKS.length][1]} ${label}`], {style: {fontWeight: '500'}}),
      ...groups.map((g) => swatch(g === control ? 'Control' : describeArm(armParts(g, s.meta)) || g, own[g], g))],
    {style: {flexWrap: 'wrap', gap: '2px 12px', alignItems: 'center'}});
  }), {style: {fontSize: '12px', gap: '2px', marginBottom: '4px', maxHeight: '96px', overflowY: 'auto'}}) : null;

  // (v) the controls against the atlas: filled once the atlas has predicted each study's
  const controlsBox = ui.divV([muted('Reading the atlas for each study’s model…')]);
  const exposureBox = ui.divV([]);
  const doseBox = ui.divV([]);
  const groupsBox = ui.divV([]);

  const heading = ui.divV([summary, ...warning ? [noteLine(warning)] : [],
    ...skipped.map((n) => noteLine(`${n} could not be read (moved, deleted or damaged): left out.`)),
    ...empty.map((e) => noteLine(`${e.name || e.id || 'A study'} has no measurements yet: left out.`))],
  'pmx-compare-head');
  const exposureSummary = computed(() => nPlaced.value === undefined ? '' :
    `${nPlaced.value} groups with the patients’ exposure`);
  const controlsSummary = computed(() => nWithin.value === undefined ? 'reading the atlas…' :
    `${nWithin.value} studies within the atlas’s range`);
  const panels = host.build(() => panelGrid([
    {id: 'growth', title: 'Tumor growth', hint: 'Median and range (smallest to largest tumor) of each group, a color ' +
      'per group and study; each study its mark', content: growth ? [axisRow,
      ...growthKey ? [growthKey] : [], sized(growthChart!, 360)] :
      [muted('No volumes to draw.')]},
    {id: 'dose', title: 'Dose and effect', hint: 'Each drug’s TGI at each dose it was given, across the studies',
      content: [doseBox]},
    {id: 'exposure', title: 'Exposure and effect', hint: 'Each treated group’s TGI against its ' +
      'multiple of the patients’ exposure; the patients’ range shaded', content: [exposureBox],
    summary: exposureSummary},
    {id: 'controls', title: 'Controls against the atlas', hint: 'Each study’s median control over the atlas’s ' +
      'prediction for its model (its 95% range, as Analyze’s gray band)', content: [controlsBox],
    summary: controlsSummary},
    {id: 'groups', title: 'Treated groups', hint: 'Every treated group: its TGI, its exposure and its result',
      content: [groupsBox]},
  ], [{id: 'growth', span: 6}, {id: 'dose', span: 6}, {id: 'groups', span: 12},
    {id: 'exposure', span: 12, folded: true}, {id: 'controls', span: 12, folded: true}], 'pmx-compare', 'compare',
  host));
  exports.push(['Customize this page', () => customizePage(panels), 'View']);
  const back: ToolLink[] = [['Back to the Overview', () => void openStage('home'),
    'Your studies, to pick others to compare']];
  view.toolbox = stageToolbox(exports, back, compareHelp(), panels);
  stageMenu(view, exports, {open: () => void openStudyDialog()}, {study: back, help: compareHelp()},
    (d) => host.own(d));
  const page = ui.divV([heading, panels], 'pmx-page');
  Object.assign(page.style, {padding: '16px', overflow: 'auto'});
  view.root.replaceChildren(page);

  // the charts drawn last, let go when drawn again or when the view closes
  let disposeCharts = () => {};
  host.own(() => disposeCharts());
  const withRunning = ui.input.bool('Include running studies', {value: false,
    tooltipText: 'A running study’s effect is not final: off, its groups are left out of this chart',
    onValueChanged: () => draw()});
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
    if (shown.length)
      drawn.push(compareChart(shown, 'TGI (%)'));
    // a running study's effect is not final: out of the dose chart unless asked for
    const running = new Set(studies.filter(({e}) => e.stage === 2).map(({label}) => label));
    const doses = doseGroups(rows, withRunning.value ? new Set() : running);
    if (doses.length)
      drawn.push(doseChart(doses, studies.length > 1));
    disposeCharts = () => drawn.forEach((d) => d.dispose());
    const exposureChart = shown.length ? drawn[0].root : null;
    const oneDose = doses.length > 0 && doses.every((d) => new Set(d.rows.map((r) => r.doseValue)).size < 2);
    doseBox.replaceChildren(...running.size ? [withRunning.root] : [],
      ...doses.length ? [drawn[drawn.length - 1].root] : [],
      ...!doses.length && running.size && doseGroups(rows).length ? [muted('Only the running studies give a ' +
        'drug at a known dose so far: include them to see it.')] :
      !doses.length ? [muted('No group gives one drug at a known dose yet. Name a group with its drug and dose ' +
        '(e.g. Trametinib 1 mg/kg), or add the dose in the study’s details, to see its effect here.')] :
        oneDose ? [muted('Each drug has one dose so far: a study with a second dose draws its curve.')] : []);

    // built under the view's host: draw() also runs after an await (a patients' exposure added), when nothing builds
    // the treated groups under their study (its model beside it), each in its study's colour, the numbers aligned; a
    // header sorts the rows within each study; where a number comes from, on hover
    const colors = compareColors(studies);
    const colorOf = new Map(studies.flatMap(({s, label}, k) =>
      groupsIn(s).map((g): [string, string] => [`${label}|${g}`, colors[k][g]])));
    const sort = signal('');
    const value = (r: CompareRow, c: string): number | string => c === 'tgi' ? r.tgi ?? -Infinity :
      c === 'day' ? r.day ?? -Infinity : c === 'x' ? r.multiple?.ratio ?? -Infinity : c === 'result' ?
        r.result ?? '' : plainArm(r.treatment);
    const items = computed(() => {
      const by = sort.value.replace(/^!/, '');
      const sign = sort.value.startsWith('!') ? -1 : 1;
      return [...new Set(rows.map((r) => r.study))].flatMap((study) => {
        const mine = rows.filter((r) => r.study === study);
        const model = mine.find((r) => r.model)?.model;
        if (by) {
          mine.sort((p, q) => {
            const [x, y] = [value(p, by), value(q, by)];
            return sign * (typeof x === 'number' && typeof y === 'number' ? x - y : String(x).localeCompare(String(y)));
          });
        }
        return [{section: model ? `${study} · ${modelName(model)}` : study}, ...mine];
      });
    });
    const tipped = (el: HTMLElement, tip: string | null | undefined) => {
      if (tip)
        ui.tooltip.bind(el, tip);
      return el;
    };
    const exposureFrom = (r: CompareRow) => r.source ? `Patients' exposure from ${r.source}${r.cv !== null ?
      ` (CV ${sig(r.cv)}%)` : ''}` : '';
    const table = host.build(() => new BasicTable<CompareRow>({items, fixedLayout: true, sort,
      keyOf: (r) => `${r.study}|${r.group}`, columns: [
        // its drugs and doses, in its study's colour; the schedule and the codes on hover
        {header: 'Group', name: 'group', sortable: true, width: '40%', render: (r) => {
          const plain = plainArm(r.treatment);
          const short = plain.split(' + ').map((p) => p.split(' · ').slice(0, 2).join(' ')).join(' + ');
          const el = swatch(short, colorOf.get(`${r.study}|${r.group}`) ?? PREDICTED_COLOR, r.group);
          return tipped(el, plain === r.treatment ? plain : `${plain} (${r.treatment})`);
        }},
        {header: 'TGI', name: 'tgi', sortable: true, align: 'right', width: '12%',
          tooltip: 'Tumor growth inhibition on the day shown', render: (r) => r.tgi !== null ? r.toCheck ?
            // read with values Monitor asks to check: its warning mark, and why on hover
            tipped(ui.span([dot('warning'), ' ', pct(r.tgi)]), `${tgiSource(r)}. ${r.toCheck} ` +
              `${r.toCheck === 1 ? 'value' : 'values'} to check in this group or its control, read as ` +
              'they are: correct them in Monitor.') :
            tipped(ui.span([pct(r.tgi)]), tgiSource(r)) : tipped(ui.span(['—']), r.day === null ? null :
              `${r.study}, day ${r.day}: no TGI, no growth to read it from (no earlier volume, or the control ` +
              'did not grow)')},
        {header: 'Day', name: 'day', sortable: true, align: 'right', width: '8%',
          render: (r) => r.day === null ? '—' : String(r.day)},
        {header: '× patients', name: 'x', sortable: true, align: 'right', width: '14%',
          tooltip: 'The group’s exposure as a multiple of the patients’', render: (r) => {
            // a combination: each drug's multiple, its name first (each drug at the AUC read in this group's mice)
            if (r.drugs?.some((d) => d.m)) {
              return tipped(ui.divV(r.drugs.map((d) => tipped(ui.divText(`${d.drug} ${d.m?.ratio !== undefined ?
                shownMultiple(d.m) : '—'}`), d.m ? d.m.ratio !== undefined ? d.m.tip : d.m.why : null))), COMBO_AUC);
            }
            const m = r.multiple;
            return m?.ratio !== undefined ? tipped(ui.span([shownMultiple(m)]), [m.tip, exposureFrom(r)]
              .filter(Boolean).join('. ')) : tipped(muted('—'), m?.why ?? (r.key ? 'No patients’ exposure yet: add it ' +
              'under Exposure and effect' : 'A combination: not placed against one drug’s exposure'));
          }},
        {header: 'Result', name: 'result', sortable: true, width: '26%', render: (r) => r.result ?
          ui.divH([dot(r.result === 'Slows growth' ? 'success' : r.result === 'Grows faster' ? 'error' : 'default'),
            ui.span([r.result])], {style: {alignItems: 'center', gap: '6px'}}) : muted('—')},
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
    exposureBox.replaceChildren(...sentence ? [ui.divText(sentence)] : [],
      ...exposureChart ? [exposureChart] : [muted('No treated group has both a TGI and a multiple of the ' +
        'patients’ exposure yet.')], ...hints, muted(compareFooter(rows)));
    groupsBox.replaceChildren(table.root);
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
  const studyColors = compareColors(studies);
  read.forEach((r, i) => {
    const control = controlOf(groupsIn(r.s), r.s.meta);
    // each control in its study's colour, as the growth chart draws it
    const color = control ? studyColors[i][control] : PREDICTED_COLOR;
    if (r.atlas && r.g.lastDay !== null) {
      for (const b of controlsTo(r.atlas.controls, r.g.lastDay)) {
        // the atlas in the predicted grey, as Analyze draws its predicted control (dashed: told apart by study)
        rows.push({series: `${r.label}: atlas`, color: PREDICTED_COLOR, day: b.day, lower: b.q025, upper: b.q975,
          median: b.median});
      }
    }
    const byDay = new Map<number, number[]>();
    for (const x of trustedRows(r.s).filter((x) => control !== null && x.group === control))
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
