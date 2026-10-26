/* The "Tolerability" tab of the Analyze view: whether the treatments were well tolerated, from the animals' body
   weight (pmxAnalyzeStudy's tolerability and weights, Crux stats.tumor_body_weight): the figures, each group's body
   weight change from day 0 (as weighed, or without the tumor's mass), the weight-loss grades (G0–G3, CTCAE v5 as
   adapted for mice), each group's nadir against the vehicle's, the animals removed early and the NCI limits for a
   treatment. Without a body weight column, how to add one. This file only shows it. */
import * as DG from 'datagrok-api/dg';
import {badge, BasicTable, BoolInput, computed, div, divH, divV, dot, span, Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal, StatValue} from '@datagrok-libraries/u2';
import {FuncCallForm} from '@datagrok-libraries/u2/src/dg/index.js';

import {AnalysisFrames, BASES, DENSITY_TAG, DensityUsed, NCI} from './analyze';
import {BandChart} from './band-chart';
import {block, cardGrid, dayText, infoNote, noteLine, num, pct, signed, sized, statCard, swatch, wrapCells}
  from './ui';

type Row = {[k: string]: any};

/** The rows of a table as objects, an empty cell as null (an integer column holds its own null value). */
const rowsOf = (df: DG.DataFrame): Row[] => Array.from({length: df.rowCount}, (_, i) =>
  Object.fromEntries(df.columns.toList().map((c) => [c.name, c.isNone(i) ? null : c.get(i)])));

/** Why a group's animals left early, in words: "1 dead, 1 weight, 2 tumor". */
const reasons = (g: Row) => [[g.removed_death, 'dead'], [g.removed_weight, 'weight'], [g.removed_volume, 'tumor'],
  [g.removed_other, 'other']].filter(([n]) => n > 0).map(([n, why]) => `${n} ${why}`).join(', ');

/** The NCI limits for a treatment, in words. */
const NCI_RULE = `a mean loss over ${NCI.nadir}% at the nadir (as weighed), or more than ${NCI.deaths}% of the ` +
  'animals dead or removed for their weight';

/** The CTCAE v5 weight-loss grades (NCI 2017; "Weight loss" under Investigations): the badge, its colour and the loss
 *  from day 0 it stands for. */
const GRADES: [string, 'success' | 'default' | 'warning' | 'error', string][] = [['G0', 'success', 'below 5%'],
  ['G1', 'default', '5–10%'], ['G2', 'warning', '10–20%'], ['G3', 'error', '20% or more, the usual humane endpoint']];
const CTCAE_URL = 'https://ctep.cancer.gov/protocoldevelopment/electronic_applications/docs/' +
  'CTCAE_v5_Quick_Reference_8.5x11.pdf';

const NADIR = 'the group’s lowest mean weight change from day 0, on its worst day';

/** A group against the NCI limits, as weighed: over them (the criteria it meets, and in `also` what else is worth a
 *  look), within but worth a look (a mean loss in G2, an animal dead, removed for its weight or past 20%, animals gone
 *  for an unrecorded reason) or clear. */
export type NciReading = {kind: 'error' | 'warning' | 'success', text: string, also: string, criteria: string,
  nadir: number | null};

/** Each group against the NCI limits, as weighed, by group; empty without body weight. */
export function nciReadings(r: AnalysisFrames): Map<string, NciReading> {
  const control = r.analysis.get('control', 0);
  return new Map(rowsOf(r.tolerability).filter((g) => g.basis === BASES.measured).map((g) =>
    [g.group, nciReading(g, g.group === control)]));
}

function nciReading(g: Row, vehicle: boolean): NciReading {
  const loss = g.mean_nadir === null ? 0 : Math.max(0, -g.mean_nadir);
  const lost = g.removed_death + g.removed_weight;
  const animals = (k: number) => `${k} animal${k === 1 ? '' : 's'}`;
  const share = `${lost} of ${g.n} animals dead or removed for their weight (${pct(lost / g.n)}`;
  const over: string[] = [];
  const watch: string[] = [];
  // which limit a group is over is Crux's reading (over_nadir, over_deaths); this only words it
  if (g.over_nadir)
    over.push(`a mean loss of ${Math.round(loss)}% at the nadir (limit ${NCI.nadir}%)`);
  else if (loss >= 10)
    watch.push(`a mean loss of ${Math.round(loss)}% at the nadir (G2)`);
  if (g.over_deaths)
    over.push(`${share}; limit ${NCI.deaths}%)`);
  else if (lost > 0)
    watch.push(`${share})`);
  if (g.grade_3 > 0 && lost === 0)
    watch.push(`${animals(g.grade_3)} lost 20% or more (G3)`);
  const unknown = `${animals(g.removed_other)} left early for an unrecorded reason (a Removal reason column tells ` +
    'deaths)';
  if (g.removed_other > 0)
    watch.push(unknown);
  const criteria = over.join('; and ') || 'see the table';
  const base = {criteria, nadir: g.mean_nadir, also: ''};
  if (g.toxic) {
    return {...base, kind: 'error', text: `over the limits: ${criteria}.` +
      (vehicle ? ' The tumor itself costs weight: read the treated groups against it.' : ''),
    also: watch.length ? `Also to watch: ${watch.join('; ')}.` : ''};
  }
  return watch.length ? {...base, kind: 'warning', text: `within the limits, to watch: ${watch.join('; ')}.`} :
    {...base, kind: 'success', text: 'within the limits, no criterion met.'};
}

/** The tumor's density the corrected weights use, in words: the line's typical one (its studies and their range),
 *  the one given, or this study's from its weighed tumors (by group, when the groups differ). */
function densityText(r: AnalysisFrames): string {
  const used = JSON.parse(r.tolerability.getTag(DENSITY_TAG) ?? '{}') as DensityUsed;
  const gs = rowsOf(r.tolerability).filter((g) => g.basis === BASES.tumor_corrected);
  const pairs = gs.reduce((t, g) => t + g.density_pairs, 0);
  const ignored = gs.reduce((t, g) => t + g.density_ignored, 0);
  const d = (x: number) => x.toFixed(2);
  const left = ignored ? ` ${ignored} tumor weight${ignored > 1 ? 's' : ''} not used: a small tumor, or a weight out ` +
    'of line with its volume.' : '';
  if (pairs === 0) {
    return `its volume times ${used.declared ? `the density given, ${d(used.density)} mg per mm³.` :
      `${used.line}’s typical density, ${d(used.density)} mg per mm³ (${used.papers} studies that weighed their ` +
      `tumors; ${d(used.range[0])} to ${d(used.range[1])} between labs).`}${left}`;
  }
  const ds = gs.map((g) => g.density as number);
  const [lo, hi] = [Math.min(...ds), Math.max(...ds)];
  return `its volume times this study’s density, ${d(lo) === d(hi) ? d(lo) : `${d(lo)} to ${d(hi)}`} mg per mm³ ` +
    `${d(lo) === d(hi) ? '' : 'by group '}(${pairs} tumors weighed once excised), and on the day a tumor was excised ` +
    `its own weight.${left}`;
}

/** Secondary text: smaller and muted. */
const muted = (el: HTMLElement) => {
  Object.assign(el.style, {fontSize: '12px', color: 'var(--dg-text-color-muted, #68727f)', lineHeight: '1.45'});
  return el;
};

export function tolerabilityTab(form: FuncCallForm, result: ReadonlySignal<AnalysisFrames | undefined>,
  colorsOf: (r: AnalysisFrames) => {[group: string]: string}): {content: HTMLElement[], own: {dispose(): void}[]} {
  const corrected = new BoolInput({label: 'Tumor-corrected', switch: true, value: false,
    tooltipText: 'Take the tumor’s mass off each weighing'});
  const basis = computed(() => corrected.value.value ? BASES.tumor_corrected : BASES.measured);
  const has = computed(() => (result.value?.tolerability.rowCount ?? 0) > 0);
  // the control first, then the treated groups in their order, in the chosen reading
  const groups = computed((): Row[] => {
    const r = result.value;
    if (!r || !has.value)
      return [];
    const all = rowsOf(r.tolerability).filter((g) => g.basis === basis.value);
    const order = [r.analysis.get('control', 0), ...r.comparisons.col('arm')!.toList()];
    return order.map((name) => all.find((g) => g.group === name)).filter((g): g is Row => !!g);
  });
  const treated = computed(() => groups.value.slice(1));
  const lost = (g: Row) => g.grade_1 + g.grade_2 + g.grade_3;
  const removed = (g: Row) => g.removed_death + g.removed_weight + g.removed_volume + g.removed_other;

  // the study's answers first, what each figure covers in its tooltip
  const card = (label: string, value: ReadonlySignal<StatValue | undefined>, tip: () => string,
    format?: (x: StatValue) => string) => statCard(form, label, value, tip, format);
  const worst = computed(() => treated.value.reduce<Row | undefined>((w, g) =>
    g.mean_nadir !== null && (w === undefined || g.mean_nadir < w.mean_nadir) ? g : w, undefined));
  const losing10 = computed(() => treated.value.length === 0 ? undefined :
    `${treated.value.reduce((t, g) => t + g.grade_2 + g.grade_3, 0)} of ` +
    `${treated.value.reduce((t, g) => t + g.n, 0)}`);
  const toxic = computed(() => treated.value.filter((g) => g.toxic));
  const perGroup = (f: (g: Row) => string) => groups.value.map((g) => `${g.group}: ${f(g)}`).join('; ');
  const cards = cardGrid([
    card('Mean nadir', computed(() => worst.value?.mean_nadir ?? undefined), () => worst.value ?
      `Each group's lowest mean weight change from day 0 (its day): ${perGroup((g) =>
        g.mean_nadir === null ? '—' : `${signed(g.mean_nadir)} (day ${Math.round(g.mean_nadir_day)})`)}` : '',
    (x) => signed(x as number)),
    card('Treated animals losing 10% or more', losing10, () => 'Animals whose weight fell 10% or more below day 0 at ' +
      `some point: ${perGroup((g) => `${g.grade_2 + g.grade_3} of ${g.n}`)}`),
    card('Animals removed early', computed(() => has.value ? groups.value.reduce((t, g) => t + removed(g), 0) :
      undefined), () => 'Animals weighed for the last time before the rest of their group: ' +
      `${perGroup((g) => reasons(g) || 'none')}`),
    card('Treated groups over the NCI limits', computed(() => !has.value ? undefined :
      `${toxic.value.length} of ${treated.value.length}`), () => {
      return toxic.value.length ? `Over the NCI limits: ${toxic.value.map((g) => g.group).join(', ')}` :
        'No treated group is over the NCI limits';
    }),
  ], 'pmx-tolerability-cards', 150);

  // each group's mean change from day 0 and its range, in the chosen reading
  const chart = new BandChart(computed(() => {
    const r = result.value;
    if (!r || !has.value)
      return undefined;
    const w = r.weights;
    const b = basis.value;
    const df = w.clone(DG.BitSet.create(w.rowCount, (i) => w.get('basis', i) === b));
    df.col('group')!.meta.colors.setCategorical(colorsOf(r));
    return df;
  }), {x: 'day', lower: 'min', upper: 'max', line: 'mean', split: 'group', xTitle: 'Day after randomization',
    yTitle: 'Body weight change', markers: true});
  const toggleRow = divH([corrected.root]);
  Object.assign(toggleRow.style, {justifyContent: 'flex-end', alignItems: 'center'});
  // what the weight includes: the median control tumor on its last day, in grams
  const chartNote = span('', 'pmx-weight-note');
  Object.assign(chartNote.style, {display: 'block', fontSize: '12px', color: 'var(--dg-text-color-muted, #68727f)',
    lineHeight: '1.45'});
  form.effect(() => {
    const r = result.value;
    const on = corrected.value.value;
    if (!r || !has.value)
      return;
    const control = r.analysis.get('control', 0);
    const last = rowsOf(r.volumes).filter((v) => v.group === control).sort((a, b) => b.day - a.day)[0];
    const density = rowsOf(r.tolerability).find((g) => g.group === control)?.density ?? 1;
    const tumor = last ? ` On day ${last.day} the median control tumor weighs about ` +
      `${num(last.median * density / 1000, 2)} g.` : '';
    chartNote.textContent = (on ? `Each weighing without the tumor’s mass: ${densityText(r)}` :
      'Body weight as weighed, tumor included: switch on Tumor-corrected to see the animal’s own weight.') + tumor;
  });

  // the weight-loss grades of each group, how many recovered and the animals removed early
  const gradesTable = div([]);
  form.effect(() => {
    const r = result.value;
    const gs = groups.value;
    if (!r || !gs.length)
      return;
    const colors = colorsOf(r);
    const removedCell = (g: Row) => {
      const el = span(removed(g) ? reasons(g) : '0');
      if (g.removed_animals)
        Tooltip.bind(el, () => String(g.removed_animals));
      return el;
    };
    // the group's lowest mean change and, for a treated group, how far below the vehicle's (percentage points)
    const vehicle = gs[0];
    const nadirCell = (g: Row) => {
      if (g.mean_nadir === null)
        return '—';
      const d = g.mean_nadir - vehicle.mean_nadir;
      const vs = g === vehicle || vehicle.mean_nadir === null ? '' :
        ` (${d > 0.5 ? '+' : ''}${signed(d).replace('%', '')} vs vehicle)`;
      const el = span(`${signed(g.mean_nadir)}${vs}`);
      Tooltip.bind(el, () => `Lowest mean change on day ${Math.round(g.mean_nadir_day)}; the animals’ worst ` +
        `losses: ${g.median_nadir_day === null ? 'none of 5% or more' : dayText(g.median_nadir_day).toLowerCase() +
        ' (median)'}`);
      return el;
    };
    // each group against the NCI limits, read as weighed whichever reading the table shows
    const readings = nciReadings(r);
    const nciCell = (g: Row) => {
      const kind = readings.get(g.group)!.kind;
      const el = divH([dot(kind), span(kind === 'error' ? 'Over' : 'Within')]);
      Object.assign(el.style, {alignItems: 'center', gap: 'var(--dg-space-s)'});
      return el;
    };
    const columns: [string, (g: Row) => string | HTMLElement][] = [
      ['Group', (g) => swatch(g.group, colors[g.group])],
      ['Weighed', (g) => String(g.n)],
      ...GRADES.map(([grade], k): [string, (g: Row) => string] => [grade, (g) => String(g[`grade_${k}`])]),
      ['Mean nadir', nadirCell],
      ['Recovered', (g) => lost(g) ? `${g.recovered} of ${lost(g)}` : '—'],
      ['Removed', removedCell],
      ['NCI limits', nciCell],
    ];
    const table = new BasicTable<Row>({items: gs, keyOf: (g) => g.group, columns: columns.map(([header, render], j) =>
      ({header, render, align: j === 0 || j >= columns.length - 2 ? 'left' : 'center'}))});
    table.root.style.width = '100%';
    wrapCells(table.root);
    // the grade headers as the legend's badges (u2's table takes text headers), what the nadir is on hover
    const head = table.root.querySelector('.u2-table-header')?.children;
    GRADES.forEach(([grade, variant, range], k) => {
      const cell = head?.[2 + k] as HTMLElement | undefined;
      cell?.replaceChildren(badge(grade, {variant}));
      if (cell)
        Tooltip.bind(cell, () => `${grade}: largest loss from day 0 ${range} (CTCAE v5)`);
    });
    const nadirHead = head?.[2 + GRADES.length] as HTMLElement | undefined;
    if (nadirHead)
      Tooltip.bind(nadirHead, () => NADIR);
    gradesTable.replaceChildren(table.root);

    // the groups to look at, with the criteria they meet (an over group's other warnings after it, muted); the rest in
    // one line
    const flagged = gs.filter((g) => readings.get(g.group)!.kind !== 'success');
    const clear = gs.filter((g) => readings.get(g.group)!.kind === 'success').map((g) => g.group);
    nciSummary.replaceChildren(...flagged.map((g) => {
      const x = readings.get(g.group)!;
      const line = noteLine(`${g.group}: ${x.text}`, x.kind);
      if (x.also)
        (line.lastElementChild as HTMLElement).append(' ', muted(span(x.also)));
      return line;
    }), ...clear.length ? [noteLine(flagged.length ? `${clear.join(', ')}: within, no criterion met.` :
      'Every group within the NCI limits: no criterion met.', 'success')] : []);
  });
  // the grades' legend: the same badges, each one's loss, and where the grades come from
  const ctcae = document.createElement('a');
  ctcae.href = CTCAE_URL;
  ctcae.target = '_blank';
  ctcae.rel = 'noopener';
  ctcae.textContent = 'CTCAE v5';
  const gradesLegend = divH([...GRADES.map(([grade, variant, range]) => {
    const item = divH([badge(grade, {variant}), muted(span(range))]);
    Object.assign(item.style, {alignItems: 'center', gap: 'var(--dg-space-xs, 4px)'});
    return item;
  }), muted(span(''))]);
  (gradesLegend.lastElementChild as HTMLElement).append('Weight-loss grades of ', ctcae,
    ' as adapted for mice, each animal by its largest loss from day 0.');
  Object.assign(gradesLegend.style, {flexWrap: 'wrap', alignItems: 'center', gap: 'var(--dg-space-s) var(--dg-space-m)',
    marginTop: 'var(--dg-space-m)'});
  const nciSummary = divV([]);
  nciSummary.style.gap = 'var(--dg-space-xs, 4px)';
  const nciBlock = divV([muted(span(`NCI limits for a treatment (Plowman et al. 1997): ${NCI_RULE}.`)), nciSummary],
    'pmx-nci');
  Object.assign(nciBlock.style, {gap: 'var(--dg-space-s)', marginTop: 'var(--dg-space-m)'});
  const terms = muted(span(`Mean nadir: ${NADIR}; for a treated group, its difference from the vehicle’s in ` +
    'percentage points. Recovered: back within 5% of day 0 at its last weighing. Removed: weighed for the last time ' +
    'before the rest of its group (hover for which animals and when).'));
  Object.assign(terms.style, {display: 'block', marginTop: 'var(--dg-space-m)'});

  const body = divV([cards,
    block('Body weight change from day 0', 'Mean and range (lowest to highest animal) of each group', toggleRow,
      sized(chart, 340), chartNote),
    block('Weight loss', 'The animals of each group by their largest loss from day 0', gradesTable, gradesLegend,
      nciBlock, terms),
  ], 'pmx-tolerability');
  body.style.gap = 'var(--dg-space-l)';
  const missing = infoNote('<b>No body weight in this study.</b> Add a column with each animal’s body weight (g), ' +
    'weighed with the tumor on the rows of its measurements, choose it in <b>Body weight</b> and run. This tab then ' +
    'shows each group’s weight change from day 0, the weight-loss grades and the animals removed early.');
  form.effect(() => {
    body.hidden = !has.value;
    missing.hidden = has.value;
  });
  return {content: [missing, body], own: [corrected, chart]};
}
