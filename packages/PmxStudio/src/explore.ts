/* Analyze's Your data tab: the study as Datagrok's own tables, for the questions the other tabs do not ask. The tabs
   answer the study's question with fixed figures; here the animals (one row each) and the measurements (one row
   each) open as table views, each with a few viewers to start from, and from there the platform's whole kit: any
   viewer, filters, clustering, layouts. The two tables are linked by animal: what is selected or filtered among the
   animals is among their measurements. The cards are the Overview's, so a door reads as a door. */
import '../css/hub.css';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {div, span} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

import type {AnalysisFrames} from './analyze';
import {finite} from './engine';

/** The study table and its columns as the analysis read them. */
export interface StudyInputs {study: DG.DataFrame, animal: string, group: string, day: string, volume: string,
  bodyWeight: string | null}

/** What a call reads, taken before it runs. */
export function studyInputsOf(call: DG.FuncCall): StudyInputs | null {
  const study = call.inputs['study'];
  const name = (p: string) => {
    const c = call.inputs[p];
    return c instanceof DG.Column ? c.name : typeof c === 'string' && c ? c : null;
  };
  const [animal, group, day, volume] = ['animal', 'group', 'day', 'volume'].map(name);
  return study instanceof DG.DataFrame && animal && group && day && volume ?
    {study, animal, group, day, volume, bodyWeight: name('bodyWeight')} : null;
}

export const COL = {animal: 'Animal', group: 'Group', v0: 'Volume at randomization', rate: 'Growth rate',
  tte: 'Day to endpoint', reached: 'Endpoint reached', mrecist: 'mRECIST', effect: 'Smaller than its twin',
  effectLow: 'Smaller than its twin, low', effectHigh: 'Smaller than its twin, high', response: 'Response',
  day: 'Day', volume: 'Tumor volume', weight: 'Body weight'};
export const RESPONSE = {responder: 'Responder', like: 'As its twin', larger: 'Larger than its twin',
  control: 'Control'};
const RESPONSE_COLORS = {[RESPONSE.responder]: '#2ca02c', [RESPONSE.like]: '#9aa5b1',
  [RESPONSE.larger]: '#d62728', [RESPONSE.control]: '#4d5261'};

/** The group colours every chart of the analysis uses, by group. */
export type GroupColors = {[group: string]: string};


/** A column's values as numbers, an empty cell NaN (an integer column's reads as INT_NULL, which is finite). */
export const numbersOf = (c: DG.Column): Float64Array =>
  Float64Array.from({length: c.length}, (_, i) => c.isNone(i) ? NaN : Number(c.get(i)));

/** The animals, one row each: its growth, its endpoint and its response, under names a reader knows. */
export function animalsTable(r: AnalysisFrames, s: StudyInputs, colors: GroupColors): DG.DataFrame {
  const control = r.analysis.get('control', 0) as string;
  const a = r.animals;
  const e = r.effects;
  const effectRow = new Map<string, number>();
  for (let i = 0; i < e.rowCount; i++)
    effectRow.set(String(e.get('animal', i)), i);
  // the volume at randomization: each animal's measurement nearest day 0
  const v0 = new Map<string, {d: number, v: number}>();
  const [an, dy, vol] = [s.animal, s.day, s.volume].map((c) => s.study.col(c)!);
  for (let i = 0; i < s.study.rowCount; i++) {
    const [id, d, v] = [String(an.get(i)), finite(dy.get(i)), finite(vol.get(i))];
    if (d === null || v === null)
      continue;
    const best = v0.get(id);
    if (!best || Math.abs(d) < Math.abs(best.d))
      v0.set(id, {d, v});
  }
  const twinned = !!r.analysis.get('twinned', 0);
  const rows = Array.from({length: a.rowCount}, (_, i) => {
    const id = String(a.get('animal', i));
    const group = a.get('group', i) as string;
    const j = effectRow.get(id);
    const eff = (k: string) => j === undefined ? null : finite(e.get(k, j));
    const [lo, hi] = [eff('effect_lower'), eff('effect_upper')];
    const response = group === control ? RESPONSE.control : !twinned || lo === null || hi === null ? null :
      lo > 0 ? RESPONSE.responder : hi < 0 ? RESPONSE.larger : RESPONSE.like;
    return {[COL.animal]: id, [COL.group]: group, [COL.v0]: v0.get(id)?.v ?? null,
      [COL.rate]: finite(a.get('slope_per_day', i)), [COL.tte]: finite(a.get('time_to_endpoint', i)),
      [COL.reached]: !!a.get('endpoint_reached', i), [COL.mrecist]: a.get('mrecist', i) ?? null,
      [COL.effect]: eff('effect'), [COL.effectLow]: lo, [COL.effectHigh]: hi, [COL.response]: response};
  });
  const df = DG.DataFrame.fromObjects(rows)!;
  df.name = `Animals · ${s.study.name}`;
  groupsIn(df, r, colors);
  meta(df, COL.v0, 'mm³', '0.0');
  meta(df, COL.rate, 'ln V/d', '0.000');
  meta(df, COL.tte, 'd', '0');
  for (const c of [COL.effect, COL.effectLow, COL.effectHigh])
    meta(df, c, '%', '0');
  df.col(COL.response)?.meta.colors.setCategorical(RESPONSE_COLORS);
  return df;
}

/** The measurements, one row each, with the animal's response: each tumor's curve, in its group's colour. */
export function measurementsTable(s: StudyInputs, animals: DG.DataFrame, r: AnalysisFrames,
  colors: GroupColors): DG.DataFrame {
  const responseOf = new Map<string, string | null>();
  for (let i = 0; i < animals.rowCount; i++)
    responseOf.set(animals.get(COL.animal, i), animals.get(COL.response, i));
  const [an, gr, dy, vol] = [s.animal, s.group, s.day, s.volume].map((c) => s.study.col(c)!);
  const bw = s.bodyWeight ? s.study.col(s.bodyWeight) : null;
  const groupOf = new Map<string, string>();
  const rows = Array.from({length: s.study.rowCount}, (_, i) => {
    const [id, g] = [String(an.get(i)), String(gr.get(i))];
    groupOf.set(id, g);
    return {[COL.animal]: id, [COL.group]: g, [COL.day]: finite(dy.get(i)),
      [COL.volume]: finite(vol.get(i)), ...(bw ? {[COL.weight]: finite(bw.get(i))} : {}),
      [COL.response]: responseOf.get(id) ?? null};
  });
  const df = DG.DataFrame.fromObjects(rows)!;
  df.name = `Measurements · ${s.study.name}`;
  groupsIn(df, r, colors);
  meta(df, COL.day, 'd', '0');
  meta(df, COL.volume, 'mm³', '0.0');
  meta(df, COL.weight, 'g', '0.0');
  df.col(COL.response)?.meta.colors.setCategorical(RESPONSE_COLORS);
  // each animal's curve in its group's colour
  df.col(COL.animal)!.meta.colors.setCategorical(Object.fromEntries([...groupOf].map(([id, g]) =>
    [id, colors[g] ?? '#9aa5b1'])));
  return df;
}

/** The groups in the analysis's colours and order, the control first. */
function groupsIn(df: DG.DataFrame, r: AnalysisFrames, colors: GroupColors): void {
  const g = df.col(COL.group)!;
  g.meta.colors.setCategorical(colors);
  g.setCategoryOrder([r.analysis.get('control', 0), ...r.comparisons.col('arm')!.toList()]);
}

function meta(df: DG.DataFrame, col: string, units: string, format: string): void {
  const c = df.col(col);
  if (c) {
    c.meta.units = units;
    c.meta.format = format;
  }
}

/** The animals' view: the table, their growth rates by group compared with the control, and their growth rate
 *  against their volume at randomization. */
function animalsView(df: DG.DataFrame, control: string, parent: DG.ViewBase): DG.TableView {
  const tv = grok.shell.addTableView(df);
  const box = DG.Viewer.fromType(DG.VIEWER.BOX_PLOT, df, {valueColumnName: COL.rate,
    categoryColumnNames: [COL.group], showGroupComparison: true, controlComparisons: true, controlGroup: control,
    showPValue: false});
  const scatter = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, df, {xColumnName: COL.v0, yColumnName: COL.rate,
    colorColumnName: COL.group});
  const right = tv.dockManager.dock(box, DG.DOCK_TYPE.RIGHT, null, 'Growth rate by group', 0.55);
  tv.dockManager.dock(scatter, DG.DOCK_TYPE.DOWN, right, 'Growth rate against volume at randomization', 0.5);
  backTo(tv, parent);
  return tv;
}

/** The measurements' view: each tumor's curve on a log axis, in its group's colour, beside the table. */
function measurementsView(df: DG.DataFrame, parent: DG.ViewBase): DG.TableView {
  const tv = grok.shell.addTableView(df);
  const curves = DG.Viewer.fromType(DG.VIEWER.LINE_CHART, df, {xColumnName: COL.day, yColumnNames: [COL.volume],
    splitColumnNames: [COL.animal], yAxisType: 'logarithmic', legendVisibility: 'Never'});
  tv.dockManager.dock(curves, DG.DOCK_TYPE.RIGHT, null, 'Tumor volume of each animal', 0.65);
  backTo(tv, parent);
  return tv;
}

/** The crumb back to Analyze over a Your data view; the grid's coloured columns in coloured text, not filled cells
 *  (an animal's every row in its group's colour fills the grid). */
function backTo(tv: DG.TableView, parent: DG.ViewBase): void {
  const host = ViewHost.of(tv);
  host.build(() => host.ribbon({breadcrumbs: {parent}}));
  for (const c of [COL.animal, COL.group, COL.response]) {
    const gc = tv.grid.col(c);
    if (gc)
      gc.isTextColorCoded = true;
  }
}

/** The two views of one result, built once and brought back while open. */
const opened = new WeakMap<AnalysisFrames, {animals: DG.TableView, measurements: DG.TableView}>();

/** Opens (or brings back) the Your data view `which` of result `r`. */
export function explore(which: 'animals' | 'measurements', r: AnalysisFrames, s: StudyInputs, colors: GroupColors,
  parent: DG.ViewBase): void {
  const isOpen = (v: DG.TableView) => [...grok.shell.tableViews].some((x) => x.dataFrame === v.dataFrame);
  let views = opened.get(r);
  if (!views || !isOpen(views.animals) || !isOpen(views.measurements)) {
    for (const v of views ? [views.animals, views.measurements] : []) {
      if (isOpen(v))
        v.close();
    }
    const animals = animalsTable(r, s, colors);
    const measurements = measurementsTable(s, animals, r, colors);
    // what is selected or filtered among the animals is among their measurements
    grok.data.linkTables(animals, measurements, [COL.animal], [COL.animal],
      [DG.SYNC_TYPE.SELECTION_TO_SELECTION, DG.SYNC_TYPE.FILTER_TO_FILTER]);
    views = {measurements: measurementsView(measurements, parent),
      animals: animalsView(animals, r.analysis.get('control', 0), parent)};
    opened.set(r, views);
  }
  // after the click that asked for it: the click, reaching the Analyze view, makes Analyze current again
  const view = views[which];
  setTimeout(() => grok.shell.v = view);
}

/** A card as the Overview draws one: an icon, a name, what it opens. */
function doorCard(name: string, what: string, iconPaths: string, open: () => void): HTMLButtonElement {
  const b = document.createElement('button');
  b.type = 'button';
  b.className = 'pmx-card';
  const ico = span('', 'pmx-ico');
  ico.innerHTML = `<svg width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="currentColor" ` +
    `stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${iconPaths}</svg>`;
  const title = span(name + ' ', 'pmx-name');
  const go = span('→', 'pmx-go');
  go.setAttribute('aria-hidden', 'true');
  title.append(go);
  b.append(ico, title, span(what, 'pmx-what'));
  b.addEventListener('click', open);
  return b;
}

const ICON = {
  animals: '<circle cx="7" cy="16" r="2.2"/><circle cx="12" cy="9" r="2.2"/><circle cx="17" cy="14" r="2.2"/>' +
    '<circle cx="16" cy="5.5" r="1.6"/><path d="M3 21h18"/>',
  measurements: '<path d="M4 4v16h16"/><path d="M7 15c3-1 5-3 7-6s3-3 5-4"/><path d="M7 17c3 0 6-2 8-4s3-2 4-3"/>',
};

/** The Your data tab's content: the two doors, and what can be done beyond them. */
export function exploreContent(open: (which: 'animals' | 'measurements') => void): HTMLElement[] {
  const cards = div([
    doorCard('Animals', 'One row per animal: growth rate, day to endpoint, mRECIST and response against its twin. ' +
      'Opens with each group compared with the control.', ICON.animals, () => open('animals')),
    doorCard('Measurements', 'One row per measurement: the curve of each tumor. The animals you select or filter ' +
      'are the ones drawn.', ICON.measurements, () => open('measurements')),
  ], 'pmx-explore-cards');
  const tips = div([
    span('There, as with any table in Datagrok:'),
    tip('Add any chart', 'from the Viewers panel on the left: scatter, histogram, trellis, heat map…'),
    tip('Group the animals', 'with ML › Cluster, or ML › Reduce Dimensionality for a map of the animals that grow ' +
      'alike.'),
    tip('Keep your arrangement', 'from the Layouts panel on the left, to open your next study the same way.'),
  ], 'pmx-explore-tips');
  return [div([cards, tips], 'pmx-explore')];
}

function tip(what: string, how: string): HTMLElement {
  const line = div([], 'pmx-explore-tip');
  const b = document.createElement('b');
  b.textContent = what;
  line.append(b, ` ${how}`);
  return line;
}
