/* The spider plot of Analyze: each animal's tumor as a % change from its volume at randomization, day by day, one
   line per animal in its group's colour; a group's legend entry shows or hides all its animals. echarts, as the band
   chart: the platform's line chart colours lines by its split, and here the split (the animal) is not the colour
   (the group). */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Control, signal} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';
import * as echarts from 'echarts/core';

import {initQuietly, legendLayout} from './band-chart';
import {closer, esc} from './sheet';
import {EXCLUDED, ticked} from './study-file';
import type {StudyInputs} from './explore';

/** The animals selected (marked on the chart, the others dimmed) and what a click on one does. */
export interface Picking {animals: ReadonlySignal<ReadonlySet<string>>, pick(animal: string, e?: MouseEvent): void}

/** The animal an echarts event is on: its point's, else (a click on a line between points: triggerLineEvent gives no
 *  point) its line's, each animal's line a series with the animal as its id (QA 0.1.3: only the 4 px markers of a
 *  growth curve answered a click). A line's click names its series by index only (its tooltip by id too): `ids`, the
 *  chart's series ids by index. */
export function animalOf(p: any, ids: (string | undefined)[] = []): string | undefined {
  if (p.data?.animal)
    return p.data.animal;
  const id = p.seriesType === 'line' ? p.seriesId ?? ids[p.seriesIndex] : undefined;
  return typeof id === 'string' ? id : undefined;
}

/** Where an echarts click on an animal goes: to picking.pick, with the keys held. */
function pickOn(chart: echarts.ECharts, picking?: Picking): void {
  if (picking) {
    chart.on('click', (p: any) => {
      const series = (chart.getOption() as {series?: {id?: string}[]}).series ?? [];
      const animal = animalOf(p, series.map((s) => s.id));
      if (animal)
        picking.pick(animal, p.event?.event);
    });
  }
}

/** One animal's line: its days and its % change from its volume at randomization (the last measured on or before
 *  day 0, else the first after it: the randomization sheet's own rule). */
export interface SpiderLine {animal: string, group: string, days: number[], change: number[]}

/** The study's animals as spider lines, its excluded measurements left out; an animal whose volume at randomization
 *  is not positive has no line. */
export function spiderLines(s: StudyInputs): SpiderLine[] {
  const df = s.study;
  const [a, g, d, v] = [s.animal, s.group, s.day, s.volume].map((c) => df.col(c)!);
  const out = df.col(EXCLUDED);
  const by = new Map<string, {group: string, points: [number, number][]}>();
  for (let i = 0; i < df.rowCount; i++) {
    if (ticked(out, i) || v.isNone(i) || d.isNone(i))
      continue;
    const animal = String(a.get(i));
    const line = by.get(animal) ?? {group: String(g.get(i)), points: []};
    line.points.push([d.get(i), v.get(i)]);
    by.set(animal, line);
  }
  return [...by.entries()].flatMap(([animal, {group, points}]) => {
    points.sort((p, q) => p[0] - q[0]);
    const v0 = points.reduce((at, p) => closer(p[0], at[0]) ? p : at)[1];
    return v0 > 0 ? [{animal, group, days: points.map((p) => p[0]),
      change: points.map((p) => Math.round(1000 * (p[1] / v0 - 1)) / 10)}] : [];
  });
}

export class SpiderChart extends Control {
  private readonly chart: echarts.ECharts;

  /** `lines` and the groups' colours, redrawn as they change. */
  constructor(data: ReadonlySignal<{lines: SpiderLine[], colors: {[group: string]: string}} | undefined>,
    picking?: Picking) {
    super();
    this.root.classList.add('pmx-spider-chart');
    this.root.style.position = 'relative';
    const host = ui.div([], {style: {position: 'absolute', inset: '0'}});
    this.root.append(host);
    this.chart = initQuietly(host);
    // (resized: drawn again, the legend's lines depend on the width)
    const width = signal(0);
    this.subscribe(ui.onSizeChanged(host), () => {
      this.chart.resize();
      width.value = this.chart.getWidth();
    });
    this.own(() => this.chart.dispose());
    pickOn(this.chart, picking);
    this.effect(() => {
      const x = data.value;
      const on = picking?.animals.value ?? new Set<string>();
      if (!x) {
        this.chart.clear();
        return;
      }
      const groups = [...new Set(x.lines.map((l) => l.group))];
      const legend = legendLayout(groups, width.value || this.chart.getWidth());
      const dim = (animal: string) => on.size > 0 && !on.has(animal);
      const fallback = DG.Color.toHtml(DG.Color.getCategoricalColor(0));
      // each animal is a series named after its group: the legend toggles a group's animals together
      this.chart.setOption({
        animation: false,
        grid: {left: 56, right: 16, top: 8 + 24 * legend.lines, bottom: 40},
        // the legend in full colour whatever is dimmed
        // wrapped over the lines it needs, every name whole, the plot below it (as the band chart: a paged legend
        // cut names and hid groups, QA 0.1.3 round 3 designer 19)
        legend: {type: legend.type, top: 4, left: 20, right: 20, data: groups, itemStyle: {opacity: 1},
          lineStyle: {opacity: 1}},
        tooltip: {trigger: 'item', formatter: (p: any) => `${esc(animalOf(p) ?? '')} (${esc(p.seriesName)})` +
          (p.data ? `<br>Day ${p.data.value[0]}: ${p.data.value[1] > 0 ? '+' : ''}${p.data.value[1]}%` : '')},
        xAxis: {type: 'value', name: 'Day after randomization', nameLocation: 'middle', nameGap: 26},
        yAxis: {type: 'value', name: 'Change from randomization (%)', nameLocation: 'middle', nameGap: 42},
        series: x.lines.map((l) => ({type: 'line', id: l.animal, name: l.group, color: x.colors[l.group] ?? fallback,
          symbolSize: 4, triggerLineEvent: true, z: dim(l.animal) ? 1 : 2,
          lineStyle: {width: on.has(l.animal) ? 2.5 : 1.2, opacity: dim(l.animal) ? 0.12 : 0.85},
          itemStyle: {opacity: dim(l.animal) ? 0.12 : 1}, emphasis: {focus: 'series'},
          data: l.days.map((day, i) => ({value: [day, l.change[i]], animal: l.animal}))})),
        // the lines anew, the rest merged: a group hidden from the legend stays hidden when an animal is picked
      }, {replaceMerge: ['series']});
    });
  }
}

/** A waterfall bar: an animal, its group and its value. */
export interface Bar {animal: string, group: string, value: number}

/** Every animal as a thin bar in its group's colour, sorted (`order`: desc, the largest first), the animal in the
 *  tooltip; a group's legend entry shows or hides its bars; `onMenu`: an animal's right click. `compact`: a context
 *  panel's few pixels, no legend nor axis names (the tooltip names the animal and its value). */
export class WaterfallChart extends Control {
  private readonly chart: echarts.ECharts;

  constructor(
    data: ReadonlySignal<{bars: Bar[], colors: {[group: string]: string}, title: string, order: 'asc' | 'desc',
      groups: string[]} | undefined>,
    onMenu?: (animal: string, e: MouseEvent) => void, picking?: Picking, o: {compact?: boolean} = {}) {
    super();
    this.root.classList.add('pmx-waterfall-chart');
    this.root.style.position = 'relative';
    const host = ui.div([], {style: {position: 'absolute', inset: '0'}});
    this.root.append(host);
    this.chart = initQuietly(host);
    const width = signal(0);
    this.subscribe(ui.onSizeChanged(host), () => {
      this.chart.resize();
      width.value = this.chart.getWidth();
    });
    this.own(() => this.chart.dispose());
    pickOn(this.chart, picking);
    if (onMenu) {
      this.chart.on('contextmenu', (p: any) => {
        p.event?.event?.preventDefault();
        if (p.data?.animal)
          onMenu(p.data.animal, p.event.event);
      });
    }
    this.effect(() => {
      const x = data.value;
      const on = picking?.animals.value ?? new Set<string>();
      if (!x) {
        this.chart.clear();
        return;
      }
      const sorted = [...x.bars].sort((a, b) => x.order === 'desc' ? b.value - a.value : a.value - b.value);
      const fmt = (v: number) => (Math.abs(v) >= 10 ? v.toFixed(0) : v.toPrecision(2)).replace('-', '−');
      const legend = legendLayout(x.groups, width.value || this.chart.getWidth());
      // a series per group, stacked on one category per animal: each bar in its group's colour, the legend toggling
      this.chart.setOption({
        animation: false,
        grid: o.compact ? {left: 4, right: 4, top: 4, bottom: 4} : {left: 56, right: 16, top: 8 + 24 * legend.lines,
          bottom: 12},
        legend: o.compact ? {show: false} : {type: legend.type, top: 4, left: 20, right: 20, data: x.groups},
        tooltip: {trigger: 'item', formatter: (p: any) => o.compact ? `${esc(p.data.animal)}: ${fmt(p.data.value)}` :
          `${esc(p.data.animal)} (${esc(p.seriesName)})<br>${esc(x.title)}: ${fmt(p.data.value)}`},
        xAxis: {type: 'category', data: sorted.map((b) => b.animal), axisLabel: {show: false},
          axisTick: {show: false}},
        yAxis: o.compact ? {type: 'value', axisLabel: {show: false}, splitLine: {show: false}} :
          {type: 'value', name: x.title, nameLocation: 'middle', nameGap: 42},
        series: x.groups.map((g) => ({type: 'bar', name: g, stack: 'all', color: x.colors[g], barCategoryGap: '20%',
          data: sorted.map((b) => b.group === g ? {value: b.value, animal: b.animal,
            itemStyle: {opacity: on.size > 0 && !on.has(b.animal) ? 0.2 : 1}} : null)})),
      }, {replaceMerge: ['series']});
    });
  }
}
