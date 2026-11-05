/* The spider plot of Analyze: each animal's tumor as a % change from its volume at randomization, day by day, one
   line per animal in its group's colour; a group's legend entry shows or hides all its animals. echarts, as the band
   chart: the platform's line chart colours lines by its split, and here the split (the animal) is not the colour
   (the group). */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Control} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';
import * as echarts from 'echarts/core';

import {initQuietly} from './band-chart';
import {EXCLUDED, ticked} from './study-file';
import type {StudyInputs} from './explore';

/** The animals selected (marked on the chart, the others dimmed) and what a click on one does. */
export interface Picking {animals: ReadonlySignal<ReadonlySet<string>>, pick(animal: string, e?: MouseEvent): void}

/** Where an echarts click on an animal goes: to picking.pick, with the keys held. */
function pickOn(chart: echarts.ECharts, picking?: Picking): void {
  if (picking)
    chart.on('click', (p: any) => p.data?.animal && picking.pick(p.data.animal, p.event?.event));
}

/** One animal's line: its days and its % change from its first volume (day 0, randomization, when measured). */
export interface SpiderLine {animal: string, group: string, days: number[], change: number[]}

/** The study's animals as spider lines, its excluded measurements left out; an animal whose first volume is not
 *  positive has no line. */
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
    const v0 = points[0][1];
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
    const size = ui.onSizeChanged(host).subscribe(() => this.chart.resize());
    this.own(() => {
      size.unsubscribe();
      this.chart.dispose();
    });
    pickOn(this.chart, picking);
    this.effect(() => {
      const x = data.value;
      const on = picking?.animals.value ?? new Set<string>();
      if (!x) {
        this.chart.clear();
        return;
      }
      const dim = (animal: string) => on.size > 0 && !on.has(animal);
      const fallback = DG.Color.toHtml(DG.Color.getCategoricalColor(0));
      // each animal is a series named after its group: the legend toggles a group's animals together
      this.chart.setOption({
        animation: false,
        grid: {left: 56, right: 16, top: 32, bottom: 40},
        // the legend in full colour whatever is dimmed
        legend: {top: 0, data: [...new Set(x.lines.map((l) => l.group))], itemStyle: {opacity: 1},
          lineStyle: {opacity: 1}},
        tooltip: {trigger: 'item', formatter: (p: any) => `${p.data.animal} (${p.seriesName})<br>Day ` +
          `${p.data.value[0]}: ${p.data.value[1] > 0 ? '+' : ''}${p.data.value[1]}%`},
        xAxis: {type: 'value', name: 'Day after randomization', nameLocation: 'middle', nameGap: 26},
        yAxis: {type: 'value', name: 'Change from randomization (%)', nameLocation: 'middle', nameGap: 42},
        series: x.lines.map((l) => ({type: 'line', name: l.group, color: x.colors[l.group] ?? fallback,
          symbolSize: 4, triggerLineEvent: true, z: dim(l.animal) ? 1 : 2,
          lineStyle: {width: on.has(l.animal) ? 2.5 : 1.2, opacity: dim(l.animal) ? 0.12 : 0.85},
          itemStyle: {opacity: dim(l.animal) ? 0.12 : 1}, emphasis: {focus: 'series'},
          data: l.days.map((day, i) => ({value: [day, l.change[i]], animal: l.animal}))})),
      }, true);
    });
  }
}

/** A waterfall bar: an animal, its group and its value. */
export interface Bar {animal: string, group: string, value: number}

/** Every animal as a thin bar in its group's colour, sorted (`order`: desc, the largest first), the animal in the
 *  tooltip; a group's legend entry shows or hides its bars; `onMenu`: an animal's right click. */
export class WaterfallChart extends Control {
  private readonly chart: echarts.ECharts;

  constructor(
    data: ReadonlySignal<{bars: Bar[], colors: {[group: string]: string}, title: string, order: 'asc' | 'desc',
      groups: string[]} | undefined>,
    onMenu?: (animal: string, e: MouseEvent) => void, picking?: Picking) {
    super();
    this.root.classList.add('pmx-waterfall-chart');
    this.root.style.position = 'relative';
    const host = ui.div([], {style: {position: 'absolute', inset: '0'}});
    this.root.append(host);
    this.chart = initQuietly(host);
    const size = ui.onSizeChanged(host).subscribe(() => this.chart.resize());
    this.own(() => {
      size.unsubscribe();
      this.chart.dispose();
    });
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
      const fmt = (v: number) => Math.abs(v) >= 10 ? v.toFixed(0) : v.toPrecision(2);
      // a series per group, stacked on one category per animal: each bar in its group's colour, the legend toggling
      this.chart.setOption({
        animation: false,
        grid: {left: 56, right: 16, top: 32, bottom: 12},
        legend: {top: 0, data: x.groups},
        tooltip: {trigger: 'item', formatter: (p: any) => `${p.data.animal} (${p.seriesName})<br>${x.title}: ` +
          `${fmt(p.data.value)}`},
        xAxis: {type: 'category', data: sorted.map((b) => b.animal), axisLabel: {show: false},
          axisTick: {show: false}},
        yAxis: {type: 'value', name: x.title, nameLocation: 'middle', nameGap: 42},
        series: x.groups.map((g) => ({type: 'bar', name: g, stack: 'all', color: x.colors[g], barCategoryGap: '20%',
          data: sorted.map((b) => b.group === g ? {value: b.value, animal: b.animal,
            itemStyle: {opacity: on.size > 0 && !on.has(b.animal) ? 0.2 : 1}} : null)})),
      }, true);
    });
  }
}
