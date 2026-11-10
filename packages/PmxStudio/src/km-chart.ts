/* Analyze's Kaplan-Meier: the share of each group's animals still under the endpoint, day by day, in the groups'
   colours, from day 0 to the day the study is followed to, the chart always its full width. One control, under the
   plot and as wide: a slider from the last day measured (as measured) to twice it; dragged right, the axis widens
   with it and what lies past the last day measured is shaded as projected. echarts, as the band chart. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Control, Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal, Signal} from '@datagrok-libraries/u2';
import * as echarts from 'echarts/core';

import {initQuietly} from './band-chart';
import {KM, Survival} from './survival';

const GRID = {left: 48, right: 16, top: 32, bottom: 40};
/** the slider's row under the chart */
const ROW = 26;

export class KmChart extends Control {
  private readonly chart: echarts.ECharts;

  /** `survival` drawn as it changes; `horizon`: the day followed to (null: as measured), set by the slider. */
  constructor(survival: ReadonlySignal<Survival | undefined>, horizon: Signal<number | null>) {
    super();
    this.root.classList.add('pmx-km-chart');
    this.root.style.position = 'relative';
    const host = ui.div([], {style: {position: 'absolute', inset: `0 0 ${ROW}px 0`}});
    const slider = document.createElement('input');
    slider.type = 'range';
    slider.step = '1';
    slider.className = 'pmx-km-slider pmx-noprint';
    slider.setAttribute('aria-label', 'Followed to day');
    Tooltip.bind(slider, 'Drag right to carry each animal still under the endpoint on at the rate of its last three ' +
      'measurements, to that day: a projection, not a result. Back to the left end: as measured.');
    const measured = ui.div([], 'pmx-km-measured');
    const row = ui.div([ui.span(['Followed to']), slider, measured], 'pmx-km-row pmx-noprint');
    Object.assign(row.style, {left: `${GRID.left}px`, right: `${GRID.right}px`, height: `${ROW}px`});
    this.root.append(host, row);
    this.chart = initQuietly(host);
    let last = 0;
    slider.addEventListener('input', () => {
      const v = Number(slider.value);
      horizon.value = v > last ? v : null;
    });
    const size = ui.onSizeChanged(host).subscribe(() => this.chart.resize());
    this.own(() => {
      size.unsubscribe();
      this.chart.dispose();
    });
    this.effect(() => {
      const k = survival.value;
      row.hidden = !k;
      if (!k) {
        this.chart.clear();
        return;
      }
      last = k.lastDay;
      slider.min = String(last);
      slider.max = String(2 * last);
      // a horizon set for another study, or before more days were measured, is outside this slider: as measured
      const h = horizon.peek();
      if (h !== null && (h <= last || h > 2 * last))
        horizon.value = null;
      if (horizon.peek() === null)
        slider.value = String(last);
      measured.textContent = k.horizon === undefined ? `day ${last}, as measured` : `day ${k.horizon}, projected` +
        ` from day ${last}`;
      measured.classList.toggle('is-projected', k.horizon !== undefined);
      const df = k.steps;
      const [day, share, group] = [KM.day, KM.share, KM.group].map((c) => df.col(c)!);
      const by = new Map<string, [number, number][]>();
      for (let i = 0; i < df.rowCount; i++)
        by.set(group.get(i), [...by.get(group.get(i)) ?? [], [day.get(i), share.get(i)]]);
      const groups = k.groups.filter((g) => by.has(g));
      const projected = k.horizon !== undefined ? [[{xAxis: last, name: 'Projected'}, {xAxis: k.horizon}]] : [];
      this.chart.setOption({
        animation: false,
        grid: GRID,
        legend: {top: 0, data: groups},
        tooltip: {trigger: 'axis', valueFormatter: (v: number) => `${Math.round(v)}%`},
        xAxis: {type: 'value', min: 0, max: k.horizon ?? last, name: 'Day after randomization', nameLocation: 'middle',
          nameGap: 24},
        yAxis: {type: 'value', min: 0, max: 100, name: 'Under the endpoint (%)', nameLocation: 'middle',
          nameGap: 34},
        series: groups.map((g, i) => ({type: 'line', name: g, showSymbol: false, color: k.colors[g] ??
          DG.Color.toHtml(DG.Color.getCategoricalColor(i)), lineStyle: {width: 2}, data: by.get(g),
        // the measured end, and what lies past it shaded as projected (on the first series: drawn once)
        ...(i === 0 ? {markLine: {silent: true, symbol: 'none', lineStyle: {color: '#9aa5b1', type: 'dashed'},
          label: {formatter: `last measured, day ${last}`, position: 'insideEndTop', color: '#717581', fontSize: 10},
          data: [{xAxis: last}]},
        markArea: {silent: true, itemStyle: {color: 'rgba(183, 121, 31, 0.16)'},
          label: {color: '#b7791f', fontSize: 10, position: 'insideTop'}, data: projected}} : {})})),
      }, true);
    });
  }
}
