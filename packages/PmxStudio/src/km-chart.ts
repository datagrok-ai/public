/* Analyze's Kaplan-Meier: the share of each group's animals still under the endpoint, day by day, in the groups'
   colours, on an axis fixed from day 0 to twice the last day measured, so following further never rescales it. One
   control: under the axis, from the last day measured to the axis' end, a slider whose handle sits under the day the
   study is followed to; past the last day measured, the chart is shaded as projected. echarts, as the band chart. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Control, Tooltip} from '@datagrok-libraries/u2';
import type {ReadonlySignal, Signal} from '@datagrok-libraries/u2';
import * as echarts from 'echarts/core';

import {initQuietly} from './band-chart';
import {KM, Survival} from './survival';

const GRID = {left: 48, right: 16, top: 32, bottom: 62};
/** half the width of the range input's handle: the input overhangs the axis by it, so the handle's centre is the day */
const THUMB = 8;

export class KmChart extends Control {
  private readonly chart: echarts.ECharts;

  /** `survival` drawn as it changes; `horizon`: the day followed to (null: as measured), set by the slider. */
  constructor(survival: ReadonlySignal<Survival | undefined>, horizon: Signal<number | null>) {
    super();
    this.root.classList.add('pmx-km-chart');
    this.root.style.position = 'relative';
    const host = ui.div([], {style: {position: 'absolute', inset: '0'}});
    const slider = document.createElement('input');
    slider.type = 'range';
    slider.step = '1';
    slider.className = 'pmx-km-slider pmx-noprint';
    slider.setAttribute('aria-label', 'Followed to day');
    Tooltip.bind(slider, 'Drag right to carry each animal still under the endpoint on at the rate of its last three ' +
      'measurements, to that day: a projection, not a result. Back to the left end: as measured.');
    const measured = ui.div([], 'pmx-km-measured');
    this.root.append(host, slider, measured);
    this.chart = initQuietly(host);
    let last = 0;
    // the slider under the axis' second half: its ends at the last day measured and at the axis' end
    const place = () => {
      if (!last)
        return;
      const [x0, x1] = [last, 2 * last].map((d) => this.chart.convertToPixel({xAxisIndex: 0}, d) as number);
      Object.assign(slider.style, {left: `${x0 - THUMB}px`, width: `${x1 - x0 + 2 * THUMB}px`});
      Object.assign(measured.style, {left: `${GRID.left}px`, width: `${x0 - GRID.left - THUMB - 6}px`});
    };
    slider.addEventListener('input', () => {
      const v = Number(slider.value);
      horizon.value = v > last ? v : null;
    });
    const size = ui.onSizeChanged(host).subscribe(() => {
      this.chart.resize();
      place();
    });
    this.own(() => {
      size.unsubscribe();
      this.chart.dispose();
    });
    this.effect(() => {
      const k = survival.value;
      slider.hidden = measured.hidden = !k;
      if (!k) {
        this.chart.clear();
        return;
      }
      last = k.lastDay;
      slider.min = String(last);
      slider.max = String(2 * last);
      if (horizon.peek() === null)
        slider.value = String(last);
      measured.textContent = k.horizon === undefined ? `Measured to day ${last}` : `Measured to day ${last} · ` +
        `projected to day ${k.horizon}`;
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
        xAxis: {type: 'value', min: 0, max: 2 * last, name: 'Day after randomization', nameLocation: 'middle',
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
      place();
    });
  }
}
