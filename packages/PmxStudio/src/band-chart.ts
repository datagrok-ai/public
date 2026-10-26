/* The band chart of the Plan and Analyze views: the area between a lower and an upper column along X (a predicted
   range, an observed one), with a line through it, one band per category of the split column. A control of this
   package, not a platform viewer: the line chart cannot fill between two columns (its whiskers and moving-average
   band are other things), and a general viewer for it is the platform's call. echarts draws it; each band is a
   custom series drawing one polygon, since a stack on two value axes stacks X, not Y. */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Control} from '@datagrok-libraries/u2';
import type {ReadonlySignal} from '@datagrok-libraries/u2';
import * as echarts from 'echarts/core';
import {CustomChart, LineChart, ScatterChart} from 'echarts/charts';
import {GridComponent, LegendComponent, MarkLineComponent, TooltipComponent} from 'echarts/components';
import {CanvasRenderer} from 'echarts/renderers';
import type {EChartsOption} from 'echarts';

echarts.use([CustomChart, LineChart, ScatterChart, GridComponent, LegendComponent, MarkLineComponent,
  TooltipComponent, CanvasRenderer]);

export interface BandChartOptions {
  x: string;
  lower: string;
  upper: string;
  /** a line through each band (a median) */
  line?: string;
  /** one band per category, in the column's colour coding */
  split?: string;
  xTitle?: string;
  yTitle?: string;
  /** a logarithmic Y axis: the view's Linear / Log switch */
  log?: ReadonlySignal<boolean>;
  /** markers on the line */
  markers?: boolean;
}

/** One category's points, sorted by X; `row` is one of its rows, where its colour is read. */
export interface Band {name: string, x: number[], lower: number[], upper: number[], line: number[] | null, row: number}

/** A constant formula line of the table: `${X} = 19` vertical, `${any other column} = 1500` horizontal. */
export interface RefLine {value: number, vertical: boolean, title: string, color: string, dashed: boolean}

const DEFAULT_FORMAT = '#0.###';
const TEXT = '#4d5261';

/** The chart over the table `table` gives (none: an empty chart), redrawn when it or the Log switch changes. */
export class BandChart extends Control {
  /** What the last draw showed. */
  bands: Band[] = [];
  private readonly chart: echarts.ECharts;

  constructor(table: ReadonlySignal<DG.DataFrame | undefined>, private readonly options: BandChartOptions) {
    super();
    this.root.classList.add('pmx-band-chart');
    this.root.style.position = 'relative';
    const host = ui.div([], {style: {position: 'absolute', inset: '0'}});
    this.root.append(host);
    this.chart = initQuietly(host);
    const size = ui.onSizeChanged(host).subscribe(() => this.chart.resize());
    this.own(() => {
      size.unsubscribe();
      this.chart.dispose();
    });
    this.effect(() => this.draw(table.value, options.log?.value ?? false));
  }

  private draw(df: DG.DataFrame | undefined, log: boolean): void {
    const o = this.options;
    const [x, lower, upper] = [o.x, o.lower, o.upper].map((n) => df?.col(n) ?? null);
    if (!df || !x || !lower || !upper) {
      this.bands = [];
      this.chart.clear();
      return;
    }
    const line = o.line ? df.col(o.line) : null;
    const split = o.split ? df.col(o.split) : null;
    this.bands = collectBands(df, x, lower, upper, line, split);
    const y = line ?? upper;
    const fmt = {xUnits: x.meta.units ?? '', yUnits: y.meta.units ?? '', x: x.meta.format || DEFAULT_FORMAT,
      y: y.meta.format || DEFAULT_FORMAT, xTitle: o.xTitle || x.name};
    const xTitle = axisTitle(fmt.xTitle, fmt.xUnits);
    const yTitle = axisTitle(o.yTitle || (line ? line.name : `${lower.name} – ${upper.name}`), fmt.yUnits, log);
    // the split column's colour coding when it has one, else the categorical palette in band order
    const coded = split !== null && split.meta.colors.getType() !== DG.COLOR_CODING_TYPE.OFF ? split : null;
    const colors = this.bands.map((b, k) =>
      DG.Color.toHtml(coded ? DG.Color.getRowColor(coded, b.row) : DG.Color.getCategoricalColor(k)));
    this.chart.setOption(chartOption(this.bands, colors, refLines(df.meta.formulaLines.items, x.name),
      {xTitle, yTitle, log, markers: o.markers ?? false, fmt}), true);
  }
}

/** echarts.init without its warning about a host not sized yet (it is resized when it is). */
function initQuietly(host: HTMLElement): echarts.ECharts {
  const warn = console.warn;
  console.warn = () => {};
  try {
    return echarts.init(host);
  } finally {
    console.warn = warn;
  }
}

type Look = {xTitle: string, yTitle: string, log: boolean, markers: boolean,
  fmt: {xUnits: string, yUnits: string, x: string, y: string, xTitle: string}};

function chartOption(bands: Band[], colors: string[], refs: RefLine[], look: Look): EChartsOption {
  const series: object[] = [];
  bands.forEach((b, k) => {
    const color = colors[k];
    const pairs = (ys: number[]) => b.x.map((x, i) => [x, ys[i]]);
    // the polygon, upper edge forward and lower edge back, drawn once on the first item; every item's lower and
    // upper still set the Y extent
    const outline = [...pairs(b.upper), ...pairs(b.lower).reverse()];
    series.push({name: b.name, type: 'custom', silent: true, clip: true, itemStyle: {color},
      data: b.x.map((x, i) => [x, b.lower[i], b.upper[i]]), encode: {x: 0, y: [1, 2]},
      renderItem: (params: any, api: any) => params.dataIndex !== 0 ? null : {type: 'polygon',
        shape: {points: outline.map((p) => api.coord(p))}, style: {fill: color, opacity: 0.25}}});
    if (b.line) {
      series.push({name: b.name, type: 'line', data: pairs(b.line), symbol: look.markers ? 'circle' : 'none',
        lineStyle: {color, width: 2}, itemStyle: {color}});
    }
  });
  if (refs.length && bands.length) {
    // the reference lines ride an empty series; its radius-0 points keep a horizontal one inside the Y axis
    series.push({type: 'scatter', symbolSize: 0, silent: true, tooltip: {show: false},
      data: refs.filter((r) => !r.vertical).map((r) => [bands[0].x[0], r.value]),
      markLine: {silent: true, symbol: 'none', animation: false, data: refs.map((r) => ({
        [r.vertical ? 'xAxis' : 'yAxis']: r.value, name: r.title,
        lineStyle: {color: r.color, type: r.dashed ? 'dashed' : 'solid', width: 1},
        // a vertical line's title above the plot, a horizontal one's along it at its start
        label: {formatter: r.title, position: r.vertical ? 'end' : 'insideStartTop', color: TEXT,
          ...(r.vertical ? {align: 'right'} : {})},
      }))}});
  }
  const axisName = {nameLocation: 'middle', nameTextStyle: {color: TEXT}} as const;
  return {
    animation: false,
    // containLabel holds the tick labels, not the axis names: room for the Y name, the X name and a vertical
    // line's title
    grid: {left: 44, right: 20, top: (bands.length > 1 ? 32 : 12) + (refs.some((r) => r.vertical) ? 18 : 0),
      bottom: 36, containLabel: true},
    legend: {show: bands.length > 1, type: 'scroll', top: 4},
    tooltip: {trigger: 'axis', formatter: (params: any) => tooltip(params as any[], bands, look.fmt)},
    xAxis: {type: 'value', scale: true, name: look.xTitle, nameGap: 26, ...axisName},
    // on a log axis the minor ticks are labelled too (2, 5 per decade): a decade alone reads as a linear axis
    yAxis: {type: look.log ? 'log' : 'value', scale: true, name: look.yTitle, nameGap: 52, ...axisName,
      ...(look.log ? {minorTick: {show: true, splitNumber: 2},
        minorSplitLine: {show: true, lineStyle: {opacity: 0.4}}} : {})},
    series,
  } as EChartsOption;
}

/** X, then lower / line / upper of every band at that point (each series at its own index: bands of different
 *  categories need not share their X values). */
function tooltip(params: any[], bands: Band[], fmt: Look['fmt']): string {
  if (!params.length)
    return '';
  const rows = new Set<string>();
  const u = fmt.yUnits ? ` ${fmt.yUnits}` : '';
  const y = (v: number) => DG.format(v, fmt.y);
  for (const p of params) {
    const b = bands.find((b) => b.name === p.seriesName);
    const i = p.dataIndex;
    if (!b || i >= b.x.length)
      continue;
    const range = `${y(b.lower[i])} – ${y(b.upper[i])}${u}`;
    rows.add(`${b.name === '' ? '' : `<b>${b.name}</b>: `}${b.line ? `${y(b.line[i])}${u} (${range})` : range}`);
  }
  const xu = fmt.xUnits ? ` ${fmt.xUnits}` : '';
  return [`${fmt.xTitle}: ${DG.format(params[0].data[0], fmt.x)}${xu}`, ...rows].join('<br>');
}

/** "Tumor volume (mm³)", "Tumor volume (mm³, log scale)": the units once, when given and not already in the title,
 *  and the scale when it is not linear. */
export function axisTitle(title: string, units: string, log = false): string {
  const notes = [units && !title.includes(units) ? units : '', log ? 'log scale' : ''].filter((n) => n);
  return notes.length ? `${title} (${notes.join(', ')})` : title;
}

/** The table's formula lines the chart draws: the visible constant ones. Other formulas need the platform's
 *  evaluator and are skipped; `DG.FormulaLinesHelper.getMeta` names the column but does not split off the constant,
 *  hence the pattern. */
export function refLines(items: DG.FormulaLine[], xName: string): RefLine[] {
  const out: RefLine[] = [];
  for (const f of items) {
    const m = /^\s*\$\{(.+?)\}\s*=\s*(-?[\d.]+(?:[eE][-+]?\d+)?)\s*$/.exec(f.formula ?? '');
    if (!m || f.visible === false || (f.type && f.type !== 'line'))
      continue;
    out.push({value: Number(m[2]), vertical: m[1] === xName, title: f.title ?? '', color: f.color || '#9aa5b1',
      dashed: (f.style ?? 'dashed') !== 'solid'});
  }
  return out;
}

/** The filtered rows with X, lower and upper present, grouped by the split category and sorted by X. */
export function collectBands(df: DG.DataFrame, x: DG.Column, lower: DG.Column, upper: DG.Column,
  line: DG.Column | null, split: DG.Column | null): Band[] {
  const byName = new Map<string, number[]>();
  for (let i = df.filter.findNext(-1, true); i !== -1; i = df.filter.findNext(i, true)) {
    if (x.isNone(i) || lower.isNone(i) || upper.isNone(i))
      continue;
    const name = split ? (split.isNone(i) ? '' : String(split.get(i))) : '';
    let rows = byName.get(name);
    if (!rows)
      byName.set(name, rows = []);
    rows.push(i);
  }
  return [...byName].map(([name, rows]) => {
    rows.sort((a, b) => x.get(a) - x.get(b));
    const values = (c: DG.Column) => rows.map((i) => c.get(i) as number);
    return {name, row: rows[0], x: values(x), lower: values(lower), upper: values(upper),
      line: line ? rows.map((i) => line.isNone(i) ? NaN : line.get(i)) : null};
  });
}
