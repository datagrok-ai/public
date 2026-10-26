/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Observable} from 'rxjs';
import type {EChartsOption} from 'echarts';

import {EChartViewer} from '../echart/echart-viewer';
import {ERROR_CLASS, MessageHandler} from '../../utils/utils';
import {bandChartStatus} from './band-chart-status';

type AxisType = 'linear' | 'logarithmic';

/** One category's points, sorted by X. */
export interface Band {name: string, x: number[], lower: number[], upper: number[], line: number[] | null,
  /** a row of the band: the split column's colour coding is read at it */
  row: number}

const NUMERICAL = {columnTypeFilter: DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME};

/** The area between a lower and an upper column along X (a prediction or confidence band), with an optional
 * line through it, one band per category of the split column. echarts has no area between two series (a stack
 * on two value axes stacks X), so each band is a custom series drawing one polygon in axis coordinates. */
@grok.decorators.viewer({
  name: 'Band chart',
  description: 'Shows the area between a lower and an upper value along X (a prediction or confidence band), with an optional line through it, one band per category',
  icon: 'icons/band-chart-viewer.svg',
})
export class BandChartViewer extends EChartViewer {
  xColumnName: string;
  lowerColumnName: string;
  upperColumnName: string;
  lineColumnName: string;
  splitColumnName: string;
  yAxisType: AxisType;
  xAxisTitle: string;
  yAxisTitle: string;
  bandOpacity: number;
  showMarkers: boolean;
  /** What the last render drew, for the status. */
  bands: Band[] = [];
  /** The table's constant formula lines the last render drew. */
  refLines: RefLine[] = [];

  get onRendered(): Observable<void> {return this.signals.rendered;}

  get isRenderPending(): boolean {return this.signals.pending;}

  get renderError(): string | null {return MessageHandler._getMessage(this.root);}

  getWidgetStatus(): DG.IWidgetStatus {return bandChartStatus(this);}

  constructor() {
    super();
    this.xColumnName = this.string('xColumnName', null, {...NUMERICAL, category: 'Data'});
    this.lowerColumnName = this.string('lowerColumnName', null, {...NUMERICAL, category: 'Data', description: 'The lower edge of the band'});
    this.upperColumnName = this.string('upperColumnName', null, {...NUMERICAL, category: 'Data', description: 'The upper edge of the band'});
    this.lineColumnName = this.string('lineColumnName', null, {...NUMERICAL, category: 'Data', description: 'A line drawn through the band (a median, a fit); none by default'});
    this.splitColumnName = this.string('splitColumnName', null, {columnTypeFilter: DG.COLUMN_TYPE_FILTER.CATEGORICAL, category: 'Data', description: 'One band per category'});
    this.yAxisType = <AxisType> this.string('yAxisType', 'linear', {choices: ['linear', 'logarithmic'], category: 'Axes'});
    this.xAxisTitle = this.string('xAxisTitle', '', {category: 'Axes', description: 'Empty: the X column name. The column units follow in brackets'});
    this.yAxisTitle = this.string('yAxisTitle', '', {category: 'Axes', description: 'Empty: the line column name, else the lower and upper names. The units follow in brackets'});
    this.bandOpacity = this.float('bandOpacity', 0.25, {min: 0, max: 1, category: 'Style'});
    this.showMarkers = this.bool('showMarkers', false, {category: 'Style', description: 'Markers on the line'});
  }

  onTableAttached(): void {
    const columns = this.dataFrame.columns;
    const nums = Array.from(columns.numericalNoDateTime).map((c) => c.name);
    const missing = (name: string) => !name || !columns.contains(name);
    // a table with no band declared: the first three numerical columns, so the viewer shows something
    if (missing(this.xColumnName) && missing(this.lowerColumnName) && missing(this.upperColumnName) && nums.length >= 3)
      [this.xColumnName, this.lowerColumnName, this.upperColumnName] = nums;
    // resubscribe drops the base's size subscription too: the resize is ours again, as in Radar
    this.resubscribe(() => {
      this.addSelectionOrDataSubs();
      this.subs.push(this.signals.debounce(this.dataFrame.filter.onChanged, 50).subscribe(() => this.render()));
      // formula lines and units are table and column tags
      this.subs.push(this.signals.debounce(this.dataFrame.onMetadataChanged, 50).subscribe(() => this.render()));
      this.subs.push(this.signals.debounce(ui.onSizeChanged(this.root), 50).subscribe(() => {
        const release = this.signals.hold();
        requestAnimationFrame(() => {
          this.chart.resize();
          release();
          this.signals.settle();
        });
      }));
    });
    this.render();
  }

  onPropertyChanged(p: DG.Property | null): void {
    if (p?.name === 'table')
      this.updateTable();
    this.render();
  }

  render(): void {
    this.signals.render(() => this.draw());
  }

  private draw(): void {
    if (!this.dataFrame)
      return;
    const df = this.dataFrame;
    const [x, lower, upper] = [this.xColumnName, this.lowerColumnName, this.upperColumnName].map((n) => n ? df.col(n) : null);
    if (!x || !lower || !upper) {
      this.bands = [];
      this.chart.clear();
      MessageHandler._showMessage(this.root, 'The Band chart needs an X, a lower and an upper numerical column.', ERROR_CLASS);
      return;
    }
    MessageHandler._removeMessage(this.root, ERROR_CLASS);
    const line = this.lineColumnName ? df.col(this.lineColumnName) : null;
    const split = this.splitColumnName ? df.col(this.splitColumnName) : null;
    this.bands = collectBands(df, x, lower, upper, line, split);
    this.split = split;
    this.refLines = refLines(df.meta.formulaLines.items, x.name);
    this.xUnits = x.meta.units ?? '';
    this.yUnits = (line ?? upper).meta.units ?? '';
    // the tooltip's numbers as their columns' format says, when they say it
    this.xFormat = x.meta.format || DEFAULT_FORMAT;
    this.yFormat = (line ?? upper).meta.format || DEFAULT_FORMAT;
    const xTitle = axisTitle(this.xAxisTitle || x.name, this.xUnits);
    const yTitle = axisTitle(this.yAxisTitle || (line ? line.name : `${lower.name} – ${upper.name}`), this.yUnits,
      this.yAxisType === 'logarithmic');
    this.chart.setOption(this.chartOption(xTitle, yTitle), true);
  }

  private split: DG.Column | null = null;
  private xUnits = '';
  private yUnits = '';
  private xFormat = DEFAULT_FORMAT;
  private yFormat = DEFAULT_FORMAT;

  private chartOption(xTitle: string, yTitle: string): EChartsOption {
    const log = this.yAxisType === 'logarithmic';
    const series: object[] = [];
    // the split column's colour coding when it has one (as the platform's viewers read it), else the palette
    const coded = this.split !== null && this.split.meta.colors.getType() !== DG.COLOR_CODING_TYPE.OFF ? this.split : null;
    this.bands.forEach((b, k) => {
      const color = DG.Color.toHtml(coded ? DG.Color.getRowColor(coded, b.row) : DG.Color.getCategoricalColor(k));
      const pairs = (ys: number[]) => b.x.map((x, i) => [x, ys[i]]);
      // the polygon upper edge forward, lower edge back; drawn once, on the first item, while every item's
      // lower and upper still set the Y extent
      const outline = [...pairs(b.upper), ...pairs(b.lower).reverse()];
      series.push({name: b.name, type: 'custom', silent: true, clip: true, itemStyle: {color},
        data: b.x.map((x, i) => [x, b.lower[i], b.upper[i]]), encode: {x: 0, y: [1, 2]},
        renderItem: (params: any, api: any) => params.dataIndex !== 0 ? null : {type: 'polygon',
          shape: {points: outline.map((p) => api.coord(p))}, style: {fill: color, opacity: this.bandOpacity}}});
      if (b.line) {
        series.push({name: b.name, type: 'line', data: pairs(b.line), symbol: this.showMarkers ? 'circle' : 'none',
          lineStyle: {color, width: 2}, itemStyle: {color}});
      }
    });
    if (this.refLines.length && this.bands.length) {
      // the reference lines ride an empty series; their points (radius 0) keep a horizontal one inside the Y axis
      const x0 = this.bands[0].x[0];
      series.push({type: 'scatter', symbolSize: 0, silent: true, tooltip: {show: false},
        data: this.refLines.filter((r) => !r.vertical).map((r) => [x0, r.value]),
        markLine: {silent: true, symbol: 'none', animation: false, data: this.refLines.map((r) => ({
          [r.vertical ? 'xAxis' : 'yAxis']: r.value, name: r.title,
          lineStyle: {color: r.color, type: r.dashed ? 'dashed' : 'solid', width: 1},
          // a vertical line's title above the plot, upright; a horizontal one's along it, at its start
          label: {formatter: r.title, position: r.vertical ? 'end' : 'insideStartTop', color: '#4d5261',
            ...(r.vertical ? {align: 'right'} : {})},
        }))}});
    }
    const axisName = {nameLocation: 'middle', nameTextStyle: {color: '#4d5261'}} as const;
    return {
      animation: false,
      // containLabel holds the tick labels, not the axis names: room for the Y name on the left, the X name
      // below and a vertical line's title above
      grid: {left: 44, right: 20, top: (this.bands.length > 1 ? 32 : 12) + (this.refLines.some((r) => r.vertical) ? 18 : 0),
        bottom: 36, containLabel: true},
      legend: {show: this.bands.length > 1, type: 'scroll', top: 4},
      tooltip: {trigger: 'axis', formatter: (params: any) => this.tooltip(params as any[])},
      xAxis: {type: 'value', scale: true, name: xTitle, nameGap: 26, ...axisName},
      // on a log axis the minor ticks carry labels too (2, 5 per decade): a decade alone reads as a linear axis
      yAxis: {type: log ? 'log' : 'value', scale: true, name: yTitle, nameGap: 52, ...axisName,
        ...(log ? {minorTick: {show: true, splitNumber: 2}, minorSplitLine: {show: true, lineStyle: {opacity: 0.4}}} : {})},
      series,
    };
  }

  /** X, then lower / line / upper of every band at that point. */
  private tooltip(params: any[]): string {
    if (!params.length)
      return '';
    const rows = new Set<string>();
    const u = this.yUnits ? ` ${this.yUnits}` : '';
    for (const p of params) {
      // each series' own index: bands of different categories need not share their X values
      const i = p.dataIndex;
      const b = this.bands.find((b) => b.name === p.seriesName);
      if (!b || i >= b.x.length)
        continue;
      const y = (v: number) => DG.format(v, this.yFormat);
      const range = `${y(b.lower[i])} – ${y(b.upper[i])}${u}`;
      const value = b.line ? `${y(b.line[i])}${u} (${range})` : range;
      rows.add(`${b.name === '' ? '' : `<b>${b.name}</b>: `}${value}`);
    }
    const xu = this.xUnits ? ` ${this.xUnits}` : '';
    return [`${this.xAxisTitle || this.xColumnName}: ${DG.format(params[0].data[0], this.xFormat)}${xu}`, ...rows]
      .join('<br>');
  }
}

const DEFAULT_FORMAT = '#0.###';

/** "Tumour volume (mm³)", "Tumour volume (mm³, log scale)": the units once, when given and not already in
 * the title, and the scale when it is not linear. */
export function axisTitle(title: string, units: string, log = false): string {
  const notes = [units && !title.includes(units) ? units : '', log ? 'log scale' : ''].filter((n) => n);
  return notes.length ? `${title} (${notes.join(', ')})` : title;
}

/** A constant formula line of the table: `${X} = 19` vertical, `${any other column} = 1500` horizontal. */
export interface RefLine {value: number, vertical: boolean, title: string, color: string, dashed: boolean}

/** The table's formula lines this chart can draw: the visible constant ones. Other formulas need the
 * platform's evaluator and are skipped; `DG.FormulaLinesHelper.getMeta` names the column but does not split
 * off the constant, hence the pattern. */
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
