/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {HIGHLIGHT_WIDTH, LINE_MAX_WIDTH, LINE_MIN_WIDTH, MAXIMUM_COLUMN_NUMBER, MAXIMUM_ROW_NUMBER, MAXIMUM_SERIES_NUMBER, MOUSE_OVER_GROUP_COLOR, RadarIndicator} from './constants';
import {StringUtils} from '@datagrok-libraries/utils/src/string-utils';
import {EChartViewer} from '../echart/echart-viewer';
import {LegendHelper, VISIBILITY_MODE, VisibilityMode} from '../../utils/legend-utils';
import {ERROR_CLASS, MessageHandler} from '../../utils/utils';
import _ from 'lodash';

import '../../../css/radar-viewer.css';

type MinimalIndicator = '1' | '5' | '10' | '25';
type MaximumIndicator = '75' | '90' | '95' | '99';
type Normalization = 'Column' | 'Global';
const WARNING_CLASS = 'radar-warning';

// Based on this example: https://echarts.apache.org/examples/en/editor.html?c=radar
@grok.decorators.viewer({
  name: 'Radar',
  description: 'Creates a radar viewer',
  icon: 'icons/radar-viewer.svg',
})
export class RadarViewer extends EChartViewer {
  min: MinimalIndicator;
  max: MaximumIndicator;
  showCurrentRow: boolean;
  showMouseOverRow: boolean;
  showMouseOverRowGroup: boolean;
  showTooltip: boolean;
  showMin: boolean;
  showMax: boolean;
  showValues: boolean;
  normalization: Normalization;
  colorColumnName: string;
  backgroundMinColor: number;
  backgroundMaxColor: number;
  currentRowColor: number;
  mouseOverRowColor: number;
  lineColor: number;
  valuesColumnNames: string[];
  minValues: {[colName: string]: number};
  maxValues: {[colName: string]: number};
  columns: DG.Column[] = [];
  title: string;
  legendVisibility: VisibilityMode;
  legendHelper: LegendHelper = new LegendHelper();

  private static _canvas: HTMLCanvasElement | null = null;
  private static _ctx: CanvasRenderingContext2D | null = null;

  constructor() {
    super();
    this.title = this.string('title', 'Radar');
    this.min = <MinimalIndicator> this.string('min', '5', {choices: ['1', '5', '10', '25'],
      description: 'Minimum percentile value (indicated as dark blue area)'});
    this.max = <MaximumIndicator> this.string('max', '95', {choices: ['75', '90', '95', '99'],
      description: 'Maximum percentile value (indicated as light blue area)'});
    this.showCurrentRow = this.bool('showCurrentRow', true, {description: 'Highlights the current row', category: 'Selection'});
    this.showMouseOverRow = this.bool('showMouseOverRow', true, {category: 'Selection'});
    this.showMouseOverRowGroup = this.bool('showMouseOverRowGroup', true, {category: 'Selection'});
    this.showTooltip = this.bool('showTooltip', true);
    this.colorColumnName = this.string('colorColumnName', null, {category: 'Color'});
    this.backgroundMinColor = this.int('backgroundMinColor', 0xFFBB845D);
    this.backgroundMaxColor = this.int('backgroundMaxColor', 0xFFE7CDCD);
    this.currentRowColor = this.int('currentRowColor', 0xFF00FF00);
    this.mouseOverRowColor = this.int('mouseOverRowColor', 0xAAAAAA);
    this.lineColor = this.int('lineColor', 0xADD8E6);
    this.showMin = this.bool('showMin', false);
    this.showMax = this.bool('showMax', false);
    this.showValues = this.bool('showValues', false);
    this.normalization = <Normalization> this.string('normalization', 'Column', {choices: ['Column', 'Global'],
      category: 'Value', description: 'Column: scales each axis to its own range; Global: one scale across all columns'});
    this.valuesColumnNames = this.addProperty('valuesColumnNames', DG.TYPE.COLUMN_LIST, null,
      {columnTypeFilter: DG.COLUMN_TYPE_FILTER.NUMERICAL_NO_DATE_TIME, category: 'Value', groupWith: 'minValues, maxValues'});
    this.minValues = this.addProperty('minValues', DG.TYPE.MAP, null,
      {category: 'Value', description: 'Axis minimum per column; defaults to 0, or below the minimum for negative columns'}) ?? {};
    this.maxValues = this.addProperty('maxValues', DG.TYPE.MAP, null,
      {category: 'Value', description: 'Axis maximum per column; defaults to the column maximum'}) ?? {};
    this.legendVisibility = <VisibilityMode> this.string('legendVisibility', VISIBILITY_MODE.AUTO,
      {choices: Object.values(VISIBILITY_MODE)});
    this.legendHelper.onCategoriesChanged = () => this.render();

    this.option = {
      animation: false,
      silent: false,
      legend: {
        show: true,
      },
      color: ['#b0d7ff', '#4fbcf7', '#4287cc'],
      radar: {
        axisName: {
          backgroundColor: 'transparent',
          fontSize: '13px',
          fontFamily: 'Roboto',
          padding: [3, 5],
          color: '#4d5261',
        },
        radius: '60%',
        indicator: [],
      },
      tooltip: {
        show: false,
      },
      series: [{
        type: 'radar',
        data: [],
      }, {
        type: 'radar',
        data: [],
      }, {
        type: 'radar',
        data: [],
      }],
    };
  }

  private static getCanvasContext(): CanvasRenderingContext2D {
    if (!RadarViewer._canvas) {
      RadarViewer._canvas = ui.canvas();
      RadarViewer._ctx = RadarViewer._canvas.getContext('2d')!;
    }
    return RadarViewer._ctx!;
  }

  init() {
    this.columns = this.getColumns();
    this.option.radar.indicator = this.createRadarIndicators();

    this.updateMin();
    this.updateMax();

    const color = DG.Color.toHtml(this.showCurrentRow ? this.currentRowColor : this.lineColor);
    const currentRow = Math.max(this.dataFrame.currentRowIdx, 0);
    this.updateRow(color, currentRow);

    this.setupChartEvents();

    this.helpUrl = '/help/visualize/viewers/radar.md';
  }

  highlightRowIfEnabled(params: any) {
    if (!this.showMouseOverRow)
      return;


    const optionCopy = _.cloneDeep(this.option);
    const series = optionCopy.series[2].data.find((series: any) => series.name === params.name);
    if (series) {
      series.lineStyle.width = HIGHLIGHT_WIDTH;
      series.itemStyle.color = DG.Color.toHtml(this.mouseOverRowColor);
    }
    this.chart.setOption(optionCopy);
  }

  setupChartEvents() {
    this.chart.on('mouseover', (params: any) => {
      ui.tooltip.showRowGroup(this.dataFrame, (i) => {
        const currentRow = Math.max(this.dataFrame.currentRowIdx, 0);
        if (i === currentRow)
          return true;
        return false;
      }, params.event.event.x, params.event.event.y);

      const tooltipText = this.getTooltip(params);
      if (tooltipText)
        ui.tooltip.root.innerText = tooltipText;

      this.highlightRowIfEnabled(params);
    });

    this.chart.on('mouseout', () => {
      ui.tooltip.hide();
      this.chart.setOption(this.option);
    });

    this.chart.on('click', (params: any) => {
      const idx = parseInt(params.name.replace(/\D/g, ''), 10) - 1;
      if (!isNaN(idx) && idx >= 0) {
        this.dataFrame.currentRowIdx = idx;
        this.render();
      }
    });
  }

  getTooltip(params: any): string | null {
    const idx = parseInt(params.name.replace(/\D/g, ''), 10) - 1;
    if (params.componentType === 'series') {
      if (params.seriesIndex === 2) {
        const rows: string[] = [];
        for (let i = 0; i < this.columns.length; ++i) {
          const colName = this.columns[i].name;
          rows[i] = `${colName} : ${this.dataFrame.get(colName, idx)}`;
        }
        return rows.join('\n');
      }
    }
    return null;
  }

  onTableAttached() {
    this.init();
    this.root.appendChild(this.legendHelper.legendDiv);
    this.updateLegend();
    this.filter = this.dataFrame.filter;
    this.valuesColumnNames = Array.from(this.dataFrame.columns.numericalNoDateTime)
      .map((c: DG.Column) => c.name).slice(0, MAXIMUM_COLUMN_NUMBER);
    this.resubscribe(() => this.addSubs());
    this.render();
  }

  addSubs() {
    this.subs.push(this.dataFrame.onCurrentRowChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.onMouseOverRowChanged.subscribe((_) => {
      if (this.showMouseOverRow)
        this.render();
    }));
    this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.onColumnsRemoved.subscribe((data) => {
      const columnNamesToRemove = data.columns.map((column: DG.Column) => column.name);
      this.valuesColumnNames = this.valuesColumnNames.filter((columnName) => !columnNamesToRemove.includes(columnName));
      this.render();
    }));
    this.subs.push(this.dataFrame.onValuesChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.onMetadataChanged.subscribe((ev) => {
      if (ev?.args?.key?.includes('color-coding'))
        this.refreshLegendOnColorCodingChange();
      this.render();
    }));

    // Color edits made from the grid header fire this event rather than onMetadataChanged.
    this.subs.push(grok.events.onEvent('d4-grid-color-coding-changed').subscribe((_) => {
      this.refreshLegendOnColorCodingChange();
      this.render();
    }));
    this.subs.push(this.dataFrame.onMouseOverRowGroupChanged.subscribe((_) => {
      if (!this.showMouseOverRowGroup)
        return;
      const func = this.dataFrame.rows.mouseOverRowFunc;
      if (func) {
        const indexes = this.dataFrame.rows.where(func);
        this.render(Array.from(indexes));
      }
    }));
    this.subs.push(
      DG.debounce(ui.onSizeChanged(this.root), 50).subscribe((_) => {
        requestAnimationFrame(() => {
          this.chart?.resize();
          this.render();
        });
      }),
    );
  }

  public override onPropertyChanged(property: DG.Property) {
    if (property.name === 'table')
      this.updateTable();
    if (property.name === 'colorColumnName' || property.name === 'legendVisibility')
      this.updateLegend();
    this.render();
  }

  // The column colors derive from: the linked source column, or the column itself when not linked.
  get colorSourceColumn(): DG.Column | null {
    const column = this.dataFrame?.col(this.colorColumnName);
    if (!column)
      return null;
    const linkedName = column.getTag(DG.Tags.ColorCodingLinkedColumnName);
    return linkedName ? (this.dataFrame.col(linkedName) ?? column) : column;
  }

  updateLegend(): void {
    const legendColumn = this.colorSourceColumn;
    if (!legendColumn) {
      this.legendHelper.hide();
      return;
    }
    this.legendHelper.update(legendColumn);
    this.legendHelper.switchVisibility(this.legendVisibility, legendColumn);
  }

  refreshLegendOnColorCodingChange(): void {
    if (this.colorSourceColumn !== this.legendHelper.column)
      setTimeout(() => this.updateLegend(), 0);
  }

  getSeriesData(indexes?: number[]): void {
    this.clearData([0, 1, 2]);
    this.columns = this.getColumns();
    this.option.radar.indicator = this.createRadarIndicators();

    this.option.series[2].data = this.createSeriesData(indexes);

    if (this.filter.trueCount > MAXIMUM_ROW_NUMBER)
      MessageHandler._showMessage(this.root, 'Only first 1000 shown', WARNING_CLASS);

    if (this.showMin)
      this.updateMin();

    if (this.showMax)
      this.updateMax();

    this.updateCurrentRow();
    this.updateMouseOverRow();
    this.option.legend.show = false;
    this.option.silent = !this.showTooltip;
  }

  createSeriesData(filter?: number[]): any[] {
    const seriesData = [];
    const colorSourceColumn = this.colorSourceColumn;
    const selectedCategories = this.legendHelper.selectedCategories;

    for (let i = 0; i < this.filter.length && seriesData.length < MAXIMUM_ROW_NUMBER; i++) {
      if (!this.filter.get(i)) continue;

      if (selectedCategories && colorSourceColumn) {
        const category = colorSourceColumn.get(i);
        if (!selectedCategories.includes(category))
          continue;
      }

      const value = this.columns.map((c, colIdx) => this.cellValue(c, colIdx, i));

      const color = colorSourceColumn ? DG.Color.getRowColor(colorSourceColumn, i) : this.lineColor;

      seriesData.push({
        value: value,
        name: `row ${i + 1}`,
        symbol: 'none',
        lineStyle: {
          width: this.filter.trueCount > MAXIMUM_SERIES_NUMBER ? LINE_MIN_WIDTH : LINE_MAX_WIDTH,
          opacity: 0.8,
        },
        itemStyle: {
          color: filter && filter.includes(i) ? MOUSE_OVER_GROUP_COLOR : DG.Color.toHtml(color),
        },
        label: {
          show: this.showValues,
          formatter: (params: any) => StringUtils.formatNumber(params.value) as string,
        },
      });
    }

    return seriesData;
  }

  calculateRadarLabelWidths(n: number, padding: number = 25) {
    const {clientWidth: canvasWidth, clientHeight: canvasHeight} = this.root;
    const radiusPercent = parseInt(this.option.radar.radius) / 100;
    const centerX = canvasWidth / 2;
    const radius = radiusPercent * Math.min(canvasWidth, canvasHeight) / 2;

    const labels = Array.from({length: n}, (_, i) => {
      const theta = (2 * Math.PI * i) / n - Math.PI / 2;
      const x = centerX + radius * Math.cos(theta);
      let maxWidth: number;

      if (Math.abs(Math.cos(theta)) < 0.1)
        maxWidth = canvasWidth - 2 * padding;
      else if (Math.cos(theta) > 0)
        maxWidth = canvasWidth - x - padding;
      else
        maxWidth = x - padding;
      maxWidth = Math.max(0, maxWidth);
      return maxWidth;
    });

    return labels;
  }

  formatLabel(text: string, maxWidth: number): string {
    const ctx = RadarViewer.getCanvasContext();
    const {fontSize, fontFamily} = this.option.radar;
    ctx.font = `${fontSize} ${fontFamily}`;

    const ellipsis = '…';
    const ellipsisWidth = ctx.measureText(ellipsis).width;

    if (ctx.measureText(text).width <= maxWidth) return text;

    let left = 0;
    let right = text.length;
    let truncated = '';

    while (left < right) {
      const mid = Math.floor((left + right) / 2);
      const substr = text.slice(0, mid);
      const width = ctx.measureText(substr).width + ellipsisWidth;

      if (width <= maxWidth) {
        truncated = substr + ellipsis;
        left = mid + 1;
      } else
        right = mid;
    }

    return truncated || ellipsis;
  }

  updateCurrentRow(): void {
    const currentRowIdx = this.dataFrame.currentRowIdx;
    if (currentRowIdx < 0) return;
    const currentIn = this.filter.get(currentRowIdx);
    if (currentIn) {
      const color = DG.Color.toHtml(this.showCurrentRow ? this.currentRowColor : this.lineColor);
      this.updateRow(color, currentRowIdx);
    }
  }

  updateMouseOverRow(): void {
    const mouseOverIn = this.filter.get(this.dataFrame.mouseOverRowIdx);
    if (mouseOverIn && this.showMouseOverRow) {
      const color = DG.Color.toHtml(this.mouseOverRowColor);
      const currentRow = this.dataFrame.mouseOverRowIdx;
      if (currentRow !== -1)
        this.updateRow(color, currentRow);
    }
  }

  createRadarIndicators(): RadarIndicator[] {
    const indicators = this.columns.map((c) => this.createRadarIndicator(c));
    if (this.normalization === 'Global') {
      const mins = indicators.map((i) => i.min).filter((v) => isFinite(v));
      const maxs = indicators.map((i) => i.max).filter((v) => isFinite(v));
      if (mins.length > 0) {
        const min = Math.min(...mins);
        const max = Math.max(...maxs);
        for (const indicator of indicators) {
          indicator.min = min;
          indicator.max = max;
        }
      }
    }
    return indicators;
  }

  createRadarIndicator(c: DG.Column): Required<RadarIndicator> {
    return {
      name: c.name,
      min: this.minValues?.[c.name] ?? (c.min < 0 ? c.min + c.min * 0.1 : 0),
      max: this.maxValues?.[c.name] ?? c.max,
    };
  }

  private cellValue(c: DG.Column, colIdx: number, rowIdx: number): number | null {
    return c.isNone(rowIdx) ? null : this.clampToIndicator(colIdx, Number(c.get(rowIdx)));
  }

  private clampToIndicator(colIdx: number, value: number): number {
    const indicator = this.option.radar.indicator[colIdx];
    return indicator ? Math.min(Math.max(value, indicator.min), indicator.max) : value;
  }

  updateMin() {
    this.option.series[0].data[0] = {
      value: this.getQuantile(this.columns, this.getOptions(true).look.min / 100),
      name: `min percentile`,
      areaStyle: {
        color: DG.Color.toHtml(this.backgroundMinColor),
        opacity: 0.4,
      },
      lineStyle: {
        opacity: 0,
      },
      emphasis: {
        disabled: true,
      },
      symbolSize: 0,
    };
    this.option.color[0] = DG.Color.toHtml(this.backgroundMinColor);
  }

  updateMax() {
    this.option.series[1].data[0] = {
      value: this.getQuantile(this.columns, this.getOptions(true).look.max / 100),
      name: `max percentile`,
      areaStyle: {
        color: DG.Color.toHtml(this.backgroundMaxColor),
        opacity: 0.4,
      },
      lineStyle: {
        opacity: 0,
      },
      emphasis: {
        disabled: true,
      },
      symbolSize: 0,
    };
    this.option.color[1] = DG.Color.toHtml(this.backgroundMaxColor);
  }

  updateRow(color: string, currentRow: number) {
    this.option.series[2].data.push({
      value: this.columns.map((c, colIdx) => this.cellValue(c, colIdx, currentRow)),
      name: `row ${currentRow + 1}`,
      lineStyle: {
        width: 2,
        color: color,
      },
      symbolSize: 6,
      itemStyle: {
        color: color,
      },
      label: {
        show: this.showValues,
        formatter: function(params: any) {
          return StringUtils.formatNumber(params.value) as string;
        },
      },
    });
  }

  clearData(indexes: number[]) {
    for (let i = 0; i < indexes.length; ++i)
      this.option.series[indexes[i]].data = [];
  }

  getColumns() : DG.Column<any>[] {
    const columns: DG.Column<any>[] = [];
    const numericalColumns: DG.Column<any>[] = Array.from(this.dataFrame.columns.numericalNoDateTime);
    if (this.valuesColumnNames?.length > 0) {
      const selectedColumns = this.dataFrame.columns.byNames(this.valuesColumnNames);
      for (let i = 0; i < selectedColumns.length; ++i) {
        if (numericalColumns.includes(selectedColumns[i]))
          columns.push(selectedColumns[i]);
      }
    }
    return columns;
  }

  _testColumns(): boolean {
    const columnSet = new Set(this.dataFrame.columns.names());
    return this.valuesColumnNames.every((colName) => columnSet.has(colName));
  }

  render(indexes?: number[]) {
    if (!this.dataFrame)
      return;

    if (!this._testColumns() || this.valuesColumnNames.length === 0) {
      MessageHandler._showMessage(this.root, 'The Radar viewer requires a minimum of 1 numerical column.', ERROR_CLASS);
      return;
    }
    MessageHandler._removeMessage(this.root, WARNING_CLASS);
    MessageHandler._removeMessage(this.root, ERROR_CLASS);
    this.getSeriesData(indexes!);

    const radarLabelWidths = this.calculateRadarLabelWidths(this.columns.length);
    this.option.radar.axisName.formatter = (param: string) => {
      const idx = this.columns.findIndex((c) => c.name === param);
      if (idx === -1) return param;
      return this.formatLabel(param, radarLabelWidths[idx]);
    };

    this.chart.setOption(this.option, false, true);
  }

  detach() {
    for (const sub of this.subs)
      sub.unsubscribe();
    super.detach();
  }

  /* Going to be replaced with perc in stats */
  getQuantile(columns: DG.Column<any>[], percent: number): number[] {
    return columns.map((column, colIdx) => {
      const validSorted: number[] = [];
      for (let i = 0; i < column.length; i++) {
        if (!column.isNone(i))
          validSorted.push(Number(column.get(i)));
      }
      validSorted.sort((a, b) => a - b);
      const idx = Math.floor(percent * (validSorted.length - 1));
      return this.clampToIndicator(colIdx, validSorted[idx]);
    });
  }
}
