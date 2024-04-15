import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import * as echarts from 'echarts';
import {MAXIMUM_SERIES_NUMBER, option} from './constants';
import {StringUtils} from '@datagrok-libraries/utils/src/string-utils';
import { EChartViewer } from '../echart/echart-viewer';

import '../../../css/radar-viewer.css';

type MinimalIndicator = '1' | '5' | '10' | '25';
type MaximumIndicator = '75' | '90' | '95' | '99';
const ERROR_CLASS = 'd4-viewer-error';
const WARNING_CLASS = 'radar-warning';

// Based on this example: https://echarts.apache.org/examples/en/editor.html?c=radar
@grok.decorators.viewer({
  name: 'Radar',
  description: 'Creates a radar viewer',
  icon: 'icons/radar-viewer.svg',
})
export class RadarViewer extends EChartViewer {
  get type(): string {return 'RadarViewer';}

  chart: echarts.ECharts;
  min: MinimalIndicator;
  max: MaximumIndicator;
  showOnlyCurrentRow: boolean;
  showTooltip: boolean;
  showMin: boolean;
  showMax: boolean;
  showValues: boolean;
  backgroundMinColor: number;
  backgroundMaxColor: number;
  currentRowColor: number;
  valuesColumnNames: string[];
  columns: DG.Column[] = [];
  title: string;
  resizeScheduled: boolean = false;

  constructor() {
    super();
    this.rowSource = 'CurrentRow';
    this.title = this.string('title', 'Radar');
    this.min = <MinimalIndicator> this.string('min', '5', { choices: ['1', '5', '10', '25'],
      description: 'Minimum percentile value (indicated as dark blue area)' });
    this.max = <MaximumIndicator> this.string('max', '95', { choices: ['75', '90', '95', '99'],
      description: 'Maximum percentile value (indicated as light blue area)' });
    this.showOnlyCurrentRow = this.bool('showOnlyCurrentRow', true, {description: 'Hides max and min values'});
    this.showTooltip = this.bool('showTooltip', true);
    this.backgroundMinColor = this.int('backgroundMinColor', 0xFFB0D7FF);
    this.backgroundMaxColor = this.int('backgroundMaxColor', 0xFFBCE2F5);
    this.currentRowColor = this.int('currentRowColor', 0xFFBAFFBA);
    this.showMin = this.bool('showMin', false);
    this.showMax = this.bool('showMax', false);
    this.showValues = this.bool('showValues', true);
    this.valuesColumnNames = this.addProperty('valuesColumnNames', DG.TYPE.COLUMN_LIST, null,
      {columnTypeFilter: DG.TYPE.NUMERICAL});
    
    const chartDiv = ui.div([], { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
    this.root.appendChild(chartDiv);
    this.chart = echarts.init(chartDiv);
    this.subs.push(ui.onSizeChanged(chartDiv).subscribe((_) => {
      if (!this.resizeScheduled) {
        this.resizeScheduled = true;
        requestAnimationFrame(() => {
          this.chart.resize();
          this.resizeScheduled = false;
        });
      }
    }));
  }

  init() {
    option.radar.indicator = [];
    const columnNames: string[] = [];
    for (const column of this.dataFrame.columns.numerical)
      columnNames.push(column.name);

    this.columns = this.getColumns();
    for (const c of this.columns) {
      const minimalVal = c.min < 0 ? (c.min + c.min * 0.1) : 0;
      if (c.type === 'datetime')
        option.radar.indicator.push({name: c.name, max: this.getYearFromDate(c.max), min: this.getYearFromDate(c.min)});
      else
        option.radar.indicator.push({name: c.name, max: c.max, min: minimalVal});
    }
    this.updateMin();
    this.updateMax();
    this.updateRow();

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
    });
    this.chart.on('mouseout', () => ui.tooltip.hide());
    this.helpUrl = 'https://datagrok.ai/help/visualize/viewers/radar';
  }

  getTooltip(params: any): string | null {
    if (params.componentType === 'series') {
      if (params.seriesIndex === 2) {
        const rows: string[] = [];
        for (let i = 0; i < this.columns.length; ++i)
          rows[i] = `${this.columns[i].name} : ${params.data.value[i]}`;
        return rows.join('\n');
      }
    }
    return null;
  }

  onTableAttached() {
    this.init();
    this.valuesColumnNames = Array.from(this.dataFrame.columns.numerical)
      .map((c: DG.Column) => c.name);
    this.subs.push(this.dataFrame.onCurrentRowChanged.subscribe((_) => this.render()));
    this.render();
  }

  public override onPropertyChanged(property: DG.Property) {
    this.render();
  }

  getSeriesData(): void {
    option.radar.indicator = [];
    this.clearData([0, 1, 2]);
    this.columns = this.getColumns();
    for (const c of this.columns) {
      const minimalVal = c.min < 0 ? (c.min + c.min * 0.1) : 0;
      if (c.type === 'datetime')
        option.radar.indicator.push({ name: c.name, max: this.getYearFromDate(c.max), min: this.getYearFromDate(c.min) });
      else
        option.radar.indicator.push({ name: c.name, max: c.max, min: minimalVal });
    }
  
    if (this.showMin)
      this.updateMin();
  
    if (this.showMax)
      this.updateMax();

    if (this.showOnlyCurrentRow)
      this.updateRow();

    const currentRowIndex = this.dataFrame.currentRowIdx ?? 0;
    const currentIn = this.filter.get(currentRowIndex);
    
    if (!currentIn || !this.showOnlyCurrentRow) {
      const filter = this.formulaFilter ? this.filter.clone().invert() : this.filter;
      const seriesData = [];
      let colorIndex = 0;

      for (let i = 0; i < filter.length && seriesData.length < MAXIMUM_SERIES_NUMBER; i++) {
        if (!filter.get(i)) continue;
        if (i === currentRowIndex && currentIn) continue;
        
        const value = this.columns.map((c) => {
          if (c.type === 'datetime')
            return this.getDate(c, c.getRawData()[i]);
          const value = Number(c.get(i));
          return value !== -2147483648 ? value : 0;
        });

        const colors = DG.Color.categoricalPalette;
        const color = colors[colorIndex % colors.length];
        ++colorIndex;

        seriesData.push({
          value: value,
          name: `row ${i + 1}`,
          itemStyle: { color: DG.Color.toHtml(color) },
          label: { show: this.showValues, formatter: (params: any) => StringUtils.formatNumber(params.value) as string },
        });
      }
      
      option.series[2].data = seriesData;
      if (currentIn && !this.showOnlyCurrentRow)
        this.updateRow();
    }

    option.legend.show = !(this.filter.length > 1);
    option.silent = !this.showTooltip;
  }

  getDate(c: DG.Column, value: number) {
    const date = this.getYearFromDate(value);
    const isRight = date >= this.getYearFromDate(c.min) && date <= this.getYearFromDate(c.max);
    return isRight ? date : this.getYearFromDate(c.min);
  }

  getYearFromDate(value: number) {
    return new Date(Math.floor(value) / 1000).getFullYear();
  }

  updateMin() {
    option.series[0].data[0] = {
      value: this.getQuantile(this.columns, this.getOptions(true).look.min / 100),
      name: `min percentile`,
      areaStyle: {
        color: DG.Color.toHtml(this.backgroundMinColor),
        opacity: 0.4,
      },
      lineStyle: {
        opacity: 0,
      },
      symbolSize: 0,
    };
    option.color[0] = DG.Color.toHtml(this.backgroundMinColor);
  }

  updateMax() {
    option.series[1].data[0] = {
      value: this.getQuantile(this.columns, this.getOptions(true).look.max / 100),
      name: `max percentile`,
      areaStyle: {
        color: DG.Color.toHtml(this.backgroundMaxColor),
        opacity: 0.4,
      },
      lineStyle: {
        opacity: 0,
      },
      symbolSize: 0,
    };
    option.color[1] = DG.Color.toHtml(this.backgroundMaxColor);
  }

  updateRow() {
    const currentRow = Math.max(this.dataFrame.currentRowIdx, 0);
    option.series[2].data.unshift({
      value: this.columns.map((c) => {
        if (c.type === 'datetime')
          return this.getDate(c, c.getRawData()[currentRow]);
        const value = Number(c.get(currentRow));
        return value != -2147483648 ? value : 0;
      }),
      name: `row ${currentRow + 1}`,
      lineStyle: {
        width: 2,
        type: 'dashed',
        color: DG.Color.toHtml(this.currentRowColor),
      },
      symbolSize: 6,
      itemStyle: {
        color: DG.Color.toHtml(this.currentRowColor),
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
      option.series[indexes[i]].data = [];
  }

  getColumns() : DG.Column<any>[] {
    let columns: DG.Column<any>[] = [];
    const numericalColumns: DG.Column<any>[] = Array.from(this.dataFrame.columns.numerical);
    if (this.valuesColumnNames?.length > 0) {
      const selectedColumns = this.dataFrame.columns.byNames(this.valuesColumnNames);
      for (let i = 0; i < selectedColumns.length; ++i) {
        if (numericalColumns.includes(selectedColumns[i]))
          columns.push(selectedColumns[i]);
      }
    }
    else
      columns = numericalColumns.slice(0, 20);
    return columns;
  }

  _showMessage(msg: string, className: string) {this.root.appendChild(ui.divText(msg, className));}

  _removeMessage(className: string) {
    const divTextElement = this.root.getElementsByClassName(className)[0];
    if (divTextElement)
      this.root.removeChild(divTextElement);
  }

  render() {
    if (this.valuesColumnNames.length < 1) {
      this._showMessage('The Radar viewer requires a minimum of 1 numerical column.', ERROR_CLASS);
      return;
    }

    this._removeMessage(WARNING_CLASS);
    this._removeMessage(ERROR_CLASS);
    this.getSeriesData();
    this.chart.setOption(option);

    if ((!this.showOnlyCurrentRow && this.filter.trueCount > MAXIMUM_SERIES_NUMBER) ||
      (this.showOnlyCurrentRow && option.series[2].data.length > 1))
      this._showMessage('Only first 100 rows are shown', WARNING_CLASS);
  }

  /* Going to be replaced with perc in stats */
  getQuantile(columns: DG.Column<any>[], percent: number) {
    const result = [];
    for (const c of columns) {
      const datetime = c.getRawData().map((value: number) => this.getDate(c, value));
      const values = c.type === 'datetime' ? datetime : Array.from(c.values());
      const isValidValue = (value: number) => typeof value === 'bigint' 
        ? value !== BigInt('-2147483648')
        : value !== -2147483648;
      const sortedValues = values.filter(isValidValue).sort((a, b) => Number(a) - Number(b));

      const idx = Math.floor(percent * (sortedValues.length - 1));
      const value = sortedValues[idx];
      result.push(Number(value));
    }
    return result;
  }
}
