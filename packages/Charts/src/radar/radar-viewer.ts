import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as echarts from 'echarts';
import {option} from '../radar/constants';

type MinimalIndicator = '1' | '5' | '10' | '25';
type MaximumIndicator = '75' | '90' | '95' | '99';

// Based on this example: https://echarts.apache.org/examples/en/editor.html?c=radar
export class RadarViewer extends DG.JsViewer {
  myChart: echarts.ECharts;
  min: MinimalIndicator;
  max: MaximumIndicator;
  showCurrentRow: boolean;
  showTooltip: boolean;
  showAllRows: boolean;
  backgroundMinColor: number;
  backgroundMaxColor: number;
  showMin: boolean;
  showMax: boolean;
  showValues: boolean;
  columnNames: string[] | null;

  constructor() {
    super();
    this.min = <MinimalIndicator>this.string('min', '5', { choices: ['1', '5', '10', '25'], description: 'Minimum percentile value (indicated as dark blue area)' });
    this.max = <MaximumIndicator>this.string('max', '95', { choices: ['75', '90', '95', '99'], description: 'Maximum percentile value (indicated as light blue area)' });
    this.showCurrentRow = this.bool('showCurrentRow', false);
    this.showTooltip = this.bool('showTooltip', true);
    this.showAllRows = this.bool('showAllRows', false);
    this.backgroundMinColor = this.int('backgroundMinColor', 0xFFB0D7FF);
    this.backgroundMaxColor = this.int('backgroundMaxColor', 0xFFBCE2F5);
    this.showMin = this.bool('showMin', true);
    this.showMax = this.bool('showMax', true);
    this.showValues = this.bool('showValues', true);
    this.columnNames = this.stringList('columnNames', null);
  
    const chartDiv = ui.div([], { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
    this.root.appendChild(chartDiv);
    this.myChart = echarts.init(chartDiv);
    this.subs.push(ui.onSizeChanged(chartDiv).subscribe((_) => this.myChart.resize()));
  }

  init() {
    option.radar.indicator = [];
    const columns = this.getColumns();
    for (const c of columns) {
      option.radar.indicator.push({name: c.name, max: c.max});
    }
    this.updateMin();
    this.updateMax();
    this.updateRow();
    //@ts-ignore
    this.myChart.on('mouseover', function(params) {
      if (params.componentType === 'series') {
        if (params.seriesIndex === 2) {
          let divs: HTMLElement[] = [];
          for (let i = 0; i < columns.length; ++i) {
            divs[i] = ui.divText(`${columns[i].name} : ${params.data.value[i]}`);
          } 
          ui.tooltip.show(ui.div(divs), params.event.event.x,  params.event.event.y);
        }
      }
    });
    this.myChart.on('mouseout', () => ui.tooltip.hide());
    this.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/master/help/visualize/viewers/radar-viewer.md';
  }

  onTableAttached() {
    this.init();
    this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.onCurrentRowChanged.subscribe((_) => {
      this.updateRow();
      this.myChart.setOption(option);
    }));
    this.subs.push(this.dataFrame.onColumnsRemoved.subscribe((_) => {
      this.init();
      this.myChart.setOption(option);
    }));
    this.render();
  }

  public override onPropertyChanged(property: DG.Property): void {
    const columns = this.getColumns();
    super.onPropertyChanged(property);
    switch (property.name) {
      case 'min':
        this.updateMin();
        break;
      case 'max':
        this.updateMax();
        break;
      case 'showMin':
        if (this.showMin != true) {
          this.clearData([0]);
        } else {
          this.clearData([1, 2]);
          this.init();
        }
        break;
      case 'showMax':
        if (this.showMax != true) {
          this.clearData([1]);
        } else {
          this.clearData([0, 2]);
          this.init();
        }
        break;
      case 'showCurrentRow':
        if (this.showCurrentRow === true) {
          this.clearData([0, 1]);
        } else {
          this.clearData([2]);
          this.init();
        }
        break;
      case 'showAllRows':
        if (this.showAllRows === true) {
          this.clearData([0, 1, 2]);
          let data = option.series[2].data; 
          for (let i = 0; i < this.dataFrame.rowCount; i++) {
            data.push({
              name: `row ${i}`,
              value: columns.map((c) => c.get(i)),
            });
          }
        } else {
          this.clearData([2]);
          this.init();
        }
        break;
      case 'showTooltip':
        if (this.showTooltip === false) {
          option['silent'] = true;
        } else {
          option['silent'] = false;
        }
        break;
      case 'backgroundMinColor':
        this.clearData([0, 1, 2]);
        this.init();
        break;
      case 'backgroundMaxColor':
        this.clearData([0, 1, 2]);
        this.init();
        break;
      case 'showValues':
        if (this.showValues === false) {
          option.series[2].data = [];
          option.series[2].data.push({
            value: columns.map((c) => c.get(this.dataFrame.currentRowIdx)),
            name: `row ${this.dataFrame.currentRowIdx + 1}`,
            lineStyle: {
              width: 2,
              type: 'dashed',
            },
            label: {
              show: false,
            },
            symbolSize: 6,
          });
        } else {
          this.clearData([0, 1, 2]);
          this.init();
        }
        break;
      case 'columnNames':
        this.init();
        break;
    }

    this.render();
  }

  updateMin() {
    const columns = this.getColumns();
    option.series[0].data[0] = {
      value: this.getQuantile(columns, this.getOptions(true).look.min / 100),
      name: `${this.getOptions(true).look.min}th percentile`,
      areaStyle: {
        color: DG.Color.toHtml(this.backgroundMinColor),
        opacity: 0.4
      },
      lineStyle: {
        opacity: 0
      },
      symbolSize: 0,
    }
  }

  updateMax() {
    const columns = this.getColumns();
    option.series[1].data[0] = {
      value: this.getQuantile(columns, this.getOptions(true).look.max / 100),
      name: `${this.getOptions(true).look.max}th percentile`,
      areaStyle: {
        color: DG.Color.toHtml(this.backgroundMaxColor),
        opacity: 0.4
      },
      lineStyle: {
        opacity: 0
      },
      symbolSize: 0,
    }
  }

  updateRow() {
    const columns = this.getColumns();
    option.series[2].data[0] = {
      value: columns.map((c) => c.get(this.dataFrame.currentRowIdx)),
      name: `row ${this.dataFrame.currentRowIdx + 1}`,
      lineStyle: {
        width: 2,
        type: 'dashed',
      },
      symbolSize: 6,
      label: {
        show: true,
        formatter: function (params: any) {
          return params.value as string;
        }
      }
    };
  }

  clearData(indexes: number[]) {
    for (let i = 0; i < indexes.length; ++i) {
      option.series[indexes[i]].data = [];
    }
  }

  getColumns() : DG.Column<any>[] {
    let columns: DG.Column<any>[] = [];
    if (this.columnNames === null) {
      columns = Array.from(this.dataFrame.columns.numerical);
      for (let i = 0; i < columns.length; ++i) {
        if (columns[i].type === 'datetime') {
          columns.splice(i, 1);
        }
      }
    } else {
      columns = this.dataFrame.columns.byNames(this.columnNames);
    }
    return columns;
  }

  render() {
    this.myChart.setOption(option);
  }

  /* Going to be replaced with perc in stats */
  getQuantile(columns: DG.Column<any>[], percent: number) {
    let result = [];
    for (const c of columns) {
      let idx = Math.floor(percent * c.length);
      let sortedIndexes = Array.from(c.getSortedOrder());
      result.push(c.get(sortedIndexes[idx]));
    }
    return result;
  }
}
