import * as echarts from 'echarts';
import { EChartViewer } from './echart-viewer';

export class TimelinesViewer extends EChartViewer {
  constructor() {
    super();
    
    this.initCommonProperties();
    this.subjectColumnName = this.string('subjectColumnName', 'USUBJID');
    this.startColumnName = this.string('startColumnName', 'AESTDY');
    this.endColumnName = this.string('endColumnName', 'AEENDY');
    this.colorByColumnName = this.string('colorByColumnName', 'EVENT');

    this.data = [];

    this.option = {
      tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' }
      },
      grid: {
        left: '3%',
        right: '4%',
        top: '4%',
        bottom: '3%',
        containLabel: true
      },
      xAxis: {
        type: 'value',
        min: value => value.min - 1,
        max: value => value.max + 1
      },
      yAxis: {
        type: 'category',
        axisTick: { show: false },
        axisLine: { show: false },
      },
      dataZoom: [
        { type: 'inside', xAxisIndex: [1, 2]},
        { type: 'inside', yAxisIndex: 0}
      ],
      series: [
        {
          type: 'custom',
          progressive: 0,   // Disable progressive rendering
          encode: {
            x: [1, 2], 
            y: 0,
            tooltip: 3
          }
        }
      ]
    };

    this.onPropertyChanged(null, false);
  }

  onTableAttached() {
    this.columns = this.dataFrame.columns.byNames([
      this.subjectColumnName, this.startColumnName,
      this.endColumnName, this.colorByColumnName
    ]).map((col, ind) => (col === null) ? this.dataFrame.columns.byIndex(ind) : col);

    this.colorMap = this.columns[3].categories.reduce((colorMap, c, i) => {
      colorMap[c] = DG.Color.toRgb(DG.Color.getCategoricalColor(i));
      return colorMap;
    }, {});

    this.option.series[0].renderItem = (params, api) => {
      const categoryIndex = api.value(0);
      const start = api.coord([api.value(1), categoryIndex]);
      const end = api.coord([api.value(2), categoryIndex]);
      const height = api.size([0, 1])[1] * 0.4;
  
      const rectShape = echarts.graphic.clipRectByRect({
        x: start[0],
        y: start[1] - height / 2,
        width: end[0] - start[0],
        height: height
      }, {
        x: params.coordSys.x,
        y: params.coordSys.y,
        width: params.coordSys.width,
        height: params.coordSys.height
      });
  
      return rectShape && {
        type: 'rect',
        transition: ['shape'],
        shape: rectShape,
        style: { fill: this.colorMap[api.value(3)] }
      };
    };

    super.onTableAttached();
  }

  getSeriesData() {
    this.data.length = 0;
    let max = this.columns[2].max;
    let min = this.columns[1].min;

    for (let i = 0; i < this.dataFrame.rowCount; i++) {
      let row = [];
      for (let j = 0; j < this.columns.length; j++) {
        row.push((this.columns[j].type === 'datetime') ?
          new Date(`${this.columns[j].get(i)}`) : (this.columns[j].isNone(i)) ?
          (j === 1) ? min : max : this.columns[j].get(i));
      }
      this.data.push(row);
    }

    this.chart.on('click', params => this.dataFrame.selection.handleClick(
      i => this.columns[0].get(i) === params.data[0], params.event));
    return this.data;
  }
  
  render() {  
    this.option.series[0].data = this.getSeriesData();
    this.chart.setOption(this.option);
  }
}
