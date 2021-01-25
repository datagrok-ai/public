import * as echarts from 'echarts';
import { EChartViewer } from './echart-viewer';

export class TimelinesViewer extends EChartViewer {
  constructor() {
    super();
    
    this.initCommonProperties();
    this.subjectColumnName = this.string('subjectColumnName', 'USUBJID');
    this.startColumnName = this.string('startColumnName', 'AESTDY');
    this.endColumnName = this.string('endColumnName', 'AEENDY');
    this.colorByColumnName = this.string('colorByColumnName', 'SEX');

    this.data = [];

    this.option = {
      tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' }
      },
      grid: {
        left: '3%',
        right: '4%',
        bottom: '3%',
        containLabel: true
      },
      xAxis: {
        type: 'time',
        boundaryGap: ['5%', '5%'],
        axisLabel: { formatter: '{yyyy}-{MM}-{dd}' }
      },
      yAxis: {
        type: 'category',
        axisTick: { show: false },
        axisLine: { show: false },
      },
      series: [
        {
          type: 'custom',
          renderItem: (params, api) => {
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
              style: api.style()  // { fill: 'green' }
            };
          },
          encode: {
            x: [1, 2], 
            y: 0,
            tooltip: [3, 4]
          }
        },
      ]
    };

    this.onPropertyChanged(null, false);
  }

  getSeriesData() {
    this.data.length = 0;

    let columns = this.dataFrame.columns.byNames([
      this.subjectColumnName, this.startColumnName,
      this.endColumnName, this.colorByColumnName
    ]);

    for (let i = 0; i < this.dataFrame.rowCount; i++) {
      let row = [];
      for (let j = 0; j < columns.length; j++) {
        row.push((columns[j].type === 'datetime') ?
          new Date(`${columns[j].get(i)}`) : columns[j].get(i));
      }
      this.data.push(row);
    }
    return this.data;
  }
  
  render() {  
    this.option.series[0].data = this.getSeriesData(); 
    this.chart.setOption(this.option);
  }
}
