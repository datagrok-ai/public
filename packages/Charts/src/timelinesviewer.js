import * as echarts from 'echarts';
import { EChartViewer } from './echart-viewer';

export class TimelinesViewer extends EChartViewer {
  constructor() {
    super();
    
    this.initCommonProperties();

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
              style: api.style()
            };
          },
          encode: {
            x: [1, 2], 
            y: 0,
            tooltip: [3, 4, 5]
          }
        },
      ]
    };

    this.onPropertyChanged(null, false);
  }

  getSeriesData() {
    return [
      ['Proj1', '2019-10-10', '2019-10-14', 320, 120, 220],
      ['Proj2', '2019-10-12', '2019-10-16', 302, 132, 182],
      ['Proj3', '2019-10-11', '2019-10-15', 301, 101, 191],
      ['Proj4', '2019-10-14', '2019-10-17', 334, 134, 234],
      ['Proj5', '2019-10-13', '2019-10-16', 390, 90, 290],
      ['Proj6', '2019-10-12', '2019-10-13', 330, 230, 330],
      ['Proj7', '2019-10-11', '2019-10-17', 320, 210, 310],
    ];
  }
  
  render() {  
    this.option.series[0].data = this.getSeriesData(); 
    this.chart.setOption(this.option);
  }
}
