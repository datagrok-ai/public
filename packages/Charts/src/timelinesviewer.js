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
    this.markerSize = this.int('markerSize', 6);
    this.lineWidth = this.int('lineWidth', 3);

    this.data = [];

    this.option = {
      tooltip: {
        trigger: 'axis',
        axisPointer: { type: 'shadow' }
      },
      grid: {
        left: '3%',
        right: '4%',
        top: '2%',
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
      {
        type: 'inside',
        xAxisIndex: [1, 2],
        start: 0,
        end: 100
      },
      {
        start: 0,
        end: 100,
        height: 10,
        bottom: '1%'
      },
      { 
        type: 'inside',
        yAxisIndex: 0,
        start: 0,
        end: 100
      },
      { 
        type: 'slider',
        yAxisIndex: 0,
        start: 0,
        end: 100,
        width: 10,
      }
      ],
      series: [
        {
          type: 'custom',
          progressive: 0,   // Disable progressive rendering
          animation: false,
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
        y: start[1] - this.lineWidth / 2,
        width: end[0] - start[0],
        height: this.lineWidth
      }, {
        x: params.coordSys.x,
        y: params.coordSys.y,
        width: params.coordSys.width,
        height: params.coordSys.height
      });
  
      return {
        type: 'group',
        children: [{
          type: 'rect',
          transition: ['shape'],
          shape: rectShape,
          style: { fill: isNaN(api.value(3)) ? 'blue' : this.colorMap[api.value(3)] }
        },
        {
          type: 'circle',
          shape: {
            cx: start[0], cy: end[1], r: this.markerSize / 2
          },
          style: { fill: 'darkblue' }
        }, {
          type: 'circle',
          shape: {
            cx: end[0], cy: end[1], r: this.markerSize / 2
          },
          style: { fill: 'coral' }
        }
        ]
      };
    };

    super.onTableAttached();
  }

  getSeriesData() {
    this.data.length = 0;
    let tempObj = {};
    let getTime = (i, j) => this.columns[j].type === 'datetime' ?
      new Date(`${this.columns[j].get(i)}`) : this.columns[j].isNone(i) ?
      null : this.columns[j].get(i);

    for (let i = 0; i < this.dataFrame.rowCount; i++) {
      let id = this.columns[0].get(i);
      let start = getTime(i, 1);
      let end = getTime(i, 2);
      let event = this.columns[3].get(i);
      let key = `${id}-${start}-${end}`;
      if (tempObj.hasOwnProperty(key)) {
        tempObj[key][3].push(event);
      } else {
        tempObj[key] = [id, start, end, [event,]];
      }
    }
    this.data = Object.values(tempObj);

    this.chart.on('click', params => this.dataFrame.selection.handleClick(
      i => this.columns[0].get(i) === params.data[0], params.event));
    return this.data;
  }
  
  render() {  
    this.option.series[0].data = this.getSeriesData();
    this.chart.setOption(this.option);
  }
}
