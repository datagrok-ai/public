import {EChartViewer} from './echart-viewer';
import 'echarts-gl';


export class SurfacePlot extends EChartViewer {
  XColumnName: string;
  YColumnName: string;
  ZColumnName: string;
  colsDict: {[name: string]: number[];} = {};
  colsNames: string[] = [];
  zip = (a: number[], b: number[], c: number[]) => a.map((k, i) => [k, b[i], c[i]]);

  constructor() {
    super();
    console.log('constructor');

    this.XColumnName = this.string('XColumnName', '');
    this.YColumnName = this.string('YColumnName', '');
    this.ZColumnName = this.string('ZColumnName', '');

    this.option = {
      tooltip: {},
      backgroundColor: '#fff',
      visualMap: {
        show: false,
        dimension: 2,
        inRange: {
          color: [
            '#313695',
            '#4575b4',
            '#74add1',
            '#abd9e9',
            '#e0f3f8',
            '#ffffbf',
            '#fee090',
            '#fdae61',
            '#f46d43',
            '#d73027',
            '#a50026'
          ]
        }
      },
      xAxis3D: {
        type: 'value'
      },
      yAxis3D: {
        type: 'value'
      },
      zAxis3D: {
        type: 'value'
      },
      grid3D: {
        viewControl: {
          // projection: 'orthographic'
        }
      },
      series: [
        {
          type: 'surface',
          wireframe: {
            // show: false
          }
        }
      ]
    };
  }

  render() {
    console.log('render');
    
    const cols = this.dataFrame.columns;
    const num = Array.from(cols.numerical);
    num.forEach(c => this.colsDict[c.name] = c.toList());
    this.colsNames = Object.keys(this.colsDict);

    const X_arr = this.colsDict[this.XColumnName];
    const Y_arr = this.colsDict[this.YColumnName];
    const Z_arr = this.colsDict[this.ZColumnName];

    this.option.visualMap.min = Math.min(...Z_arr);
    this.option.visualMap.max = Math.max(...Z_arr);
    this.option.series[0].data = this.zip(X_arr, Y_arr, Z_arr);
    this.chart.setOption(this.option);
  }
}