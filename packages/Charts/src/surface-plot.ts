import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {EChartViewer} from './echart-viewer';
import 'echarts-gl';


export class SurfacePlot extends EChartViewer {
  X: string;
  Y: string;
  Z: string;
  XArr: [number[], number, number] = [[], 0, 0];
  YArr: [number[], number, number] = [[], 0, 0];
  ZArr: [number[], number, number] = [[], 0, 0];
  colsDict: {[name: string]: [number[], number, number]} = {};
  //colsNames: string[] = [];
  zip = (a: number[], b: number[], c: number[]) => a.map((k, i) => [k, b[i], c[i]]);

  constructor() {
    super();

    this.X = this.string('XColumnName', null);
    this.Y = this.string('YColumnName', null);
    this.Z = this.string('ZColumnName', null);

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


  onTableAttached() {
    const cols = this.dataFrame.columns;
    const num = Array.from(cols.numerical);
    if (num.length < 3) grok.shell.error(`Error: insufficient number of numerical columns: expected 3+, got ${num.length}`);
    num.forEach(c => this.colsDict[c.name] = [c.toList(), c.min, c.max]);
    //this.colsNames = Object.keys(this.colsDict);

    this.XArr = this.colsDict[Object.keys(this.colsDict)[0]];
    this.YArr = this.colsDict[Object.keys(this.colsDict)[1]];
    this.ZArr = this.colsDict[Object.keys(this.colsDict)[2]];

    this.render();
  }


  onPropertyChanged(property: DG.Property) {
    switch (property.name) {
      case 'XColumnName':
        this.XArr = this.colsDict[property.get(this)];
        break;
      case 'YColumnName':
        this.YArr = this.colsDict[property.get(this)];
        break;
      case 'ZColumnName':
        this.ZArr = this.colsDict[property.get(this)];
        break;
    }
    this.render();
  }


  render() {
    this.option.visualMap.min = this.ZArr[1];
    this.option.visualMap.max = this.ZArr[2];
    this.option.series[0].data = this.zip(this.XArr[0], this.YArr[0], this.ZArr[0]);
    this.chart.setOption(this.option);
  }
}