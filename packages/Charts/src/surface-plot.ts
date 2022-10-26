import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {EChartViewer} from './echart-viewer';
import 'echarts-gl';


export class SurfacePlot extends EChartViewer {
  X: string;
  Y: string;
  Z: string;
  projection: string;
  bkgcolor: string | number;
  visualMapComponent: boolean;
  grid: boolean;
  axisLabel: boolean;
  XArr: [number[], number, number] = [[], 0, 0];
  YArr: [number[], number, number] = [[], 0, 0];
  ZArr: [number[], number, number] = [[], 0, 0];
  colsDict: {[name: string]: [number[], number, number]} = {};
  zip = (a: number[], b: number[], c: number[]) => {
    const arr = a.map((k, i) => [k, b[i], c[i]]);
    return arr.sort((a, b) => a[1] - b[1] || a[0] - b[0]);
  }


  constructor() {
    super();

    this.X = this.string('XColumnName', null);
    this.Y = this.string('YColumnName', null);
    this.Z = this.string('ZColumnName', null);
    this.projection = this.string('projection', 'perspective', {choices: ['perspective', 'orthographic']});
    this.bkgcolor = this.int('backgroundColor', 0xFFF);
    this.visualMapComponent = this.bool('legendVisualMapComponent', false);
    this.grid = this.bool('axisGrid', true);
    this.axisLabel = this.bool('axisLabel', true);

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
        type: 'category'
      },
      yAxis3D: {
        type: 'category'
      },
      zAxis3D: {
        type: 'category'
      },
      grid3D: {
        show: true,
        viewControl: {
          projection: 'perspective'
        },
        axisLabel: {
          show: true
        }
      },
      series: [
        {
          type: 'surface'
        }
      ]
    };
  }


  onTableAttached() {
    const num = Array.from(this.dataFrame.columns);
    if (num.length < 3) grok.shell.error(`Error: insufficient number of columns: expected 3+, got ${num.length}`);
    num.forEach(c => this.colsDict[c.name] = [c.toList(), c.min, c.max]);

    this.X = Object.keys(this.colsDict)[0];
    this.Y = Object.keys(this.colsDict)[1];
    this.Z = Object.keys(this.colsDict)[2];

    this.XArr = this.colsDict[this.X];
    this.YArr = this.colsDict[this.Y];
    this.ZArr = this.colsDict[this.Z];

    this.render();
  }


  onPropertyChanged(property: DG.Property) {
    const newVal = property.get(this);

    if (property.name.endsWith('ColumnName')) {
      const col = this.colsDict[newVal];
      switch (property.name) {
        case 'XColumnName':
          this.X = newVal;
          this.XArr = col;
          break;
        case 'YColumnName':
          this.Y = newVal;
          this.YArr = col;
          break;
        case 'ZColumnName':
          this.Z = newVal;
          this.ZArr = col;
          break;
      }
    } else
      switch (property.name) {
        case 'projection':
          this.projection = newVal;
          break;
        case 'backgroundColor':
          this.bkgcolor = DG.Color.toHtml(newVal);
          break;
        case 'legendVisualMapComponent':
          this.visualMapComponent = newVal;
          break;
        case 'axisGrid':
          this.grid = newVal;
          break;
        case 'axisLabel':
          this.axisLabel = newVal;
          break;
      }

    this.render();
  }


  render() {
    this.option.visualMap.min = this.ZArr[1];
    this.option.visualMap.max = this.ZArr[2];
    this.option.grid3D.viewControl.projection = this.projection;
    this.option.backgroundColor = this.bkgcolor;
    this.option.visualMap.show = this.visualMapComponent; 
    this.option.grid3D.show = this.grid;
    this.option.grid3D.axisLabel.show = this.axisLabel;
    this.option.xAxis3D.name = this.X;
    this.option.yAxis3D.name = this.Y;
    this.option.zAxis3D.name = this.Z;
    this.option.series[0].data = this.zip(this.XArr[0], this.YArr[0], this.ZArr[0]);
    this.chart.setOption(this.option);
  }
}