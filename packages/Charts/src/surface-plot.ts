import * as DG from 'datagrok-api/dg';
//import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {EChartViewer} from './echart-viewer';
import 'echarts-gl';

type AxisArray = {data: number[], min: number, max: number, type: string}

export class SurfacePlot extends EChartViewer {
  X: string;
  Y: string;
  Z: string;
  projection: string;
  bkgcolor: string | number;
  visualMapComponent: boolean;
  grid: boolean;
  axisLabel: boolean;
  wireframe: boolean;
  XArr: AxisArray = {data: [], min: 0, max: 0, type: ''};
  YArr: AxisArray = {data: [], min: 0, max: 0, type: ''};
  ZArr: AxisArray = {data: [], min: 0, max: 0, type: ''};
  rawData: number[][] = [[]];
  colsDict: {[name: string]: {data: any[], min: number, max: number, type: string}} = {};

  zip = (a: any[], b: any[], c: any[]) => a.map((k, i) => [k, b[i], c[i]]);
  sort = (a: any[], b: any[]) => {
    const x = this.switch(a, b, 0, this.XArr.type);
    const y = this.switch(a, b, 1, this.YArr.type);
    return y || x;
  };
  filter = () => {
    const ind = Array.from(this.dataFrame.filter.getSelectedIndexes());
    return ind.map((id) => this.rawData[id]);
  };
  switch = (a: any[], b: any[], n: number, type: string) => {
    if (a[n] === undefined || b[n] === undefined) return 0;
    let res;
    switch (type) {
    case 'category':
      res = a[n].localeCompare(b[n]);
      break;
    case 'time':
      res = a[n].getTime() - b[n].getTime();
      break;
    default:
      res = a[n] - b[n];
      break;
    }
    return res;
  };

  constructor() {
    super();
    this.initCommonProperties();

    this.X = this.string('XColumnName', null);
    this.Y = this.string('YColumnName', null);
    this.Z = this.string('ZColumnName', null);
    this.projection = this.string('projection', 'perspective', {choices: ['perspective', 'orthographic']});
    this.bkgcolor = this.int('backgroundColor', 0xFFF);
    this.visualMapComponent = this.bool('legendVisualMapComponent', false);
    this.grid = this.bool('axisGrid', true);
    this.axisLabel = this.bool('axisLabel', true);
    this.wireframe = this.bool('wireframe', true);

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
            '#a50026',
          ],
        },
      },
      xAxis3D: {
        type: 'value',
        min: 'dataMin',
        max: 'dataMax',
      },
      yAxis3D: {
        type: 'value',
        min: 'dataMin',
        max: 'dataMax',
      },
      zAxis3D: {
        type: 'value',
        min: 'dataMin',
        max: 'dataMax',
      },
      grid3D: {
        show: true,
        viewControl: {
          projection: 'perspective',
        },
        axisLabel: {
          show: true,
        },
      },
      series: [
        {
          type: 'surface',
          wireframe: {
            show: true,
          },
        },
      ],
    };
  }

  onTableAttached() {
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render(false)));
    this.subs.push(DG.debounce(this.dataFrame.onDataChanged, 50).subscribe((_) => this.render()));

    const num = Array.from(this.dataFrame.columns);
    if (num.length < 3) grok.shell.error(`Error: insufficient number of columns: expected 3+, got ${num.length}`);
    num.forEach((c) => {
      let type: string; let isTime = false;
      switch (c.type) {
      case 'string':
        type = 'category';
        break;
      case 'datetime':
        type = 'time';
        isTime = true;
        break;
      default:
        type = 'value';
        break;
      }
      if (isTime)
        this.colsDict[c.name] = {data: c.toList().map((val) => new Date(val)), min: c.min, max: c.max, type: type};
      else this.colsDict[c.name] = {data: c.toList(), min: c.min, max: c.max, type: type};
    });

    this.X = Object.keys(this.colsDict)[0];
    this.Y = Object.keys(this.colsDict)[1];
    this.Z = Object.keys(this.colsDict)[2];

    this.XArr = this.colsDict[this.X];
    this.YArr = this.colsDict[this.Y];
    this.ZArr = this.colsDict[this.Z];

    this.render(true, false);
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
      this.render(true);
    } else {
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
      case 'wireframe':
        this.wireframe = newVal;
        break;
      }
      this.render();
    }
  }

  render(computeData=false, filter=true) {
    this.option.grid3D.viewControl.projection = this.projection;
    this.option.backgroundColor = this.bkgcolor;
    this.option.visualMap.show = this.visualMapComponent;
    this.option.grid3D.show = this.grid;
    this.option.grid3D.axisLabel.show = this.axisLabel;
    this.option.series[0].wireframe.show = this.wireframe;

    if (computeData) {
      this.option.visualMap.min = this.ZArr.min;
      this.option.visualMap.max = this.ZArr.max;
      this.option.xAxis3D.name = this.X;
      this.option.yAxis3D.name = this.Y;
      this.option.zAxis3D.name = this.Z;
      this.option.xAxis3D.type = this.XArr.type;
      this.option.yAxis3D.type = this.YArr.type;
      this.option.zAxis3D.type = this.ZArr.type;
      this.rawData = this.zip(this.XArr.data, this.YArr.data, this.ZArr.data);
    }

    if (filter) this.option.series[0].data = this.filter().sort(this.sort);
    else this.option.series[0].data = this.rawData.sort(this.sort);

    this.chart.setOption(this.option);
  }
}
