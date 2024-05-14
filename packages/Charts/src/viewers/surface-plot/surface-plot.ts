import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {EChartViewer} from '../echart/echart-viewer';
import 'echarts-gl';

type AxisArray = {data: number[], min: number, max: number, type: string}

/** Represents a surface plot viewer */
@grok.decorators.viewer({
  name: 'Surface plot',
  description: 'Creates a surface plot viewer',
  icon: 'icons/surfaceplot-viewer.svg',
})
export class SurfacePlot extends EChartViewer {
  XColumnName: string;
  YColumnName: string;
  ZColumnName: string;
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
  plotFilter = () => {
    const ind = Array.from(this.filter.getSelectedIndexes());
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
    this.XColumnName = this.string('XColumnName', null);
    this.YColumnName = this.string('YColumnName', null);
    this.ZColumnName = this.string('ZColumnName', null);
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
        axisLabel: {
          formatter: formatter,
        },
      },
      yAxis3D: {
        type: 'value',
        min: 'dataMin',
        max: 'dataMax',
        axisLabel: {
          formatter: formatter,
        },
      },
      zAxis3D: {
        type: 'value',
        min: 'dataMin',
        max: 'dataMax',
        axisLabel: {
          formatter: formatter,
        },
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
          connectNulls: false,
        },
      ],
    };

    this.chart.on('click', (params: any) => {
      this.dataFrame.selection.handleClick((i) => {
        return this.XArr.data[i] === params.data[0] &&
          this.YArr.data[i] === params.data[1] &&
          this.ZArr.data[i] === params.data[2];
      }, params.event.event);
    });
  }

  onTableAttached() {
    if (!this.dataFrame.columns.length) {
      this._showErrorMessage('Table is empty');
      return;
    }
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe(() => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.onDataChanged, 50).subscribe(() => this.render()));
    this.colsDict = {};
    const num = Array.from(this.dataFrame.columns);
    num.forEach((c) => {
      let type: string;
      let isTime = false;
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
    let indx = [0, 1, 2];
    if (num.length === 1) {
      grok.shell.error('Error: insufficient number of columns: expected 3+, got 1');
      indx = [0, 0, 0];
    } else if (num.length === 2) {
      grok.shell.warning('Error: insufficient number of columns: expected 3+, got 2');
      indx = [0, 1, 1];
    }

    this.XColumnName = Object.keys(this.colsDict)[indx[0]];
    this.YColumnName = Object.keys(this.colsDict)[indx[1]];
    this.ZColumnName = Object.keys(this.colsDict)[indx[2]];
    this.XArr = this.colsDict[this.XColumnName];
    this.YArr = this.colsDict[this.YColumnName];
    this.ZArr = this.colsDict[this.ZColumnName];
    this.render(true, false);
  }

  onPropertyChanged(property: DG.Property) {
    const newVal = property.get(this);
    if (property.name.endsWith('ColumnName')) {
      const col = this.colsDict[newVal];
      switch (property.name) {
      case 'XColumnName':
        this.XColumnName = newVal;
        this.XArr = col;
        this.option.xAxis3D.axisLabel.formatter = col.type === 'time' ? undefined : formatter;
        break;
      case 'YColumnName':
        this.YColumnName = newVal;
        this.YArr = col;
        this.option.yAxis3D.axisLabel.formatter = col.type === 'time' ? undefined : formatter;
        break;
      case 'ZColumnName':
        this.ZColumnName = newVal;
        this.ZArr = col;
        this.option.zAxis3D.axisLabel.formatter = col.type === 'time' ? undefined : formatter;
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
      case 'table':
        this.updateTable();
        this.onTableAttached();
        return;
      }
      this.render();
    }
  }

  _testColumns() {
    return this.dataFrame.columns.toList().length >= 3;
  }

  _showErrorMessage(msg: string) {this.root.appendChild(ui.divText(msg, 'd4-viewer-error'));}

  render(computeData=false, filter=true) {
    // if (!this._testColumns()) {
    //   this._showErrorMessage('The Surface Plot viewer requires a minimum of 3 columns.');
    //   return;
    // }

    this.option.grid3D.viewControl.projection = this.projection;
    this.option.backgroundColor = this.bkgcolor;
    this.option.visualMap.show = this.visualMapComponent;
    this.option.grid3D.show = this.grid;
    this.option.grid3D.axisLabel.show = this.axisLabel;
    this.option.series[0].wireframe.show = this.wireframe;

    if (computeData) {
      this.option.visualMap.min = this.ZArr.min;
      this.option.visualMap.max = this.ZArr.max;
      this.option.xAxis3D.name = this.XColumnName;
      this.option.yAxis3D.name = this.YColumnName;
      this.option.zAxis3D.name = this.ZColumnName;
      this.option.xAxis3D.type = this.XArr.type;
      this.option.yAxis3D.type = this.YArr.type;
      this.option.zAxis3D.type = this.ZArr.type;
      const s = Math.round(Math.sqrt(this.XArr.data.length));
      this.option.series[0].dataShape = [s, s];
      this.rawData = this.zip(this.XArr.data, this.YArr.data, this.ZArr.data);
    }

    if (filter) this.option.series[0].data = this.plotFilter().sort(this.sort);
    else this.option.series[0].data = this.rawData.sort(this.sort);
    this.chart.resize();
    this.chart.setOption(this.option);
  }
}

const formatter = (val: any) => {
  if (typeof val === 'number' && !Number.isInteger(val))
    return val.toFixed(2);
  return val;
};
