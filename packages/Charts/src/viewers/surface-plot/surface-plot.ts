import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {EChartViewer} from '../echart/echart-viewer';
import 'echarts-gl';


type AxisArray = {data: number[], min: number, max: number, type: string}

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

  generateX = (start: number, stop: number, rows: number) => {
    const arr: number[] = [];
    const num = Math.ceil(Math.sqrt(rows));
    const step = (stop - start) / (num - 1);
    for (let i = 0; i < num; i++) arr.push(start + (step * i));
    return Array(...Array(num)).map(() => arr).flat();
  };

  generateY = (start: number, stop: number, rows: number) => {
    return this.generateX(start, stop, rows).sort((a, b) => a - b);
  };

  zip = (a: any[], b: any[], c: any[]) => a.map((k, i) => [k, b[i], c[i]]);

  filter = () => {
    const ind = Array.from(this.dataFrame.filter.getSelectedIndexes());
    return this.rawData.map((v, i) => ind.includes(i) ? v : [null, null, null]);
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
    this.axisLabel = this.bool('axisLabel', true);
    this.wireframe = this.bool('wireframe', true);
    this.grid = this.bool('axisGrid', true);

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
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe(() => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe(() => this.render(false)));
    this.subs.push(DG.debounce(this.dataFrame.onDataChanged, 50).subscribe(() => this.render()));

    const cols = Array.from(this.dataFrame.columns);
    const numCols = Array.from(this.dataFrame.columns.numerical);
    if (cols.length < 3)
      grok.shell.error(`Error: insufficient number of columns: expected 3+, got ${cols.length}`);
    if (numCols.length < 2)
      grok.shell.error(`Error: insufficient number of numerical columns: expected 2+, got ${numCols.length}`);

    for (const c of cols) {
      let type: string; let isTime = false;
      switch (c.type) {
      case ('string' || 'object' || 'bool' || 'dataframe'):
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
      else
        this.colsDict[c.name] = {data: c.toList(), min: c.min, max: c.max, type: type};
    };

    this.XColumnName = numCols[0].name;
    this.YColumnName = numCols[1].name;
    this.ZColumnName = numCols.length > 2 ? numCols[2].name : Array.from(this.dataFrame.columns.categorical)[0].name;

    this.XArr = {data: this.generateX(this.colsDict[this.XColumnName].min,
      this.colsDict[this.XColumnName].max, this.dataFrame.rowCount),
    min: this.colsDict[this.XColumnName].min,
    max: this.colsDict[this.XColumnName].max,
    type: 'value'};

    this.YArr = {data: this.generateY(this.colsDict[this.YColumnName].min,
      this.colsDict[this.YColumnName].max, this.dataFrame.rowCount),
    min: this.colsDict[this.YColumnName].min,
    max: this.colsDict[this.YColumnName].max,
    type: 'value'};

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
        if (col.type !== 'value') {
          grok.shell.warning(`${newVal} column is not numerical`);
          this.XArr = col;
        } else {
          this.XArr = {data: this.generateX(col.min, col.max, this.dataFrame.rowCount),
            min: col.min, max: col.max, type: 'value'};
        }
        break;
      case 'YColumnName':
        this.YColumnName = newVal;
        if (col.type !== 'value') {
          grok.shell.warning(`${newVal} column is not numerical`);
          this.YArr = col;
        } else {
          this.YArr = {data: this.generateY(col.min, col.max, this.dataFrame.rowCount),
            min: col.min, max: col.max, type: 'value'};
        }
        break;
      case 'ZColumnName':
        this.ZColumnName = newVal;
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
      this.option.xAxis3D.name = this.XColumnName;
      this.option.yAxis3D.name = this.YColumnName;
      this.option.zAxis3D.name = this.ZColumnName;
      this.option.xAxis3D.type = this.XArr.type;
      this.option.yAxis3D.type = this.YArr.type;
      this.option.zAxis3D.type = this.ZArr.type;
      this.rawData = this.zip(this.XArr.data, this.YArr.data, this.ZArr.data);
    }

    if (filter) this.option.series[0].data = this.filter();
    else this.option.series[0].data = this.rawData;

    this.chart.setOption(this.option);
  }
}
