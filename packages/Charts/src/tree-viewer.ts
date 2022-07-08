import * as DG from 'datagrok-api/dg';
import { EChartViewer, Utils } from './echart-viewer';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class TreeViewer extends EChartViewer {
  layout: layoutType;
  orient: orientation;
  expandAndCollapse: boolean;
  edgeShape: edgeShape;
  symbol: symbolType;
  symbolSize: number;
  hierarchyColumnNames: string[];
  animation: boolean;
  sizeColumnName: string;
  sizeAggrType: string;
  aggregations: string[] = Object.values(DG.AGG).filter((f) => f !== DG.AGG.KEY && f !== DG.AGG.PIVOT);

  constructor() {
    super();

    this.initCommonProperties();
    this.animation = this.bool('animation', true);
    this.layout = <layoutType>this.string('layout', 'orthogonal', { choices: ['orthogonal', 'radial'] });
    this.orient = <orientation>this.string('orient', 'LR', { choices: ['LR', 'RL', 'TB', 'BT'] });
    this.expandAndCollapse = this.bool('expandAndCollapse', true);
    this.animationDuration = this.int('animationDuration', 750);
    this.edgeShape = <edgeShape>this.string('edgeShape', 'curve', { choices: ['curve', 'polyline'] });
    this.symbol = <symbolType>this.string('symbol', 'emptyCircle', { choices: [
      'circle', 'emptyCircle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none',
    ] });
    this.symbolSize = this.int('symbolSize', 7);
    this.sizeColumnName = this.string('sizeColumnName');
    this.sizeAggrType = this.string('sizeAggrType', DG.AGG.AVG, { choices: this.aggregations });
    this.hierarchyColumnNames = this.addProperty('hierarchyColumnNames', DG.TYPE.COLUMN_LIST);

    this.option = {
      tooltip: {
        trigger: 'item',
        triggerOn: 'mousemove',
      },
      series: [
        {
          type: 'tree',

          label: {
            position: 'left',
            verticalAlign: 'middle',
            align: 'right',
            fontSize: 9,
          },

          leaves: {
            label: {
              position: 'right',
              verticalAlign: 'middle',
              align: 'left',
            },
          },
        },
      ],
    };

    this.onPropertyChanged(null);
  }

  initChartEventListeners() {
    this.chart.on('click', (params: {[key: string]: any}) => this.dataFrame.selection.handleClick((i) => {
      if (params.componentType !== 'series')
        return false;
      if (params.data.path === null)
        return true;
      else {
        const categories: string[] = params.data.path.split(' | ');
        let isMatch = true;
        categories.forEach((category, idx) => {
          const col = this.dataFrame.col(this.hierarchyColumnNames[idx]);
          isMatch = isMatch && category === (col!.type === DG.COLUMN_TYPE.BOOL ? col!.get(i).toString() : col!.get(i));
        });
        return isMatch;
      }
    }, params.event!.event, true));
  }

  onPropertyChanged(p: DG.Property | null, render: boolean = true) {
    if (p?.name === 'hierarchyColumnNames' || p?.name === 'sizeColumnName' ||
        p?.name === 'sizeAggrType')
      this.render();
    else
      super.onPropertyChanged(p, render);
  }

  onTableAttached() {
    const categoricalColumns = [...this.dataFrame.columns.categorical].sort((col1, col2) =>
      col1.categories.length - col2.categories.length);

    if (categoricalColumns.length < 1) {
      return;
    }

    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      this.hierarchyColumnNames = categoricalColumns.slice(0, 3).map((col) => col.name);

    super.onTableAttached();
    this.initChartEventListeners();
  }

  _mapToRange(x: number, min1: number, max1: number, min2: number, max2: number) {
    const range1 = max1 - min1;
    const range2 = max2 - min2;
    return (((x - min1) * range2) / range1) + min2;
  }

  getSeriesData() {
    const aggregations = [];
    if (this.sizeColumnName)
      aggregations.push({ type: <DG.AggregationType>this.sizeAggrType, columnName: this.sizeColumnName, propertyName: 'size' });

    return [Utils.toTree(this.dataFrame, this.hierarchyColumnNames, this.dataFrame.filter, null, aggregations)];
  }

  render() {
    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      return;

    this.option.series[0].data = this.getSeriesData();

    this.option.series[0]['symbolSize'] = this.sizeColumnName ?
      (value: number, params: {[key: string]: any}) => params.data.path ? this._mapToRange(
      params.data.size, this.option.series[0].data[0]['size-meta']['min'],
      this.option.series[0].data[0]['size-meta']['max'], 5, 20) : this.symbolSize : this.symbolSize;

    this.chart.setOption(this.option);
  }
}

type layoutType = 'orthogonal' | 'radial';
type orientation = 'LR' | 'RL' | 'TB' | 'BT';
type edgeShape = 'curve' | 'polyline';
type symbolType = 'circle' | 'emptyCircle' | 'rect' | 'roundRect' | 'triangle' | 'diamond' | 'pin' | 'arrow' | 'none';
