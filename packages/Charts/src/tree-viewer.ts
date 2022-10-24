import * as DG from 'datagrok-api/dg';
import { EChartViewer } from './echart-viewer';
import { TreeUtils, treeDataType } from './utils/tree-utils';
import * as utils from './utils/utils';

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
  initialTreeDepth: number;
  sizeColumnName: string;
  sizeAggrType: DG.AggregationType;
  symbolSizeRange: [number, number] = [5, 20];
  aggregations: string[] = Object.values(DG.AGG).filter((f) => f !== DG.AGG.KEY && f !== DG.AGG.PIVOT);
  colorColumnName: string;
  colorAggrType: DG.AggregationType;
  selectionColor = DG.Color.toRgb(DG.Color.selectedRows);

  constructor() {
    super();

    this.initCommonProperties();
    this.animation = this.bool('animation', true);
    this.layout = <layoutType>this.string('layout', 'orthogonal', { choices: ['orthogonal', 'radial'] });
    this.orient = <orientation>this.string('orient', 'LR', { choices: ['LR', 'RL', 'TB', 'BT'] });
    this.expandAndCollapse = this.bool('expandAndCollapse', true);
    this.initialTreeDepth = this.int('initialTreeDepth', 2, { min: -1, max: 5 });
    this.edgeShape = <edgeShape>this.string('edgeShape', 'curve', { choices: ['curve', 'polyline'] });
    this.symbol = <symbolType>this.string('symbol', 'emptyCircle', { choices: [
      'circle', 'emptyCircle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none',
    ] });
    this.symbolSize = this.int('symbolSize', 7);
    this.sizeColumnName = this.string('sizeColumnName');
    this.sizeAggrType = <DG.AggregationType>this.string('sizeAggrType', DG.AGG.AVG, { choices: this.aggregations });
    this.colorColumnName = this.string('colorColumnName');
    this.colorAggrType = <DG.AggregationType>this.string('colorAggrType', DG.AGG.AVG, { choices: this.aggregations });
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
    this.animationDuration = 750;

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
        p?.name === 'sizeAggrType' || p?.name === 'colorColumnName' || p?.name === 'colorAggrType')
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
    this.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/Charts/README.md#treeviewer';
  }

  colorCodeTree(data: treeDataType): void {
    const min = data['color-meta']['min'];
    const max = data['color-meta']['max'];
    const setItemStyle = (data: treeDataType) => {
      const nodeColor = DG.Color.toRgb(DG.Color.scaleColor(data.color, min, max));
      data.itemStyle = data.itemStyle ?? {};
      data.itemStyle.color = data.itemStyle.color === this.selectionColor ? this.selectionColor : nodeColor;
      if (data.children)
        data.children.forEach((child) => setItemStyle(child));
    }
    setItemStyle(data);
  }

  getSeriesData() {
    const aggregations = [];
    if (this.sizeColumnName)
      aggregations.push({ type: <DG.AggregationType>this.sizeAggrType, columnName: this.sizeColumnName, propertyName: 'size' });
    if (this.colorColumnName)
      aggregations.push({ type: <DG.AggregationType>this.colorAggrType, columnName: this.colorColumnName, propertyName: 'color' });

    return [TreeUtils.toTree(this.dataFrame, this.hierarchyColumnNames, this.dataFrame.filter, null, aggregations)];
  }

  render() {
    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      return;

    this.option.series[0].data = this.getSeriesData();

    this.option.series[0]['symbolSize'] = this.sizeColumnName ?
      (value: number, params: {[key: string]: any}) => utils.data.mapToRange(
      params.data.size, this.option.series[0].data[0]['size-meta']['min'],
      this.option.series[0].data[0]['size-meta']['max'], ...this.symbolSizeRange) : this.symbolSize;
    if (this.colorColumnName)
      this.colorCodeTree(this.option.series[0].data[0]);

    this.chart.setOption(this.option);
  }
}

type layoutType = 'orthogonal' | 'radial';
type orientation = 'LR' | 'RL' | 'TB' | 'BT';
type edgeShape = 'curve' | 'polyline';
type symbolType = 'circle' | 'emptyCircle' | 'rect' | 'roundRect' | 'triangle' | 'diamond' | 'pin' | 'arrow' | 'none';
