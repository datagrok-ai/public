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

  constructor() {
    super();

    this.initCommonProperties();
    this.animation = this.bool('animation', false);
    this.layout = <layoutType>this.string('layout', 'orthogonal', { choices: ['orthogonal', 'radial'] });
    this.orient = <orientation>this.string('orient', 'LR', { choices: ['LR', 'RL', 'TB', 'BT'] });
    this.expandAndCollapse = this.bool('expandAndCollapse', false);
    this.animationDuration = this.int('animationDuration', 750);
    this.edgeShape = <edgeShape>this.string('edgeShape', 'curve', { choices: ['curve', 'polyline'] });
    this.symbol = <symbolType>this.string('symbol', 'emptyCircle', { choices: [
      'circle', 'emptyCircle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none',
    ] });
    this.symbolSize = this.int('symbolSize', 7);
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
    }, params.event!.event));
  }

  onPropertyChanged(p: DG.Property | null, render: boolean = true) {
    if (p?.name === 'hierarchyColumnNames')
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

  getSeriesData() {
    return [Utils.toTree(this.dataFrame, this.hierarchyColumnNames, this.dataFrame.filter)];
  }

  render() {
    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      return;

    super.render();
  }
}

type layoutType = 'orthogonal' | 'radial';
type orientation = 'LR' | 'RL' | 'TB' | 'BT';
type edgeShape = 'curve' | 'polyline';
type symbolType = 'circle' | 'emptyCircle' | 'rect' | 'roundRect' | 'triangle' | 'diamond' | 'pin' | 'arrow' | 'none';
