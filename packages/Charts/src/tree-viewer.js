import * as echarts from 'echarts';
import { EChartViewer, Utils } from './echart-viewer';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class TreeViewer extends EChartViewer {

  constructor() {
    super();
    
    this.initCommonProperties();
    this.layout = this.string('layout', 'orthogonal', { choices: ['orthogonal', 'radial'] });
    this.orient = this.string('orient', 'LR', { choices: ['LR', 'RL', 'TB', 'BT'] });
    this.expandAndCollapse = this.bool('expandAndCollapse', true);
    this.animationDuration = this.int('animationDuration', 750);
    this.edgeShape = this.string('curve', '', { choices: ['curve', 'polyline'] });
    this.symbol = this.string('symbol', 'emptyCircle', { choices: ['circle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none'] });
    this.symbolSize = this.int('symbolSize', 7);

    this.option = {
      tooltip: {
        trigger: 'item',
        triggerOn: 'mousemove'
      },
      series: [
        {
          type: 'tree',

          label: {
            position: 'left',
            verticalAlign: 'middle',
            align: 'right',
            fontSize: 9
          },

          leaves: {
            label: {
              position: 'right',
              verticalAlign: 'middle',
              align: 'left'
            }
          },
        }
      ]
    };

    this.onPropertyChanged(null);
  }

  getSeriesData() {
    return [ Utils.toTree(this.dataFrame, ['sex', 'race', 'dis_pop'], this.dataFrame.filter) ];
  }
}
