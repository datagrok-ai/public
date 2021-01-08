import * as echarts from 'echarts';
import { EChartViewer, Utils } from './echart-viewer';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class TreeViewer extends EChartViewer {

  constructor() {
    super();
    
    this.layout = this.string('layout', 'orthogonal', { choices: ['orthogonal', 'radial'] });
    this.orient = this.string('orient', 'LR', { choices: ['LR', 'RL', 'TB', 'BT'] });
    this.expandAndCollapse = this.bool('expandAndCollapse', true);
    this.animationDuration = this.int('animationDuration', 750);
    this.edgeShape = this.string('curve', '', { choices: ['curve', 'polyline'] });
    this.symbol = this.string('symbol', 'emptyCircle', { choices: ['circle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none'] });

    this.top = this.string('top', '10%');

    this.option = {
      tooltip: {
        trigger: 'item',
        triggerOn: 'mousemove'
      },
      series: [
        {
          type: 'tree',

          top: '1%',
          left: '7%',
          bottom: '1%',
          right: '20%',

          symbolSize: 7,

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

          expandAndCollapse: true,
          animationDuration: 550,
          animationDurationUpdate: 750
        }
      ]
    };
  }

  onTableAttached() {
    this.subs.push(this.dataFrame.selection.onChanged.subscribe((_) => this.render()));
    this.subs.push(this.dataFrame.filter.onChanged.subscribe((_) => this.render()));

    this.render();
  }

  onPropertyChanged(p) {
    if (p !== null) {
      this.option.series[0][p.name] = p.get(this);
      this.myChart.setOption(this.option);
    }
  }

  render() {
    let data = Utils.toHierarchy(this.dataFrame, ['sex', 'race', 'dis_pop'], this.dataFrame.filter);
    this.option.series[0].data = [ data ];
    this.myChart.setOption(this.option);
  }
}
