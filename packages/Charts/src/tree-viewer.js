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
    this.edgeShape = this.string('edgeShape', 'curve', { choices: ['curve', 'polyline'] });
    this.symbol = this.string('symbol', 'emptyCircle', { choices: [
      'circle', 'emptyCircle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none',
    ] });
    this.symbolSize = this.int('symbolSize', 7);

    this.legendColors = ['#beb0de', '#deb2b0', '#6495ed', '#b0c4de', '#ffc0cb'];

    this.option = {
      tooltip: {
        trigger: 'item',
        triggerOn: 'mousemove',
      },
      legend: {
        top: '2%',
        left: '2%',
        orient: 'vertical',
        data: [],
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

  setForestOpts() {
    const treeOptions = this.option.series[0];
    this.option.series.length = 0;

    for (const i of this.dataFrame.filter.getSelectedIndexes()) {
      const tree = this.parseNwk(this.newickCol.get(i));
      const name = `tree-${i}`;
      const color = this.legendColors[i % this.legendColors.length];

      this.option.series.push({
        ...treeOptions,
        name: name,
        data: [tree],

        orient: ['LR', 'TB', 'RL', 'BT'][i % 4],
        top: ['10%', '60%', '10%', '10%'][i % 4],
        left: ['60%', '0%', '10%', '0%'][i % 4],
        bottom: ['10%', '10%', '10%', '50%'][i % 4],
        right: ['10%', '10%', '60%', '10%'][i % 4],

        lineStyle: { color },
        emphasis: { focus: 'descendant' },
      });

      this.option.legend.data.push({
        name: name,
        icon: 'roundRect',
        itemStyle: { color },
      });
    }
  }

  parseNwk(s) {
    return { name: 'A', children: [
      { name: 'B', value: 123 }, { name: 'C', children: [{ name: 'D', value: 456 }, { name: 'E', value: 789 },
      ] }] };
  }

  updateSeriesData() {
    this.newickCol = this.dataFrame.columns.bySemType('newick');
    if (this.newickCol)
      this.setForestOpts();
    else
      this.option.series[0].data = [Utils.toTree(this.dataFrame, ['sex', 'race', 'dis_pop'], this.dataFrame.filter)];
  }

  onPropertyChanged(p, render = true) {
    const properties = p !== null ? [p] : this.props.getProperties();
    const isForest = this.option.series.length > 1;

    for (const p of properties) {
      for (const chart of this.option.series) {
        if (isForest && p.name === 'orient') continue;
        chart[p.name] = p.get(this);
      }
    }

    if (render)
      this.chart.setOption(this.option);
  }

  render() {
    this.updateSeriesData();
    this.chart.setOption(this.option);
  }
}
