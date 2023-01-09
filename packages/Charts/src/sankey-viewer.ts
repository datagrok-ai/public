import { EChartViewer } from './echart-viewer';
import { TreeUtils } from './utils/tree-utils';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class SankeyViewer extends EChartViewer {
  constructor() {
    super();

    this.initCommonProperties();

    this.option = {
      tooltip: {
        trigger: 'item',
        triggerOn: 'mousemove',
      },
      series: [
        {
          type: 'sankey',
          emphasis: {
            focus: 'adjacency',
          },
          nodeAlign: 'left',
          lineStyle: {
            color: 'source',
            curveness: 0.5,
          },
        },
      ]};

    this.onPropertyChanged(null, false);
  }

  render() {
    const fromCol = this.dataFrame.getCol('source');
    const toCol = this.dataFrame.getCol('target');
    const nodes = [];
    for (const name of new Set(fromCol.categories.concat(toCol.categories)))
      nodes.push({name: name});

    this.option.series[0].data = nodes;
    this.option.series[0].links = TreeUtils.mapRowsToObjects(this.dataFrame, ['source', 'target', 'value']);

    this.chart.setOption(this.option);
  }
}
