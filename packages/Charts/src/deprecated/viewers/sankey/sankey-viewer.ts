import { EChartViewer } from '../../../viewers/echart/echart-viewer';
import { TreeUtils } from '../../../utils/tree-utils';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class SankeyViewer extends EChartViewer {
  constructor() {
    super();

    this.initCommonProperties();

    this.option = {
      tooltip: {},
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

  initChartEventListeners() {
    const dataFrameSourceColumn = this.dataFrame.getCol('source');
    const dataFrameTargetColumn = this.dataFrame.getCol('target');

    this.chart.on('click', { dataType: 'node' }, (params: any) => {
      this.dataFrame.selection.handleClick((i) => {
        return dataFrameSourceColumn.get(i) === params.data.name ||
          dataFrameTargetColumn.get(i) === params.data.name;
      }, params.event.event);
    });

    this.chart.on('click', { dataType: 'edge' }, (params: any) => {
      this.dataFrame.selection.handleClick((i) => {
        return dataFrameSourceColumn.get(i) === params.data.source &&
          dataFrameTargetColumn.get(i) === params.data.target;
      }, params.event.event);
    });
  }

  onTableAttached() {
    super.onTableAttached();
    this.initChartEventListeners();
  }

  getNodes() {
    const nodes = [];

    const dataFrameSourceColumn = this.dataFrame.getCol('source');
    const dataFrameTargetColumn = this.dataFrame.getCol('target');
    const filteredIndexList = this.dataFrame.filter.getSelectedIndexes();

    const sourceList: Array<string> = new Array<string>(filteredIndexList.length);
    const targetList: Array<string> = new Array<string>(filteredIndexList.length);

    for (let i = 0; i < filteredIndexList.length; i++) {
      sourceList[i] = dataFrameSourceColumn.get(filteredIndexList[i]);
      targetList[i] = dataFrameTargetColumn.get(filteredIndexList[i]);
    }

    for (const name of new Set(sourceList.concat(targetList)))
      nodes.push({name: name});

    return nodes;
  }

  render() {
    this.option.series[0].data = this.getNodes();
    this.option.series[0].links = TreeUtils.mapRowsToObjects(this.dataFrame, ['source', 'target', 'value']);

    this.chart.setOption(this.option);
  }
}
