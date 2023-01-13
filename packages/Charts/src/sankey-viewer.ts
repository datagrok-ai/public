import { EChartViewer } from './echart-viewer';
import { TreeUtils } from './utils/tree-utils';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class SankeyViewer extends EChartViewer {
  chartSourceColumnValues: string[] = [];
  chartTargetColumnValues: string[] = [];
  chartValueColumnValues: number[] = [];

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

    this.dataFrame.onRowsFiltered.subscribe((_) => {
      this.refreshColumnsOnFilter();
    });
  }

  onTableAttached() {
    this.chartSourceColumnValues = this.dataFrame.getCol('source').toList();
    this.chartTargetColumnValues = this.dataFrame.getCol('target').toList();
    this.chartValueColumnValues = this.dataFrame.getCol('value').toList();

    super.onTableAttached();
    this.initChartEventListeners();
  }

  refreshColumnsOnFilter() {
    const dataFrameSourceColumn = this.dataFrame.getCol('source');
    const dataFrameTargetColumn = this.dataFrame.getCol('target');
    const dataFrameValueColumn = this.dataFrame.getCol('value');
    const filteredIndexList = this.dataFrame.filter.getSelectedIndexes();

    const sourceList: Array<string> = new Array<string>(filteredIndexList.length);
    const targetList: Array<string> = new Array<string>(filteredIndexList.length);
    const valueList: Array<number> = new Array<number>(filteredIndexList.length);

    for (let i = 0; i < filteredIndexList.length; i++) {
      sourceList[i] = dataFrameSourceColumn.get(filteredIndexList[i]);
      targetList[i] = dataFrameTargetColumn.get(filteredIndexList[i]);
      valueList[i] = dataFrameValueColumn.get(filteredIndexList[i]);
    }

    this.chartSourceColumnValues = sourceList;
    this.chartTargetColumnValues = targetList;
    this.chartValueColumnValues = valueList;
  }

  render() {
    const nodes = [];
    for (const name of new Set(this.chartSourceColumnValues.concat(this.chartTargetColumnValues)))
      nodes.push({name: name});

    this.option.series[0].data = nodes;
    this.option.series[0].links = TreeUtils.mapRowsToObjects(
      [this.chartSourceColumnValues, this.chartTargetColumnValues, this.chartValueColumnValues],
      ['source', 'target', 'value']);

    this.chart.setOption(this.option);
  }
}
