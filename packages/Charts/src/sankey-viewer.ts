import * as DG from 'datagrok-api/dg';

import { EChartViewer } from './echart-viewer';
import { TreeUtils } from './utils/tree-utils';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class SankeyViewer extends EChartViewer {
  sourceColumn: DG.Column;
  targetColumn: DG.Column;

  constructor() {
    super();

    this.sourceColumn = DG.Column.fromList('string', 'source', []);
    this.targetColumn = DG.Column.fromList('string', 'target', []);

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
    this.chart.on('click', { dataType: 'node' }, (params: any) => {
      this.dataFrame.selection.handleClick((i) => {
        return this.sourceColumn.get(i) === params.data.name ||
          this.targetColumn.get(i) === params.data.name;
      }, params.event.event);
    });

    this.chart.on('click', { dataType: 'edge' }, (params: any) => {
      this.dataFrame.selection.handleClick((i) => {
        return this.sourceColumn.get(i) === params.data.source &&
          this.targetColumn.get(i) === params.data.target;
      }, params.event.event);
    });

    this.dataFrame.onRowsFiltering.subscribe((_) => {
      this.refreshFilteredRows();
    });
  }

  onTableAttached() {
    this.sourceColumn = this.dataFrame.getCol('source');
    this.targetColumn = this.dataFrame.getCol('target');

    super.onTableAttached();
    this.initChartEventListeners();
  }

  refreshFilteredRows() {
    const dataframeSourceColumn = this.dataFrame.getCol('source');
    const dataframeTargetColumn = this.dataFrame.getCol('target');
    const filteredIndexList = this.dataFrame.filter.getSelectedIndexes();

    const sourceList: Array<string> = new Array<string>(filteredIndexList.length);
    const targetList: Array<string> = new Array<string>(filteredIndexList.length);

    for (let i = 0; i < filteredIndexList.length; i++) {
      sourceList[i] = dataframeSourceColumn.get(filteredIndexList[i]);
      targetList[i] = dataframeTargetColumn.get(filteredIndexList[i]);
    }

    this.sourceColumn = DG.Column.fromList('string', 'source', sourceList);
    this.targetColumn = DG.Column.fromList('string', 'target', targetList);
  }

  render() {
    const nodes = [];
    for (const name of new Set(this.sourceColumn.categories.concat(this.targetColumn.categories)))
      nodes.push({name: name});

    this.option.series[0].data = nodes;
    this.option.series[0].links = TreeUtils.mapRowsToObjects(this.dataFrame, ['source', 'target', 'value']);

    this.chart.setOption(this.option);
  }
}
