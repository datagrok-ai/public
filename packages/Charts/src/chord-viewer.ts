import * as DG from 'datagrok-api/dg';

import { EChartViewer } from './echart-viewer';
import { TreeUtils } from './utils/tree-utils';


export class ChordViewer extends EChartViewer {
  chartSourceColumn: DG.Column;
  chartTargetColumn: DG.Column;

  constructor() {
    super();

    this.chartSourceColumn = DG.Column.fromList('string', 'source', []);
    this.chartTargetColumn = DG.Column.fromList('string', 'target', []);

    this.top = this.string('top', '50px');
    this.left = this.string('left', '100px');
    this.bottom = this.string('bottom', '50px');
    this.right = this.string('right', '100px');

    this.animationDuration = this.int('animationDuration', 500);
    this.animationDurationUpdate = this.int('animationDurationUpdate', 750);

    this.option = {
      tooltip: {},
      series: [
        {
          type: 'graph',
          layout: 'circular',
          circular: {
            rotateLabel: true,
          },
          label: {
            position: 'right',
            formatter: '{b}',
          },
          roam: true,
          focusNodeAdjacency: true,
          lineStyle: {
            color: 'source',
            curveness: 0.3,
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

    this.dataFrame.onRowsFiltering.subscribe((_) => {
      this.refreshColumnsOnFilter();
    });
  }

  onTableAttached() {
    this.chartSourceColumn = this.dataFrame.getCol('source');
    this.chartTargetColumn = this.dataFrame.getCol('target');

    super.onTableAttached();
    this.initChartEventListeners();
  }

  refreshColumnsOnFilter() {
    const dataframeSourceColumn = this.dataFrame.getCol('source');
    const dataframeTargetColumn = this.dataFrame.getCol('target');
    const filteredIndexList = this.dataFrame.filter.getSelectedIndexes();

    const sourceList: Array<string> = new Array<string>(filteredIndexList.length);
    const targetList: Array<string> = new Array<string>(filteredIndexList.length);

    for (let i = 0; i < filteredIndexList.length; i++) {
      sourceList[i] = dataframeSourceColumn.get(filteredIndexList[i]);
      targetList[i] = dataframeTargetColumn.get(filteredIndexList[i]);
    }

    this.chartSourceColumn = DG.Column.fromList('string', 'source', sourceList);
    this.chartTargetColumn = DG.Column.fromList('string', 'target', targetList);
  }

  render() {
    const nodes = [];

    const categories = Array.from(new Set(this.chartSourceColumn.categories.concat(this.chartTargetColumn.categories)));
    const map: { [key: string]: any } = {};
    categories.forEach((cat, ind) => map[cat] = {id: ind, value: 0});
    const rowCount = this.chartSourceColumn.length;
    for (let i = 0; i < rowCount; i++) {
      map[this.chartSourceColumn.get(i)]['value']++;
      map[this.chartTargetColumn.get(i)]['value']++;
    }

    const min = 1; const max = rowCount * 2;
    const minSize = 5; const maxSize = 150;
    const scale = (n: number) => maxSize*(n - min)/(max - min) + minSize;

    for (const name of categories) {
      nodes.push({
        name: name,
        value: map[name]['value'],
        symbolSize: scale(map[name]['value']),
        itemStyle: {
          color: DG.Color.toRgb(DG.Color.getCategoricalColor(map[name]['id'])),
        },
        label: { show: map[name]['value'] > min, color: 'inherit', fontSize: 12 },
      });
    };

    this.option.series[0].data = nodes;
    this.option.series[0].links = TreeUtils.mapRowsToObjects(this.dataFrame, ['source', 'target', 'value']);

    this.chart.setOption(this.option);
  }
}
