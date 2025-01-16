import * as DG from 'datagrok-api/dg';

import { EChartViewer } from '../../../viewers/echart/echart-viewer';
import { TreeUtils } from '../../../utils/tree-utils';

export class ChordViewer extends EChartViewer {
  constructor() {
    super();

    this.top = this.string('top', '50px');
    this.bottom = this.string('bottom', '50px');

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
  }

  onTableAttached() {
    this.addSelectionOrDataSubs();
    this.initChartEventListeners();
    this.render();
  }

  getNodes() {
    const nodes = [];

    const dataFrameSourceColumn = this.dataFrame.getCol('source');
    const dataFrameTargetColumn = this.dataFrame.getCol('target');
    const filteredIndexList = this.filter.getSelectedIndexes();

    const sourceList: Array<string> = new Array<string>(filteredIndexList.length);
    const targetList: Array<string> = new Array<string>(filteredIndexList.length);

    for (let i = 0; i < filteredIndexList.length; i++) {
      sourceList[i] = dataFrameSourceColumn.get(filteredIndexList[i]);
      targetList[i] = dataFrameTargetColumn.get(filteredIndexList[i]);
    }

    const categories = Array.from(new Set(sourceList.concat(targetList)));
    const map: { [key: string]: any } = {};
    categories.forEach((cat, ind) => map[cat] = {id: ind, value: 0});
    const rowCount = filteredIndexList.length;
    for (let i = 0; i < rowCount; i++) {
      map[sourceList[i]]['value']++;
      map[targetList[i]]['value']++;
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

    return nodes;
  }

  render() {
    this.option.series[0].data = this.getNodes();
    this.option.series[0].links = TreeUtils.mapRowsToObjects(this.dataFrame, ['source', 'target', 'value']);

    this.chart.setOption(this.option);
  }
}
