import * as DG from 'datagrok-api/dg';
import { EChartViewer } from './echart-viewer';
import { TreeUtils } from './utils/tree-utils';


export class ChordViewer extends EChartViewer {
  constructor() {
    super();

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

  render() {
    const fromCol = this.dataFrame.getCol('source');
    const toCol = this.dataFrame.getCol('target');
    const nodes = [];

    const categories = Array.from(new Set(fromCol.categories.concat(toCol.categories)));
    const map: { [key: string]: any } = {};
    categories.forEach((cat, ind) => map[cat] = {id: ind, value: 0});
    const rowCount = this.dataFrame.rowCount;
    for (let i = 0; i < rowCount; i++) {
      map[fromCol.get(i)]['value']++;
      map[toCol.get(i)]['value']++;
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
