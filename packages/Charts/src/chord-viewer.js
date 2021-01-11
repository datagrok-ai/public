import * as echarts from 'echarts';
import { EChartViewer, Utils } from './echart-viewer';

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
            rotateLabel: true
          },
          label: {
            position: 'right',
            formatter: '{b}'
          },
          roam: true,
          focusNodeAdjacency: true,
          lineStyle: {
            color: 'source',
            curveness: 0.3
          }
        }
      ]};

    this.onPropertyChanged(null, false);
  }

  render() {
    let fromCol = this.dataFrame.getCol('source');
    let toCol = this.dataFrame.getCol('target');
    let nodes = [];

    let categories = Array.from(new Set(fromCol.categories.concat(toCol.categories)));
    let map = {};
    categories.forEach((cat, ind) => map[cat] = {id: ind, value: 0});
    let rowCount = this.dataFrame.rowCount;
    for (let i = 0; i < rowCount; i++) {
      map[fromCol.get(i)]['value']++;
      map[toCol.get(i)]['value']++;
    }

    let min = 1, max = rowCount * 2;
    let minSize = 5, maxSize = 150;
    let scale = (n) => maxSize*(n - min)/(max - min) + minSize;

    for (let name of categories) {
      nodes.push({
        name: name,
        value: map[name]['value'],
        symbolSize: scale(map[name]['value']),
        itemStyle: {
          color: DG.Color.toRgb(DG.Color.getCategoricalColor(map[name]['id']))
        },
        label: { show: map[name]['value'] > min, color: 'inherit', fontSize: 12 },
      });
    };

    this.option.series[0].data = nodes;
    this.option.series[0].links = Utils.mapRowsToObjects(this.dataFrame, ['source', 'target', 'value']);

    this.chart.setOption(this.option);
  }
}
