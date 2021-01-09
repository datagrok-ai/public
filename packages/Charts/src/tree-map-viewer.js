import * as echarts from 'echarts';
import { EChartViewer, Utils } from './echart-viewer';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class TreeMapViewer extends EChartViewer {

  constructor() {
    super();
    this.initCommonProperties();

    this.option = {
      tooltip: {
        formatter: function (info) {
          var value = info.value;
          var treePathInfo = info.treePathInfo;
          var treePath = [];

          for (var i = 1; i < treePathInfo.length; i++) {
            treePath.push(treePathInfo[i].name);
          }

          return '<div class="tooltip-title">' + treePath.join('/') + '</div>';
        }
      },
      series: [
        {
          type: 'treemap',
          levels: this.getLevelOption(),
        }
      ]
    };

    this.onPropertyChanged(null);
    this.chart.on('click', params => {
      console.log(params);
      let pattern = Utils.pathToPattern(this.getColumnNames(), params.data.path);
      this.dataFrame.rows.match(pattern).select();
    });
  }

  getColumnNames() {
    return ['sex', 'race', 'dis_pop'];
  }

  getSeriesData() {
    return Utils.toForest(
      this.dataFrame,
      this.getColumnNames(),
      this.dataFrame.filter);
  }

  getLevelOption() {
    return [
      {
        itemStyle: {
          borderWidth: 0,
          gapWidth: 5
        }
      },
      {
        itemStyle: {
          gapWidth: 1
        }
      },
      {
        colorSaturation: [0.35, 0.5],
        itemStyle: {
          gapWidth: 1,
          borderColorSaturation: 0.6
        }
      }
    ];
  }

}
