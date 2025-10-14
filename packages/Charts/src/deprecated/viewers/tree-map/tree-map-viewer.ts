import { EChartViewer } from '../../../viewers/echart/echart-viewer';
import { TreeUtils } from '../../../utils/tree-utils';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class TreeMapViewer extends EChartViewer {
  constructor() {
    super();
    this.initCommonProperties();

    this.option = {
      tooltip: {
        formatter: function(info: { [key: string]: any }) {
          const treePathInfo = info.treePathInfo;
          const treePath = [];

          for (let i = 1; i < treePathInfo.length; i++)
            treePath.push(treePathInfo[i].name);


          return '<div class="tooltip-title">' + treePath.join('/') + '</div>';
        },
      },
      series: [
        {
          type: 'treemap',
          levels: this.getLevelOption(),
        },
      ],
    };

    this.onPropertyChanged(null);
    this.chart.on('click', (params: { [key: string]: any }) => {
      console.log(params);
      const pattern = TreeUtils.pathToPattern(this.getColumnNames(), params.data.path);
      this.dataFrame.rows.match(pattern).select();
    });
  }

  getColumnNames() {
    return ['sex', 'race', 'dis_pop'];
  }

  getSeriesData() {
    return TreeUtils.toForest(
      this.dataFrame,
      this.getColumnNames(),
      this.filter,
    );
  }

  getLevelOption() {
    return [
      {
        itemStyle: {
          borderWidth: 0,
          gapWidth: 5,
        },
      },
      {
        itemStyle: {
          gapWidth: 1,
        },
      },
      {
        colorSaturation: [0.35, 0.5],
        itemStyle: {
          gapWidth: 1,
          borderColorSaturation: 0.6,
        },
      },
    ];
  }
}
