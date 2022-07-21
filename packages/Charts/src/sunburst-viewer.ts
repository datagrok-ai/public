import { EChartViewer } from './echart-viewer';
import { TreeUtils } from './utils/tree-utils';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
export class SunburstViewer extends EChartViewer {
  constructor() {
    super();
    this.initCommonProperties();

    this.option = {
      series: [
        {
          type: 'sunburst',
          label: {
            rotate: 'radial',
          },
        },
      ],
    };

    this.onPropertyChanged(null);
  }

  getSeriesData() {
    return TreeUtils.toForest(
      this.dataFrame,
      ['sex', 'race', 'dis_pop'],
      this.dataFrame.filter);
  }
}
