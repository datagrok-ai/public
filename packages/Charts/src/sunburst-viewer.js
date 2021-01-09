import * as echarts from 'echarts';
import { EChartViewer, Utils } from './echart-viewer';

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
            rotate: 'radial'
          }
        }
      ]
    };

    this.onPropertyChanged(null);
  }

  getSeriesData() {
    return Utils.toForest(
      this.dataFrame,
      ['sex', 'race', 'dis_pop'],
      this.dataFrame.filter);
  }
}
