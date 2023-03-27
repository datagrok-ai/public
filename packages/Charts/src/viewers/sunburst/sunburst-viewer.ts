import * as DG from 'datagrok-api/dg';

import { EChartViewer } from '../echart/echart-viewer';
import { TreeUtils } from '../../utils/tree-utils';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic

/** Represents a sunburst viewer */
export class SunburstViewer extends EChartViewer {
  hierarchyColumnNames: string[];
  hierarchyLevel: number;

  constructor() {
    super();
    this.initCommonProperties();

    this.hierarchyColumnNames = this.addProperty('hierarchyColumnNames', DG.TYPE.COLUMN_LIST);
    this.hierarchyLevel = 3;

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

  onPropertyChanged(p: DG.Property | null, render: boolean = true) {
    if (p?.name === 'hierarchyColumnNames')
      this.render();
    else
      super.onPropertyChanged(p, render);
  }

  onTableAttached() {
    const categoricalColumns = [...this.dataFrame.columns.categorical].sort((col1, col2) =>
      col1.categories.length - col2.categories.length);

    if (categoricalColumns.length < 1)
      return;

    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      this.hierarchyColumnNames = categoricalColumns.slice(0, this.hierarchyLevel).map((col) => col.name);

    super.onTableAttached();
  }

  getSeriesData() {
    return TreeUtils.toForest(this.dataFrame, this.hierarchyColumnNames, this.dataFrame.filter);
  }

  render() {
    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      return;

    this.option.series[0].data = this.getSeriesData();

    this.chart.setOption(this.option);
  }
}
