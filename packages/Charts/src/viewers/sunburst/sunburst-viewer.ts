import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {EChartViewer} from '../echart/echart-viewer';
import {TreeUtils, treeDataType} from '../../utils/tree-utils';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic

/** Represents a sunburst viewer */
@grok.decorators.viewer({
  name: 'Sunburst',
  description: 'Creates a sunburst viewer',
  icon: 'icons/sunburst-viewer.svg',
})
export class SunburstViewer extends EChartViewer {
  hierarchyColumnNames: string[];
  hierarchyLevel: number;

  constructor() {
    super();
    this.initCommonProperties();
    this.initEventListeners();

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

  initEventListeners(): void {
    this.chart.on('click', (params: any) => {
      if (params.event.event.ctrlKey) {
        const path: string[] = params.data.path.split('|').map((str: string) => str.trim());
        this.dataFrame.selection.handleClick((i) => {
          if (!this.dataFrame.filter.get(i))
            return false;
          for (let j = 0; j < path.length; j++) {
            if (this.dataFrame.getCol(this.hierarchyColumnNames[j]).get(i).toString() !== path[j])
              return false;
          }
          return true;
        }, params.event.event);
      }
    });
  }

  onPropertyChanged(p: DG.Property | null, render: boolean = true): void {
    if (p?.name === 'hierarchyColumnNames')
      this.render();
    else
      super.onPropertyChanged(p, render);
  }

  onTableAttached(): void {
    const categoricalColumns = [...this.dataFrame.columns.categorical].sort((col1, col2) =>
      col1.categories.length - col2.categories.length);

    if (categoricalColumns.length < 1)
      return;

    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      this.hierarchyColumnNames = categoricalColumns.slice(0, this.hierarchyLevel).map((col) => col.name);

    super.onTableAttached();
  }

  getSeriesData(): treeDataType[] | undefined {
    return TreeUtils.toForest(this.dataFrame, this.hierarchyColumnNames, this.dataFrame.filter);
  }

  render(): void {
    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      return;

    this.option.series[0].data = this.getSeriesData();

    this.chart.setOption(this.option);
  }
}
