import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {EChartViewer} from '../echart/echart-viewer';
import {TreeUtils, treeDataType} from '../../utils/tree-utils';
import { StringUtils } from '@datagrok-libraries/utils/src/string-utils';

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
      animation: false,
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

  isCanvasEmpty(ctx: any, x: any, y: any) {
    const pixel = ctx.getImageData(x, y, 1, 1).data;
    return pixel[3] === 0;
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
    this.chart.on('mouseover', (params: any) => {
      const divs: HTMLElement[] = [];
      const { hierarchyColumnNames, dataFrame } = this;
      for (const columnName of hierarchyColumnNames) {
        const column = dataFrame.columns.byName(columnName);
        const idx = Array.from(column.values()).map((item) => String(item)).indexOf(params.name);
        if (idx !== -1) {
          for (let j = 0; j < dataFrame.columns.length; ++j) {
            const columnAtIndex = dataFrame.columns.byIndex(j);
            const value = columnAtIndex.get(idx);
            const formattedValue = typeof value === 'string' ? value : StringUtils.formatNumber(value);
            divs[j] = ui.divText(`${columnAtIndex.name} : ${formattedValue}`);
          }
        }
      }
      ui.tooltip.show(ui.div(divs), params.event.event.x, params.event.event.y);
    });      
    this.chart.on('mouseout', () => ui.tooltip.hide());
    this.chart.getDom().ondblclick = (event: MouseEvent) => {
      const canvas = this.chart.getDom().querySelector('canvas');
      const rect = canvas!.getBoundingClientRect();
      const scaleX = canvas!.width / rect.width;
      const scaleY = canvas!.height / rect.height;
      const clickX = (event.clientX - rect.left) * scaleX;
      const clickY = (event.clientY - rect.top) * scaleY;
      if (this.isCanvasEmpty(canvas!.getContext('2d'), clickX, clickY)) {
        this.render();
      }
    };
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
    
    this.dataFrame.onMetadataChanged.subscribe((_) => {this.render()});
    super.onTableAttached();
  }

  getSeriesData(): Promise<treeDataType[] | undefined> {
    return TreeUtils.toForest(this.dataFrame, this.hierarchyColumnNames, this.dataFrame.filter);
  }

  render(): void {
    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      return;
  
    this.getSeriesData().then((seriesData: treeDataType[] | undefined) => {
      if (seriesData) {
        this.option.series[0].data = seriesData;
        this.chart.setOption(this.option);
      }
    }).catch((error) => {
      console.error(error);
    });
  }
  
}
