import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {EChartViewer} from '../echart/echart-viewer';
import {TreeUtils, treeDataType} from '../../utils/tree-utils';
import { delay } from '@datagrok-libraries/utils/src/test';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic

type onClickOptions = 'Select' | 'Filter';

/** Represents a sunburst viewer */
@grok.decorators.viewer({
  name: 'Sunburst',
  description: 'Creates a sunburst viewer',
  icon: 'icons/sunburst-viewer.svg',
})

export class SunburstViewer extends EChartViewer {
  hierarchyColumnNames: string[];
  hierarchyLevel: number;
  onClick: onClickOptions;

  constructor() {
    super();
    this.initCommonProperties();
    this.initEventListeners();

    this.hierarchyColumnNames = this.addProperty('hierarchyColumnNames', DG.TYPE.COLUMN_LIST);
    this.hierarchyLevel = 3;
    this.onClick = <onClickOptions> this.string('onClick', 'Select', { choices: ['Select', 'Filter'] });

    this.option = {
      animation: false,
      series: [
        {
          type: 'sunburst',
          nodeClick: false,
          label: {
            rotate: 'radial',
            fontSize: 10,
          }
        },
      ],
    };

    this.onPropertyChanged(null);
  }

  isCanvasEmpty(ctx: any, x: any, y: any) {
    const pixel = ctx.getImageData(x, y, 1, 1).data;
    return pixel[3] === 0;
  }

  handleDataframeSelection(path: string[], event: any) {
    this.dataFrame.selection.handleClick((i) => {
      if (!this.filter.get(i))
        return false;
      for (let j = 0; j < path.length; j++) {
        if (this.dataFrame.getCol(this.hierarchyColumnNames[j]).get(i).toString() !== path[j])
          return false;
      }
      return true;
    }, event);
  }

  createQueryMatcher(path: string[]): DG.RowMatcher {
    const conditions = path.map((value, i) => {
      return `${this.hierarchyColumnNames[i]} = ${value}`;
    }).join(' and ');
    return this.dataFrame.rows.match(conditions);
  }

  removeFiltering() {
    if (this.dataFrame.filter.trueCount !== this.dataFrame.rowCount) {
      this.dataFrame.filter.setAll(true);
    }
  }

  initEventListeners(): void {
    this.chart.on('click', (params: any) => {
      const selectedSectors: string[] = [];
      if (!params.data.path)
        return;
      const path: string[] = params.data.path.split('|').map((str: string) => str.trim());
      const pathString: string = path.join('|');
      if (this.onClick === 'Filter') {
        this.createQueryMatcher(path).filter();
        return;
      }
      const isSectorSelected = selectedSectors.includes(pathString);
      if (params.event.event.shiftKey || params.event.event.ctrlKey || params.event.event.metaKey) {
        if (!isSectorSelected) {
          selectedSectors.push(pathString);
          this.handleDataframeSelection(path, params.event.event);
        }
      } else if ((params.event.event.shiftKey && params.event.event.ctrlKey) || 
                (params.event.event.shiftKey && params.event.event.metaKey)) {
        if (isSectorSelected) {
          const index = selectedSectors.indexOf(pathString);
          selectedSectors.splice(index, 1);
          this.handleDataframeSelection(path, params.event.event);
        }
      } else {
        this.handleDataframeSelection(path, params.event.event);
      }
    });
    this.chart.on('mouseover', (params: any) => {
      const path: string[] = params.data.path.split('|').map((str: string) => str.trim());
      const matchDf = this.createQueryMatcher(path).toDataFrame();
      const matchCount = matchDf.rowCount;
      ui.tooltip.showRowGroup(this.dataFrame, (i) => {
        const { hierarchyColumnNames, dataFrame } = this;
        for (let j = 0; j < hierarchyColumnNames.length; ++j) {
          const column = dataFrame.getCol(hierarchyColumnNames[j]);
          const format = column.getTag(DG.TAGS.FORMAT);
          if (format) {
            const number = format.indexOf('.');
            const len = format.length - number - 1;
            if ((column.get(i)).toFixed(len) === params.name)
              return true;
          }
          if (column.get(i).toString() === params.name) {
            return true;
          }
        }
        return false;
      }, params.event.event.x, params.event.event.y);
      ui.tooltip.root.innerText = `${matchCount}\n${params.name}`;
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
      this.removeFiltering();
    };
  }

  onContextMenuHandler(menu: DG.Menu): void {
    menu.item('Reset View', () => {
      this.render();
      this.removeFiltering();
    });
  }

  onPropertyChanged(p: DG.Property | null, render: boolean = true): void {
    if (p?.name === 'hierarchyColumnNames')
      this.render();
    if (p?.name === 'table') {
      this.updateTable();
      this.onTableAttached(true);
    }
    else
      super.onPropertyChanged(p, render);
  }

  onTableAttached(propertyChanged?: boolean): void {
    let categoricalColumns = [...this.dataFrame.columns.categorical].sort((col1, col2) =>
      col1.categories.length - col2.categories.length);
    categoricalColumns = categoricalColumns.filter((col: DG.Column) => col.stats.missingValueCount != col.length);

    if (categoricalColumns.length < 1)
      return;

    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0 || propertyChanged)
      this.hierarchyColumnNames = categoricalColumns.slice(0, this.hierarchyLevel).map((col) => col.name);
    
    this.subs.push(this.dataFrame.onMetadataChanged.subscribe((_) => {this.render()}));
    this.subs.push(this.onContextMenu.subscribe(this.onContextMenuHandler.bind(this)));
    this.addSelectionOrDataSubs();
    this.render();
  }

  getSeriesData(): treeDataType[] | undefined {
    return TreeUtils.toForest(this.dataFrame, this.hierarchyColumnNames, this.filter);
  }

  formatLabel(params: any) {
    //@ts-ignore
    const ItemAreaInfoArray = this.chart.getModel().getSeriesByIndex(0).getData()._itemLayouts.slice(1);
    const getCurrentItemIndex = params.seriesIndex;
    const ItemLayoutInfo = ItemAreaInfoArray.find((item: any, index: number) => {
        if (getCurrentItemIndex === index) {
            return item;
        }
    });
    const r = ItemLayoutInfo.r;
    const startAngle = ItemLayoutInfo.startAngle;
    const endAngle = ItemLayoutInfo.endAngle;
    const cx = ItemLayoutInfo.cx;
    const cy = ItemLayoutInfo.cy;
    const width = Math.abs(cx + r * Math.cos(startAngle) - (cx + r * Math.cos(endAngle)));
    const height = Math.abs(cy + r / 1.5 * Math.sin(startAngle) - (cy + r / 1.5 * Math.sin(endAngle)));

    const averageCharWidth = 10 * 0.6;
    const averageCharHeight = 10 * 1.2;
    const maxWidthCharacters = Math.floor(width / averageCharWidth);
    const maxHeightCharacters = Math.floor(height / averageCharHeight);

    const maxLength = maxWidthCharacters;
    let name = params.name;
    let lines = [name];

    if (name.length > maxLength) {
      lines = [];
      let remainingHeight = maxHeightCharacters;
      while (name.length > 0 && remainingHeight > 0) {
        let line = name.substring(0, maxLength);
        const lastSpaceIndex = line.lastIndexOf(' ');
        if (lastSpaceIndex !== -1) {
          line = line.substring(0, lastSpaceIndex);
        }
        line = line.trimRight();
        lines.push(line);  
        remainingHeight -= 1;
        name = name.substring(line.length).trim();
      }
      if (name.length > 0) {
        lines[lines.length - 1] += '...';
      }
    }
    return lines.join('\n');
  }

  async handleStructures(data: treeDataType[] | undefined) {
    for (const entry of data!) {
      const name = entry.name;
      const isSmiles = await grok.functions.call('Chem:isSmiles', {s: name});
      if (isSmiles && entry.semType === 'Molecule') {
        const imageContainer = await grok.functions.call('Chem:drawMolecule', {
          'molStr': name, 'w': 70, 'h': 80, 'popupMenu': false
        });
        const image = imageContainer.querySelector(".chem-canvas");
        await delay(5);
        const img = new Image();
        img.src = image.toDataURL('image/png');
        entry.label = {
          show: true,
          formatter: '{b}',
          color: 'rgba(0,0,0,0)',
          height: '80',
          width: '70',
          backgroundColor: {
            image: img.src,
          },
        }
      } 
      if (entry.children) {
        await this.handleStructures(entry.children);
      }
    }
    return data;
  }

  render() {
    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      return;

    this.handleStructures(this.getSeriesData()).then((data) => {
      this.option.series[0].data = data;
      this.option.series[0].label.formatter = (params: any) => this.formatLabel(params);
      this.chart.setOption(this.option);
    });
  }

  detach() {
    for (const sub of this.subs)
      sub.unsubscribe();
    super.detach();
  }
  
}
