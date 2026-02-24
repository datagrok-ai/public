/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import { EChartViewer } from '../echart/echart-viewer';
import { TreeUtils, TreeDataType } from '../../utils/tree-utils';
import * as echarts from 'echarts';
import { fromEvent } from 'rxjs';
import { debounceTime } from 'rxjs/operators';
import _ from 'lodash';
import { ERROR_CLASS, MessageHandler } from '../../utils/utils';

/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic

type onClickOptions = 'Select' | 'Filter';
const CATEGORIES_NUMBER = 500;
const MAXIMUM_COLUMN_NUMBER = 20;
let sunburstId = 0;
const rowSourceMap: Record<onClickOptions, string> = {
  Select: 'Filtered',
  Filter: 'All',
};

/** Represents a sunburst viewer */
@grok.decorators.viewer({
  name: 'Sunburst',
  description: 'Creates a sunburst viewer',
  icon: 'icons/sunburst-viewer.svg',
})

export class SunburstViewer extends EChartViewer {
  private renderQueue: Promise<void> = Promise.resolve();
  hierarchyColumnNames: string[];
  eligibleHierarchyNames!: string[];
  hierarchyLevel: number;
  onClick: onClickOptions;
  selectedOptions: string[] = ['Selected', 'SelectedOrCurrent', 'FilteredSelected'];
  inheritFromGrid: boolean;
  title: string;
  sunburstVersion: number | null = null;
  currentVersion: number | null = null;
  includeNulls: boolean;
  private moleculeRenderQueue: Promise<void> = Promise.resolve();
  private latestRenderToken = 0;
  viewerFilter: DG.BitSet | null = null;
  constructor() {
    super();
    this.initEventListeners();

    this.title = this.string('title', 'Sunburst', {category: 'Description'});
    this.hierarchyColumnNames = this.addProperty('hierarchyColumnNames', DG.TYPE.COLUMN_LIST, null,
      {columnTypeFilter: DG.TYPE.CATEGORICAL});
    this.hierarchyLevel = 3;
    this.onClick = <onClickOptions> this.string('onClick', 'Select', { choices: ['Select', 'Filter']});
    this.inheritFromGrid = this.bool('inheritFromGrid', true, { category: 'Color' });
    this.includeNulls = this.bool('includeNulls', true, {category: 'Value'});

    this.option = {
      animation: false,
      silent: false,
      series: [
        {
          type: 'sunburst',
          radius: '98%',
          nodeClick: false,
          emphasis: {
            focus: 'series',
          },
          label: {
            rotate: 'radial',
            fontSize: 10,
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

  applySelectionFilter(bitset: DG.BitSet, path: string[], event: any) {
    bitset.handleClick((index: number) => {
      if (!this.filter.get(index) && this.rowSource !== 'Selected')
        return false;

      return path.every((segment, j) => {
        const columnValue = this.dataFrame.getCol(this.eligibleHierarchyNames[j]).get(index);
        return columnValue == null ? segment === ' ' : columnValue.toString() === segment;
      });
    }, event);
  }

  handleDataframeFiltering(path: string[], event: MouseEvent) {
    if (this.viewerFilter === null)
      this.viewerFilter = DG.BitSet.create(this.dataFrame.rowCount);

    this.viewerFilter.handleClick((index: number) => {
      return path.every((segment, j) => {
        const columnValue = this.dataFrame.getCol(this.eligibleHierarchyNames[j]).get(index);
        return columnValue == null ? segment === ' ' : columnValue.toString() === segment;
      });
    }, event);

    if (this.viewerFilter.trueCount === 0)
      this.viewerFilter = null;

    this.dataFrame.rows.requestFilter();
  }

  initEventListeners(): void {
    if (!this.chart) return;

    let selectedSectors: string[] = [];
    const handleChartClick = (params: any) => {
      const path = params.treePathInfo.slice(1).map((obj: any) => obj.name);
      const pathString = path.join('|');
      const isSectorSelected = selectedSectors.includes(pathString);
      const event = params.event.event;
      const isMultiSelect = event.shiftKey || event.ctrlKey || event.metaKey;
      const isMultiDeselect = (event.shiftKey && event.ctrlKey) || (event.shiftKey && event.metaKey);
      if (isMultiSelect && !isSectorSelected)
        selectedSectors.push(pathString);
      else if (isMultiDeselect && isSectorSelected)
        selectedSectors = selectedSectors.filter((sector) => sector !== pathString);

      if (this.onClick === 'Filter') {
        this.handleDataframeFiltering(path, event);
        return;
      } else
        this.applySelectionFilter(this.dataFrame.selection, path, event);
    };

    const handleChartMouseover = async (params: any) => {
      const { x, y } = params.event.event;
      const { name, value, data } = params;
      const displayName = name || 'Nulls';
      const tooltipDiv = ui.div();

      ui.tooltip.show(tooltipDiv, x + 10, y);

      if (data.semType !== DG.SEMTYPE.MOLECULE) {
        tooltipDiv.innerText = `${value}\n${displayName}`;
        return;
      }

      const image = await TreeUtils.getMoleculeImage(name, 150, 100);
      tooltipDiv.appendChild(ui.divText(`${value}\n`));
      if (name)
        tooltipDiv.appendChild(image);
      else
        tooltipDiv.appendChild(ui.divText(displayName));
    };

    const handleCanvasDblClick = (event: MouseEvent) => {
      const canvas = this.chart?.getDom().querySelector('canvas');
      if (!canvas) return;

      const { left, top, width, height } = canvas.getBoundingClientRect();
      const scaleX = canvas.width / width;
      const scaleY = canvas.height / height;
      const clickX = (event.clientX - left) * scaleX;
      const clickY = (event.clientY - top) * scaleY;

      if (this.isCanvasEmpty(canvas.getContext('2d'), clickX, clickY)) {
        this.viewerFilter = null;
        this.dataFrame.rows.requestFilter();
      }
    };

    this.chart.on('click', handleChartClick);
    this.chart.on('mouseout', () => ui.tooltip.hide());
    this.chart.getDom().ondblclick = handleCanvasDblClick;

    this.subs.push(ui.onSizeChanged(this.root).subscribe((_) => {
      requestAnimationFrame(() => this.chart?.resize());
    }));

    fromEvent(this.chart, 'mouseover')
      .pipe(debounceTime(100))
      .subscribe((params: any) => handleChartMouseover(params));
  }

  onContextMenuHandler(menu: DG.Menu): void {
    menu.item('Reset View', () => {
      this.viewerFilter = null;
      this.dataFrame.rows.requestFilter();
    });
  }

  onPropertyChanged(p: DG.Property | null, render: boolean = true): void {
    if (!p) return;
    switch (p.name) {
    case 'table':
      this.updateTable();
      this.onTableAttached(true);
      break;

    case 'onClick':
      this.rowSource = rowSourceMap[this.onClick as onClickOptions] || this.rowSource;
      break;

    default:
      this.render();
      break;
    }
  }

  addSubs() {
    if (!this.dataFrame)
      return;
    this.subs.push(this.dataFrame.onMetadataChanged.subscribe((_) => this.render()));
    this.subs.push(grok.events.onEvent('d4-grid-color-coding-changed').subscribe(() => {
      if (this.inheritFromGrid)
        this.render();
    }));
    this.subs.push(this.dataFrame.onValuesChanged.subscribe((_) => this.render()));
    this.subs.push(grok.events.onEvent('d4-current-viewer-changed').subscribe((args) => {
      const {viewer} = args.args;
      if (viewer instanceof SunburstViewer)
        this.currentVersion = viewer.sunburstVersion;
    }));
    this.subs.push(grok.events.onEvent('d4-drag-drop').subscribe((args) => {
      if (this.sunburstVersion != this.currentVersion) return;
      const grid = (args.args.dragObject.grid as DG.Grid);
      const gridOrder: Int32Array = new Int32Array(grid.getRowOrder().buffer);
      const names = this.hierarchyColumnNames;
      this.hierarchyColumnNames = Array.from(gridOrder)
        .map((index) => grid.table.row(index).get('name'))
        .filter((columnName) => names.includes(columnName!));
      this.render();
    }));
    this.subs.push(this.onContextMenu.subscribe(this.onContextMenuHandler.bind(this)));
    this.subs.push(this.dataFrame.onColumnsRemoved.subscribe((data) => {
      const columnNamesToRemove = data.columns.map((column: DG.Column) => column.name);
      this.hierarchyColumnNames = this.hierarchyColumnNames.filter((columnName) =>
        !columnNamesToRemove.includes(columnName));
      this.render();
    }));
    this.subs.push(this.dataFrame.onRowsFiltering.subscribe(() => {
      if (this.viewerFilter)
        this.dataFrame.filter.and(this.viewerFilter);
    }));
    this.subs.push(DG.debounce(grok.events.onResetFilterRequest, 10).subscribe(() => {
      this.viewerFilter = null;
    }));
  }

  onTableAttached(propertyChanged?: boolean): void {
    let categoricalColumns = [...this.dataFrame.columns.categorical].sort((col1, col2) =>
      col1.categories.length - col2.categories.length);
    categoricalColumns = categoricalColumns.filter((col: DG.Column) => col.stats.missingValueCount != col.length &&
      !col.name.startsWith('~'));

    if (categoricalColumns.length < 1)
      return;

    this.hierarchyColumnNames = categoricalColumns.slice(0, this.hierarchyLevel).map((col) => col.name);
    this.sunburstVersion = sunburstId;
    sunburstId++;

    this.addSubs();
    this.render();
  }

  _showMessage(msg: string, className: string) {
    const errorDiv = ui.divText(msg, className);
    errorDiv.style.textAlign = 'center';
    this.root.appendChild(errorDiv);
  }

  async getSeriesData(): Promise<TreeDataType[] | undefined> {
    const rowSource = this.selectedOptions.includes(this.rowSource!);
    return await TreeUtils.toForest(this.dataFrame, this.eligibleHierarchyNames, this.filter,
      this.includeNulls, rowSource, this.inheritFromGrid);
  }

  async renderMolecule(params: any, width: number, height: number) {
    const image = await TreeUtils.getMoleculeImage(params.name, width, height);
    const img = new Image();
    img.src = image!.toDataURL('image/png');
    params.data.label = {
      show: true,
      color: 'rgba(0,0,0,0)',
      height: height.toString(),
      width: width.toString(),
      backgroundColor: {
        image: img.src,
      },
    };
  }

  async renderMoleculeQueued(params: any, width: number, height: number): Promise<void> {
    this.moleculeRenderQueue = this.moleculeRenderQueue.then(() =>
      this.renderMolecule(params, width, height),
    );
    await this.moleculeRenderQueue;
  }

  formatLabel(params: any): string {
    //@ts-ignore
    const ItemAreaInfoArray = this.chart.getModel().getSeriesByIndex(0).getData()._itemLayouts.slice(1);
    const getCurrentItemIndex = params.dataIndex - 1;
    const ItemLayoutInfo = ItemAreaInfoArray.find((item: any, index: number) => {
      if (getCurrentItemIndex === index)
        return item;
    });

    const { r, r0, startAngle, endAngle } = ItemLayoutInfo;
    const { width, height } = this.calculateRingDimensions(r0, r, startAngle, endAngle);

    if (params.data.semType === 'Molecule') {
      const minImageWidth = 70;
      const minImageHeight = 80;

      if (width >= minImageWidth && height >= minImageHeight) {
        const scaleByWidth = width / minImageWidth;
        const scaleByHeight = height / minImageHeight;
        const scale = Math.min(scaleByWidth, scaleByHeight);

        const renderWidth = Math.max(minImageWidth, minImageWidth * scale);
        const renderHeight = Math.max(minImageHeight, minImageHeight * scale);

        this.renderMoleculeQueued(params, renderWidth, renderHeight);
        return ' ';
      } else {
        if (params.data.label)
          delete params.data.label;
        return ' ';
      }
    }

    const averageCharWidth = 5;
    const averageCharHeight = 10;
    const maxWidthCharacters = Math.floor(width / averageCharWidth);
    const maxHeightCharacters = Math.floor(height / averageCharHeight);
    const maxLength = maxWidthCharacters;

    if (width < averageCharWidth || height < averageCharHeight)
      return ' ';

    const name = params.name;
    const lines = name.split(/[\s-]/);
    let result = '';
    let remainingHeight = maxHeightCharacters;

    for (const line of lines) {
      if (line.length > maxLength || remainingHeight <= 0) {
        result = '';
        break;
      }
      if (result.length > 0) {
        result += '\n';
        remainingHeight--;
      }
      result += line;
      if (result.length >= maxLength || remainingHeight <= 0)
        break;
    }

    if (result.length < name.length)
      result += '...';

    const resultWidth = result.length * averageCharWidth;
    const resultHeight = result.split('\n').length * averageCharHeight;
    if (resultWidth > width) {
      const maxChars = Math.floor(width / averageCharWidth) - 3;
      result = result.slice(0, maxChars) + '...';
    }

    if (resultHeight > height) {
      const maxLines = Math.floor(height / averageCharHeight) - 1;
      const truncatedLines = result.split('\n').slice(0, maxLines);
      result = truncatedLines.join('\n') + '...';
    }

    if (result === '...')
      result = ' ';

    return result;
  }

  calculateRingDimensions(innerRadius: number, outerRadius: number, startAngle: number, endAngle: number) {
    const width = outerRadius - innerRadius;
    const height = Math.abs(endAngle - startAngle) * outerRadius;
    return { height, width };
  }

  render(orderedHierarchyNames?: string[]): void {
    const currentToken = ++this.latestRenderToken;

    this.renderQueue = this.renderQueue
      .then(() => this._renderWithToken(currentToken, orderedHierarchyNames));
  }

  private async _renderWithToken(token: number, orderedHierarchyNames?: string[]) {
    if (token !== this.latestRenderToken)
      return;

    await this._render(orderedHierarchyNames);
  }

  async _render(orderedHierarchyNames?: string[]) {
    if (!this.dataFrame)
      return;

    if (this.filter.trueCount >= CATEGORIES_NUMBER) {
      this.eligibleHierarchyNames = (orderedHierarchyNames ?? this.hierarchyColumnNames).filter(
        (name) => {
          const column = this.dataFrame.col(name);
          if (column)
            return column.categories.length <= CATEGORIES_NUMBER;
          return false;
        },
      );
    } else {
      const validColumnNames = new Set(this.dataFrame.columns.names());
      this.eligibleHierarchyNames = (orderedHierarchyNames ?? this.hierarchyColumnNames)
        .filter((name) => validColumnNames.has(name));
    }

    if (!this.eligibleHierarchyNames.length) {
      this._showMessage('The Sunburst viewer requires at least one categorical column with fewer than 500 unique categories', ERROR_CLASS);
      return;
    }

    MessageHandler._removeMessage(this.root, ERROR_CLASS);

    this.eligibleHierarchyNames = this.eligibleHierarchyNames.slice(0, MAXIMUM_COLUMN_NUMBER);
    const data = await this.getSeriesData();
    Object.assign(this.option.series[0], {
      data,
    });

    // Reinitialize the chart (needed in order to prevent memory leak)
    if (this.chart) {
      this.chart.clear();
      this.chart.dispose();
      this.detach();
      this.chart = null;
    }

    this.chart = echarts.init(this.root);
    this.initEventListeners();
    this.addSubs();

    this.option.series[0].label.formatter = (params: any) => this.formatLabel(params);
    this.chart.setOption(this.option, false, true);
  }

  detach() {
    for (const sub of this.subs)
      sub.unsubscribe();
    super.detach();
  }
}
