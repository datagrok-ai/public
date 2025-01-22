import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import { EChartViewer } from '../echart/echart-viewer';
import { TreeUtils, TreeDataType } from '../../utils/tree-utils';

import * as utils from '../../utils/utils';

type onClickOptions = 'Select' | 'Filter';
const rowSourceMap: Record<onClickOptions, string> = {
  Select: 'Filtered',
  Filter: 'All'
};
/// https://echarts.apache.org/examples/en/editor.html?c=tree-basic
@grok.decorators.viewer({
  name: 'Tree',
  description: 'Creates a tree viewer',
  trellisable: false,
  icon: 'icons/tree-viewer.svg',
})
export class TreeViewer extends EChartViewer {
  private renderQueue: Promise<void> = Promise.resolve();
  layout: layoutType;
  orient: orientation;
  expandAndCollapse: boolean;
  initialTreeDepth: number;
  symbolSizeRange: [number, number] = [5, 20];
  edgeShape: edgeShape;
  symbol: symbolType;
  symbolSize: number;
  fontSize: number;
  showCounts: boolean;
  sizeColumnName: string = '';
  sizeAggrType: DG.AggregationType = 'avg';
  colorColumnName: string = '';
  colorAggrType: DG.AggregationType = 'avg';
  hierarchyColumnNames: string[];
  aggregations: string[] = Object.values(DG.AGG).filter((f) => f !== DG.AGG.KEY && f !== DG.AGG.PIVOT);
  aggregationsStr: string[] = Object.values({...DG.STR_AGG, ...DG.STAT_COUNTS});
  selectionColor = DG.Color.toRgb(DG.Color.selectedRows);
  applySizeAggr: boolean = false;
  applyColorAggr: boolean = false;
  onClick: onClickOptions;
  includeNulls: boolean;
  labelRotate: number;

  constructor() {
    super();

    this.layout = <layoutType> this.string('layout', 'orthogonal', { choices: ['orthogonal', 'radial'], category: 'Style'});
    this.orient = <orientation> this.string('orient', 'LR', { choices: ['LR', 'RL', 'TB', 'BT'], category: 'Style' });
    this.expandAndCollapse = this.bool('expandAndCollapse', true, {category: 'Style'});
    this.initialTreeDepth = this.int('initialTreeDepth', 3, { min: -1, max: 5, category: 'Style'});
    this.edgeShape = <edgeShape> this.string('edgeShape', 'curve', { choices: ['curve', 'polyline'], category: 'Style' });
    this.symbol = <symbolType> this.string('symbol', 'emptyCircle', { choices: [
      'circle', 'emptyCircle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none',
    ], category: 'Style' });
    this.symbolSize = this.int('symbolSize', 7, {category: 'Style'});
    this.fontSize = this.int('fontSize', 12, {category: 'Style', max: 30});
    this.labelRotate = this.int('labelRotate', 45, {category: 'Style', max: 360});
    this.showCounts = this.bool('showCounts', false, {category: 'Style'});

    this.sizeColumnName = this.string('sizeColumnName', '', {category: 'Size'});
    this.sizeAggrType = <DG.AggregationType> this.string('sizeAggrType', DG.AGG.AVG, { choices: this.aggregations, category: 'Size' });
    
    this.colorColumnName = this.string('colorColumnName', '', {category: 'Color'});
    this.colorAggrType = <DG.AggregationType> this.string('colorAggrType', DG.AGG.AVG, { choices: this.aggregations, category: 'Color' });

    this.hierarchyColumnNames = this.addProperty('hierarchyColumnNames', DG.TYPE.COLUMN_LIST, null, {category: 'Data', columnTypeFilter: DG.TYPE.CATEGORICAL});
    this.onClick = <onClickOptions> this.string('onClick', 'Select', { choices: ['Select', 'Filter']});
    this.includeNulls = this.bool('includeNulls', true, {category: 'Value'});

    this.option = {
      series: [
        {
          type: 'tree',
          label: {
            position: 'left',
            verticalAlign: 'middle',
            align: 'right',
            rotate: this.labelRotate,
            fontSize: this.fontSize,
          },

          leaves: {
            label: {
              position: 'right',
              verticalAlign: 'middle',
              align: 'left',
              rotate: this.labelRotate,
            },
          },
        },
      ],
    };

    this.onPropertyChanged(null);
  }

  handleDataframeSelection(path: string[], event: any) {
    this.dataFrame.selection.handleClick((index: number) => {
      if (!this.filter.get(index) && this.rowSource !== 'Selected')
        return false;

      return path.every((segment, j) => {
        const columnValue = this.dataFrame.getCol(this.hierarchyColumnNames[j]).get(index);
        return (columnValue && columnValue.toString() === segment) || (!columnValue && segment === '');
      });
    }, event);
  }

  handleDataframeFiltering(path: string[], dataFrame: DG.DataFrame) {
    const filterFunction = this.buildFilterFunction(path);
    dataFrame.rows.filter(filterFunction);
  }

  buildFilterFunction(path: string[]): (row: any) => boolean {
    return (row) => {
      return path.every((expectedValue, i) => {
        const column = this.dataFrame.getCol(this.hierarchyColumnNames[i]);
        const columnValue = row.get(this.hierarchyColumnNames[i]);
        const formattedValue = columnValue ?
          (column.type !== DG.TYPE.STRING ? columnValue.toString() : columnValue) : '';
        return formattedValue === expectedValue;
      });
    };
  }

  initChartEventListeners() {
    let selectedSectors: string[] = [];
    const handleChartClick = (params: any) => {
      const path = params.treeAncestors.slice(2).map((obj: any) => obj.name);
      const pathString = path.join('|');
      const isSectorSelected = selectedSectors.includes(pathString);
      if (this.onClick === 'Filter') {
        this.handleDataframeFiltering(path, this.dataFrame);
        return;
      }

      const event = params.event.event;
      const isMultiSelect = event.shiftKey || event.ctrlKey || event.metaKey;
      const isMultiDeselect = (event.shiftKey && event.ctrlKey) || (event.shiftKey && event.metaKey);
      if (isMultiSelect && !isSectorSelected)
        selectedSectors.push(pathString);
      else if (isMultiDeselect && isSectorSelected)
        selectedSectors = selectedSectors.filter((sector) => sector !== pathString);
      this.handleDataframeSelection(path, event);
    };
    
    this.chart.on('mouseover', async (params: any) => {
      const ancestors = params.treeAncestors.filter((item: any) => item.name).map((item: any) => item.name).join('.');
      let div = ui.divV([]);
      if (params.data.semType === DG.SEMTYPE.MOLECULE) {
        const image = await TreeUtils.getMoleculeImage(params.name, 150, 100);
        const { width, height } = image;
        if (width && height) {
          const pixels = image!.getContext('2d')!.getImageData(0, 0, width, height).data;
          if (pixels.some((_, i) => i % 4 === 3 && pixels[i] !== 0))
            div.appendChild(image);
        }
      }
      div.appendChild(ui.divText(ancestors));
      div.appendChild(ui.divText(params.value, { style: { fontWeight: 'bold' } }));
      ui.tooltip.show(div, params.event.event.x, params.event.event.y);
    });
    this.chart.on('mouseout', () => ui.tooltip.hide());
    this.chart.on('click', handleChartClick);
  }

  onPropertyChanged(p: DG.Property | null, render: boolean = true) {
    if (p?.name === 'edgeShape') {
      this.getProperty('layout')?.set(this, 'orthogonal');
      this.option.series[0].layout = 'orthogonal';
      //this.option.series[0].label.rotate = 0;
      this.chart.clear();
    }
    if (p?.name === 'layout') {
      const layout: layoutType = p.get(this);
      this.option.series[0].layout = layout;
      
      const es = this.getProperty('edgeShape');
      if (es?.get(this) !== 'curve') {
        es?.set(this, 'curve');
        this.option.series[0].edgeShape = 'curve';
      }
      this.render();
      return;
    }
    if (p?.name === 'table') {
      this.updateTable();
      this.onTableAttached();
      this.render();
    }
    if (p?.name === 'onClick')
      this.rowSource = rowSourceMap[this.onClick as onClickOptions] || this.rowSource;
    if (p?.name === 'hierarchyColumnNames' || p?.name === 'sizeColumnName' ||
        p?.name === 'sizeAggrType' || p?.name === 'colorColumnName' || p?.name === 'colorAggrType' ||
        p?.name === 'fontSize' || p?.name === 'showCounts' || p?.name === 'includeNulls' || p?.name === 'orient' || p?.name === 'labelRotate') {
      if (p?.name === 'hierarchyColumnNames')
        this.chart.clear();
      if (p?.name === 'colorColumnName' || p?.name === 'colorAggrType')
        this.applyColorAggr = this.shouldApplyAggregation(this.colorColumnName, this.colorAggrType);
      if (p?.name === 'sizeColumnName' || p?.name === 'sizeAggrType')
        this.applySizeAggr = this.shouldApplyAggregation(this.sizeColumnName, this.sizeAggrType);
      if (p?.name === 'fontSize')
        this.option.series[0].label.fontSize = p.get(this);
      if (p?.name === 'orient' || p?.name === 'labelRotate') {
        if (p?.name === 'orient')
          this.option.series[0].orient = p.get(this);
        this.updateOrient();
      }
      if (p?.name === 'includeNulls')
        this.chart.clear();
      this.render();
    } else
      super.onPropertyChanged(p, render);
  }

  updateOrient() {
    // Default option structure
    let labelOptions = {
      position: 'left',
      verticalAlign: 'middle',
      align: 'right',
      rotate: this.labelRotate,
      fontSize: this.fontSize,
    };
  
    let leavesLabelOptions = {
      position: 'right',
      verticalAlign: 'middle',
      rotate: this.labelRotate,
      align: 'left',
    };
  
    // Update label options based on orientation
    switch (this.orient) {
      case 'LR':
        labelOptions = {
          position: 'left',
          verticalAlign: 'middle',
          align: 'right',
          rotate: this.labelRotate,
          fontSize: this.fontSize,
        };
        leavesLabelOptions = {
          position: 'right',
          verticalAlign: 'middle',
          rotate: this.labelRotate,
          align: 'left',
        };
        break;
  
      case 'RL':
        labelOptions = {
          position: 'right',
          verticalAlign: 'middle',
          align: 'left',
          rotate: this.labelRotate,
          fontSize: this.fontSize
        };
        leavesLabelOptions = {
          position: 'left',
          verticalAlign: 'middle',
          rotate: this.labelRotate,
          align: 'right',
        };
        break;
  
      case 'BT':
        labelOptions = {
          position: 'bottom',
          rotate: this.labelRotate,
          verticalAlign: 'middle',
          align: 'right',
          fontSize: this.fontSize
        };
        leavesLabelOptions = {
          position: 'top',
          rotate: this.labelRotate,
          verticalAlign: 'middle',
          align: 'left',
        };
        break;
  
      case 'TB':
        labelOptions = {
          position: 'top',
          rotate: -(this.labelRotate),
          verticalAlign: 'middle',
          align: 'right',
          fontSize: this.fontSize
        };
        leavesLabelOptions = {
          position: 'bottom',
          rotate: -(this.labelRotate),
          verticalAlign: 'middle',
          align: 'left',
        };
        break;

      default:
        break;
    }
  
    this.option.series[0].label = labelOptions;
    this.option.series[0].leaves.label = leavesLabelOptions;
  }  

  shouldApplyAggregation(columnName: string, aggrType: string): boolean {
    const numericalColumns = this.dataFrame.columns.byName(columnName);
    const isColumnNumerical = numericalColumns ? numericalColumns.matches('numerical') : false;
    const isAggregationApplicable = this.aggregationsStr.includes(aggrType);
    return isColumnNumerical || (!isColumnNumerical && isAggregationApplicable);
  }

  _testColumns() {
    return this.dataFrame.columns.length >= 1;
  }

  onTableAttached() {
    const categoricalColumns = [...this.dataFrame.columns.categorical].sort((col1, col2) =>
      col1.categories.length - col2.categories.length);

    if (categoricalColumns.length < 1)
      return;

    this.filter = this.dataFrame.filter;
    this.hierarchyColumnNames = categoricalColumns.slice(0, 3).map((col) => col.name);
    this.sizeColumnName = '';
    this.colorColumnName = '';

    super.onTableAttached();
    this.subs.push(this.dataFrame.onColumnsRemoved.subscribe((data) => {
      const columnNamesToRemove = data.columns.map((column: DG.Column) => column.name);
      this.hierarchyColumnNames = this.hierarchyColumnNames.filter((columnName) => !columnNamesToRemove.includes(columnName));
      this.render();
    }));
    this.initChartEventListeners();
    this.helpUrl = 'https://datagrok.ai/help/visualize/viewers/tree';
  }

  colorCodeTree(data: TreeDataType): void {
    const min = data['color-meta']['min'];
    const max = data['color-meta']['max'];
    const setItemStyle = (data: TreeDataType) => {
      const safeMin = isFinite(min) && !isNaN(min) ? min : 0;
      const safeMax = isFinite(max) && !isNaN(max) ? max : 0;

      const nodeColor = DG.Color.toRgb(DG.Color.scaleColor(data.color, safeMin, safeMax));
      data.itemStyle = data.itemStyle ?? {};
      data.itemStyle.color = data.itemStyle.color === this.selectionColor ? this.selectionColor : nodeColor;
      if (data.children)
        data.children.forEach((child) => setItemStyle(child));
    };
    setItemStyle(data);
  }

  async getSeriesData() {
    const aggregations = [];

    if (this.sizeColumnName && this.applySizeAggr)
      aggregations.push({ type: <DG.AggregationType> this.sizeAggrType,
        columnName: this.sizeColumnName, propertyName: 'size' });

    if (this.colorColumnName && this.applyColorAggr)
      aggregations.push({ type: <DG.AggregationType> this.colorAggrType,
        columnName: this.colorColumnName, propertyName: 'color' });
    return [await TreeUtils.toTree(this.dataFrame, this.hierarchyColumnNames, this.filter, null, aggregations, true, undefined, undefined, this.includeNulls)];
  }

  async renderMolecule(params: any, width: number, height: number) {
    const image = await TreeUtils.getMoleculeImage(params.name, width, height);
    const img = new Image();
    img.src = image!.toDataURL('image/png');
    params.data.label = {
      show: true,
      formatter: '{b}',
      color: 'rgba(0,0,0,0)',
      height: height.toString(),
      width: width.toString(),
      backgroundColor: {
        image: img.src,
      },
    }
  }

  formatLabel(params: any) {
    //@ts-ignore
    const ItemAreaInfoArray = this.chart.getModel().getSeriesByIndex(0).getData()._itemLayouts.slice(1);
    const getCurrentItemIndex = params.dataIndex - 1;
    const ItemLayoutInfo = ItemAreaInfoArray.find((item: any, index: number) => getCurrentItemIndex === index);

    const { x, y } = ItemLayoutInfo;
    const isVerticalOrientation = this.isVerticalOrientation();

    const sortedItems = [...ItemAreaInfoArray]
      .filter((item: any) => item && (item.x !== undefined && item.y !== undefined))
      .sort((a: any, b: any) => isVerticalOrientation ? a.y - b.y : a.x - b.x);

    let positions: number[];
    let distances: number[];

    if (isVerticalOrientation) {
      positions = sortedItems.map((item: any) => item.y);
      distances = positions.slice(1).map((y: number, index: number) => y - positions[index]);
    } else {
      positions = sortedItems.map((item: any) => item.x);
      distances = positions.slice(1).map((x: number, index: number) => x - positions[index]);
    }

    const averageDistance = distances.length ? distances.reduce((acc, val) => acc + val, 0) / distances.length : 0;
    const chartSize = isVerticalOrientation ? this.chart.getHeight() : this.chart.getWidth();
    const tolerance = Math.max(20, Math.min(averageDistance, chartSize * 0.05));

    const sortedPositions = [...new Set(positions)];
    const nodesAtLevel = sortedItems.filter((item: any) => Math.abs(isVerticalOrientation ? item.y - y : item.x - x) <= tolerance);

    const index = sortedPositions.indexOf(isVerticalOrientation ? y : x);
    const availableSpace = (index === 0)
      ? (sortedPositions[1] - sortedPositions[0])
      : (index === sortedPositions.length - 1)
      ? (chartSize - sortedPositions[index])
      : (sortedPositions[index + 1] - sortedPositions[index - 1]) / 2;

    const averageCharWidth = this.fontSize;
    const maxWidthChars = Math.floor(availableSpace / averageCharWidth);
    
    let labelText = this.showCounts ? `${params.name}: ${params.value}` : `${params.name}`;
    if (labelText.length > maxWidthChars) {
      labelText = labelText.slice(0, maxWidthChars - 3) + '...';
    }

    const labelHeight = Math.floor(chartSize / nodesAtLevel.length);
    const maxHeightChars = Math.floor(labelHeight / this.fontSize);

    if (params.data.semType === 'Molecule') {
      const minImageWidth = 70;
      const minImageHeight = 80;

      if (availableSpace >= minImageWidth && labelHeight >= minImageHeight) {
        this.renderMolecule(params, minImageWidth, minImageHeight);
        return ' ';
      }
      return ' ';
    }

    const lines = labelText.split(' ');
    let result = '';
    let remainingHeight = maxHeightChars;

    for (const line of lines) {
      if (remainingHeight <= 0) break;
      result += line + ' ';
      remainingHeight--;
    }

    if (result.trim() === '...') {
      result = '';
    }

    return result.trim();
  }

  isVerticalOrientation() {
    return this.orient === 'BT' || this.orient === 'TB';
  }

  render(): void {
    this.renderQueue = this.renderQueue
      .then(() => this._render());
  }

  async _render() {
    if (!this.dataFrame)
      return;

    if (this.hierarchyColumnNames?.some((colName) => !this.dataFrame.columns.names().includes(colName)))
      this.hierarchyColumnNames = this.hierarchyColumnNames.filter((value) => this.dataFrame.columns.names().includes(value));
    if (this.hierarchyColumnNames == null || this.hierarchyColumnNames.length === 0)
      return;

    const data = await this.getSeriesData();

    Object.assign(this.option.series[0], {
      data,
      label: Object.assign({}, this.option.series[0].label, {
        formatter: (params: any) => this.formatLabel(params),
      }),
    });    

    this.option.series[0]['symbolSize'] = this.sizeColumnName && this.applySizeAggr ?
      (value: number, params: {[key: string]: any}) => utils.data.mapToRange(
        params.data.size, this.option.series[0].data[0]['size-meta']['min'],
        this.option.series[0].data[0]['size-meta']['max'], ...this.symbolSizeRange) : this.symbolSize;
    if (this.colorColumnName && this.applyColorAggr)
      this.colorCodeTree(this.option.series[0].data[0]);

    console.log(this.option);
    this.chart.setOption(this.option);
  }
}

type layoutType = 'orthogonal' | 'radial';
type orientation = 'LR' | 'RL' | 'TB' | 'BT';
type edgeShape = 'curve' | 'polyline';
type symbolType = 'circle' | 'emptyCircle' | 'rect' | 'roundRect' | 'triangle' | 'diamond' | 'pin' | 'arrow' | 'none';
