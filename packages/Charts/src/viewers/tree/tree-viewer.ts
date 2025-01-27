import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import { EChartViewer } from '../echart/echart-viewer';
import { TreeUtils, TreeDataType } from '../../utils/tree-utils';

import * as utils from '../../utils/utils';
import * as echarts from 'echarts';

type onClickOptions = 'Select' | 'Filter';
const rowSourceMap: Record<onClickOptions, string> = {
  Select: 'Filtered',
  Filter: 'All'
};
const aggregationMap = new Map<string, string[]>([
  // Aggregation for numeric columns: int, bigint, float, qnum
  ['int', Object.values(DG.AGG)],
  ['bigint', Object.values(DG.AGG)],
  ['float', Object.values(DG.AGG)],
  ['qnum', Object.values(DG.AGG)],

  // Aggregation for string columns
  ['string', [
    ...Object.values(DG.STR_AGG),
    ...Object.values(DG.STAT_COUNTS)
  ]],

  // Aggregation for datetime columns
  ['datetime', [
    ...Object.values(DG.STAT_COUNTS), 
    DG.AGG.MIN, 
    DG.AGG.MAX, 
    DG.AGG.AVG
  ]],

  // Aggregation for virtual columns
  ['virtual', [
    DG.AGG.TOTAL_COUNT, 
    DG.AGG.MISSING_VALUE_COUNT
  ]]
]);
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
  aggregations: string[] = [
    ...Object.values(DG.AGG),
    ...Object.values(DG.STR_AGG),
    ...Object.values(DG.STAT_COUNTS)
  ].filter((f) => f !== DG.AGG.KEY && f !== DG.AGG.PIVOT);
  aggregationsStr: string[] = Object.values({...DG.STR_AGG, ...DG.STAT_COUNTS});
  selectionColor = DG.Color.toRgb(DG.Color.selectedRows);
  applySizeAggr: boolean = false;
  applyColorAggr: boolean = false;
  onClick: onClickOptions;
  includeNulls: boolean;
  labelRotate: number;

  private moleculeRenderQueue: Promise<void> = Promise.resolve();
  constructor() {
    super();

    this.layout = <layoutType> this.string('layout', 'orthogonal', { choices: ['orthogonal', 'radial'], category: 'Style'});
    this.orient = <orientation> this.string('orient', 'LR', { choices: ['LR', 'RL', 'TB', 'BT'], category: 'Style' });
    this.initialTreeDepth = this.int('initialTreeDepth', 3, { min: 0, max: 5, category: 'Style'});
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
          roam: true,
          initialTreeDepth: this.initialTreeDepth,
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
              fontSize: this.fontSize,
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
        return (columnValue !== null && columnValue.toString() === segment) || (columnValue == null && segment === ' ');
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
        const formattedValue = columnValue !== null ?
          (column.type !== DG.TYPE.STRING ? columnValue.toString() : columnValue) : '';
        return formattedValue === expectedValue;
      });
    };
  }

  initChartEventListeners() {
    let selectedSectors: string[] = [];
  
    const handleChartClick = (path: string[], event: any) => {
      const pathString = path.join('|');
      const isSectorSelected = selectedSectors.includes(pathString);
  
      if (this.onClick === 'Filter') {
        this.handleDataframeFiltering(path, this.dataFrame);
        return;
      }

      const isMultiSelect = event.shiftKey || event.ctrlKey || event.metaKey;
      if (isMultiSelect) {
        if (isSectorSelected) {
          selectedSectors = selectedSectors.filter((sector) => sector !== pathString);
        } else {
          selectedSectors.push(pathString);
        }
      }
  
      this.handleDataframeSelection(path, event);
    };
  
    const showTooltip = async (params: any) => {
      const ancestors = params.treeAncestors
        .filter((item: any) => item.name)
        .map((item: any) => item.name);
    
      const div = ui.divH([]);
      let firstAncestor = true;
    
      for (let i = 0; i < ancestors.length; i++) {
        const ancestor = ancestors[i];
        const [isSmiles, isMolBlock, isSmarts] = await Promise.all([
          grok.chem.checkSmiles(ancestor),
          grok.chem.isMolBlock(ancestor),
          grok.chem.isSmarts(ancestor),
        ]);
        const isStructure = isSmiles || isMolBlock || isSmarts;
        if (isStructure) {
          const image = await TreeUtils.getMoleculeImage(ancestor, 150, 100);
          if (image) {
            const { width, height } = image;
            const pixels = image.getContext('2d')!.getImageData(0, 0, width, height).data;
            if (pixels.some((_, i) => i % 4 === 3 && pixels[i] !== 0)) {
              div.appendChild(image);
            }
          }
        } else {
          if (!firstAncestor) {
            div.append(ui.divText('.'));
          }
          div.append(ui.divText(ancestor));
          firstAncestor = false;
        }
      }
      const resultDiv = ui.divV([div, ui.divText(params.value, { style: { fontWeight: 'bold' } })]);
      ui.tooltip.show(resultDiv, params.event.event.x, params.event.event.y);
    };    
  
    const handleZrClick = (params: any) => {
      if (!params.target) return;
  
      //@ts-ignore
      const sortedChildren = params.target.parent._children.sort((a, b) => a.id - b.id);
      const targetIndex = sortedChildren.findIndex((child: { id: number }) => child && child.id === params.target.id);
      const nextElement = sortedChildren[targetIndex - 1];
  
      if (!nextElement) return;
  
      //@ts-ignore
      const seriesModel = this.chart.getModel().getSeriesByIndex(0);
      const itemGraphics = seriesModel.getData()._graphicEls.slice(1);
      const itemIds = seriesModel.getData()._idList.slice(1);
  
      const idx = itemGraphics.findIndex((item: any) => item && nextElement && item.id === nextElement.id);
      let name = itemIds[idx];
      if (name === "") name = " ";
  
      if (name) {
        const seriesData = this.chart.getOption()?.series?.[0]?.data ?? [];
        const targetPath = this.findByPathOccurrenceUsingStrings(seriesData, name);
  
        if (targetPath) {
          handleChartClick(targetPath.split(' | '), params);
          this.paintBranchByPath(targetPath);
        }
      }
    };
  
    this.chart.on('mouseover', showTooltip);
    this.chart.on('mouseout', () => ui.tooltip.hide());
    this.chart.getZr().on('click', handleZrClick);
  }

  findByPathOccurrenceUsingStrings(arr: any[], targetOccurrence: string): string | undefined {
    const pathOccurrences: string[] = [];
  
    const flattenTree = (nodes: any[]) => {
      for (const node of nodes) {
        if (node.path) pathOccurrences.push(node.path);
        if (node.children?.length) flattenTree(node.children);
      }
    };
  
    flattenTree(arr);
  
    const occurrenceRegex = /^(.*?)(?:__ec__(\d+))?$/;
    const match = targetOccurrence.match(occurrenceRegex);
  
    if (!match) return undefined;
  
    let [_, basePath, occurrenceStr] = match;
    if (basePath === "") basePath = " ";
    const targetIndex = occurrenceStr ? parseInt(occurrenceStr, 10) : 1;
    let currentCount = 0;
  
    for (const path of pathOccurrences) {
      if (path.endsWith(basePath)) {
        currentCount++;
        if (currentCount === targetIndex) return path;
      }
    }
    return undefined;
  }
  
  paintBranchByPath(path: string): void {
    const hoverStyle = {
      lineStyle: { color: 'orange' },
      itemStyle: { color: 'orange' },
    };
  
    const updatedData = this.buildSeriesConfig(path, hoverStyle);
    if (updatedData) {
      const updatedOption = { 
        ...this.option, 
        series: [
          {
            ...this.option.series[0], 
            data: [updatedData],
          },
        ],
      };
      this.chart.setOption(updatedOption, false, true);
    }
  }
  
  buildSeriesConfig(path: string, hoverStyle: any): any {
    const originalTree = this.chart!.getOption()?.series?.[0]?.data?.[0];
    if (!originalTree) return undefined;
  
    const cloneTree = echarts.util.clone(originalTree);
  
    const applyHoverStyle = (node: any) => {
      if (node.path !== null) {
        const nodePathParts = node.path.split('|').map((p: string) => p.trim());
        const pathParts = path.split('|').map((p: string) => p.trim());
        const isMatch = nodePathParts.every((part: string, index: number) => pathParts[index] === part) || node.path.includes(path);
  
        if (isMatch) Object.assign(node, hoverStyle);
      }
  
      if (node.children) node.children.forEach(applyHoverStyle);
    };
  
    applyHoverStyle(cloneTree);
    return cloneTree;
  }

  onPropertyChanged(p: DG.Property | null, render: boolean = true) {
    if (!p) return;
  
    const { name } = p;

    switch (name) {
      case 'edgeShape':
        if (p.get(this) === 'polyline') {
          this.getProperty('layout')?.set(this, 'orthogonal');
          this.option.series[0].layout = 'orthogonal';
        }
        this.option.series[0].edgeShape = p.get(this);
        this.chart.clear();
        this.render();
        break;
  
      case 'layout':
        const layout: layoutType = p.get(this);
        this.option.series[0].layout = layout;
  
        if (layout === 'orthogonal') {
          this.option.series[0].label.rotate = this.labelRotate;
        } else {
          delete this.option.series[0].label.rotate;
          const es = this.getProperty('edgeShape');
          if (es?.get(this) !== 'curve') {
            es?.set(this, 'curve');
            this.option.series[0].edgeShape = 'curve';
          }
        }
        this.render();
        return;
  
      case 'table':
        this.updateTable();
        this.onTableAttached();
        this.render();
        break;
  
      case 'onClick':
        this.rowSource = rowSourceMap[this.onClick as onClickOptions] || this.rowSource;
        break;
  
      case 'hierarchyColumnNames':
        this.chart.clear();
        this.render();
        break;
  
      case 'colorColumnName':
      case 'colorAggrType':
        this.applyColorAggr = this.shouldApplyAggregation(this.colorColumnName, this.colorAggrType);
        this.render();
        break;
  
      case 'sizeColumnName':
      case 'sizeAggrType':
        this.applySizeAggr = this.shouldApplyAggregation(this.sizeColumnName, this.sizeAggrType);
        this.render();
        break;
  
      case 'fontSize':
        this.option.series[0].label.fontSize = p.get(this);
        this.option.series[0].leaves.label.fontSize = p.get(this);
        this.render();
        break;
  
      case 'orient':
      case 'labelRotate':
        if (name === 'orient') this.option.series[0].orient = p.get(this);
        this.updateOrient();
        this.render();
        break;  
  
      case 'initialTreeDepth':
        this.option.series[0].initialTreeDepth = p.get(this);
        this.render();

      case 'showCounts':
      case 'includeNulls':
        this.render();
        break;
  
      default:
        super.onPropertyChanged(p, render);
        break;
    }
  }  

  updateOrient() {
    const defaultLabelOptions = (position: string, align: string, rotate: number) => ({
      position,
      verticalAlign: 'middle',
      align,
      rotate,
      fontSize: this.fontSize,
    });
  
    const orientations = {
      LR: {
        label: defaultLabelOptions('left', 'right', this.labelRotate),
        leavesLabel: defaultLabelOptions('right', 'left', this.labelRotate),
      },
      RL: {
        label: defaultLabelOptions('right', 'left', this.labelRotate),
        leavesLabel: defaultLabelOptions('left', 'right', this.labelRotate),
      },
      BT: {
        label: defaultLabelOptions('bottom', 'right', this.labelRotate),
        leavesLabel: defaultLabelOptions('top', 'left', this.labelRotate),
      },
      TB: {
        label: defaultLabelOptions('top', 'right', -this.labelRotate),
        leavesLabel: defaultLabelOptions('bottom', 'left', -this.labelRotate),
      },
    };
  
    const { label, leavesLabel } = orientations[this.orient] || {};
    if (label && leavesLabel) {
      this.option.series[0].label = label;
      this.option.series[0].leaves.label = leavesLabel;
    }
  }    

  shouldApplyAggregation(columnName: string, aggrType: string): boolean {
    const column = this.dataFrame.getCol(columnName);
    return aggregationMap.get(column.type)?.includes(aggrType)!;
  }

  _testColumns() {
    return this.dataFrame.columns.length >= 1;
  }

  onTableAttached() {
    const categoricalColumns = [...this.dataFrame.columns.categorical].sort((col1, col2) =>
      col1.categories.length - col2.categories.length);

    if (categoricalColumns.length < 1)
      return;
    
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

  async renderMoleculeQueued(params: any, width: number, height: number) {
    this.moleculeRenderQueue = this.moleculeRenderQueue.then(() => 
      this.renderMolecule(params, width, height)
    );
    await this.moleculeRenderQueue;
  }

  formatLabel(params: any) {
    if (params.data.semType === 'Molecule') {
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
        
      const labelHeight = Math.floor(chartSize / nodesAtLevel.length);
      const minImageWidth = 70;
      const minImageHeight = 80;

      if (availableSpace >= minImageWidth && labelHeight >= minImageHeight) {
        this.renderMoleculeQueued(params, minImageWidth, minImageHeight);
        return ' ';
      }
      return ' ';
    }
    let labelText = this.showCounts ? `${params.name}: ${params.value}` : `${params.name}`;
    return labelText;
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
      data
    });

    this.option.series[0]['symbolSize'] = this.sizeColumnName && this.applySizeAggr ?
      (value: number, params: {[key: string]: any}) => utils.data.mapToRange(
        params.data.size, this.option.series[0].data[0]['size-meta']['min'],
        this.option.series[0].data[0]['size-meta']['max'], ...this.symbolSizeRange) : this.symbolSize;
    if (this.colorColumnName && this.applyColorAggr)
      this.colorCodeTree(this.option.series[0].data[0]);

    this.option.series[0].label.formatter = (params: any) => this.formatLabel(params);
    this.chart.setOption(this.option, false, true);
  }
}

type layoutType = 'orthogonal' | 'radial';
type orientation = 'LR' | 'RL' | 'TB' | 'BT';
type edgeShape = 'curve' | 'polyline';
type symbolType = 'circle' | 'emptyCircle' | 'rect' | 'roundRect' | 'triangle' | 'diamond' | 'pin' | 'arrow' | 'none';