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
  edgeShape: edgeShape;
  symbol: symbolType;
  symbolSize: number;
  hierarchyColumnNames: string[];
  initialTreeDepth: number;
  sizeColumnName: string;
  sizeAggrType: DG.AggregationType;
  symbolSizeRange: [number, number] = [5, 20];
  aggregations: string[] = Object.values(DG.AGG).filter((f) => f !== DG.AGG.KEY && f !== DG.AGG.PIVOT);
  aggregationsStr: string[] = Object.values({...DG.STR_AGG, ...DG.STAT_COUNTS});
  colorColumnName: string;
  colorAggrType: DG.AggregationType;
  selectionColor = DG.Color.toRgb(DG.Color.selectedRows);
  applySizeAggr: boolean = false;
  applyColorAggr: boolean = false;
  fontSize: number;
  showCounts: boolean;
  onClick: onClickOptions;
  includeNulls: boolean;

  constructor() {
    super();

    this.layout = <layoutType> this.string('layout', 'orthogonal', { choices: ['orthogonal', 'radial'] });
    this.orient = <orientation> this.string('orient', 'LR', { choices: ['LR', 'RL', 'TB', 'BT'] });
    this.expandAndCollapse = this.bool('expandAndCollapse', true);
    this.initialTreeDepth = this.int('initialTreeDepth', 3, { min: -1, max: 5 });
    this.edgeShape = <edgeShape> this.string('edgeShape', 'curve', { choices: ['curve', 'polyline'] });
    this.symbol = <symbolType> this.string('symbol', 'emptyCircle', { choices: [
      'circle', 'emptyCircle', 'rect', 'roundRect', 'triangle', 'diamond', 'pin', 'arrow', 'none',
    ] });
    this.symbolSize = this.int('symbolSize', 7);
    this.sizeColumnName = this.string('sizeColumnName');
    this.sizeAggrType = <DG.AggregationType> this.string('sizeAggrType', DG.AGG.AVG, { choices: this.aggregations });
    this.colorColumnName = this.string('colorColumnName');
    this.colorAggrType = <DG.AggregationType> this.string('colorAggrType', DG.AGG.AVG, { choices: this.aggregations });
    this.hierarchyColumnNames = this.addProperty('hierarchyColumnNames', DG.TYPE.COLUMN_LIST);
    this.fontSize = this.int('fontSize', 12);
    this.showCounts = this.bool('showCounts', false);
    this.onClick = <onClickOptions> this.string('onClick', 'Select', { choices: ['Select', 'Filter']});
    this.includeNulls = this.bool('includeNulls', true);

    this.option = {
      animation: false,
      silent: false,
      series: [
        {
          type: 'tree',

          label: {
            position: 'left',
            verticalAlign: 'middle',
            align: 'right',
            fontSize: this.fontSize
          },

          leaves: {
            label: {
              position: 'right',
              verticalAlign: 'middle',
              align: 'left',
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
    
    this.chart.on('mouseover', (params: any) => {
      const ancestors = params.treeAncestors.filter((item: any) => item.name).map((item: any) => item.name).join('.'); 
      const div = ui.divV([
        ui.divText(ancestors), 
        ui.divText(params.value, { style: { fontWeight: 'bold' } })
      ]);
      ui.tooltip.show(div, params.event.event.x, params.event.event.y);
    });
    this.chart.on('mouseout', () => ui.tooltip.hide());
    this.chart.on('click', handleChartClick);
  }

  onPropertyChanged(p: DG.Property | null, render: boolean = true) {
    if (p?.name === 'edgeShape') {
      this.getProperty('layout')?.set(this, 'orthogonal');
      this.option.series[0].layout = 'orthogonal';
      this.option.series[0].label.rotate = 0;
      this.chart.clear();
    }
    if (p?.name === 'layout') {
      const layout: layoutType = p.get(this);
      this.option.series[0].layout = layout;
      if (layout === 'orthogonal')
        this.option.series[0].label.rotate = 0;
      else {
        this.option.series[0].label.rotate = null;
        const es = this.getProperty('edgeShape');
        if (es?.get(this) !== 'curve') {
          es?.set(this, 'curve');
          this.option.series[0].edgeShape = 'curve';
        }
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
        p?.name === 'fontSize' || p?.name === 'showCounts' || p?.name === 'includeNulls') {
      if (p?.name === 'hierarchyColumnNames')
        this.chart.clear();
      if (p?.name === 'colorColumnName' || p?.name === 'colorAggrType')
        this.applyColorAggr = this.shouldApplyAggregation(this.colorColumnName, this.colorAggrType);
      if (p?.name === 'sizeColumnName' || p?.name === 'sizeAggrType')
        this.applySizeAggr = this.shouldApplyAggregation(this.sizeColumnName, this.sizeAggrType);
      if (p?.name === 'fontSize')
        this.option.series[0].label.fontSize = p.get(this);
      if (p?.name === 'showCounts')
        this.option.series[0].label.formatter = p.get(this) ? '{b}: {c}' : '{b}';
      this.render();
    } else
      super.onPropertyChanged(p, render);
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
      const nodeColor = DG.Color.toRgb(DG.Color.scaleColor(data.color, min, max));
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
    // need to add heuristic to render only in case there is enough place for this
    if (params.data.semType === 'Molecule') {
      const minImageWidth = 70;
      const minImageHeight = 80;
      this.renderMolecule(params, minImageWidth, minImageHeight);
      return ' ';
    }
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
      label: { formatter: (params: any) => this.formatLabel(params) },
    });

    this.option.series[0]['symbolSize'] = this.sizeColumnName && this.applySizeAggr ?
      (value: number, params: {[key: string]: any}) => utils.data.mapToRange(
        params.data.size, this.option.series[0].data[0]['size-meta']['min'],
        this.option.series[0].data[0]['size-meta']['max'], ...this.symbolSizeRange) : this.symbolSize;
    if (this.colorColumnName && this.applyColorAggr)
      this.colorCodeTree(this.option.series[0].data[0]);

    this.chart.setOption(this.option);
  }
}

type layoutType = 'orthogonal' | 'radial';
type orientation = 'LR' | 'RL' | 'TB' | 'BT';
type edgeShape = 'curve' | 'polyline';
type symbolType = 'circle' | 'emptyCircle' | 'rect' | 'roundRect' | 'triangle' | 'diamond' | 'pin' | 'arrow' | 'none';
