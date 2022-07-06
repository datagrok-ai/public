import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as echarts from 'echarts';


export class Utils {
  static toTree(dataFrame: DG.DataFrame, splitByColumnNames: string[], rowMask: DG.BitSet,
    visitNode: ((arg0: treeDataType) => void) | null = null, aggregations:
      aggregationInfo[] = [], linkSelection: boolean = true): treeDataType {
    const data: treeDataType = {
      name: 'All',
      value: 0,
      path: null,
      children: [],
    };

    const builder = dataFrame
      .groupBy(splitByColumnNames)
      .count()
      .whereRowMask(rowMask);

    for (const aggregation of aggregations) {
      data[aggregation.propertyName] = 0; //TODO: accept custom initial values
      builder.add(aggregation.type, aggregation.columnName, aggregation.propertyName);
    }

    const aggregated = builder.aggregate();

    if (linkSelection) {
      grok.data.linkTables(dataFrame, aggregated, splitByColumnNames,
        splitByColumnNames, [DG.SYNC_TYPE.SELECTION_TO_SELECTION], true);
    }

    const countCol = aggregated.columns.byName('count');
    const columns = aggregated.columns.byNames(splitByColumnNames);
    const propNames = aggregations.map((a) => a.propertyName);
    const aggrColumns = aggregated.columns.byNames(propNames);
    const parentNodes: (treeDataType | null)[] = columns.map((_) => null);

    const selectedPaths: string[] = [];
    const selectedNodeStyle = { color: DG.Color.toRgb(DG.Color.selectedRows) };

    const markSelectedNodes = (node: treeDataType): boolean => {
      if (selectedPaths.includes(node.path!)) {
        node.itemStyle = selectedNodeStyle;
        return true;
      }
      if (node.children && node.children.length > 0) {
        let parentSelected = true;
        for (const child of node.children)
          parentSelected = markSelectedNodes(child) && parentSelected;

        if (parentSelected) {
          node.itemStyle = selectedNodeStyle;
          return true;
        }
      }
      return false;
    };

    function aggregateParentNodes(node: treeDataType): void {
      // grok.functions.call(`${DG.AGG}(node.children.map((child) => child[prop]))`)
      if (!node.children)
        return;
      node.children.forEach(aggregateParentNodes);
      for (const prop of propNames)
        node[prop] = node.children.reduce((sum, child) => sum += child[prop], 0) / node.children.length;
    };

    for (let i = 0; i < aggregated.rowCount; i++) {
      const idx = i === 0 ? 0 : columns.findIndex((col) => col.get(i) !== col.get(i - 1));
      const value = countCol.get(i);
      const aggrValues = aggrColumns.reduce((obj, col) => (obj[col.name] = col.get(i), obj), <{ [key: string]: number }>{});
      if (aggregated.selection.get(i))
        selectedPaths.push(columns.map((col) => col.getString(i)).join(' | '));

      for (let colIdx = idx; colIdx < columns.length; colIdx++) {
        const parentNode = colIdx === 0 ? data : parentNodes[colIdx - 1];
        const name = columns[colIdx].getString(i);
        const node: treeDataType = {
          name: name,
          path: parentNode?.path == null ? name : parentNode.path + ' | ' + name,
          value: 0,
        };
        if (colIdx === columns.length - 1)
          propNames.forEach((prop) => node[prop] = aggrValues[prop]);

        parentNodes[colIdx] = node;

        if (!parentNode!.children)
          parentNode!.children = [];
        parentNode!.children.push(node);
        if (visitNode !== null)
          visitNode(node);
      }

      for (let i = 0; i < parentNodes.length; i++)
        parentNodes[i]!.value += value;
      data.value += value;
    }

    aggregateParentNodes(data);
    console.log(JSON.stringify(data));
    markSelectedNodes(data);

    return data;
  }

  static toForest(dataFrame: DG.DataFrame, splitByColumnNames: string[], rowMask: DG.BitSet) {
    const tree = Utils.toTree(dataFrame, splitByColumnNames, rowMask, (node) => node.value = 10);
    return tree.children;
  }

  static mapRowsToObjects(dataFrame: DG.DataFrame, columnNames: string[],
    objectKeys: string[] | null = null): {[key: string]: any}[] {
    const columns = dataFrame.columns.byNames(columnNames);
    if (objectKeys === null)
      objectKeys = columnNames;

    const result = [];
    for (let i = 0; i < dataFrame.rowCount; i++) {
      const object: {[key: string]: any} = {};
      for (let j = 0; j < columns.length; j++)
        object[objectKeys[j]] = columns[j].get(i);
      result.push(object);
    }
    return result;
  }

  /**
   * @param {String[]} columnNames
   * @param {String} path - pipe-separated values
   */
  static pathToPattern(columnNames: string[], path: string): {[key: string]: string} {
    const values = path.split(' | ');
    const pattern: {[key: string]: string} = {};
    for (let i = 0; i < columnNames.length; i++)
      pattern[columnNames[i]] = values[i];
    return pattern;
  }
}

type treeDataType = { name: string, value: number, path: null | string, children?: treeDataType[], itemStyle?: { color?: string }, [prop: string]: any };
type aggregationInfo = { type: DG.AggregationType, columnName: string, propertyName: string };


export class EChartViewer extends DG.JsViewer {
  chart: echarts.ECharts;
  option: any;

  top?: string;
  left?: string;
  bottom?: string;
  right?: string;
  animationDuration?: number;
  animationDurationUpdate?: number;
  tableName?: string;

  constructor() {
    super();

    //common properties
    this.tableName = this.string('table', null, { fieldName: 'tableName', category: 'Data', editor: 'table' });

    const chartDiv = ui.div([], { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
    this.root.appendChild(chartDiv);
    this.chart = echarts.init(chartDiv);
    this.subs.push(ui.onSizeChanged(chartDiv).subscribe((_) => this.chart.resize()));
  }

  initCommonProperties() {
    this.top = this.string('top', '5px');
    this.left = this.string('left', '5px');
    this.bottom = this.string('bottom', '5px');
    this.right = this.string('right', '5px');

    this.animationDuration = this.int('animationDuration', 500);
    this.animationDurationUpdate = this.int('animationDurationUpdate', 750);
  }

  onTableAttached() {
    this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));
    this.subs.push(DG.debounce(this.dataFrame.onDataChanged, 50).subscribe((_) => this.render()));

    this.render();
  }

  prepareOption() {}

  onPropertyChanged(p: DG.Property | null, render: boolean = true) {
    const properties = p !== null ? [p] : this.props.getProperties();

    for (const p of properties)
      this.option.series[0][p.name] = p.get(this);

    if (render)
      this.chart.setOption(this.option);
  }

  render() {
    this.option.series[0].data = this.getSeriesData();
    this.chart.setOption(this.option);
  }

  getSeriesData() {
    throw new Error('Method not implemented.');
  }
}
