import * as echarts from 'echarts';


export class Utils {
  /** @param {DataFrame} dataFrame
   * @param {String[]} splitByColumnNames
   * @param {DG.BitSet} rowMask
   * @param {Function} visitNode
   * @return {Object} */
  static toTree(dataFrame, splitByColumnNames, rowMask, visitNode = null) {
    const data = {
      name: 'All',
      value: 0,
      path: null,
      children: [],
    };

    const aggregated = dataFrame
      .groupBy(splitByColumnNames)
      .count()
      .whereRowMask(rowMask)
      .aggregate();

    const countCol = aggregated.columns.byName('count');
    const columns = aggregated.columns.byNames(splitByColumnNames);
    const parentNodes = columns.map((_) => null);

    for (let i = 0; i < aggregated.rowCount; i++) {
      const idx = i === 0 ? 0 : columns.findIndex((col) => col.get(i) !== col.get(i - 1));
      const value = countCol.get(i);

      for (let colIdx = idx; colIdx < columns.length; colIdx++) {
        const parentNode = colIdx === 0 ? data : parentNodes[colIdx - 1];
        const name = columns[colIdx].getString(i);
        const node = {
          name: name,
          path: parentNode.path == null ? name : parentNode.path + ' | ' + name,
          value: 0,
        };
        parentNodes[colIdx] = node;

        if (!parentNode.children)
          parentNode.children = [];
        parentNode.children.push(node);
        if (visitNode !== null)
          visitNode(node);
      }

      for (let i = 0; i < parentNodes.length; i++)
        parentNodes[i].value += value;
      data.value += value;
    }

    console.log(JSON.stringify(data));

    return data;
  }

  static toForest(dataFrame, splitByColumnNames, rowMask) {
    const tree = Utils.toTree(dataFrame, splitByColumnNames, rowMask, (node) => node.value = 10);
    return tree.children;
  }

  /**
   * @param {DataFrame} dataFrame
   * @param {String[]} columnNames
   * @param {String[]} objectKeys
   * @return {Object[]}
   */
  static mapRowsToObjects(dataFrame, columnNames, objectKeys = null) {
    const columns = dataFrame.columns.byNames(columnNames);
    if (objectKeys === null)
      objectKeys = columnNames;

    const result = [];
    for (let i = 0; i < dataFrame.rowCount; i++) {
      const object = {};
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
  static pathToPattern(columnNames, path) {
    const values = path.split(' | ');
    const pattern = {};
    for (let i = 0; i < columnNames.length; i++)
      pattern[columnNames[i]] = values[i];
    return pattern;
  }
}


export class EChartViewer extends DG.JsViewer {
  constructor() {
    super();
    const chartDiv = ui.div(null, { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
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

  onPropertyChanged(p, render = true) {
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
}
