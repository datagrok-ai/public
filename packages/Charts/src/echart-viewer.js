import * as echarts from "echarts";
import {wu} from "wu";

export class Utils {

  /** @param {DataFrame} dataFrame
   * @param {String[]} splitByColumnNames
   * @param {DG.BitSet} rowMask
   * @param {Function} visitNode
   * @returns {Object} */
  static toTree(dataFrame, splitByColumnNames, rowMask, visitNode = null) {
    let data = {
      name: 'All',
      value: 0,
      path: null,
      children: []
    };

    let aggregated = dataFrame
      .groupBy(splitByColumnNames)
      .count()
      .whereRowMask(rowMask)
      .aggregate();

    let countCol = aggregated.columns.byName('count');
    let columns = aggregated.columns.byNames(splitByColumnNames);
    let parentNodes = columns.map(_ => null);

    for (let i = 0; i < aggregated.rowCount; i++) {
      let idx = i === 0 ? 0 : columns.findIndex((col) => col.get(i) !== col.get(i - 1));
      let value = countCol.get(i);

      for (let colIdx = idx; colIdx < columns.length; colIdx++) {
        let parentNode = colIdx === 0 ? data : parentNodes[colIdx - 1];
        let name = columns[colIdx].getString(i);
        let node = {
          name: name,
          path: parentNode.path == null ? name : parentNode.path + ' | ' + name,
          value: 0
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
    let tree = Utils.toTree(dataFrame, splitByColumnNames, rowMask, (node) => node.value = 10);
    return tree.children;
  }

  /**
   * @param {DataFrame} dataFrame
   * @param {String[]} columnNames
   * @param {String[]} objectKeys
   * @returns {Object[]}
   */
  static mapRowsToObjects(dataFrame, columnNames, objectKeys = null) {
    let columns = dataFrame.columns.byNames(columnNames);
    if (objectKeys === null)
      objectKeys = columnNames;

    let result = [];
    for (let i = 0; i < dataFrame.rowCount; i++) {
      let object = {};
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
    let values = path.split(' | ');
    let pattern = {};
    for (let i = 0; i < columnNames.length; i++)
      pattern[columnNames[i]] = values[i];
    return pattern;
  }
}


export class EChartViewer extends DG.JsViewer {
  constructor() {
    super();
    let chartDiv = ui.div(null, { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
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
    let properties = p !== null ? [p] : this.props.getProperties();

    for (let p of properties)
      this.option.series[0][p.name] = p.get(this);

    if (render)
      this.chart.setOption(this.option);
  }

  render() {
    this.option.series[0].data = this.getSeriesData();
    this.chart.setOption(this.option);
  }
}