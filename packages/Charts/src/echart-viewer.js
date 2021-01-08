import * as echarts from "echarts";

export class Utils {

  /** @param {DataFrame} dataFrame
   * @param {String[]} splitByColumnNames
   * @returns {Object} */
  static toHierarchy(dataFrame, splitByColumnNames) {
    let data = {
      name: 'All',
      children: []
    };

    let aggregated = dataFrame.groupBy(splitByColumnNames).aggregate();
    let columns = aggregated.columns.byNames(splitByColumnNames);
    let parentNodes = columns.map(_ => null);

    for (let i = 0; i < aggregated.rowCount; i++) {
      let idx = i === 0 ? 0 : columns.findIndex((col) => col.get(i) !== col.get(i - 1));

      for (let colIdx = idx; colIdx < columns.length; colIdx++) {
        let parentNode = colIdx === 0 ? data : parentNodes[colIdx - 1];
        let node = { name: columns[colIdx].getString(i) };
        parentNodes[colIdx] = node;
        if (!parentNode.children)
          parentNode.children = [];
        parentNode.children.push(node);
      }
    }

    return data;
  }
}

export class EChartViewer extends DG.JsViewer {
  constructor() {
    super();
    let chartDiv = ui.div(null, { style: { position: 'absolute', left: '0', right: '0', top: '0', bottom: '0'}} );
    this.root.appendChild(chartDiv);
    this.myChart = echarts.init(chartDiv);
    this.subs.push(ui.onSizeChanged(chartDiv).subscribe((_) => this.myChart.resize()));
  }
}