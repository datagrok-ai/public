import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import * as utils from './utils';


export class TreeUtils {
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
      data[aggregation.propertyName] = utils.data.aggToStat(dataFrame, aggregation.columnName, aggregation.type);
      data[`${aggregation.propertyName}-meta`] = {};
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

    function aggregateParentNodes(): void {
      const paths: {[key: string]: {[key: string]: number}} = {};
      for (let i = 1; i < columns.length; i++) {
        const builder = dataFrame
          .groupBy(splitByColumnNames.slice(0, -i))
          .whereRowMask(rowMask);
        for (const aggregation of aggregations)
          builder.add(aggregation.type, aggregation.columnName, aggregation.propertyName);
        const df = builder.aggregate();
        const rowCount = df.rowCount;
        for (let i = 0; i < rowCount; i++) {
          let path = '';
          const props: {[key: string]: number} = {};
          for (const column of df.columns) {
            if (propNames.includes(column.name))
              props[column.name] = column.get(i);
            else
              path = (path ? path + ' | ' : '') + column.getString(i);
          }
          paths[path] = props;
        }
      }

      function updatePropMeta(node: treeDataType) {
        for (const prop of propNames) {
          if (!node.path) {
            data[`${prop}-meta`] = { min: Infinity, max: -Infinity };
            continue;
          }
          node[prop] = node[prop] ?? paths[node.path][prop];
          if (!data[`${prop}-meta`])
            continue;
          data[`${prop}-meta`].min = Math.min(data[`${prop}-meta`].min, node[prop]);
          data[`${prop}-meta`].max = Math.max(data[`${prop}-meta`].max, node[prop]);
        }
        node.children?.forEach(updatePropMeta);
      }

      updatePropMeta(data);
    };

    for (let i = 0; i < aggregated.rowCount; i++) {
      const idx = i === 0 ? 0 : columns.findIndex((col) => col.get(i) !== col.get(i - 1));
      const value = countCol.get(i);
      const aggrValues = aggrColumns.reduce((obj, col) =>
        (obj[col.name] = col.get(i), obj), <{ [key: string]: number }>{});
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

    if (aggregations.length > 0)
      aggregateParentNodes();

    markSelectedNodes(data);

    return data;
  }

  static toForest(dataFrame: DG.DataFrame, splitByColumnNames: string[], rowMask: DG.BitSet) {
    const tree = TreeUtils.toTree(dataFrame, splitByColumnNames, rowMask, (node) => node.value = 10);
    return tree.children;
  }

  static mapRowsToObjects(dataFrame: DG.DataFrame, columnNames: string[],
    objectKeys: string[] | null = null): {[key: string]: any}[] {
    const columns = dataFrame.columns.byNames(columnNames);
    if (objectKeys === null)
      objectKeys = columnNames;

    const result = [];
    const rowIndexes = dataFrame.filter.getSelectedIndexes();
    for (let i = 0; i < rowIndexes.length; i++) {
      const object: {[key: string]: any} = {};
      for (let j = 0; j < columns.length; j++)
        object[objectKeys[j]] = columns[j].get(rowIndexes[i]);
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

export type treeDataType = { name: string, value: number, path: null | string, children?: treeDataType[],
  itemStyle?: { color?: string }, [prop: string]: any };
export type aggregationInfo = { type: DG.AggregationType, columnName: string, propertyName: string };
