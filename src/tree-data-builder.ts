import * as d3 from 'd3';
import { HierarchyNode } from 'd3';
import { Column } from 'datagrok-api/dg';

export type TreeData = HierarchyNode<TreeItem>

export interface TreeItem {
    parent: string;
    id: string;
    value: number;
}

export class TreeDataBuilder {
    private treeData!: TreeData;

    constructor(private rootId: string) {
        this.clearTreeData();
    }

    getTreeData() {
        return this.treeData;
    }

    buildTreeData(
        categoryColumns: Column[],
        valueColumn: "" | Column,
        selectedRows: Int32Array
    ) {
        this.clearTreeData();

        if (!categoryColumns.length) {
            return;
        }

        const table = new Map();
        table.set(this.getKey('', this.rootId), 0);
        const [lastColumn] = categoryColumns.slice(-1);
        for (const index of selectedRows) {
            let prevPrevColumnValue = null;
            let prevColumnValue = this.rootId;
            for (const column of categoryColumns) {
                let columnValue = column.get(index);
                let isLast = column === lastColumn;
                if (columnValue === '') {
                    if (prevPrevColumnValue == null) {
                        break;
                    }
                    isLast = true;
                    columnValue = prevColumnValue;
                    prevColumnValue = prevPrevColumnValue;
                }
                const key = this.getKey(prevColumnValue, columnValue);
                const value = table.get(key) || 0;
                const valueIncrement = isLast ? (valueColumn ? valueColumn.get(index) : 1) : 0;
                table.set(key, value + valueIncrement);
                if (isLast) {
                    break;
                }
                prevPrevColumnValue = prevColumnValue;
                prevColumnValue = columnValue;
            }
        }

        const list: TreeItem[] = [];
        table.forEach((value, k) => list.push({ ...JSON.parse(k), value }));
        // console.error(list);

        this.treeData = this.treeItemsToTreeData(list);
    }

    private clearTreeData() {
        this.treeData = this.treeItemsToTreeData([{ id: this.rootId, parent: '', value: 0 }]);
    }

    private getKey(parent: string, id: string) {
        return JSON.stringify({ parent, id })
    }

    private treeItemsToTreeData(list: TreeItem[]) {
        return d3.stratify<TreeItem>()
            .id(d => d.id)
            .parentId(d => d.parent)
            (list)
            .sum(d => d.value);
    }
}
