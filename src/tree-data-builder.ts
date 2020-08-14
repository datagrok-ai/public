import * as d3 from 'd3';
import { HierarchyNode } from 'd3';
import * as DG from 'datagrok-api/dg';

export type TreeData = HierarchyNode<TreeItem>
export type AltTreeData = HierarchyNode<Branch>

export interface TreeItem {
    parent: string;
    id: string;
    value: number;
}

export type ColumnValue = any;
export type RowId = number;
export class Branch {
    readonly value: ColumnValue;
    readonly children: Map<ColumnValue, Branch>;
    readonly leafIds: RowId[] = [];
    
    constructor(value: ColumnValue) {
        this.value = value;
        this.children = new Map();
    }
    
    toArray(): (ColumnValue | ColumnValue[])[] {
        const childrenArray = Array.from(this.children.values()).map(x => x.toArray());
        return [ this.value, childrenArray ];
    }
}

export class AlternativeTreeDataBuilder {
    buildTreeData(
        categoryColumns: DG.Column[],
        valueColumn: "" | DG.Column,
        selectedRows: Int32Array
    ) {
        const root = new Branch("root");
        
        const lastColumn = categoryColumns[categoryColumns.length - 1]; 
        
        for(const rowId of selectedRows) {
            let currentBranch = root;
            for (const column of categoryColumns) {
                const value = column.get(rowId);
                if (value === null || value === undefined) {
                    currentBranch.leafIds.push(rowId);
                    break;
                }

                const existingBranch = currentBranch.children.get(value);
                if (existingBranch) {
                    currentBranch = existingBranch;
                }
                else {
                    const parent = currentBranch;
                    currentBranch = new Branch(value);
                    parent.children.set(value, currentBranch);
                }
                
                if (lastColumn === column) {
                    currentBranch.leafIds.push(rowId);
                }
            }
        }
        
        console.log(root);
        
        const h1 = d3
            .hierarchy(root, x => Array.from(x.children.values()))
            .sum(x => x.leafIds.length);
        
        console.log(h1);
        
        return h1;
    }
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
        categoryColumns: DG.Column[],
        valueColumn: "" | DG.Column,
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
