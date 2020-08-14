import * as d3 from 'd3';
import {HierarchyNode} from 'd3';
import * as DG from 'datagrok-api/dg';

export type TreeData = HierarchyNode<Branch>

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
}

export class TreeDataBuilder {
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

        return d3
            .hierarchy(root, x => Array.from(x.children.values()))
            .sum(x => x.leafIds.length);
    }
}
