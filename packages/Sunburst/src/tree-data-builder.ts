import * as d3 from 'd3';
import {HierarchyNode} from 'd3';
import * as DG from 'datagrok-api/dg';

export type TreeData = HierarchyNode<Branch>
export type AggregationType = 'count' | 'sum' | 'avg';

export type ColumnCategory = any;
export type ColumnValue = number;
export type RowId = number;
export class Branch {
    readonly category: ColumnCategory;
    value: ColumnValue = 0;
    readonly children: Map<ColumnCategory, Branch>;
    readonly leafIds: RowId[] = [];

    constructor(category: ColumnCategory) {
        this.category = category;
        this.children = new Map();
    }
}

type Aggregator = (x: Branch) => number;

export class TreeDataBuilder {
    buildTreeData(
        categoryColumns: DG.Column[],
        valueColumn: DG.Column | undefined,
        selectedRows: Int32Array,
        aggregationType: AggregationType
    ) {
        const root = new Branch("root");

        const lastCategoryColumn = categoryColumns[categoryColumns.length - 1];

        for(const rowId of selectedRows) {
            let currentBranch = root;
            for (const categoryColumn of categoryColumns) {
                const category = categoryColumn.get(rowId);
                if (category === null || category === undefined) {
                    break;
                }

                const existingBranch = currentBranch.children.get(category);
                if (existingBranch) {
                    currentBranch = existingBranch;
                }
                else {
                    const parent = currentBranch;
                    currentBranch = new Branch(category);
                    parent.children.set(category, currentBranch);
                }

                if (lastCategoryColumn === categoryColumn) {
                }
            }
            currentBranch.leafIds.push(rowId);
            if (valueColumn) {
                currentBranch.value = valueColumn.get(rowId);
            }
        }

        return d3
            .hierarchy(root, x => Array.from(x.children.values()))
            .sum(this.getAggregator(aggregationType));
    }

    private getAggregator(aggregationType: AggregationType): Aggregator {
        switch (aggregationType) {
            case 'sum':
                return x => x.value;
            case 'avg':
                return x => x.leafIds.length
                    ? (x.value / x.leafIds.length)
                    : 0;
            case 'count':
            default:
        }
        return x => x.leafIds.length;
    }
}
