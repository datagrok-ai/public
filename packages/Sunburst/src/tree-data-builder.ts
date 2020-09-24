import * as d3 from 'd3';
import {HierarchyNode} from 'd3';
import { BitSet, Column, Row, RowList } from 'datagrok-api/dg';

export type TreeData = HierarchyNode<Branch>

export type ColumnCategory = any;
export type ColumnValue = number;
export type RowId = number;

export class Branch {
    readonly category: ColumnCategory;
    valueSum: ColumnValue = 0;
    valueAvg: ColumnValue = 0;
    readonly children: Map<ColumnCategory, Branch>;
    readonly leafRowIds: RowId[] = [];
    selectedRowsNumber: number = 0;
    branchRowsNumber: number = 0;

    constructor(category: ColumnCategory) {
        this.category = category;
        this.children = new Map();
    }

    public traverseTree(mapFn: (branch: Branch) => void, currentBranch?: Branch): void {
        if (currentBranch == null) {
            this.traverseTree(mapFn, this);
            return;
        }
        if (currentBranch.children && currentBranch.children.size) {
            for (var branch of currentBranch.children.values()) {
                this.traverseTree(mapFn, branch);
            }
        }
        mapFn(currentBranch);
    }
}

export class TreeDataBuilder {
    buildTreeData(
        categoryColumns: Column[],
        valueColumn: Column | undefined,
        rowNumber: number,
        selection: BitSet
    ) {
        const root = new Branch("root");

        for (let rowId = 0; rowId < rowNumber; rowId++) {
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
            }
            currentBranch.leafRowIds.push(rowId);
            currentBranch.branchRowsNumber++;
            if (selection.get(rowId)) {
                currentBranch.selectedRowsNumber++;
            }
            if (valueColumn) {
                currentBranch.valueSum += parseFloat(valueColumn.get(rowId));
            }
        }

        // Recalculate selectedChildrenNumber and branchChildrenNumber for all the nodes
        root.traverseTree(branch => {
            if (!branch || !branch.children || !branch.children.size) {
                branch.valueAvg = this.getAvgValue(branch);
                return;
            }
            var selectedChildrenSum = 0;
            var branchChildrenSum = 0;
            var valueSum = 0;
            for (const childBranch of branch.children.values()) {
                selectedChildrenSum += childBranch.selectedRowsNumber;
                branchChildrenSum += childBranch.branchRowsNumber;
                valueSum += childBranch.valueSum;
            }
            branch.branchRowsNumber = branchChildrenSum;
            branch.selectedRowsNumber = selectedChildrenSum;
            branch.valueSum = valueSum;
            branch.valueAvg = this.getAvgValue(branch);
        });

        return d3
            .hierarchy(root, x => Array.from(x.children.values()))
            .sum(x => x.leafRowIds.length);
    }

    private getAvgValue(branch: Branch): number {
        return branch.branchRowsNumber
            ? (branch.valueSum / branch.branchRowsNumber)
            : 0;
    }
}
