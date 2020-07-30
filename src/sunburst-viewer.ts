import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import { d3sunburst } from './sunburst';
import { Property } from 'datagrok-api/dg';
import { TreeDataBuilder } from './tree-data-builder';

export class Sunburst2Viewer extends DG.JsViewer {
    private level1ColumnName: string;
    private level2ColumnName: string;
    private level3ColumnName: string;
    private valueColumnName: string;

    private containerId = 'sunburst-id-123456';
    private containerDiv?: HTMLDivElement;

    private treeDataBuilder = new TreeDataBuilder(this.containerId);

    constructor() {
        super();

        // properties
        this.level1ColumnName = this.string('level1ColumnName');
        this.level2ColumnName = this.string('level2ColumnName');
        this.level3ColumnName = this.string('level3ColumnName');
        this.valueColumnName = this.string('valueColumnName');
    }

    init() {
        this.containerDiv = ui.div([], 'd4-viewer-host');
        this.containerDiv.setAttribute("id", this.containerId);
        this.root.appendChild(this.containerDiv);
    }

    onTableAttached() {
        this.init();

        // this.level1ColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.TEXT);
        // this.level2ColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.TEXT);
        // this.level3ColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.TEXT);
        // this.valueColumnName = this.dataFrame.columns.bySemType(DG.SEMTYPE.TEXT);

        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()) as any);
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()) as any);

        this.render();
    }

    onPropertyChanged(prop: Property) {
        this.render();
    }

    onSizeChanged(w: number, h: number) {
        // if (w !== 0 && h !== 0)
        //     this.map.invalidateSize();
    }

    detach() {
        this.subs.forEach((sub) => sub.unsubscribe());
    }

    render() {
        this.buildTreeData();

        this.containerDiv!.innerHTML = '';
        d3sunburst(this.containerId, this.treeDataBuilder.getTreeData()!);
    }

    private buildTreeData() {
        const selectedRows = this.dataFrame.filter.getSelectedIndexes();

        const categoryColumns = [];
        this.level1ColumnName && categoryColumns.push(this.dataFrame.getCol(this.level1ColumnName));
        this.level2ColumnName && categoryColumns.push(this.dataFrame.getCol(this.level2ColumnName));
        this.level3ColumnName && categoryColumns.push(this.dataFrame.getCol(this.level3ColumnName));
        const valueColumn = this.valueColumnName && this.dataFrame.getCol(this.valueColumnName);

        this.treeDataBuilder.buildTreeData(categoryColumns, valueColumn, selectedRows);
    }
}
