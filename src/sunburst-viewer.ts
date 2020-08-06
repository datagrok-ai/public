import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import { COLUMN_TYPE, Property } from "datagrok-api/dg";
import { d3sunburst } from './sunburst';
import { TreeDataBuilder } from './tree-data-builder';
import { Column } from 'datagrok-api/src/dataframe';

export class Sunburst2Viewer extends DG.JsViewer {

    private containerId = 'sunburst-id-' + Math.random().toString(36).slice(2);
    private chartDiv!: HTMLDivElement;
    private selectorDiv!: HTMLDivElement;
    private selectors: HTMLSelectElement[] = [];
    private valueSelector!: HTMLSelectElement;

    private treeDataBuilder = new TreeDataBuilder(this.containerId);

    constructor() {
        super();
    }

    init() {
        this.selectorDiv = ui.div([ui.span(["Categories:"] as any)], 'sunburst-selectors-container');
        this.root.appendChild(this.selectorDiv);
        this.addSelector(true);

        const valueContainer = ui.div([ui.span(["Value:"] as any)], 'sunburst-value-container');
        this.root.appendChild(valueContainer);
        this.valueSelector = this.createNumberColumnSelector(true);
        valueContainer.appendChild(this.valueSelector);

        this.chartDiv = ui.div([], 'sunburst-chart-container');
        this.chartDiv.setAttribute("id", this.containerId);
        this.root.appendChild(this.chartDiv);
    }

    addSelector(setDefault = false) {
        const selectorIndex = this.selectors.length;
        const selectorName = 'sunburst-selector-' + selectorIndex;
        const selector = this.createStringColumnSelector(setDefault);
        selector.name = selectorName;
        selector.onchange = (event) => {
            this.render();
            if (selector.selectedIndex != 0 && selectorIndex == this.selectors.length - 1) {
                this.addSelector();
            }
        }
        this.selectors.push(selector);
        selectorIndex && this.selectorDiv.appendChild(ui.span(['>' as any]));
        this.selectorDiv.appendChild(selector);
    }

    private createStringColumnSelector(setDefault: boolean) {
        const columnNames = this.getColumnNames([COLUMN_TYPE.STRING]);
        const defaultColumnName = setDefault && columnNames.length ? columnNames[0] : '';
        return this.createSelector(columnNames, defaultColumnName);
    }

    private createNumberColumnSelector(setDefault: boolean) {
        const columnNames = this.getColumnNames([COLUMN_TYPE.INT, COLUMN_TYPE.FLOAT]);
        const defaultColumnName = setDefault && columnNames.length ? columnNames[0] : '';
        return this.createSelector(columnNames, defaultColumnName);
    }

    private getColumnNames(type: COLUMN_TYPE[]) {
        return this.dataFrame.columns.toList().filter(c => type.some(t => t === c.type)).map(c => c.name);
    }

    createSelector(columnNames: string[], selectedName = ''): HTMLSelectElement {
        const select = document.createElement('select');
        select.className = 'sunburst-selector';
        select.add(this.createSelectOption());
        for (const columnName of columnNames) {
            select.add(this.createSelectOption(columnName, columnName, columnName === selectedName));
        }
        return select;
    }

    createSelectOption(text: string = "", value?: string, selected = false): HTMLOptionElement {
        const option = document.createElement('option');
        option.innerText = text;
        option.value = value || text;
        option.selected = selected;
        return option;
    }

    getSelectedColumnNames(): string[] {
        return this.selectors
            .map(selector => {
                const selectedOptions = selector.selectedOptions;
                if (!selectedOptions.length) {
                    return '';
                }
                return selectedOptions.item(0)!.value!;
            })
            .filter(s => !!s);
    }

    onTableAttached() {
        this.init();

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

        this.chartDiv.innerHTML = '';
        const width = this.root.parentElement!.offsetWidth;
        const height = this.root.parentElement!.offsetHeight;
        const radius = Math.min(width, height) / 2 * 0.9;
        console.error({root: this.root, width, height, radius});

        d3sunburst(this.chartDiv, this.treeDataBuilder.getTreeData()!, radius);
    }

    private buildTreeData() {
        const selectedRows = this.dataFrame.filter.getSelectedIndexes();

        const categoryColumns: Column[] = this.getSelectedColumnNames()
            .map(columnName => this.dataFrame.getCol(columnName));

        this.treeDataBuilder.buildTreeData(categoryColumns, '', selectedRows);
    }
}
