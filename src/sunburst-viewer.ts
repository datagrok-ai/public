import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import {SunburstRenderer} from './sunburst';
import {TreeDataBuilder} from './tree-data-builder';

export class SunburstViewer extends DG.JsViewer {

    private containerId = 'sunburst-id-' + Math.random().toString(36).slice(2);
    private chartDiv!: HTMLDivElement;
    private selectorDiv!: HTMLDivElement;
    private selectors: HTMLSelectElement[] = [];
    //private valueSelector!: HTMLSelectElement;

    private colors?: string[];

    constructor() {
        super();
    }

    init() {
        this.selectorDiv = ui.div([ui.span(["Categories:"] as any)], 'sunburst-selectors-container');
        this.root.appendChild(this.selectorDiv);
        this.addSelector(true);
        this.addSelector(false);

        //const valueContainer = ui.div([ui.span(["Value:"] as any)], 'sunburst-value-container');
        //this.root.appendChild(valueContainer);
        //this.valueSelector = this.createNumberColumnSelector(true);
        //valueContainer.appendChild(this.valueSelector);

        this.chartDiv = ui.div([], 'sunburst-chart-container');
        this.chartDiv.setAttribute("id", this.containerId);
        this.root.appendChild(this.chartDiv);
    }

    onTableAttached() {
        this.init();

        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));

        this.render();
    }

    detach() {
        this.subs.forEach((sub) => sub.unsubscribe());
    }

    render() {
        const data = this.buildTreeData();

        this.chartDiv.innerHTML = '';
        const width = this.root.parentElement!.offsetWidth;
        const height = this.root.parentElement!.offsetHeight;
        const radius = Math.min(width, height) / 2 * 0.9;

        const renderer = new SunburstRenderer(radius, this.getColors(), this.clickHandler);
        renderer.render(this.chartDiv, data);
    }

    private buildTreeData() {
        const selectedRows = this.dataFrame.filter.getSelectedIndexes();

        const categoryColumns: DG.Column[] = this.getSelectedColumnNames()
            .map(columnName => this.dataFrame.getCol(columnName));


        const alt = new TreeDataBuilder();
        return alt.buildTreeData(categoryColumns, '', selectedRows);
    }

    private clickHandler = (selectedRowIds: number[]) => {
        const selection = this.dataFrame.selection;
        selection.setAll(false);
        for (const rowId of selectedRowIds) {
            selection.set(rowId, true);
        }
    }

    getColors(): string[] {
        if (this.colors) {
            return this.colors;
        }
        this.colors = DG.Color.categoricalPalette.map(DG.Color.toRgb);
        return this.colors;
    }

    // UI COLUMN SELECTION STUFF

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
        const columnNames = this.getColumnNames([DG.COLUMN_TYPE.STRING]);
        const defaultColumnName = setDefault && columnNames.length ? columnNames[0] : '';
        return this.createSelector(columnNames, defaultColumnName);
    }

    private createNumberColumnSelector(setDefault: boolean) {
        const columnNames = this.getColumnNames([DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT]);
        const defaultColumnName = setDefault && columnNames.length ? columnNames[0] : '';
        return this.createSelector(columnNames, defaultColumnName);
    }

    private getColumnNames(type: DG.COLUMN_TYPE[]) {
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
}
