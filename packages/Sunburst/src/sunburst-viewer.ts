import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SunburstRenderer} from './sunburst-renderer';
import { Branch, TreeData, TreeDataBuilder } from './tree-data-builder';
import { HierarchyNode } from 'd3-hierarchy';
import { BitSet } from 'datagrok-api/src/dataframe';

interface BitSetEx extends BitSet {
    set(i: number, x: boolean, notify?: boolean): BitSet;
    fireChanged(): void;
}

export class SunburstViewer extends DG.JsViewer {

    private chartDiv!: HTMLDivElement;
    private selectorDiv!: HTMLDivElement;
    private selectors: HTMLSelectElement[] = [];

    private renderer: SunburstRenderer;

    private colors?: string[];
    private stringColumnNames?: string[];
    private initialized = false;

    private treeData?: TreeData;

    private valueColumnName = this.string('valueColumnName');
    private categoryColumnsSerialized = this.string('categoryColumnsSerialized');

    constructor() {
        super();
        // FIXME: colors will be dynamical in the future
        this.renderer = new SunburstRenderer(this.getColors(), this.clickHandler);
    }

    init(): void {
        const label = ui.span(['Categories:'] as any);
        label.className = 'label';
        this.selectorDiv = ui.div([label], 'sunburst-selectors-container');
        this.root.appendChild(this.selectorDiv);
        this.addSelector(true);
        this.addSelector(false);

        this.chartDiv = ui.div([], 'sunburst-chart-container');
        this.root.appendChild(this.chartDiv);
        setTimeout(() => {
            this.initialized = true;
        }, 0);
    }

    onTableAttached(): void {
        this.init();

        this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe((_) => this.render()));
        this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe((_) => this.render()));

        this.loadSelectedColumnNames();
        this.render();
    }

    detach(): void {
        this.treeData = undefined;
        this.subs.forEach((sub) => sub.unsubscribe());
    }

    onPropertyChanged(property: DG.Property): void {
        super.onPropertyChanged(property);
        if (this.initialized && (!property || property.name !== 'categoryColumnsSerialized')) {
            this.render();
        }
    }

    onSizeChanged(width: number, height: number): void {
        this.render(false);
    }

    private render(dataChanged = true): void {
        if (dataChanged) {
            const selectedColumnNames = this.getSelectedColumnNames();
            this.saveSelectedColumnNames(selectedColumnNames);
            this.treeData = this.buildTreeData(selectedColumnNames);
        }

        const selectorsHeight = this.selectorDiv.getBoundingClientRect().height;
        const height = this.root.parentElement!.offsetHeight - selectorsHeight;
        const width = this.root.parentElement!.offsetWidth;

        if (this.treeData) {
            this.renderer.render(this.chartDiv, this.treeData, width, height);
        }
    }

    private buildTreeData(selectedColumnNames: string[]): HierarchyNode<Branch> {
        const selection = this.dataFrame.filter;
        const rowCount = this.dataFrame.rowCount;

        const categoryColumns: DG.Column[] = selectedColumnNames
            .map(columnName => this.dataFrame.getCol(columnName));
        const valueColumn = this.valueColumnName ? this.dataFrame.getCol(this.valueColumnName) : undefined;

        const alt = new TreeDataBuilder();
        return alt.buildTreeData(categoryColumns, valueColumn, rowCount, selection);
    }

    private clickHandler = (selectedRowIds: number[]) => {
        const selection = this.dataFrame.selection as BitSetEx;
        selection.setAll(false);
        for (const rowId of selectedRowIds) {
            selection.set(rowId, true, false);
        }
        selection.fireChanged();
    }

    private getColors(): string[] {
        if (this.colors) {
            return this.colors;
        }
        this.colors = DG.Color.categoricalPalette.map(DG.Color.toRgb);
        return this.colors;
    }

    // UI COLUMN SELECTION STUFF

    private addSelector(setDefault = false): void {
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
        const span = ui.span([selector as any]);
        span.className = 'selector';
        this.selectorDiv.appendChild(span);
    }

    private createStringColumnSelector(setDefault: boolean): HTMLSelectElement {
        const columnNames = this.getStringColumnNames();
        const defaultColumnName = setDefault && columnNames.length ? columnNames[0] : '';
        return this.createSelector(columnNames, defaultColumnName);
    }

    private getStringColumnNames(): string[] {
        if (this.stringColumnNames) {
            return this.stringColumnNames;
        }
        this.stringColumnNames = this.getColumnNames([DG.COLUMN_TYPE.STRING]);
        return this.stringColumnNames;
    }

    private getColumnNames(type: DG.COLUMN_TYPE[]): string[] {
        return this.dataFrame.columns.toList().filter(c => type.some(t => t === c.type)).map(c => c.name);
    }

    private createSelector(columnNames: string[], selectedName = ''): HTMLSelectElement {
        const select = document.createElement('select');
        select.className = 'sunburst-selector';
        select.add(this.createSelectOption());
        for (const columnName of columnNames) {
            select.add(this.createSelectOption(columnName, columnName, columnName === selectedName));
        }
        return select;
    }

    private createSelectOption(text: string = '', value?: string, selected = false): HTMLOptionElement {
        const option = document.createElement('option');
        option.innerText = text;
        option.value = value || text;
        option.selected = selected;
        return option;
    }

    private getSelectedColumnNames(): string[] {
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

    private saveSelectedColumnNames(selectedColumnNames: string[]): void {
        (this.getProperty('categoryColumnsSerialized') as any).set('categoryColumnsSerialized', JSON.stringify(selectedColumnNames));
    }

    private loadSelectedColumnNames(): void {
        const selectedColumnNames = this.deserializeSelectedColumnNames();
        if (selectedColumnNames == null || !selectedColumnNames.length) {
            return;
        }
        const columnNames = this.getStringColumnNames();
        let selector: HTMLSelectElement | undefined = undefined;
        for (let i = 0; i < selectedColumnNames.length; i++) {
            selector = this.selectors[i];
            if (!selector) {
                this.addSelector();
                selector = this.selectors[i];
            }
            const selectedColumnName = selectedColumnNames[i];
            const index = columnNames.indexOf(selectedColumnName);
            if (index === -1) {
                selector.value = '';
            } else {
                selector.value = selectedColumnName;
            }
        }
        if (selector) {
            // Add the empty selector if needed
            setTimeout(() => {
                selector!.dispatchEvent(new Event('change'));
            }, 0) ;
        }
    }

    private deserializeSelectedColumnNames(): string[] | undefined {
        const categoryColumnsSerialized = this.getProperty('categoryColumnsSerialized').get(this);
        try {
            const selectedColumnNames = JSON.parse(categoryColumnsSerialized);
            if (Array.isArray(selectedColumnNames)) {
                return selectedColumnNames
            }
        } catch (e) {
        }
    }
}
