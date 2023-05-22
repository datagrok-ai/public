import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import Fuse from 'fuse.js';
import { COLUMN_TYPE, InputBase } from 'datagrok-api/dg';

export class FuzzyFilter extends DG.Filter {
    inputElement: InputBase<string>;
    get input(): string {
        return this.inputElement.value;
    }
    
    bitset: DG.BitSet | undefined;
    lastColumnVersion: number | undefined;
    lastCol: DG.Column | null = null;
    lastInput: String | undefined;

    async applyFilter() {
        if (this.lastCol == this.column && this.lastColumnVersion == this.column?.version 
            && this.lastInput == this.input && this.bitset != undefined) {
                this.dataFrame?.filter?.and(this.bitset);
                return;
            }

        let fuse = new Fuse(this.column?.toList() ?? []);
        let filterResult = fuse.search(this.input)
        this.bitset = DG.BitSet.create(this.dataFrame?.rowCount ?? 0, (i: number) => false)
        for (let res of filterResult)
            this.bitset.set(res.refIndex, true, false)

        this.lastCol = this.column;
        this.lastColumnVersion = this.column?.version;
        this.lastInput = this.input;
        this.dataFrame?.filter?.and(this.bitset);
    }

    get filterSummary(): string {
        return this.input;
    }

    constructor() {
        super();
        this.inputElement = ui.textInput('input', '', () => {
            this.dataFrame?.rows?.requestFilter();
        });
        this.root = ui.divV([this.inputElement]);
    }

    attach(dataFrame: DG.DataFrame): void {
        super.attach(dataFrame);
        this.column ??= dataFrame.columns.byName('Questions');
        this.column ??= dataFrame.columns.toList().find((el) => el.type == COLUMN_TYPE.STRING) ?? null
    }

    saveState(): any {
        const state = super.saveState();
        state.type = 'PowerPack:fuzzyFilter';
        state.input = this.input;
        return state;
      }
    
    applyState(state: any): void {
        super.applyState(state);
        if (state.input) {
        this.inputElement.value = state.input
        }

        if (state.columnName) {
        this.column = this.dataFrame?.columns.byName(state.columnName) ?? null
        } 
    }
}