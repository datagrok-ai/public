import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import Fuse from 'fuse.js';
import { COLUMN_TYPE, InputBase } from 'datagrok-api/dg';
import { isEmpty } from 'rxjs/operators';
import { timer } from 'rxjs';

export class FuzzyFilter extends DG.Filter {
    inputElement: InputBase<string>;
    get input(): string {
        return this.inputElement.value;
    }
    
    get caption() {
        return this.column?.name ?? 'fuzzy filter'
    }
    bitset: DG.BitSet | undefined;
    lastColumnVersion: number | undefined;
    lastCol: DG.Column | null = null;
    lastInput: String | undefined;

    async applyFilter() {
        if (this.input == null || this.input == undefined || this.input == '')
            return

        if (this.lastCol == this.column && this.lastColumnVersion == this.column?.version 
            && this.lastInput == this.input && this.bitset != undefined) {
                this.dataFrame?.filter?.and(this.bitset);
                return;
            }

        let fuse = new Fuse(this.column?.toList() ?? [], 
            {includeScore: true, isCaseSensitive: false, ignoreLocation: true,
                 useExtendedSearch: true, includeMatches: true});

        let filterResult = fuse.search(this.input)
        this.bitset = DG.BitSet.create(this.dataFrame?.rowCount ?? 0, (i: number) => false)
        for (let res of filterResult)
            if (res.score != undefined && res.score < 0.5)
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
            setTimeout(() => this.dataFrame?.rows?.requestFilter(), 500)
        });
        this.inputElement.input.setAttribute('placeholder', 'search')
        this.root = ui.divV([this.inputElement.input]);
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