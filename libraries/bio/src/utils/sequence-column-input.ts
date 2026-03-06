import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export interface ISequenceColumnInput {
    get root(): HTMLElement;
    get value(): DG.Column | null;
    set value(col: DG.Column | null);
    get inputBase(): DG.InputBase<DG.Column | null> | undefined | null;
}

class SequenceColumnInputBase implements ISequenceColumnInput {
    private readonly _inputBase: DG.InputBase<DG.Column | null>;
    constructor(inp: DG.InputBase<DG.Column | null>) {
        this._inputBase = inp;
    }
    get root(): HTMLElement { return this._inputBase.root; }
    get value(): DG.Column | null { return this._inputBase.value; }
    set value(col: DG.Column | null) { this._inputBase.value = col; }
    get inputBase(): DG.InputBase<DG.Column | null> { return this._inputBase; }
}

export async function createSequenceColumnInput(name: string, options?: ui.input.IColumnInputInitOptions<DG.Column>): Promise<ISequenceColumnInput> {
    const func = DG.Func.find({package: 'Bio', name: 'sequenceColumnInput'})[0];
    if (!func)
        return new SequenceColumnInputBase(ui.input.column(name, options));
    return await func.apply({name: name, options: options}) as ISequenceColumnInput;
}