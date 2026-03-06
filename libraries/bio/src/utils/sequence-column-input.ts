import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

export interface ISequenceColumnInput {
    get root(): HTMLElement;
    get value(): DG.Column | null;
    set value(col: DG.Column | null);
    get inputBase(): DG.InputBase<DG.Column | null> | undefined | null;
}

export async function createSequenceColumnInput(name: string, options: ui.input.IColumnInputInitOptions<DG.Column>): Promise<ISequenceColumnInput> {
    const func = DG.Func.find({package: 'Bio', name: 'sequenceColumnInput'})[0];
    if (!func)
        return ui.input.column(name, options) as unknown as ISequenceColumnInput;
    return await func.apply({name: name, options: options}) as ISequenceColumnInput;
}