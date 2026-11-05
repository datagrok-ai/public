import * as DG from 'datagrok-api/dg';
import {ChoiceInput} from '@datagrok-libraries/u2';
import {columnInput} from '@datagrok-libraries/u2/dg';

export function moleculeColumnInput(table: DG.DataFrame): ChoiceInput {
  const input = columnInput('Molecules', table,
    {filter: (c) => c.semType === DG.SEMTYPE.MOLECULE && !c.name.startsWith('~')});
  input.reset(table.columns.bySemType(DG.SEMTYPE.MOLECULE)?.name ?? null);
  return input;
}
