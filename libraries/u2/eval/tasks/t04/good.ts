import * as DG from 'datagrok-api/dg';
import {ChoiceInput} from '@datagrok-libraries/u2';
import {Platform, columnInput} from '@datagrok-libraries/u2/dg';

export function moleculeColumnInput(table: DG.DataFrame): ChoiceInput {
  const isMolecule = (c: DG.Column) => c.semType === DG.SEMTYPE.MOLECULE;
  const input = columnInput('Molecules', table, {filter: isMolecule});
  input.reset(Platform.listedColumns(table, isMolecule)[0]?.name ?? null);
  return input;
}
