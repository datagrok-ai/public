/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getMacroMol} from './atomic-works';
import {getMolfilesFromSeq} from './monomer-utils';

export async function _toAtomicLevel(
  df: DG.DataFrame,
  macroMolecule: DG.Column,
  monomersLibObject: any[]
): Promise<void> {
  if (DG.Func.find({package: 'Chem', name: 'getRdKitModule'}).length === 0) {
    grok.shell.warning('Transformation to atomic level requires package "Chem" installed.');
    return;
  }
  const atomicCodes = getMolfilesFromSeq(macroMolecule, monomersLibObject);
  const result = await getMacroMol(atomicCodes!);
  const col = DG.Column.fromStrings('regenerated', result);
  col.semType = DG.SEMTYPE.MOLECULE;
  col.tags[DG.TAGS.UNITS] = 'molblock';
  df.columns.add(col, true);
  await grok.data.detectSemanticTypes(df);
}
