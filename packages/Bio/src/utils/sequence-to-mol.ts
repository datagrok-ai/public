/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {helm2mol} from './helm-to-molfile';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {checkInputColumnUI} from './check-input-column';
import {getMonomerLibHelper} from '../package';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

export async function sequenceToMolfile(df: DG.DataFrame, macroMolecule: DG.Column, nonlinear: boolean): Promise<void> {
  if (DG.Func.find({package: 'Chem', name: 'getRdKitModule'}).length === 0) {
    grok.shell.warning('Transformation to atomic level requires package "Chem" installed.');
    return;
  }
  if (nonlinear) {
    const seqUh = UnitsHandler.getOrCreate(macroMolecule);
    if (!seqUh.isHelm())
      macroMolecule = seqUh.convert(NOTATION.HELM);
    helm2mol(df, macroMolecule);
    return;
  }
  if (!checkInputColumnUI(macroMolecule, 'To Atomic Level'))
    return;
  const monomerLib: IMonomerLib = getMonomerLibHelper().getBioLib();
  const atomicLevelRes = await _toAtomicLevel(df, macroMolecule, monomerLib);
  if (atomicLevelRes.col !== null) {
    df.columns.add(atomicLevelRes.col, true);
    await grok.data.detectSemanticTypes(df);
  }

  if (atomicLevelRes.warnings && atomicLevelRes.warnings.length > 0)
    grok.shell.warning(ui.list(atomicLevelRes.warnings));
}
