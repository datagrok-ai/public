/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {helm2mol} from './helm-to-molfile/utils';
import {checkInputColumnUI} from './check-input-column';

export async function sequenceToMolfile(
  df: DG.DataFrame, macroMolecule: DG.Column, nonlinear: boolean, monomerLib: IMonomerLib
): Promise<void> {
  if (DG.Func.find({package: 'Chem', name: 'getRdKitModule'}).length === 0) {
    grok.shell.warning('Transformation to atomic level requires package "Chem" installed.');
    return;
  }
  if (nonlinear) {
    const seqSh = SeqHandler.forColumn(macroMolecule);
    if (!seqSh.isHelm())
      macroMolecule = seqSh.convert(NOTATION.HELM);
    return helm2mol(df, macroMolecule);
  }
  if (!checkInputColumnUI(macroMolecule, 'To Atomic Level'))
    return;
  const atomicLevelRes = await _toAtomicLevel(df, macroMolecule, monomerLib);
  if (atomicLevelRes.col !== null) {
    df.columns.add(atomicLevelRes.col, true);
    await grok.data.detectSemanticTypes(df);
  }

  if (atomicLevelRes.warnings && atomicLevelRes.warnings.length > 0)
    grok.shell.warning(ui.list(atomicLevelRes.warnings));
}
