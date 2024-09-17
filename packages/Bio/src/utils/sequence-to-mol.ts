import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {getSeqHelper, ToAtomicLevelRes} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {checkInputColumnUI} from './check-input-column';

export async function sequenceToMolfile(
  df: DG.DataFrame, macroMolecule: DG.Column, nonlinear: boolean, monomerLib: IMonomerLib
): Promise<ToAtomicLevelRes> {
  let rdKitModule: RDModule;
  try {
    rdKitModule = await getRdKitModule();
  } catch (ex: any) {
    grok.shell.warning('Transformation to atomic level requires package "Chem" installed.');
    return {mol: null, warnings: []};
  }

  let atomicLevelRes: ToAtomicLevelRes;
  if (nonlinear) {
    const seqHelper = await getSeqHelper();
    const seqSh = SeqHandler.forColumn(macroMolecule);

    let helmCol: DG.Column<string>;
    let seqColName!: string;
    if (seqSh.isHelm())
      helmCol = macroMolecule;
    else {
      seqColName = macroMolecule.name;
      macroMolecule.name = `__${seqColName}`;
      helmCol = seqSh.convert(NOTATION.HELM);
      helmCol.name = seqColName;
      df.columns.add(helmCol, false);
    }
    try {
      atomicLevelRes = await seqHelper.helmToAtomicLevel(helmCol, true, true);
    } finally {
      if (helmCol !== macroMolecule) {
        df.columns.remove(helmCol.name);
        macroMolecule.name = seqColName;
      }
    }
    df.columns.add(atomicLevelRes.mol!.col);
    df.columns.add(atomicLevelRes.mol!.highlightCol);
  } else {
    if (!checkInputColumnUI(macroMolecule, 'To Atomic Level'))
      return {mol: null, warnings: ['Column is not suitable']};

    atomicLevelRes = await _toAtomicLevel(df, macroMolecule, monomerLib);
    if (atomicLevelRes.mol!.col !== null) {
      df.columns.add(atomicLevelRes.mol!.col, true);
    }
  }
  if (atomicLevelRes.mol?.col) await grok.data.detectSemanticTypes(df);
  return atomicLevelRes;
}
