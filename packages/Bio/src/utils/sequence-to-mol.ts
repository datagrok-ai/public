import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getSeqHelper, ToAtomicLevelRes} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ChemTemps} from '@datagrok-libraries/chem-meta/src/consts';
import {buildMonomerHoverLink} from '@datagrok-libraries/bio/src/monomer-works/monomer-hover';

import {checkInputColumnUI} from './check-input-column';
import {getMolColName, getMolHighlightColName} from '@datagrok-libraries/bio/src/monomer-works/utils';

export async function sequenceToMolfile(
  df: DG.DataFrame, macroMolecule: DG.Column, nonlinear: boolean, highlight: boolean,
  monomerLib: IMonomerLib, rdKitModule: RDModule
): Promise<ToAtomicLevelRes> {
  let res: ToAtomicLevelRes;
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
      res = await seqHelper.helmToAtomicLevel(helmCol, true, true);
    } finally {
      if (helmCol !== macroMolecule) {
        df.columns.remove(helmCol.name);
        macroMolecule.name = seqColName;
      }
    }
  } else { // linear
    if (!checkInputColumnUI(macroMolecule, 'To Atomic Level'))
      return {mol: null, warnings: ['Column is not suitable']};

    res = await _toAtomicLevel(df, macroMolecule, monomerLib, rdKitModule);
  }


  if (res.mol) {
    const molColName = getMolColName(df, macroMolecule.name);
    const molHlColName = getMolHighlightColName(df, molColName);
    res.mol.col.name = molColName;
    df.columns.add(res.mol!.col, true);

    if (highlight) {
      res.mol.highlightCol.name = molHlColName;
      df.columns.add(res.mol.highlightCol);
      res.mol!.col.temp[ChemTemps.SUBSTRUCT_COL] = molHlColName;
    }
    buildMonomerHoverLink(macroMolecule, res.mol!.col, monomerLib, rdKitModule);
    await grok.data.detectSemanticTypes(df);
  }
  return res;
}
