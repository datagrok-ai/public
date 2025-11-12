/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ISubstruct, ISubstructProvider} from '@datagrok-libraries/chem-meta/src/types';
import {MMP_NAMES} from './mmp-constants';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getUncommonAtomsAndBonds} from '../../../utils/chem-common';

export class MMPSubstructProvider implements ISubstructProvider {
  idxsToSmiles: string[];
  df: DG.DataFrame;
  colName: string;
  rdkit: RDModule;
  color: string;

  constructor(idxsToSmiles: string[], df: DG.DataFrame, colName: string, rdKit: RDModule, color: string) {
    this.idxsToSmiles = idxsToSmiles;
    this.df = df;
    this.colName = colName;
    this.rdkit = rdKit;
    this.color = color;
  };

  getSubstruct(tableRowIdx: number | null): ISubstruct | undefined {
    if (tableRowIdx == null || tableRowIdx === -1)
      return;

    let mcsMol;
    try {
      const molecule = this.df.get(this.colName, tableRowIdx);
      const core = this.idxsToSmiles[this.df.get(MMP_NAMES.CORE_NUM, tableRowIdx)].replace('[*:1]', '[H]');
      mcsMol = this.rdkit.get_mol(core);

      const res = getUncommonAtomsAndBonds(molecule, mcsMol, this.rdkit, this.color);
      if (!res)
        return;
        //@ts-ignore
      res['highlightBondWidthMultiplier'] = 40;
        res!.alignByScaffold = mcsMol.get_molblock();
        return res;
    } finally {
      mcsMol?.delete();
    }
  }
}
