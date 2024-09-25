import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IMonomerLibBase} from '@datagrok-libraries/bio/src/types';

import {SeqHelper} from '../seq-helper';

import {_package, getMonomerLib} from '../../package';

/** Translate HELM column into molfile column and append to the dataframe */
export async function getMolColumnFromHelm(
  df: DG.DataFrame, helmCol: DG.Column<string>, chiralityEngine: boolean = true, monomerLib: IMonomerLibBase
): Promise<DG.Column<string>> {
  const seqHelper = await SeqHelper.getInstance();
  const converter = seqHelper.getHelmToMolfileConverter(monomerLib);
  const molCol = converter.convertToRdKitBeautifiedMolfileColumn(helmCol, chiralityEngine, _package.rdKitModule, monomerLib);
  molCol.semType = DG.SEMTYPE.MOLECULE;
  return molCol;
}

export async function getSmilesColumnFromHelm(
  helmCol: DG.Column<string>
): Promise<DG.Column<string>> {
  const seqHelper = await SeqHelper.getInstance();
  const monomerLib = getMonomerLib();
  const converter = seqHelper.getHelmToMolfileConverter(monomerLib);
  const smilesCol = converter.convertToSmiles(helmCol);
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  return smilesCol;
}
