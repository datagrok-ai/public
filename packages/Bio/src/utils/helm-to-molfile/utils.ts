import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {SeqHelper} from '../seq-helper';


/** Translate HELM column into molfile column and append to the dataframe */
export async function getMolColumnFromHelm(
  df: DG.DataFrame, helmCol: DG.Column<string>, chiralityEngine: boolean = true
): Promise<DG.Column<string>> {
  const seqHelper = await SeqHelper.getInstance();
  const converter = seqHelper.getHelmToMolfileConverter(df, helmCol);
  const molCol = converter.convertToRdKitBeautifiedMolfileColumn(chiralityEngine, _package.rdKitModule);
  molCol.semType = DG.SEMTYPE.MOLECULE;
  return molCol;
}

export async function getSmilesColumnFromHelm(
  df: DG.DataFrame, helmCol: DG.Column<string>
): Promise<DG.Column<string>> {
  const seqHelper = await SeqHelper.getInstance();
  const converter = seqHelper.getHelmToMolfileConverter(df, helmCol);
  const smilesCol = converter.convertToSmiles(_package.rdKitModule);
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  return smilesCol;
}
