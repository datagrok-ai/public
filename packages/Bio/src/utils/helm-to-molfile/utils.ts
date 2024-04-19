/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {HelmToMolfileConverter} from './converter';

/** Translate HELM column into molfile column and append to the dataframe */
export async function helm2mol(df: DG.DataFrame, helmCol: DG.Column<string>): Promise<void> {
  const molCol = await getMolColumnFromHelm(df, helmCol);
  df.columns.add(molCol, true);
  await grok.data.detectSemanticTypes(df);
}


/** Translate HELM column into molfile column and append to the dataframe */
export async function getMolColumnFromHelm(
  df: DG.DataFrame, helmCol: DG.Column<string>, chiralityEngine?: boolean
): Promise<DG.Column<string>> {
  const converter = new HelmToMolfileConverter(helmCol, df);
  const molCol = await converter.convertToRdKitBeautifiedMolfileColumn(chiralityEngine);
  molCol.semType = DG.SEMTYPE.MOLECULE;
  return molCol;
}

export async function getSmilesColumnFromHelm(
  df: DG.DataFrame, helmCol: DG.Column<string>
): Promise<DG.Column<string>> {
  const converter = new HelmToMolfileConverter(helmCol, df);
  const smilesCol = await converter.convertToSmiles();
  smilesCol.semType = DG.SEMTYPE.MOLECULE;
  return smilesCol;
}
