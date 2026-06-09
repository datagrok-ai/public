import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {MAX_MCS_ROW_COUNT} from '../constants';
import * as chemCommonRdKit from './chem-common-rdkit';

export async function getMCS(
  molecules: DG.Column<string>, exactAtomSearch: boolean, exactBondSearch: boolean,
): Promise<string> {
  if (molecules.length > MAX_MCS_ROW_COUNT) {
    grok.shell.warning(`Too many rows, maximum for MCS is ${MAX_MCS_ROW_COUNT}`);
    return '';
  }


  const rdkitService = await chemCommonRdKit.getRdKitService();
  const mcs = await rdkitService.getMCS(molecules.toList(), exactAtomSearch, exactBondSearch);
  return mcs ?? '';
}
