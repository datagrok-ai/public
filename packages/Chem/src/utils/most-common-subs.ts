import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {MAX_MCS_ROW_COUNT} from '../constants';
import {getRdKitService} from './chem-common-rdkit';


export async function getMCS(molecules: DG.Column<string>, exactAtomSearch: boolean, exactBondSearch: boolean): Promise<string> {
  if (molecules.length > MAX_MCS_ROW_COUNT) {
    grok.shell.warning(`Too many rows, maximum for MCS is ${MAX_MCS_ROW_COUNT}`);
    return '';
  }
  return await (await getRdKitService()).getMCS(molecules.toList(), exactAtomSearch, exactBondSearch);
}

