import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {MAX_MCS_ROW_COUNT} from '../constants';
import * as chemCommonRdKit from './chem-common-rdkit';
import {chemBeginCriticalSection, chemEndCriticalSection} from './chem-common';

export async function getMCS(
  molecules: DG.Column<string>, exactAtomSearch: boolean, exactBondSearch: boolean,
): Promise<string> {
  if (molecules.length > MAX_MCS_ROW_COUNT) {
    grok.shell.warning(`Too many rows, maximum for MCS is ${MAX_MCS_ROW_COUNT}`);
    return '';
  }

  // MCS runs a single blocking call on worker 0 — the same worker R-Group uses. Take the section so the
  // two serialize instead of sharing the worker, so cancelling R-Group can't collaterally kill this call.
  await chemBeginCriticalSection();
  try {
    const rdkitService = await chemCommonRdKit.getRdKitService();
    const mcs = await rdkitService.getMCS(molecules.toList(), exactAtomSearch, exactBondSearch);
    return mcs ?? '';
  } finally {
    chemEndCriticalSection();
  }
}
