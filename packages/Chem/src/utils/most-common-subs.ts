import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {MAX_MCS_ROW_COUNT} from '../constants';
import {mcsByEngine} from '../engine/analysis-jobs';

/** FMCS's SMARTS of the molecules ('' for none), by the "Chemical engine" setting (engine/analysis-jobs.ts). */
export async function getMCS(
  molecules: DG.Column<string>, exactAtomSearch: boolean, exactBondSearch: boolean,
): Promise<string> {
  if (molecules.length > MAX_MCS_ROW_COUNT) {
    grok.shell.warning(`Too many rows, maximum for MCS is ${MAX_MCS_ROW_COUNT}`);
    return '';
  }


  const mcs = await mcsByEngine(molecules.toList(), exactAtomSearch, exactBondSearch);
  return mcs ?? '';
}
