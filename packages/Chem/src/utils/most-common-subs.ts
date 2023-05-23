import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getRdKitModule} from '../package';
import {getMolSafe} from './mol-creation_rdkit';
import {MAX_MCS_ROW_COUNT} from '../constants';


export function getMCS(molecules: DG.Column<string>, exactAtomSearch: boolean, exactBondSearch: boolean): string {
  if (molecules.length > MAX_MCS_ROW_COUNT) {
    grok.shell.warning(`Too many rows, maximum for MCS is ${MAX_MCS_ROW_COUNT}`);
    return '';
  }
  const rdkit = getRdKitModule();
  const mols: RDMol[] = [];

  for (let i = 0; i < molecules.length; i++) {
    const molString = molecules.get(i);
    const molSafe = getMolSafe(molString!, {}, rdkit);
    if (molSafe.mol !== null && !molSafe.isQMol && molSafe.mol.is_valid())
      mols.push(molSafe.mol);
    else
      molSafe.mol?.delete();
  }

  if (mols.length > 0) {
    const arr = new Uint32Array(mols.length);

    for (let i = 0; i < mols.length; i++) {
      //@ts-ignore
      arr[i] = mols[i].$$.ptr;
    }

    //as wasm works with 8 bit and 32 bit are used
    const buff = rdkit._malloc(mols.length*4);

    // >> 2 is the reduction of element number of 32 bit vs 8 bit for offset
    //@ts-ignore
    rdkit.HEAPU32.set(arr, buff >> 2);

    const smarts: string = rdkit.get_mcs(buff, mols.length, exactAtomSearch, exactBondSearch);

    rdkit._free(buff);

    for (let j = 0; j < mols.length; j++)
      mols[j].delete();

    return smarts;
  }

  return '';
}
