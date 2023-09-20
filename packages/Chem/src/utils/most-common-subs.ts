import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from '../package';
import {getMolSafe} from './mol-creation_rdkit';
import {MAX_MCS_ROW_COUNT} from '../constants';


export function getMCS(molecules: DG.Column<string>, exactAtomSearch: boolean, exactBondSearch: boolean): string|null {
  if (molecules.length > MAX_MCS_ROW_COUNT) {
    grok.shell.warning(`Too many rows, maximum for MCS is ${MAX_MCS_ROW_COUNT}`);
    return null;
  }
  const rdkit = getRdKitModule();
  let mols;
  try {
    mols = new rdkit.MolList();

    for (let i = 0; i < molecules.length; i++) {
      const molString = molecules.get(i);
      if (!molString)
        continue;
      let molSafe;
      try {
        molSafe = getMolSafe(molString!, {}, rdkit);
        if (molSafe.mol !== null && !molSafe.isQMol)
          mols.append(molSafe.mol);
      } finally {
        molSafe?.mol?.delete();
      }
    }

    let mcsSmarts: string|null = null;
    if (mols.size() > 1) {
      mcsSmarts = rdkit.get_mcs_as_smarts(mols, JSON.stringify({
        AtomCompare: exactAtomSearch ? 'Elements' : 'Any',
        BondCompare: exactBondSearch ? 'OrderExact' : 'Order',
      }));
    }

    return mcsSmarts;
  } finally {
    mols?.delete();
  }
}
