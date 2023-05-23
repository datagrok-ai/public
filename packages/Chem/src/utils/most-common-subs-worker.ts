import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getMolSafe} from './mol-creation_rdkit';
// @ts-ignore
import initRDKitModule from '../RDKit_minimal.js';
//@ts-ignore
import rdKitLibVersion from '../rdkit_lib_version';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

onmessage = async (event) => {
  const {molecules, exactAtomSearch, exactBondSearch} = event.data;
  const data: { error?: any; smarts?: string } = {};
  try {
    // Unfortunately, we cannot use the getRDKitModule() function here, because it does not work in the worker
    // context. So we have to use the initRDKitModule() function instead. it takes about 200MS
    const rdkit: RDModule = await initRDKitModule({locateFile: () => `../dist/${rdKitLibVersion}.wasm`});
    const mols: RDMol[] = [];

    for (let i = 0; i < molecules.length; i++) {
      const molString = molecules[i];
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
      const buff = rdkit._malloc(mols.length * 4);

      // >> 2 is the reduction of element number of 32 bit vs 8 bit for offset
      //@ts-ignore
      rdkit.HEAPU32.set(arr, buff >> 2);

      const smarts: string = rdkit.get_mcs(buff, mols.length, exactAtomSearch, exactBondSearch);

      rdkit._free(buff);

      for (let j = 0; j < mols.length; j++)
        mols[j].delete();
      data.smarts = smarts;
    }
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};
