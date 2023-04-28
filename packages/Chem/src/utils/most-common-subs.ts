import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {RDMol} from "@datagrok-libraries/chem-meta/src/rdkit-api";
import {getRdKitModule} from "../package";
import {getMolSafe} from './mol-creation_rdkit';


export function getMCS(molecules: string, df: DG.DataFrame, compareElements: boolean, compareBonds: boolean): string {
  let rdkit = getRdKitModule();

  let molCol = df.columns.byName(molecules);
  const nmolecules = molCol.length;
  const batchSize = 100;  // batches by 100 molecules due to memory limits
  const nbatches = molCol.length/batchSize;
  let mcs = "";
  
  for(let i = 0; i < nbatches; i++) {
    let mols: RDMol[] = [];
    const currentBatchSize = i < nbatches - 1 ? batchSize : nmolecules - i * batchSize;

    for(let j = 0; j < currentBatchSize; j++) {
      const add = batchSize * i;

      let molString = molCol.get(add + j);

      let molSafe = getMolSafe(molString, {}, rdkit);
      if(molSafe.mol !== null && !molSafe.isQMol && molSafe.mol.is_valid())
        mols.push(molSafe.mol);
      else
        molSafe.mol?.delete();
    }

    let n = mols.length;
    if(n > 0) {
      n = mcs == "" ? n : n + 1;

      let arr = new Uint32Array(n);

      for(let j = 0; j < mols.length; j++) {
        //@ts-ignore
        arr[j] = mols[j].$$.ptr;
      }

      let mcsMol = null;
      if(mcs != "") {
        mcsMol = rdkit.get_qmol(mcs);
        //@ts-ignore
        arr[n - 1] = mcsMol.$$.ptr;
      }

      let buff = rdkit._malloc(n*4);

      // >> 2 is the reduction of element number of 32 bit vs 8 bit for offset
      //@ts-ignore
      rdkit.HEAPU32.set(arr, buff >> 2);
      //rdkit.writeArrayToMemory(arr, buff);

      let smarts: string = rdkit.get_mcs(buff, n, compareElements, compareBonds);
      mcs = smarts;

      rdkit._free(buff);

      for(let j = 0; j < mols.length; j++)
        mols[j].delete();

      if(mcs != "")
        mcsMol?.delete();
    }
  }

  return mcs;
}
