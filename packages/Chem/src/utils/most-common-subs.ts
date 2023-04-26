import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {RDMol} from "@datagrok-libraries/chem-meta/src/rdkit-api";
import {getRdKitModule} from "../package";
import {getMolSafe} from './mol-creation_rdkit';


export function getMCS(molecules: string, df: DG.DataFrame, compareElements: boolean, compareBonds: boolean): string {
    let rdkit = getRdKitModule();

    let molCol = df.columns.byName(molecules);
    let n = molCol.length;

    let mols: RDMol[] = [];

    for(let i = 0; i < molCol.length; i++) {
      let molSafe = getMolSafe(molCol.get(i), {}, rdkit);
      if(molSafe.mol !== null && !molSafe.isQMol)
        mols.push(molSafe.mol);
    }

    let arr = new Uint32Array(mols.length);

    for(let i = 0; i < mols.length; i++) {
      //@ts-ignore
      arr[i] = mols[i].$$.ptr;
    }

    //@ts-ignore
    let buff = rdkit.asm.Zb(mols.length*4);

    // >> 2 is the reduction of element number of 32 bit vs 8 bit for offset
    //@ts-ignore
    rdkit.HEAPU32.set(arr, buff >> 2);

    let smarts: string = rdkit.get_mcs(buff, mols.length, compareElements, compareBonds);

    //@ts-ignore
    //rdkit.asm.gb("get_mcs", "void", ["number", "number"], [2, buff]);

    //@ts-ignore
    rdkit.asm.$b(buff);

    for(let i = 0; i < molCol.length; i++)
      mols[i].delete();

    return smarts;
  }
