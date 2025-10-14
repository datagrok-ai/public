import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';

import {getUncommonAtomsAndBonds} from '../../../utils/chem-common';


export async function getInverseSubstructuresAndAlign(cores: string[],
  from: string[], to: string[], module: RDModule):
  Promise<{
    inverse1: (ISubstruct | null)[],
    inverse2: (ISubstruct | null)[],
    fromAligned: string[],
    toAligned: string[]}> {
  const fromAligned = new Array<string>(from.length);
  const toAligned = new Array<string>(from.length);
  const res1 = new Array<(ISubstruct | null)>(from.length);
  const res2 = new Array<(ISubstruct | null)>(from.length);


  for (let i = 0; i < from.length; i++) {
    //aligning molecules
    let mol1 = null;
    let mol2 = null;
    let mcsMol = null;

    const opts = JSON.stringify({
      useCoordGen: true,
      allowRGroups: true,
      acceptFailure: false,
      alignOnly: true,
    });

    try {
      const core = cores[i].replace('[*:1]', '[H]');
      mcsMol = module.get_mol(core);
      mol1 = module.get_mol(from[i]);
      mol2 = module.get_mol(to[i]);
      mcsMol.set_new_coords();
      mol1.generate_aligned_coords(mcsMol, opts);
      mol2.generate_aligned_coords(mcsMol, opts);
      fromAligned[i] = mol1.get_molblock();
      toAligned[i] = mol2.get_molblock();
      res1[i] = getUncommonAtomsAndBonds(from[i], mcsMol, module, '#bc131f');
      //@ts-ignore
      res1[i]['highlightBondWidthMultiplier'] = 40;
      res2[i] = getUncommonAtomsAndBonds(to[i], mcsMol, module, '#49bead');
      //@ts-ignore
      res2[i]['highlightBondWidthMultiplier'] = 40;
    } catch (e: any) {
      fromAligned[i] = '';
      toAligned[i] = '';
    } finally {
      mol1?.delete();
      mol2?.delete();
      mcsMol?.delete();
    }
  }

  return {inverse1: res1, inverse2: res2, fromAligned: fromAligned, toAligned: toAligned};
}
