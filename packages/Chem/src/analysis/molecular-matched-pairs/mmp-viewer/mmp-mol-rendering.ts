import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';

import {getRdKitService, getUncommonAtomsAndBonds} from '../../../utils/chem-common-rdkit';

async function getMmpMcs(molecules1: string[], molecules2: string[]): Promise<string[]> {
  const service = await getRdKitService();
  const molecules: [string, string][] = new Array<[string, string]>(molecules1.length);
  for (let i = 0; i < molecules1.length; i++)
    molecules[i] = [molecules1[i], molecules2[i]];

  return await service.mmpGetMcs(molecules);
}

export async function getInverseSubstructuresAndAlign(from: string[], to: string[], module: RDModule):
  Promise<{
    inverse1: (ISubstruct | null)[],
    inverse2: (ISubstruct | null)[],
    fromAligned: string[],
    toAligned: string[]}> {
  const fromAligned = new Array<string>(from.length);
  const toAligned = new Array<string>(from.length);
  const res1 = new Array<(ISubstruct | null)>(from.length);
  const res2 = new Array<(ISubstruct | null)>(from.length);

  const mcs = await getMmpMcs(from, to);

  for (let i = 0; i < from.length; i++) {
    //aligning molecules
    let mcsMol = null;
    let mol1 = null;
    let mol2 = null;
    const opts = JSON.stringify({
      useCoordGen: true,
      allowRGroups: true,
      acceptFailure: false,
      alignOnly: true,
    });

    try {
      mcsMol = module.get_qmol(mcs[i]);
      mol1 = module.get_mol(from[i]);
      mol2 = module.get_mol(to[i]);
      mcsMol.set_new_coords();
      mol1.generate_aligned_coords(mcsMol, opts);
      mol2.generate_aligned_coords(mcsMol, opts);
      fromAligned[i] = mol1.get_molblock();
      toAligned[i] = mol2.get_molblock();
      res1[i] = getUncommonAtomsAndBonds(from[i], mcsMol, module, '#bc131f');
      res2[i] = getUncommonAtomsAndBonds(to[i], mcsMol, module, '#49bead');
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
