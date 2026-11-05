import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';

import {getUncommonAtomsAndBonds} from '../../../utils/chem-common';
import {byEngineModule} from '../../../engine/chem-engine';


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

  // "Chemical Engine": in Crux mode each pair on Crux, on `module` (RDKit's) where Crux fails (engine/README.md)
  for (let i = 0; i < from.length; i++) {
    [fromAligned[i], toAligned[i], res1[i], res2[i]] =
      byEngineModule('mmpAlign', module, (m, crux) => alignPair(cores[i], from[i], to[i], m, crux))!;
  }

  return {inverse1: res1, inverse2: res2, fromAligned: fromAligned, toAligned: toAligned};
}

/** A pair aligned on its core, and each molecule's atoms and bonds outside the core; '' and nulls where it fails (under
 * Crux it throws, so that RDKit aligns the pair). */
function alignPair(core: string, from: string, to: string, module: RDModule, crux: boolean):
  [string, string, ISubstruct | null, ISubstruct | null] {
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
    mcsMol = module.get_mol(core.replace('[*:1]', '[H]'));
    mol1 = module.get_mol(from);
    mol2 = module.get_mol(to);
    mcsMol.set_new_coords();
    mol1.generate_aligned_coords(mcsMol, opts);
    mol2.generate_aligned_coords(mcsMol, opts);
    const inverse1 = getUncommonAtomsAndBonds(from, mcsMol, module, '#bc131f');
    //@ts-ignore
    inverse1['highlightBondWidthMultiplier'] = 40;
    const inverse2 = getUncommonAtomsAndBonds(to, mcsMol, module, '#49bead');
    //@ts-ignore
    inverse2['highlightBondWidthMultiplier'] = 40;
    return [mol1.get_molblock(), mol2.get_molblock(), inverse1, inverse2];
  } catch (e: any) {
    if (crux)
      throw e; // RDKit aligns the pair (a trap is noted)
    return ['', '', null, null];
  } finally {
    mol1?.delete();
    mol2?.delete();
    mcsMol?.delete();
  }
}
