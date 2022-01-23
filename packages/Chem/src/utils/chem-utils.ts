import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getRdKitModule} from "../package";

export function isMolBlock(s: string) {
  return s.includes('M  END');
}

export function _convertMolecule(molecule: string, from: string, to: string): string {
  let mol;
  try {
    mol = getRdKitModule().get_mol(molecule);
    if (to === 'molblock')
      return mol.get_molblock();
    if (to === 'smiles')
      return mol.get_smiles();
    if (to === 'v3Kmolblock')
      return mol.get_v3Kmolblock();
    if (to == 'inchi')
      return mol.get_inchi();
    throw `Failed to convert molecule: unknown target unit: "${to}"`;
  } finally {
    mol?.delete();
  }
}