import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export enum MolNotation {
  Smiles = 'smiles',
  MolBlock = 'molblock',
  V3KMolBlock = 'v3Kmolblock',
  Inchi = 'inchi',
  Unknown = 'unknown',
}

export function isMolBlock(s: string) {
  return s.includes('M  END');
}

export function convertMoleculeImpl(molecule: string, from: MolNotation, to: MolNotation, module: RDModule): string {
  let mol;
  try {
    mol = module.get_mol(molecule);
    if (to === MolNotation.MolBlock)
      return mol.get_molblock();
    if (to === MolNotation.Smiles)
      return mol.get_smiles();
    if (to === MolNotation.V3KMolBlock)
      return mol.get_v3Kmolblock();
    if (to == MolNotation.Inchi)
      return mol.get_inchi();
    throw `Failed to convert molecule: unknown target unit: "${to}"`;
  } finally {
    mol?.delete();
  }
}

export function molToMolblock(molStr: string, module: RDModule): string {
  return isMolBlock(molStr)
    ? molStr
    : convertMoleculeImpl(molStr, MolNotation.Unknown, MolNotation.MolBlock, module);
}

export function molToSmiles(molStr: string, module: RDModule): string {
  return convertMoleculeImpl(molStr, MolNotation.Unknown, MolNotation.Smiles, module);
}