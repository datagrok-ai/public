import {RDModule} from '../rdkit-api';

export enum MolNotation {
  Smiles = 'smiles',
  Smarts = 'smarts',
  MolBlock = 'molblock', // molblock V2000
  V3KMolBlock = 'v3Kmolblock', // molblock V3000
  Inchi = 'inchi',
  Unknown = 'unknown',
}

export function isMolBlock(s: string) {
  return s.includes('M  END');
}

/**
 * Convert between the following notations: SMILES, SMARTS, Molfile V2000 and Molfile V3000
 *
 * @param {string} moleculeString  admissible notations: listed below
 * @param {string} sourceNotation  possible values: 'smiles', 'smarts',
 * 'molblock' (corresponding to Molfile V2000), 'V3KMolBlock' (corresponding to
 * Molfile V3000) and 'inchi'
 * @param {string} targetNotation  possible values: same as for sourceNotation
 * @return {string} the converted representation
 */
export function _convertNotation(
  moleculeString: string,
  sourceNotation: string,
  targetNotation: string,
  rdKitModule: RDModule,
): string {
  if (sourceNotation === targetNotation)
    throw new Error(`Convert molecule notation: source and target notations must differ: "${sourceNotation}"`);
  const mol = rdKitModule.get_mol(moleculeString);
  try {
    const converter =
      (targetNotation === MolNotation.MolBlock) ?
        mol.get_molblock :
        (targetNotation === MolNotation.Smiles) ?
          mol.get_smiles :
          (targetNotation === MolNotation.V3KMolBlock) ?
            mol.get_v3Kmolblock :
            (targetNotation === MolNotation.Smarts) ?
              mol.get_smarts :
              (targetNotation == MolNotation.Inchi) ?
                mol.get_inchi : null;
    if (converter)
      return converter();
    else
      throw new Error(`Failed to convert molucule notation, target notation unknown: ${targetNotation}`);
  } finally {
    mol?.delete();
  }
}

export function molToMolblock(molStr: string, module: RDModule): string {
  return isMolBlock(molStr)
    ? molStr
    : _convertNotation(molStr, MolNotation.Unknown, MolNotation.MolBlock, module);
}

export function molToSmiles(molStr: string, module: RDModule): string {
  return _convertNotation(molStr, MolNotation.Unknown, MolNotation.Smiles, module);
}
