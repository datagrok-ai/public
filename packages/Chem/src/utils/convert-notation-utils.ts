import * as DG from 'datagrok-api/dg';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import MolNotation = DG.chem.Notation;


/**
 * Convert between the following notations: SMILES, SMARTS, Molfile V2000 and Molfile V3000
 *
 * @param {string} moleculeString  admissible notations: listed below
 * @param {string} sourceNotation  possible values: 'smiles', 'smarts',
 * 'molblock' (corresponding to Molfile V2000) and 'V3KMolBlock' (corresponding to
 * Molfile V3000)
 * @param {string} targetNotation  possible values: same as for sourceNotation
 * @param {RDModule} rdKitModule
 * @return {string} the converted representation
 */
export function _convertMolNotation(
  moleculeString: string,
  sourceNotation: string,
  targetNotation: string,
  rdKitModule: RDModule,
): string {
  if (sourceNotation === targetNotation)
    throw new Error(`Convert molecule notation: source and target notations must differ: "${sourceNotation}"`);
  const mol = rdKitModule.get_mol(moleculeString);
  try {
    if (targetNotation === MolNotation.MolBlock) {
      //when converting from smiles set coordinates and rendering parameters
      if (sourceNotation === MolNotation.Smiles) {
        if (!mol.has_coords())
          mol.set_new_coords();
        mol.normalize_depiction(1);
        mol.straighten_depiction(false);
      }
      return mol.get_molblock();
    }
    if (targetNotation === MolNotation.Smiles)
      return mol.get_smiles();
    if (targetNotation === MolNotation.V3KMolBlock)
      return mol.get_v3Kmolblock();
    if (targetNotation === MolNotation.Smarts)
      return mol.get_smarts();
    // if (targetNotation === MolNotation.Inchi)
    //   return mol.get_inchi();
    throw new Error(`Failed to convert molucule notation, target notation unknown: ${targetNotation}`);
  } finally {
    mol?.delete();
  }
}

export function molToMolblock(molStr: string, module: RDModule): string {
  return DG.chem.isMolBlock(molStr) ?
    molStr : _convertMolNotation(molStr, MolNotation.Unknown, MolNotation.MolBlock, module);
}

export function molToSmiles(molStr: string, module: RDModule): string {
  return _convertMolNotation(molStr, MolNotation.Unknown, MolNotation.Smiles, module);
}
