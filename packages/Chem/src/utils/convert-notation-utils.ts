import * as DG from 'datagrok-api/dg';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import MolNotation = DG.chem.Notation;

// datagrok libraries dependencies
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';


const MALFORMED_MOL_V2000 = `
Malformed

  0  0  0  0  0  0  0  0  0  0999 V2000
M  END`;
const MALFORMED_MOL_V3000 = `
Malformed

  0  0  0  0  0  0            999 V3000
M  END`;

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
  let result = (targetNotation === MolNotation.MolBlock) ? MALFORMED_MOL_V2000 :
    (targetNotation === MolNotation.V3KMolBlock) ? MALFORMED_MOL_V3000 :
      'MALFORMED_INPUT_VALUE';
  let mol: RDMol | null = null;
  try {
    mol = rdKitModule.get_mol(moleculeString);
    if (targetNotation === MolNotation.MolBlock) {
      //when converting from smiles set coordinates and rendering parameters
      if (sourceNotation === MolNotation.Smiles) {
        if (!mol.has_coords())
          mol.set_new_coords();
        mol.normalize_depiction(1);
        mol.straighten_depiction(false);
      }
      result = mol.get_molblock();
    }
    if (targetNotation === MolNotation.Smiles)
      result = mol.get_smiles();
    if (targetNotation === MolNotation.V3KMolBlock)
      result = mol.get_v3Kmolblock();
    if (targetNotation === MolNotation.Smarts)
      result = mol.get_smarts();
    // if (targetNotation === MolNotation.Inchi)
    //   return mol.get_inchi();
  } catch (err) {
    console.error(errorToConsole(err));
  } finally {
    mol?.delete();
    return result;
  }
}

export function molToMolblock(molStr: string, module: RDModule): string {
  return DG.chem.isMolBlock(molStr) ?
    molStr : _convertMolNotation(molStr, MolNotation.Unknown, MolNotation.MolBlock, module);
}

export function molToSmiles(molStr: string, module: RDModule): string {
  return _convertMolNotation(molStr, MolNotation.Unknown, MolNotation.Smiles, module);
}
