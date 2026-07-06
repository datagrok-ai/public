/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import MolNotation = DG.chem.Notation;
// temporary polyfills for extended notations, remove these when js-api is updated. added to avoid local dependency
// @ts-ignore
MolNotation.CxSmarts = 'cxsmarts';
// @ts-ignore
MolNotation.CxSmiles = 'cxsmiles';


// enum MolNotation {
//   Smiles = 'smiles',
//   Smarts = 'smarts',
//   MolBlock = 'molblock', // molblock V2000
//   V3KMolBlock = 'v3Kmolblock', // molblock V3000
//   Unknown = 'unknown',
//   CxSmiles = 'cxsmiles', // extended smiles
//   CxSmarts = 'cxsmarts', // extended smarts
// }

// datagrok libraries dependencies
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {getMolSafe} from './mol-creation_rdkit';
import {getRdKitService} from './chem-common-rdkit';


export const MALFORMED_MOL_V2000 = `
Malformed

  0  0  0  0  0  0  0  0  0  0999 V2000
M  END`;
const MALFORMED_MOL_V3000 = `
Malformed

  0  0  0  0  0  0            999 V3000
M  END`;

/** True when the molecule has a defined absolute tetrahedral stereocenter (CIP R/S). */
export function hasDefinedStereo(mol: RDMol): boolean {
  try {
    return /\([RrSs]\)/.test(mol.get_stereo_tags());
  } catch (e) {
    return false;
  }
}

// Sets the V2000 counts-line chiral flag so absolute stereo is not loaded as racemate by OCL/Ketcher.
export function setMolBlockChiralFlag(molBlock: string, chiral: boolean): string {
  if (!chiral)
    return molBlock;
  const eol = molBlock.includes('\r\n') ? '\r\n' : '\n';
  const lines = molBlock.split(/\r?\n/);
  if (lines.length < 4 || !lines[3].includes('V2000'))
    return molBlock;
  lines[3] = lines[3].substring(0, 12) + '  1' + lines[3].substring(15);
  return lines.join(eol);
}

/**
 * Convert between the following notations: SMILES, SMARTS, Molfile V2000 and Molfile V3000
 *
 * @param {string} moleculeString  admissible notations: listed below
 * @param {string} sourceNotation  possible values: 'smiles', 'smarts',
 * 'molblock' (corresponding to Molfile V2000) and 'V3KMolBlock' (corresponding to
 * Molfile V3000)
 * @param {string} targetNotation  possible values: same as for sourceNotation
 * @param {RDModule} rdKitModule
 * @param {boolean} addHs
 * @return {string} the converted representation
 */
export function _convertMolNotation(
  moleculeString: string,
  sourceNotation: DG.chem.Notation,
  targetNotation: DG.chem.Notation,
  rdKitModule: RDModule,
  addHs: boolean = false,
): string {
  if (sourceNotation === targetNotation)
    throw new Error(`Convert molecule notation: source and target notations must differ: "${sourceNotation}"`);
  let result = (targetNotation === MolNotation.MolBlock) ? MALFORMED_MOL_V2000 :
    (targetNotation === MolNotation.V3KMolBlock) ? MALFORMED_MOL_V3000 :
      'MALFORMED_INPUT_VALUE';
  let mol: RDMol | null = null;
  try {
    mol = getMolSafe(moleculeString, {}, rdKitModule).mol;
    if (mol) {
      if (targetNotation === MolNotation.MolBlock) {
        //when converting from smiles set coordinates and rendering parameters
        if (sourceNotation === MolNotation.Smiles) {
          if (!mol.has_coords())
            mol!.set_new_coords();
          mol.normalize_depiction(1);
          mol.straighten_depiction(false);
        }
        if (addHs) {
          try {
            mol.add_hs_in_place();
          } catch (e) {}
        }
        result = setMolBlockChiralFlag(mol.get_molblock(), hasDefinedStereo(mol));
      }
      if (targetNotation === MolNotation.Smiles)
        result = mol.get_smiles();
      if (targetNotation === MolNotation.V3KMolBlock)
        result = mol.get_v3Kmolblock();
      if (targetNotation === MolNotation.Smarts)
        result = mol.get_smarts(); // @ts-ignore
      if (targetNotation === MolNotation.CxSmiles)
        result = mol.get_cxsmiles(); // @ts-ignore
      if (targetNotation === MolNotation.CxSmarts)
        result = mol.get_cxsmarts();
      // if (targetNotation === MolNotation.Inchi)
      //   return mol.get_inchi();
    }
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

export async function convertNotationForColumn(molecules: DG.Column<string>, targetNotation: MolNotation, kekulize = false): Promise<string[]> {
  return (await getRdKitService()).convertMolNotation(molecules.toList(), targetNotation, kekulize);
}
