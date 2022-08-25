import {getRdKitModule} from './chem-common-rdkit';

/**
 * Convert between any of the notations: SMILES, Molfile V2000 and Molfile V3000
 *
 * @param {string} moleculeString  admissible notations: same as possible values
 * for source/target notation parameters
 * @param {string} sourceNotation  possible values: 'smiles', 'smarts',
 * 'molblockV2000', 'molblockV3000'
 * @param {string} targetNotation  possible values: 'smiles', 'smarts',
 * 'molblockV2000', 'molblockV3000'
 * @return {string} the converted representation
 */
export async function _convertNotation(
  moleculeString: string,
  sourceNotation: string,
  targetNotation: string,
): Promise<string> {
  const admissibleNotations = new Set(['smiles', 'smarts', 'molblockV2000', 'molblockV3000']);
  if (!admissibleNotations.has(sourceNotation) ||
    !admissibleNotations.has(targetNotation) ||
    (sourceNotation == targetNotation)) {
    throw new Error(
      `The source notation ${sourceNotation} or target notation ${targetNotation} is incorrect, or both coincide`,
    );
  }
  const rdKitModule = await getRdKitModule();
  const mol = rdKitModule.get_mol(moleculeString);
  const converter =
    (targetNotation == 'molblockV2000') ?
      mol.get_molblock :
      (targetNotation == 'molblockV3000') ?
        mol.get_v3Kmolblock :
        (targetNotation == 'smarts') ?
          mol.get_smarts : mol.get_smiles;
  const result = converter();
  mol?.delete();
  return result;
}
