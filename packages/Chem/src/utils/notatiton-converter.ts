import {getRdKitModule} from './chem-common-rdkit';

/**
 * Convert between any of the notations: SMILES, Molfile V2000 and Molfile V3000
 *
 * @param {string} moleculeString  possible notations: smiles, molfile
 * @param {string} sourceNotation  possible values: 'smiles', 'molblockV2000',
 * @param {string} targetNotation  possible values: 'smiles', 'molblockV2000',
 * 'molblockV3000'
 * @return {string} the converted representation
 */
export function convertNotation(
  moleculeString: string,
  sourceNotation: string,
  targetNotation: string,
): string {
  const admissibleNotations = new Set(['smiles', 'molblockV2000', 'molblockV3000']);
  if (!admissibleNotations.has(sourceNotation) ||
    !admissibleNotations.has(targetNotation) ||
    (sourceNotation == targetNotation)) {
    throw new Error(
      `The source notation ${sourceNotation} or target notation ${targetNotation} is incorrect`,
    );
  }
  const rdKitModule = getRdKitModule();
  const mol = rdKitModule.get_mol(moleculeString);
  const converter =
    (targetNotation == 'molblockV2000') ?
      mol.get_molblock :
      (targetNotation == 'molblockV3000') ?
        mol.get_v3Kmolblock : mol.get_smiles;
  const result = converter();
  mol?.delete();
  return result;
}
