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
  const rdKitModule = getRdKitModule();
  const mol = rdKitModule.get_mol(moleculeString);
  const converter =
    (sourceNotation == 'smiles' && targetNotation == 'molblockV3000') ?
      mol.get_v3Kmolblock :
      (sourceNotation == 'molblockV3000' && targetNotation == 'molblockV2000') ?
        mol.get_molblock :
        (sourceNotation == 'molblockV2000' && targetNotation == 'molblockV3000') ?
          mol.get_v3Kmolblock :
          ( (sourceNotation == 'molblockV2000' || sourceNotation == 'molblockV3000') && targetNotation == 'smiles'
          ) ? mol.get_smiles : null;
  if (converter) {
    const result = converter();
    mol?.delete();
    return result;
  } else {
    mol?.delete();
    throw new Error(
      'The source notation ${sourceNotation} or target notation ${targetNotation} is incorrect',
    );
  }
}
