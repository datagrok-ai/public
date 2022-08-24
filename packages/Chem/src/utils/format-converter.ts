import {getRdKitModule} from './chem-common-rdkit';

/**
 * Convert between any of the formats: SMILES, Molfile V2000 and Molfile V3000
 *
 * @param {string} moleculeString  possible formats: smiles, molfile
 * @param {string} sourceFormat  possible values: 'smiles', 'v2KMolblock',
 * @param {string} targetFormat  possible values: 'smiles', 'v2KMolblock',
 * 'v3KMolblock'
 * @return {string} the converted representation
 */
export function convertFormat(
  moleculeString: string,
  sourceFormat: string,
  targetFormat: string,
): string {
  const rdKitModule = getRdKitModule();
  if (sourceFormat == 'smiles' && targetFormat == 'v3KMolblock') {
    const mol = rdKitModule.get_mol(moleculeString);
    const v3KMolblock = mol.get_v3Kmolblock();
    mol?.delete();
    return v3KMolblock;
  } else if (sourceFormat == 'v3KMolblock' && targetFormat == 'v2KMolblock') {
    const mol = rdKitModule.get_mol(moleculeString);
    const v2KMolblock = mol.get_molblock();
    mol?.delete();
    return v2KMolblock;
  } else if (sourceFormat == 'v2KMolblock' && targetFormat == 'v3KMolblock') {
    const mol = rdKitModule.get_mol(moleculeString);
    const v3KMolblock = mol.get_v3Kmolblock();
    mol?.delete();
    return v3KMolblock;
  } else if (
    (sourceFormat == 'v2KMolblock' || sourceFormat == 'v3KMolblock') &&
    targetFormat == 'smiles'
  ) {
    const mol = rdKitModule.get_mol(moleculeString);
    const smiles = mol.get_smiles();
    mol?.delete();
    return smiles;
  } else {
    throw new Error(
      'The source format ${sourceFormat} or target format ${targetFormat} is incorrect',
    );
  }
}
