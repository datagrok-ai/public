import {category, expectArray, test, before} from '@datagrok-libraries/test/src/test';

import {loadFileAsText} from './utils';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';

type TestedData = {
  X: number[],
  Y: number[],
  Z: number[],
  ATOM_TYPES: string[],
  BONDED_ATOMS: number[][],
  BOND_TYPES: number[],
}

const PARSED_X = [-4.0300, -2.5342, -1.7842, -2.6292, -0.3009, 0.7986, 0.6865, 1.9856, 3.2846, 4.5837, 4.5837, 3.2846,
  1.9856, -0.5528, -0.2190, -1.3186, -2.7520, -3.0858, -1.9862];

const PARSED_Y = [-1.5401, -1.4280, -2.7270, -3.9664, -2.9506, -1.9303, -0.4345, 0.3155, -0.4345, 0.3155,
  1.8155, 2.5655, 1.8155, 0.4104, 1.8728, 2.8931, 2.4510, 0.9886, -0.0317];

const ATOM_TYPES_V2K = ['C', 'N', 'C', 'O', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'];

const ATOM_TYPES_V3K = ['NOT [C,N]', 'N', 'C', 'O', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
  'C', 'C', 'C', 'C', 'C'];

const PARSED_BONDED_ATOMS = [
  [1, 2], [2, 3], [3, 4], [3, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12], [12, 13], [7, 14],
  [14, 15], [15, 16], [16, 17], [17, 18], [18, 19], [19, 2], [13, 8], [19, 14],
];

const PARSED_BOND_TYPES = [1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 1, 1];

const OUTPUT_V2K: TestedData = {
  X: PARSED_X,
  Y: PARSED_Y,
  Z: new Array<number>(19).fill(0),
  ATOM_TYPES: ATOM_TYPES_V2K,
  BONDED_ATOMS: PARSED_BONDED_ATOMS,
  BOND_TYPES: PARSED_BOND_TYPES,
};

const OUTPUT_V3K = Object.assign({}, OUTPUT_V2K);
OUTPUT_V3K.ATOM_TYPES = ATOM_TYPES_V3K;

category('chemical table parsing', async () => {
  let molfileV2K: string;
  let molfileV3K: string;

  before(async () => {
    molfileV2K = await loadFileAsText('tests/molfileV2000.mol');
    molfileV3K = await loadFileAsText('tests/molfileV3000.mol');
  });

  function getRoundedNumberArray(floatArray: Float32Array): number[] {
    return Array.from(floatArray)
      .map((item: number) => Math.round(item * 10000) / 10000);
  }

  function _testCoordinates(molfile: string, expectedData: TestedData): void {
    const mol = MolfileHandler.getInstance(molfile);
    const expected = [expectedData.X, expectedData.Y, expectedData.Z];
    const obtained = [
      getRoundedNumberArray(mol.x),
      getRoundedNumberArray(mol.y),
      getRoundedNumberArray(mol.z),
    ];
    expectArray(expected, obtained);
  }

  function _testAtomTypes(molfile: string, expectedData: TestedData): void {
    const mol = MolfileHandler.getInstance(molfile);
    console.log('comment:', mol.atomTypes);
    expectArray(expectedData.ATOM_TYPES, mol.atomTypes);
  }

  function _testBondedAtoms(molfile: string, expectedData: TestedData): void {
    const mol = MolfileHandler.getInstance(molfile);
    const obtained = mol.pairsOfBondedAtoms.map(
      (item: Uint16Array) => Array.from(item),
    );
    expectArray(expectedData.BONDED_ATOMS, obtained);
  }

  function _testBondTypes(molfile: string, expectedData: TestedData): void {
    const mol = MolfileHandler.getInstance(molfile);
    const obtained = Array.from(mol.bondTypes);
    expectArray(expectedData.BOND_TYPES, obtained);
  }

  test('parse coordinates V2K', async () => {
    _testCoordinates(molfileV2K, OUTPUT_V2K);
  });

  test('parse atom types V2K', async () => {
    _testAtomTypes(molfileV2K, OUTPUT_V2K);
  });

  test('parse bonded atoms V2K', async () => {
    _testBondedAtoms(molfileV2K, OUTPUT_V2K);
  });

  test('parse bond types V2K', async () => {
    _testBondTypes(molfileV2K, OUTPUT_V2K);
  });

  test('parse coordinates V3K', async () => {
    _testCoordinates(molfileV3K, OUTPUT_V3K);
  });

  test('parse atom types V3K', async () => {
    _testAtomTypes(molfileV3K, OUTPUT_V3K);
  });

  test('parse bonded atoms V3K', async () => {
    _testBondedAtoms(molfileV3K, OUTPUT_V3K);
  });

  test('parse bond types V3K', async () => {
    _testBondTypes(molfileV3K, OUTPUT_V3K);
  });
});
