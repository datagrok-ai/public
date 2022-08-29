import {category, expectArray, test} from '@datagrok-libraries/utils/src/test';

// import * as grok from 'datagrok-api/grok';
// import * as DG from 'datagrok-api/dg';

import {getRdKitModule} from '../package';
import {_convertMolNotation, MolNotation} from '../utils/convert-notation-utils';

category('convertMolNotationTest', () => {
  const molecules: { [key: string]: string[]} = {
    smiles: [
      'CN1C(=O)CN=C(C2CCCCC2)c2ccccc21',
      // 'CCCCN1C(=O)CN=C(c2ccccc12)C3CCCCC3',
      // 'CC(C)CCN1C(=O)CN=C(c2ccccc12)C3CCCCC3',
      // 'O=C1CN=C(c2ccccc2N1CC3CCCCC3)C4CCCCC4',
      // 'O=C1CN=C(c2cc(Cl)ccc2N1)C3CCCCC3',
      // 'CN1C(=O)CN=C(c2cc(Cl)ccc12)C3CCCCC3',
    ],
    smarts: [
      `[#6]-[#7]1-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-[#6](-[#6]2-[#6]-[#6]-[#6]-[#6]-[#6]-2)=[#7]-[#6]-[#6]-1=[#8]`,
    ],
    molblock: [
      `
     RDKit          2D

 19 21  0  0  0  0  0  0  0  0999 V2000
    3.1476   -2.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.6852   -2.4718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6852   -2.4718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.1476   -2.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.5898   -0.7046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.0522   -0.3709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0724   -1.4704    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6303   -2.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.1679   -3.2376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3515   -3.9342    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -4.5850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3515   -3.9342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5242   -4.8694    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  1  0
 10 11  1  0
 11 12  1  0
 12 13  1  0
 13 14  1  0
 14 15  1  0
  9 16  2  0
 16 17  1  0
  2 18  1  0
 18 19  2  0
  8  3  1  0
 15 10  1  0
 18 17  1  0
M  END`,
    ],
    v3Kmolblock: [
      `
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 19 21 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 3.147627 -2.138004 0.000000 0
M  V30 2 N 1.685235 -2.471785 0.000000 0
M  V30 3 C 0.750000 -1.299038 0.000000 0
M  V30 4 C 1.500000 0.000000 0.000000 0
M  V30 5 C 0.750000 1.299038 0.000000 0
M  V30 6 C -0.750000 1.299038 0.000000 0
M  V30 7 C -1.500000 0.000000 0.000000 0
M  V30 8 C -0.750000 -1.299038 0.000000 0
M  V30 9 C -1.685235 -2.471785 0.000000 0
M  V30 10 C -3.147627 -2.138004 0.000000 0
M  V30 11 C -3.589759 -0.704645 0.000000 0
M  V30 12 C -5.052151 -0.370863 0.000000 0
M  V30 13 C -6.072410 -1.470441 0.000000 0
M  V30 14 C -5.630278 -2.903800 0.000000 0
M  V30 15 C -4.167886 -3.237582 0.000000 0
M  V30 16 N -1.351453 -3.934177 0.000000 0
M  V30 17 C -0.000000 -4.585003 0.000000 0
M  V30 18 C 1.351453 -3.934177 0.000000 0
M  V30 19 O 2.524201 -4.869412 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 7
M  V30 7 2 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 11 12
M  V30 12 1 12 13
M  V30 13 1 13 14
M  V30 14 1 14 15
M  V30 15 2 9 16
M  V30 16 1 16 17
M  V30 17 1 2 18
M  V30 18 2 18 19
M  V30 19 1 8 3
M  V30 20 1 15 10
M  V30 21 1 18 17
M  V30 END BOND
M  V30 END CTAB
M  END`,
    ],
    inchi: [
      'InChI=1S/C16H20N2O/c1-18-14-10-6-5-9-13(14)16(17-11-15(18)19)12-7-3-2-4-8-12/h5-6,9-10,12H,2-4,7-8,11H2,1H3',
    ],
  };

  function _testConvert(srcNotation: MolNotation, tgtNotation: MolNotation) {
    const result = [];
    for (const mol of molecules.srcNotation)
      result.push(_convertMolNotation(mol, srcNotation, tgtNotation, getRdKitModule()));
    expectArray(result, molecules.tgtNotation);
  }

  test('testSmilesToMolfileV2000', async () => {
    _testConvert(MolNotation.Smiles, MolNotation.MolBlock);
  });
  test('testMolfileV2000ToV3000', async () => {
    _testConvert(MolNotation.MolBlock, MolNotation.V3KMolBlock);
  });
  test('testMolfileV3000ToSmarts', async () => {
    _testConvert(MolNotation.V3KMolBlock, MolNotation.Smarts);
  });
  test('testSmartsToInchi', async () => {
    _testConvert(MolNotation.Smarts, MolNotation.Inchi);
  });
  test('testInchiToSmiles', async () => {
    _testConvert(MolNotation.Inchi, MolNotation.Smiles);
  });
});
