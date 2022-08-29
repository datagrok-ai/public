import { category, expectArray, test, before } from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {getRdKitModule} from '../package';
import {_convertMolNotation, MolNotation} from '../utils/convert-notation-utils';

category('convert mol notation test', () => {
  before(async () => { // wait until RdKit module gets initialized
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

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
      `[#6]-[#7]1-[#6](=[#8])-[#6]-[#7]=[#6](-[#6]2-[#6]-[#6]-[#6]-[#6]-[#6]-2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-1`,
    ],
    molblock: [
`
     RDKit          2D

 19 21  0  0  0  0  0  0  0  0999 V2000
   -3.1476   -2.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6852   -2.4718    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3515   -3.9342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5242   -4.8694    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -4.5850    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3515   -3.9342    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.6852   -2.4718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1476   -2.1380    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.1679   -3.2376    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.6303   -2.9038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0724   -1.4704    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0522   -0.3709    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.5898   -0.7046    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  3  5  1  0
  5  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  1  0
  9 10  1  0
 10 11  1  0
 11 12  1  0
 12 13  1  0
  7 14  1  0
 14 15  2  0
 15 16  1  0
 16 17  2  0
 17 18  1  0
 18 19  2  0
 19  2  1  0
 13  8  1  0
 19 14  1  0
M  END
`    ],
    v3Kmolblock: [
`
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 19 21 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -3.147600 -2.138000 0.000000 0
M  V30 2 N -1.685200 -2.471800 0.000000 0
M  V30 3 C -1.351500 -3.934200 0.000000 0
M  V30 4 O -2.524200 -4.869400 0.000000 0
M  V30 5 C -0.000000 -4.585000 0.000000 0
M  V30 6 N 1.351500 -3.934200 0.000000 0
M  V30 7 C 1.685200 -2.471800 0.000000 0
M  V30 8 C 3.147600 -2.138000 0.000000 0
M  V30 9 C 4.167900 -3.237600 0.000000 0
M  V30 10 C 5.630300 -2.903800 0.000000 0
M  V30 11 C 6.072400 -1.470400 0.000000 0
M  V30 12 C 5.052200 -0.370900 0.000000 0
M  V30 13 C 3.589800 -0.704600 0.000000 0
M  V30 14 C 0.750000 -1.299000 0.000000 0
M  V30 15 C 1.500000 0.000000 0.000000 0
M  V30 16 C 0.750000 1.299000 0.000000 0
M  V30 17 C -0.750000 1.299000 0.000000 0
M  V30 18 C -1.500000 0.000000 0.000000 0
M  V30 19 C -0.750000 -1.299000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 3 5
M  V30 5 1 5 6
M  V30 6 2 6 7
M  V30 7 1 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 1 11 12
M  V30 12 1 12 13
M  V30 13 1 7 14
M  V30 14 2 14 15
M  V30 15 1 15 16
M  V30 16 2 16 17
M  V30 17 1 17 18
M  V30 18 2 18 19
M  V30 19 1 19 2
M  V30 20 1 13 8
M  V30 21 1 19 14
M  V30 END BOND
M  V30 END CTAB
M  END
`    ],
    inchi: [
      'InChI=1S/C16N2O/c1-18-14-10-6-5-9-13(14)16(17-11-15(18)19)12-7-3-2-4-8-12',
    ],
  };

  function _testConvert(srcNotation: MolNotation, tgtNotation: MolNotation) {
    const result = [];
    for (let mol of molecules[srcNotation])
      result.push(_convertMolNotation(mol, srcNotation, tgtNotation, getRdKitModule()));
    console.log("The result is");
    console.log([result.toString()]);
    console.log("The expected value is");
    console.log(molecules[tgtNotation]);
    expectArray(result, molecules[tgtNotation]);
  }

  test('testSmilesToMolfileV2000', async () => {
    _testConvert(MolNotation.Smiles, MolNotation.MolBlock);
  });
  // test('testMolfileV2000ToSmiles', async () => {
  //   _testConvert(MolNotation.MolBlock, MolNotation.Smiles);
  // });
  test('testMolfileV2000ToV3000', async () => {
    _testConvert(MolNotation.MolBlock, MolNotation.V3KMolBlock);
  });
  // test('testMolfileV3000ToV2000', async () => {
  //   _testConvert(MolNotation.V3KMolBlock, MolNotation.MolBlock);
  // });
  // test('testSmartsToMolfileV3000', async () => {
  //   _testConvert(MolNotation.Smarts, MolNotation.V3KMolBlock);
  // });
  test('testMolfileV3000ToSmarts', async () => {
    _testConvert(MolNotation.V3KMolBlock, MolNotation.Smarts);
  });
  test('testSmartsToInchi', async () => {
    _testConvert(MolNotation.Smarts, MolNotation.Inchi);
  });
  // test('testInchiToSmarts', async () => {
  //   _testConvert(MolNotation.Inchi, MolNotation.Smarts);
  // });
  // test('testSmilesToInchi', async () => {
  //   _testConvert(MolNotation.Smiles, MolNotation.Inchi);
  // });
  test('testInchiToSmiles', async () => {
    _testConvert(MolNotation.Inchi, MolNotation.Smiles);
  });
});
