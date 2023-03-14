import {category, expectArray, test, before} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {getRdKitModule} from '../package';
import {_convertMolNotation} from '../utils/convert-notation-utils';

// import {_package} from '../package-test';
import {loadFileAsText} from './utils';

category('converters', async () => {
  let molfileV2K: string;
  let molfileV3K: string;
  before(async () => { // wait until RdKit module gets initialized
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    molfileV2K = await loadFileAsText('tests/molfileV2000.mol');
    molfileV3K = await loadFileAsText('tests/molfileV3000.mol');
  });

  const molecules: { [key: string]: string[]} = {
    smiles: ['CN1C(=O)CN=C(C2CCCCC2)c2ccccc21'],
    smarts: [
      '[#6]-[#7]1-[#6](=[#8])-[#6]-[#7]=[#6](-[#6]2-[#6]-[#6]-[#6]-[#6]-[#6]-2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-1',
    ],
    molblock: [molfileV2K],
    v3Kmolblock: [
      `
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 19 21 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -4.030000 -1.540100 0.000000 0
M  V30 2 N -2.534200 -1.428000 0.000000 0
M  V30 3 C -1.784200 -2.727000 0.000000 0
M  V30 4 O -2.629200 -3.966400 0.000000 0
M  V30 5 C -0.300900 -2.950600 0.000000 0
M  V30 6 N 0.798600 -1.930300 0.000000 0
M  V30 7 C 0.686500 -0.434500 0.000000 0
M  V30 8 C 1.985600 0.315500 0.000000 0
M  V30 9 C 3.284600 -0.434500 0.000000 0
M  V30 10 C 4.583700 0.315500 0.000000 0
M  V30 11 C 4.583700 1.815500 0.000000 0
M  V30 12 C 3.284600 2.565500 0.000000 0
M  V30 13 C 1.985600 1.815500 0.000000 0
M  V30 14 C -0.552800 0.410400 0.000000 0
M  V30 15 C -0.219000 1.872800 0.000000 0
M  V30 16 C -1.318600 2.893100 0.000000 0
M  V30 17 C -2.752000 2.451000 0.000000 0
M  V30 18 C -3.085800 0.988600 0.000000 0
M  V30 19 C -1.986200 -0.031700 0.000000 0
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
`],
  };

  function _testConvert(srcNotation: DG.chem.Notation, tgtNotation: DG.chem.Notation) {
    const result = [];
    for (const mol of molecules[srcNotation])
      result.push(_convertMolNotation(mol, srcNotation, tgtNotation, getRdKitModule()));
    expectArray(result, molecules[tgtNotation]);
  }

  test('testSmilesToMolfileV2000', async () => {
    _testConvert(DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
  });
  test('testSmilesToSmarts', async () => {
    _testConvert(DG.chem.Notation.Smiles, DG.chem.Notation.Smarts);
  });
  test('testMolfileV2000ToSmiles', async () => {
    _testConvert(DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
  });
  test('testMolfileV2000ToSmarts', async () => {
    _testConvert(DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts);
  });
  test('testMolfileV2000ToV3000', async () => {
    _testConvert(DG.chem.Notation.MolBlock, DG.chem.Notation.V3KMolBlock);
  });
  test('testMolfileV3000ToSmarts', async () => {
    _testConvert(DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smarts);
  });
  test('testMolfileV3000ToSmiles', async () => {
    _testConvert(DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smiles);
  });
});
