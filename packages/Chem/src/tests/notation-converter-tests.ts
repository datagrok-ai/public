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
  let molecules: { [key: string]: string[]};
  before(async () => { // wait until RdKit module gets initialized
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    molfileV2K = await loadFileAsText('tests/molfileV2000.mol');
    molfileV3K = await loadFileAsText('tests/molfileV3000.mol');
    molecules = {
      smiles: ['CN1C(=O)CN=C(C2CCCCC2)c2ccccc21'],
      smarts: [
        '[#6]-[#7]1-[#6](=[#8])-[#6]-[#7]=[#6](-[#6]2-[#6]-[#6]-[#6]-[#6]-[#6]-2)-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-1',
      ],
      molblock: [molfileV2K],
      v3Kmolblock: [molfileV3K],
    };
  });


  function _testConvert(srcNotation: DG.chem.Notation, tgtNotation: DG.chem.Notation) {
    const result = [];
    for (const mol of molecules[srcNotation])
      result.push(_convertMolNotation(mol, srcNotation, tgtNotation, getRdKitModule()));
    expectArray(result, molecules[tgtNotation]);
  }

  test('SMILES to Molfile V2000', async () => {
    _testConvert(DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
  });

  test('SMILES to SMARTS', async () => {
    _testConvert(DG.chem.Notation.Smiles, DG.chem.Notation.Smarts);
  });
  test('Molfile V2000 to SMILES', async () => {
    _testConvert(DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
  });
  test('Molfile V2000 to SMARTS', async () => {
    _testConvert(DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts);
  });
  test('Molfile V2000 to V3000', async () => {
    _testConvert(DG.chem.Notation.MolBlock, DG.chem.Notation.V3KMolBlock);
  });
  test('Molfile V3000 to SMARTS', async () => {
    _testConvert(DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smarts);
  });
  test('Molfile V3000 to SMILES', async () => {
    _testConvert(DG.chem.Notation.V3KMolBlock, DG.chem.Notation.Smiles);
  });
});
