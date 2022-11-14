import {before, category, test, expect} from '@datagrok-libraries/utils/src/test';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {readDataframe, loadFileAsText} from './utils';

import {getSdfString} from '../utils/sdf-utils';

category('saveAsSdf', async () => {
  let inputDf: DG.DataFrame;
  let fileWithSavedSmiles: string;
  let fileWithSavedMolblock: string;

  before(async () => {
    // saveAsSdf employs RDKit for conversion
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    inputDf = await readDataframe('tests/sdf-test.csv');
    await grok.data.detectSemanticTypes(inputDf);
    fileWithSavedSmiles = await loadFileAsText('tests/sdf-test-smiles.sdf');
    fileWithSavedMolblock = await loadFileAsText('tests/sdf-test-scaffold.sdf');
    fileWithSavedSmiles = fileWithSavedSmiles.replace(/\r/g, '');
    fileWithSavedMolblock = fileWithSavedMolblock.replace(/\r/g, '');
  });

  test('saveSmilesColumn', async () => {
    const savedColumn = inputDf.col('Smiles')!;
    expect(getSdfString(inputDf, savedColumn), fileWithSavedSmiles);
  });

  test('saveMolblockColumn', async () => {
    const savedColumn = inputDf.col('Scaffold')!;
    expect(getSdfString(inputDf, savedColumn).replace(/\r/g, ''), fileWithSavedMolblock);
  });
});
