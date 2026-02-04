/* eslint-disable max-len */
import {before, category, test, expect} from '@datagrok-libraries/test/src/test';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {readDataframe, loadFileAsText} from './utils';

import {getSdfString} from '../utils/sdf-utils';
import {_importSdfString} from '../open-chem/sdf-importer';

category('save as sdf', async () => {
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
    inputDf.getCol('Smiles').semType = 'Molecule';
    fileWithSavedSmiles = await loadFileAsText('tests/sdf-test-smiles.sdf');
    fileWithSavedMolblock = await loadFileAsText('tests/sdf-test-scaffold.sdf');
    fileWithSavedSmiles = fileWithSavedSmiles.replace(/\r/g, '');
    fileWithSavedMolblock = fileWithSavedMolblock.replace(/\r/g, '');
  });

  test('save Smiles column', async () => {
    const savedColumn = inputDf.col('Smiles')!;
    const sdfString = await getSdfString(inputDf, savedColumn);
    expectSdf(sdfString.replace(/\r/g, ''), fileWithSavedSmiles);
  });

  test('save Molblock column', async () => {
    const savedColumn = inputDf.col('Scaffold')!;
    const sdfString = await getSdfString(inputDf, savedColumn);
    expectSdf(sdfString.replace(/\r/g, ''), fileWithSavedMolblock);
  });
});

function expectSdf(actual: string, expected: string) {
  const actualDF = _importSdfString(actual)[0];
  const expectedDF = _importSdfString(expected)[0];
  expect(actualDF.rowCount, expectedDF.rowCount, `Row count differs: actual ${actualDF.rowCount}, expected ${expectedDF.rowCount}`);
  expect(actualDF.columns.length, expectedDF.columns.length, `Column count differs: actual ${actualDF.columns.length}, expected ${expectedDF.columns.length}`);
  expectedDF.columns.names().forEach((col) => {
    const actualCol = actualDF.col(col)!;
    expect(actualCol != null, true, `Column ${col} is missing from actual`);

    const expectedCol = expectedDF.col(col)!;
    for (let i = 0; i < expectedDF.rowCount; i++) {
      const expectedVal = expectedCol.get(i);
      const actualVal = actualCol.get(i);
      if (typeof expectedVal === 'number' && typeof actualVal === 'number')
        expect((isNaN(actualVal) && isNaN(expectedVal)) || (actualVal as number).toFixed(2) === expectedVal.toFixed(2), true, `Column ${col}, row ${i} differs: actual ${actualVal}, expected ${expectedVal}`);
      else
        expect(actualVal, expectedVal, `Column ${col}, row ${i} differs: actual ${actualVal}, expected ${expectedVal}`);
    }
  });
}
