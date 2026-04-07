import {before, category, test, expect} from '@datagrok-libraries/test/src/test';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

// import {_package} from '../package-test';
import {loadFileAsText, loadFileAsBytes} from './utils';

import {getSdfString} from '../utils/sdf-utils';
import {_importTripos} from '../file-importers/mol2-importer';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {_package} from '../package-test';

category('mol2 to SDF', async () => {
  let mol2InputBytes: Uint8Array;
  let mol2InputDF: DG.DataFrame;
  let sdfOutputString: string;

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    mol2InputBytes = await loadFileAsBytes('tests/mol2-test.mol2');
    mol2InputDF = _importTripos(mol2InputBytes)[0];
    await grok.data.detectSemanticTypes(mol2InputDF);
    sdfOutputString = await loadFileAsText('tests/mol2-test.sdf');
    sdfOutputString = sdfOutputString.replace(/\r/g, '');
  });

  test('Save SMILES column', async () => {
    const col = mol2InputDF.col('molecules')!;
    expect(await getSdfString(mol2InputDF, col), sdfOutputString);
  });
});
