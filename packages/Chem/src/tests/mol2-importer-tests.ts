import {before, category, test, expect} from '@datagrok-libraries/utils/src/test';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';

// import {_package} from '../package-test';
import {loadFileAsText, loadFileAsBytes} from './utils';

import {getSdfString} from '../utils/sdf-utils';
import {_importTripos} from '../file-importers/mol2-importer';

category('mol2toSdf', async () => {
  let mol2InputBytes: Uint8Array;
  let mol2InputDF: DG.DataFrame;
  let sdfOutputString: string;

  before(async () => {
    mol2InputBytes = await loadFileAsBytes('tests/mol2-test.mol2');
    mol2InputDF = _importTripos(mol2InputBytes)[0];
    await grok.data.detectSemanticTypes(mol2InputDF);
    sdfOutputString = await loadFileAsText('tests/mol2-test.sdf');
    // mol2InputDF = mol2InputDF.replace(/\r/g, '');
    sdfOutputString = sdfOutputString.replace(/\r/g, '');
  });

  test('saveSmilesColumn', async () => {
    const col = mol2InputDF.col('molecules')!;
    expect(getSdfString(mol2InputDF, col), sdfOutputString);
  });
});
