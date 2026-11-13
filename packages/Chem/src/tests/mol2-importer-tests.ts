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

category('mol2 import', () => {
  const mol2 = (atomLines: string[], end: string = '\n'): Uint8Array => new TextEncoder().encode(
    '@<TRIPOS>MOLECULE\nname\n 2 1 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n' +
    atomLines.join('\n') + '\n@<TRIPOS>BOND\n     1    1    2 1' + end);
  const molblock = (bytes: Uint8Array): string => _importTripos(bytes)[0].col('molecules')!.get(0);

  test('Atoms without the charge column', async () => {
    const res = molblock(mol2(['  1 C1 1.0 2.0 3.0 1 CAN', '  2 O1 4.0 5.0 6.0 1 CAN']));
    expect(res.includes('M  V30 1 C 1 2 3 0\n'), true);
    expect(res.includes('M  V30 2 O 4 5 6 0\n'), true);
    expect(res.includes('M  V30 1 1 1 2\n'), true);
  });

  test('Eight columns are not read as a charge from the next line', async () => {
    const res = molblock(mol2(['  1 C1 1.0 2.0 3.0 1 CAN 0.0000', '  2 O1 4.0 5.0 6.0 1 CAN 0.0000']));
    expect(res.includes('CHG='), false);
    expect(res.includes('M  V30 COUNTS 2 1 0 0 0\n'), true);
  });

  test('Charge in the standard column', async () => {
    const res = molblock(mol2(['  1 C1 1.0 2.0 3.0 C.3 1 LIG 1.0000', '  2 O1 4.0 5.0 6.0 O.3 1 LIG -1.0000']));
    expect(res.includes('M  V30 1 C 1 2 3 0 CHG=1\n'), true);
    expect(res.includes('M  V30 2 O 4 5 6 0 CHG=-1\n'), true);
  });

  test('File without a trailing line break', async () => {
    const res = molblock(mol2(['  1 C1 1.0 2.0 3.0 C.3 1 LIG 0.0', '  2 O1 4.0 5.0 6.0 O.3 1 LIG 0.0'], ''));
    expect(res.includes('M  V30 1 1 1 2\n'), true);
  });

  test('Fewer records than declared throw', async () => {
    let error = '';
    try {
      molblock(mol2(['  1 C1 1.0 2.0 3.0 C.3 1 LIG 0.0']));
    }
    catch (e: any) {
      error = e.message;
    }
    expect(error.includes('fewer records'), true);
  });
});
