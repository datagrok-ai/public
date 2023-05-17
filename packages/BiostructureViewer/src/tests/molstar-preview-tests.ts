import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {
  category /*, expect*/,
  test,
  expect
} from '@datagrok-libraries/utils/src/test';
import {_packageName} from './utils';
import {previewMolstarUI} from '../viewers/molstar-viewer';

const validFileNames = ['1bdq.pdb', '1bdq.sdf', 'dc.mol2',
  '4tkx.mmcif', 'example.xyz', 'grofile.gro', 'pdbqt.pdbqt'];

category('MolstarPreview', () => {
  validFileNames.forEach((fn) => {
    test(`open${fn.substring(fn.indexOf('.'), fn.length)}`, async () => {
      let noException = true;
      const folderName: string = `System:AppData/${_packageName}/samples`;
      const file = (
        await grok.dapi.files.list(folderName, false, fn))[0];

      try {
        const {view, loadingPromise} = previewMolstarUI(file);
        grok.shell.newView('Molstar Preview', [view]);
        await loadingPromise;
      } catch (e) {
        noException = false;
      }
      expect(noException, true);
    });
  });

  // tests that opening csv through molstar causes exception. visually, errror baloon should appear
  test('negative-openCsvFile', async () => {
    let noException = true;
    const folderName: string = `System:AppData/${_packageName}/samples`;
    const file = (await grok.dapi.files.list(folderName, false, 'dock.csv'))[0];

    try {
      const {view, loadingPromise} = previewMolstarUI(file);
      grok.shell.newView('Molstar Preview', [view]);
      await loadingPromise;
    } catch (e) {
      noException = false;
    }
    expect(noException, false);
  });
});
