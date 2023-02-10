import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {test, after, before, category, expect} from '@datagrok-libraries/utils/src/test';

import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {LIB_STORAGE_NAME} from '../utils/monomer-lib';


category('monomerLibraries', () => {
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibrariesSettings: any = null;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibrariesSettings = await grok.dapi.userDataStorage.get(LIB_STORAGE_NAME, true);
  });

  after(async () => {
    await grok.dapi.userDataStorage.put(LIB_STORAGE_NAME, userLibrariesSettings, true);
  });

  test('default', async () => {
    // Clear settings to test default
    await grok.dapi.userDataStorage.put(LIB_STORAGE_NAME, {}, true);
    await monomerLibHelper.loadLibraries(true); // test defaultLib

    // Currently default monomer lib is empty
    const currentMonomerLib = monomerLibHelper.getBioLib();
    expect(currentMonomerLib.getTypes().length, 0);
  });
});
