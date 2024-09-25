import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {test, after, before, category, expect} from '@datagrok-libraries/utils/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings, setUserLibSettingsForTests
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {expectMonomerLib} from '@datagrok-libraries/bio/src/tests/monomer-lib-tests';


category('monomerLibraries', () => {
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibrariesSettings: any = null;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibrariesSettings = getUserLibSettings();
  });

  after(async () => {
    await setUserLibSettings(userLibrariesSettings);
  });

  test('default', async () => {
    // Clear settings to test default
    await setUserLibSettings({exclude: [], explicit: [], duplicateMonomerPreferences: {}});
    await monomerLibHelper.loadMonomerLib(true); // test defaultLib

    // Currently default monomer lib set is of all files at LIB_PATH (at least HELMCoreLibrary.json)
    const currentMonomerLib = monomerLibHelper.getMonomerLib();
    expect(currentMonomerLib.getPolymerTypes().length > 0, true);
  });

  test('forTests', async () => {
    await setUserLibSettingsForTests();
    await monomerLibHelper.loadMonomerLib(true); // test defaultLib

    // Currently default monomer lib set is of all files at LIB_PATH (at least HELMCoreLibrary.json)
    const currentMonomerLib = monomerLibHelper.getMonomerLib();
    // HELMCoreLibrary.json checks
    expectMonomerLib(currentMonomerLib);
  });

  test('empty', async () => {
    // exclude all monomer libraries for empty set
    const libSettings = await getUserLibSettings();
    const libFileManager = await monomerLibHelper.getFileManager();

    let libFnList = libFileManager.getValidLibraryPaths();
    if (libFnList.length === 0)
      libFnList = await libFileManager.getValidLibraryPathsAsynchronously();

    libSettings.exclude = libFnList;
    libSettings.explicit = [];
    await setUserLibSettings(libSettings);

    await monomerLibHelper.loadMonomerLib(true);
    const currentMonomerLib = monomerLibHelper.getMonomerLib();
    const polymerTypes = currentMonomerLib.getPolymerTypes();
    expect(polymerTypes.length === 0, true);
  });
});
