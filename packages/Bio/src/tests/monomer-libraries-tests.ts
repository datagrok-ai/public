import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {test, after, before, category, expect, expectObject} from '@datagrok-libraries/utils/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {expectMonomerLib} from '@datagrok-libraries/bio/src/tests/monomer-lib-tests';
import {MonomerTypes, PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {Monomer} from '@datagrok-libraries/bio/src/types/index';

import {monomerLibForTestsSummary} from '../utils/monomer-lib/consts';


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
    await monomerLibHelper.loadMonomerLibForTests(); // test defaultLib

    // Currently default monomer lib set is of all files at LIB_PATH (at least HELMCoreLibrary.json)
    const currentMonomerLib = monomerLibHelper.getMonomerLib();
    // HELMCoreLibrary.json checks
    expectMonomerLib(currentMonomerLib, monomerLibForTestsSummary);
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

  test('override', async () => {
    const overMon: Monomer = {
      symbol: 'over1',
      name: 'Test override monomer 1',
      molfile: '',
      author: 'Test Author',
      id: 0,
      rgroups: [],
      smiles: '',
      polymerType: PolymerTypes.PEPTIDE,
      monomerType: MonomerTypes.BACKBONE,
      createDate: null,
    };
    const monomerLib = monomerLibHelper.getMonomerLib();
    const absentOverrideMonomer = monomerLib.getMonomer(overMon.polymerType, overMon.symbol);
    expect(absentOverrideMonomer === null, true, `Unexpectedly found monomer '${overMon.symbol}' `);

    const overriddenMonomerLib = monomerLib.override({[overMon.polymerType]: {[overMon.symbol]: overMon}}, 'test');
    const resOverMon = overriddenMonomerLib.getMonomer(overMon.polymerType, overMon.symbol);
    if (resOverMon) resOverMon.lib = undefined; // cleanup to prevent infinite recursive comparison
    expectObject(resOverMon as any, overMon);
  });
});
