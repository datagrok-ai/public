import * as grok from 'datagrok-api/grok';
import {UserLibSettings} from './types';

// -- Monomer libraries --
export const LIB_STORAGE_NAME = 'Libraries';
export const LIB_PATH = 'System:AppData/Bio/monomer-libraries/';
const LIB_SETTINGS_FOR_TESTS: UserLibSettings =
  {explicit: ['HELMCoreLibrary.json'], exclude: [], duplicateMonomerPreferences: {}};

export const SETS_STORAGE_NAME: string = 'Monomer Sets';
export const SETS_PATH: string = 'System:AppData/Bio/monomer-sets/';

let userLibSettingsPromise: Promise<void> = Promise.resolve();

export async function getUserLibSettings(): Promise<UserLibSettings> {
  let res: UserLibSettings;
  userLibSettingsPromise = userLibSettingsPromise.then(async () => {
    const resStr: string = await grok.dapi.userDataStorage.getValue(LIB_STORAGE_NAME, 'Settings', true);
    res = resStr ? JSON.parse(resStr) : {exclude: [], explicit: [], duplicateMonomerPreferences: {}};

    // Fix empty object returned in case there is no settings stored for user
    res.exclude = res.exclude instanceof Array ? res.exclude : [];
    res.explicit = res.explicit instanceof Array ? res.explicit : [];
    res.duplicateMonomerPreferences = res.duplicateMonomerPreferences instanceof Object ?
      res.duplicateMonomerPreferences : {};
    console.debug(`Bio: getUserLibSettings()\n${JSON.stringify(res, undefined, 2)}`);
  });
  await userLibSettingsPromise;
  return res!;
}

export async function setUserLibSettings(value: UserLibSettings): Promise<void> {
  userLibSettingsPromise = userLibSettingsPromise.then(async () => {
    console.debug(`Bio: setUserLibSettings()\n${JSON.stringify(value, undefined, 2)}`);
    await grok.dapi.userDataStorage.postValue(LIB_STORAGE_NAME, 'Settings', JSON.stringify(value), true);
  });
  await userLibSettingsPromise;
}

/** Set only HELMCoreLibrary.json */
export async function setUserLibSettingsForTests(): Promise<void> {
  await setUserLibSettings(LIB_SETTINGS_FOR_TESTS);
}
