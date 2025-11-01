/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import {UserLibSettings} from './types';

// -- Monomer libraries --
const LIB_STORAGE_NAME = 'Libraries';
const LIB_SETTINGS_KEY = 'Settings';
const DUPLICATE_MONOMER_PREFERENCES_KEY = 'DuplicateMonomerPreferences';

// note: as the settings has limit of 5000 characters, we will store the settings
// in shorter JSON format at least for duplicateMonomerPreferences
type DuplicateMonomerPreferencesShortHand = {
  [polymerType: string]: {
    [libraryName: string]: string[] // list of monomer symbols to prefer from the library for the polymer type
  }
}

let userLibSettingsPromise: Promise<void> = Promise.resolve();

export async function getUserLibSettings(): Promise<UserLibSettings> {
  let res: UserLibSettings;
  userLibSettingsPromise = userLibSettingsPromise.then(async () => {
    const resStr: string | undefined = await grok.userSettings.getValue(LIB_STORAGE_NAME, LIB_SETTINGS_KEY, true);
    res = resStr ? JSON.parse(resStr) : {exclude: [], explicit: [], duplicateMonomerPreferences: {}};
    const duplicatePrefsStr = await grok.userSettings.getValue(LIB_STORAGE_NAME, DUPLICATE_MONOMER_PREFERENCES_KEY, true);
    const duplicatePrefs: DuplicateMonomerPreferencesShortHand = duplicatePrefsStr ? JSON.parse(duplicatePrefsStr) : {};
    // Fix empty object returned in case there is no settings stored for user
    res.exclude = res.exclude instanceof Array ? res.exclude : [];
    res.explicit = res.explicit instanceof Array ? res.explicit : [];
    res.duplicateMonomerPreferences = res.duplicateMonomerPreferences instanceof Object ?
      res.duplicateMonomerPreferences : {};

    // migrate from short hand if needed
    for (const polymerType in duplicatePrefs) {
      res.duplicateMonomerPreferences[polymerType] = res.duplicateMonomerPreferences[polymerType] ?? {};
      for (const libraryName in duplicatePrefs[polymerType]) {
        for (const monomerSymbol of duplicatePrefs[polymerType][libraryName])
          res.duplicateMonomerPreferences[polymerType][monomerSymbol] = libraryName;
      }
    }
    console.debug(`Bio: getUserLibSettings()\n${JSON.stringify(res, undefined, 2)}`);
  });
  await userLibSettingsPromise;
  return res!;
}

export async function setUserLibSettings(value: UserLibSettings): Promise<void> {
  userLibSettingsPromise = userLibSettingsPromise.then(async () => {
    value.duplicateMonomerPreferences = value.duplicateMonomerPreferences ?? {};
    console.debug(`Bio: setUserLibSettings()\n${JSON.stringify(value, undefined, 2)}`);
    // extract short hand for duplicateMonomerPreferences
    const duplicatePrefsShortHand: DuplicateMonomerPreferencesShortHand = {};
    for (const polymerType in value.duplicateMonomerPreferences) {
      duplicatePrefsShortHand[polymerType] = {};
      for (const monomerSymbol in value.duplicateMonomerPreferences[polymerType]) {
        const libraryName = value.duplicateMonomerPreferences[polymerType][monomerSymbol];
        duplicatePrefsShortHand[polymerType][libraryName] = duplicatePrefsShortHand[polymerType][libraryName] ?? [];
        duplicatePrefsShortHand[polymerType][libraryName].push(monomerSymbol);
      }
    }
    await grok.userSettings.add(LIB_STORAGE_NAME, DUPLICATE_MONOMER_PREFERENCES_KEY, JSON.stringify(duplicatePrefsShortHand), true);
    value.duplicateMonomerPreferences = {}; // clear to reduce size
    await grok.userSettings.add(LIB_STORAGE_NAME, LIB_SETTINGS_KEY, JSON.stringify(value), true);
  });
  await userLibSettingsPromise;
}
