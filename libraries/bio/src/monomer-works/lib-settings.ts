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

// user settings does not allow more than 5000 characters, but unfortunately we need this, so we will split
// the string into chunks and store them separately
const MAX_USER_SETTINGS_STRING_LENGTH = 4000;
const CHUNK_KEY_SUFFIX = 'nextChunk';
async function getDuplicatePrefsFull(): Promise<DuplicateMonomerPreferencesShortHand> {
// the first part of it will always live in DUPLICATE_MONOMER_PREFERENCES_KEY
  let resStr: string = await grok.userSettings.getValue(LIB_STORAGE_NAME, DUPLICATE_MONOMER_PREFERENCES_KEY, true) ?? '';
  if (!resStr)
    return {};
  // chunk will end with special suffix and just before it will have the number of the next chunk key
  while (resStr.endsWith(CHUNK_KEY_SUFFIX)) {
    const nextChunkKey = `${DUPLICATE_MONOMER_PREFERENCES_KEY}${resStr.charAt(resStr.length - CHUNK_KEY_SUFFIX.length - 1)}`;
    const nextChunkStr: string | undefined = (await grok.userSettings.getValue(LIB_STORAGE_NAME, nextChunkKey, true)) ?? '';
    resStr = resStr.slice(0, resStr.length - CHUNK_KEY_SUFFIX.length - 1) + nextChunkStr;
  }
  return JSON.parse(resStr);
}

async function setDuplicatePrefsFull(value: DuplicateMonomerPreferencesShortHand): Promise<void> {
  const valueStr: string = JSON.stringify(value);
  const chunksCount = Math.ceil(valueStr.length / MAX_USER_SETTINGS_STRING_LENGTH);
  for (let i = 0; i < chunksCount; i++) {
    const key = i === 0 ? DUPLICATE_MONOMER_PREFERENCES_KEY : `${DUPLICATE_MONOMER_PREFERENCES_KEY}${i}`;
    const nextKey = i < chunksCount - 1 ? `${i + 1}` : '';
    const chunkStr = valueStr.slice(i * MAX_USER_SETTINGS_STRING_LENGTH, (i + 1) * MAX_USER_SETTINGS_STRING_LENGTH);
    const chunkStrToStore = nextKey ? chunkStr + nextKey + CHUNK_KEY_SUFFIX : chunkStr;
    await grok.userSettings.add(LIB_STORAGE_NAME, key, chunkStrToStore, true);
  }
}

let userLibSettingsPromise: Promise<void> = Promise.resolve();

export async function getUserLibSettings(): Promise<UserLibSettings> {
  let res: UserLibSettings;
  userLibSettingsPromise = userLibSettingsPromise.then(async () => {
    const resStr: string | undefined = await grok.userSettings.getValue(LIB_STORAGE_NAME, LIB_SETTINGS_KEY, true);
    res = resStr ? JSON.parse(resStr) : {exclude: [], explicit: [], duplicateMonomerPreferences: {}};
    const duplicatePrefs: DuplicateMonomerPreferencesShortHand = await getDuplicatePrefsFull();
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
    // console.debug(`Bio: getUserLibSettings()\n${JSON.stringify(res, undefined, 2)}`);
  });
  await userLibSettingsPromise;
  return res!;
}

export async function setUserLibSettings(value: UserLibSettings): Promise<void> {
  userLibSettingsPromise = userLibSettingsPromise.then(async () => {
    value.duplicateMonomerPreferences = value.duplicateMonomerPreferences ?? {};
    // console.debug(`Bio: setUserLibSettings()\n${JSON.stringify(value, undefined, 2)}`);
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
    await setDuplicatePrefsFull(duplicatePrefsShortHand);
    //await grok.userSettings.add(LIB_STORAGE_NAME, DUPLICATE_MONOMER_PREFERENCES_KEY, JSON.stringify(duplicatePrefsShortHand), true);
    value.duplicateMonomerPreferences = {}; // clear to reduce size
    await grok.userSettings.add(LIB_STORAGE_NAME, LIB_SETTINGS_KEY, JSON.stringify(value), true);
  });
  await userLibSettingsPromise;
}
