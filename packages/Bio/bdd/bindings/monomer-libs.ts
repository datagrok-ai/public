/* The monomer library as Bio holds it: the selection of libraries in the user's settings, what
   the loaded library knows, and the library and collection files on the server that a feature
   creates and removes. Everything is read through Bio's own helper (`Bio:getMonomerLibHelper`). */
import {expect, Page} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {atFeatureEnd} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const LIB_STORAGE = 'Libraries';
const LIB_SETTINGS = 'Settings';

/** Every library on the server selected: the settings Bio keeps in user storage with nothing
 * excluded, and the library reloaded from them. */
async function selectAllLibraries(page: Page): Promise<void> {
  await page.evaluate(async ([storage, key]) => {
    const settings = JSON.parse((await grok.userSettings.getValue(storage, key, true)) || '{}');
    settings.exclude = [];
    settings.explicit = [];
    await grok.userSettings.add(storage, key, JSON.stringify(settings), true);
    const helper = await grok.functions.call('Bio:getMonomerLibHelper', {});
    await helper.loadMonomerLib(true);
  }, [LIB_STORAGE, LIB_SETTINGS]);
}

export const allLibrariesSelected = Given('all monomer libraries are selected', async (page: Page) => {
  await selectAllLibraries(page);
  atFeatureEnd(page, () => selectAllLibraries(page));
}, {tier: 'api', description: 'nothing excluded in the user\'s library settings and the library reloaded — the state the feature starts from and leaves behind, whatever a scenario toggled'});

export const noLibraryOnServer = Given('no {string} monomer library is on the server', async (page: Page, name: string) => {
  await page.evaluate(async (n) => {
    const helper = await grok.functions.call('Bio:getMonomerLibHelper', {});
    for (const provider of await helper.getProviders())
      if ((await provider.listLibraries()).includes(n))
        await provider.deleteLibrary(n);
  }, name);
}, {tier: 'api', description: 'deletes the library file left by an earlier run through its provider'});

export const noCollectionOnServer = Given('no {string} monomer collection is on the server', async (page: Page, name: string) => {
  await page.evaluate(async (n) => {
    const helper = await grok.functions.call('Bio:getMonomerLibHelper', {});
    await helper.deleteMonomerCollection(n);
  }, name);
}, {tier: 'api', description: 'deletes the collection file left by an earlier run'});

/** The sources the loaded library's monomers come from (a file name for the files provider). */
function sources(page: Page): Promise<string[]> {
  return page.evaluate(async () => {
    const lib = (await grok.functions.call('Bio:getMonomerLibHelper', {})).getMonomerLib();
    const out = new Set<string>();
    for (const pt of lib.getPolymerTypes())
      for (const symbol of lib.getMonomerSymbolsByType(pt))
        out.add(String(lib.getMonomer(pt, symbol)?.lib?.source ?? ''));
    out.delete('');
    return [...out].sort();
  });
}

export const loadedFrom = Then('the monomer library should be loaded from {string}', async (page: Page, name: string) => {
  await expect.poll(() => sources(page), {timeout: 30000, message: 'the sources of the loaded monomers'}).toContain(name);
}, {description: 'some monomer of the loaded library comes from that source; polls up to 30 s for a reload in flight'});

export const notLoadedFrom = Then('the monomer library should not be loaded from {string}', async (page: Page, name: string) => {
  await expect.poll(() => sources(page), {timeout: 30000, message: 'the sources of the loaded monomers'}).not.toContain(name);
});

function knows(page: Page, polymerType: string, symbol: string): Promise<boolean> {
  return page.evaluate(async ([pt, s]) => {
    const lib = (await grok.functions.call('Bio:getMonomerLibHelper', {})).getMonomerLib();
    return lib.getMonomer(pt, s) != null;
  }, [polymerType, symbol]);
}

export const knownMonomer = Then('{string} should be a known {string} monomer', async (page: Page, symbol: string, polymerType: string) => {
  await expect.poll(() => knows(page, polymerType, symbol), {timeout: 30000, message: `the library knows ${polymerType} monomer "${symbol}"`}).toBe(true);
}, {description: 'the loaded library resolves the symbol for the polymer type (PEPTIDE, RNA, CHEM)'});

export const unknownMonomer = Then('{string} should not be a known {string} monomer', async (page: Page, symbol: string, polymerType: string) => {
  await expect.poll(() => knows(page, polymerType, symbol), {timeout: 30000, message: `the library knows ${polymerType} monomer "${symbol}"`}).toBe(false);
});

export const collectionHolds = Then('the {string} monomer collection should hold monomers {string}', async (page: Page, name: string, list: string) => {
  const symbols = await page.evaluate(async (n) => {
    const helper = await grok.functions.call('Bio:getMonomerLibHelper', {});
    return (await helper.readMonomerCollection(n)).monomerSymbols;
  }, name);
  expect(symbols, `the monomers of collection "${name}" on the server`).toEqual(list.split(',').map((s) => s.trim()));
}, {description: 'the collection file on the server, its symbols in order, comma-separated'});

export const noSuchCollection = Then('there should be no {string} monomer collection on the server', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(async (n) => {
    const helper = await grok.functions.call('Bio:getMonomerLibHelper', {});
    return (await helper.listMonomerCollections()).map((f: string) => f.replace(/\.json$/i, ''));
  }, name), {message: 'the collections on the server'}).not.toContain(name.replace(/\.json$/i, ''));
});
