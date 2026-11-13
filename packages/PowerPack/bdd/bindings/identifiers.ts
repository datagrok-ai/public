import {Page} from '@playwright/test';
import {Given, Then, When, element} from '@datagrok-libraries/bdd';
import {atFeatureEnd, escapeRegExp, expect} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

// every identifiers configuration of the stand is in this one file
const CONFIGS = 'System:AppData/PowerPack/db-explorer/configs.json';

element('identifiers editor', {selector: '.db-explorer-editor',
  description: 'the view Configure Identifiers opens: the Identifiers, Joins and other sections, and Save'});
element('add identifier icon', {selector: '.entry-points-section [aria-label="Add a new identifier"]',
  description: 'the plus under the list of identifiers'});
element('identifier preview', {selector: '.d4-dialog .db-explorer-entrypoint-preview-card',
  description: 'the card in the Add Identifier dialog: one row of the chosen table, read from the database'});
element('configuration text', {selector: '.d4-dialog textarea',
  description: 'the text area of the Import Configuration dialog, which takes a configuration as JSON'});
element('catalog search icon', {selector: '.d4-dialog [name="input-host-Catalog"] .fa-search',
  description: 'the magnifier beside the Catalog field of the Configure Identifiers dialog, shown once the catalogs are listed'});

const covered = (config: any): string => config.catalogName ? `${config.catalogName}.${config.schemaName}` : String(config.schemaName);

/** The "Namespace:Name" of a connection named in full or by its name alone: the namespace of a connection a
 * person made is that person's, so a feature that names it would run on one account's stand only. */
function qualifiedName(page: Page, connection: string): Promise<string> {
  if (connection.includes(':'))
    return Promise.resolve(connection);
  return page.evaluate(async (name) => {
    let found = (await grok.dapi.connections.list({pageSize: 5000})).filter((c: any) => c.friendlyName === name || c.name === name);
    if (found.length > 1) {
      const me = (await grok.dapi.users.current()).id;
      found = found.filter((c: any) => c.author?.id === me);
    }
    if (found.length !== 1)
      throw new Error(`${found.length} connections named "${name}"`);
    return String(found[0].nqName);
  }, connection);
}

function savedConfigs(page: Page, connection: string): Promise<any[]> {
  return page.evaluate(async ([path, nqName]) => {
    if (!await grok.dapi.files.exists(path))
      return [];
    return JSON.parse(await grok.dapi.files.readAsText(path)).filter((c: any) => c.nqName === nqName);
  }, [CONFIGS, connection] as [string, string]);
}

/** Takes the connection's configurations out of the file and returns them. A file left with none goes, and its
 * folder with it: a stand nobody configured identifiers on has neither, and Remove leaves an empty list behind. */
function takeConfigs(page: Page, connection: string): Promise<any[]> {
  return page.evaluate(async ([path, nqName]) => {
    if (!await grok.dapi.files.exists(path))
      return [];
    const all = JSON.parse(await grok.dapi.files.readAsText(path));
    const taken = all.filter((c: any) => c.nqName === nqName);
    const kept = all.filter((c: any) => c.nqName !== nqName);
    if (kept.length > 0) {
      if (taken.length > 0)
        await grok.dapi.files.writeAsText(path, JSON.stringify(kept, null, 2));
      return taken;
    }
    await grok.dapi.files.delete(path);
    const folder = path.slice(0, path.lastIndexOf('/'));
    if ((await grok.dapi.files.list(folder, false)).length === 0)
      await grok.dapi.files.delete(folder);
    return taken;
  }, [CONFIGS, connection] as [string, string]);
}

function putConfigs(page: Page, configs: any[]): Promise<void> {
  return page.evaluate(async ([path, back]) => {
    if (back.length === 0)
      return;
    const all = await grok.dapi.files.exists(path) ? JSON.parse(await grok.dapi.files.readAsText(path)) : [];
    await grok.dapi.files.writeAsText(path, JSON.stringify([...all, ...back], null, 2));
  }, [CONFIGS, configs] as [string, any[]]);
}

export const noIdentifiersConfig = Given('no identifiers configuration of the {string} connection is on the server', async (page: Page, connection: string) => {
  const nqName = await qualifiedName(page, connection);
  const setAside = await takeConfigs(page, nqName);
  await expect.poll(async () => (await savedConfigs(page, nqName)).map(covered),
    {message: `the identifiers configurations of ${nqName} still on the server`}).toEqual([]);
  atFeatureEnd(page, async () => {
    await takeConfigs(page, nqName);
    await putConfigs(page, setAside);
    await expect.poll(async () => (await savedConfigs(page, nqName)).map(covered).sort(),
      {message: `the identifiers configurations of ${nqName} after the feature`}).toEqual(setAside.map(covered).sort());
  });
}, {tier: 'api', description: 'the configurations saved for the connection ("System:Datagrok", or a name alone for one a person made) are set aside now; when the feature ends the ones it made go and those come back. The page keeps the handlers it registered until it reloads'});

export const savedIdentifiers = Then('the saved identifiers of {string} in the {string} connection should be {string}', async (page: Page, target: string, connection: string, list: string) => {
  const expected = list.split(',').map((t) => t.trim()).filter(Boolean).sort();
  const nqName = await qualifiedName(page, connection);
  await expect.poll(async () => {
    const config = (await savedConfigs(page, nqName)).find((c) => covered(c) === target);
    return config ? Object.keys(config.entryPoints ?? {}).sort() : [`no configuration for ${target}`];
  }, {message: `the semantic types the configuration of ${target} defines, as saved`}).toEqual(expected);
}, {tier: 'api', description: 'the semantic types (comma-separated) of the configuration saved for "catalog.schema" or a schema: all of them, no other'});

export const identifiersConfigCovers = Then('the identifiers configurations of the {string} connection on the server should cover {string}', async (page: Page, connection: string, targets: string) => {
  const expected = targets.split(',').map((t) => t.trim()).filter(Boolean).sort();
  const nqName = await qualifiedName(page, connection);
  await expect.poll(async () => (await savedConfigs(page, nqName)).map(covered).sort(),
    {message: `what the identifiers configurations of ${nqName} cover, as saved`}).toEqual(expected);
}, {tier: 'api', description: 'the saved configurations of the connection, each as "catalog.schema" or the schema alone (comma-separated): all of them, no other. A configuration of the connection\'s own catalog names no catalog'});

export const findCatalog = When('user finds the catalog {string} by typing {string}', async (page: Page, catalog: string, part: string) => {
  const search = page.locator('.d4-menu-popup input').first();
  await search.waitFor({state: 'visible'});
  await search.pressSequentially(part);
  await page.locator('.d4-menu-popup .d4-menu-item-label:visible').filter({hasText: new RegExp(`^${escapeRegExp(catalog)}$`)}).click();
}, {tier: 'ui', description: 'in the list the catalog search icon opens: the part typed into its Search box, then a click on the catalog left in the list'});
