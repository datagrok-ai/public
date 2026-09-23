/* The MPO Profiles app: what the mpo domain table holds, and the gestures the profile editor needs
   (its screen parts are registered in elements.ts). A profile a feature saves is removed when the feature ends, with the pMPO
   model file a data-driven save writes beside it. */
import {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, expect, gestures, locate} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const PMPO_FOLDER = 'System:AppData/EDA/pmpo';

/** The name the pMPO model file gets (`generateMpoFileName` in the statistics library). */
const pmpoFileOf = (name: string): string =>
  `${PMPO_FOLDER}/${name.trim().replace(/[\s/\\]+/g, '-').replace(/^-|-$/g, '') || 'profile'}.json`;

const namesOf = (list: string): string[] => list.split(',').map((n) => n.trim()).filter(Boolean);

async function deleteProfiles(page: Page, names: string[]): Promise<void> {
  await page.evaluate(async ([profiles, files]) => {
    const table = grok.dapi.domains.table('mpo.profile');
    for (const row of await table.query().top(1000))
      if (profiles.includes(row.name))
        await table.delete(row.id);
    for (const file of files)
      if (await grok.dapi.files.exists(file))
        await grok.dapi.files.delete(file);
  }, [names, names.map(pmpoFileOf)] as [string[], string[]]);
}

/** The ids of the current user's profiles when the feature started: at its end every profile of
 * theirs not among them goes, whatever name it was saved under (a regression that saves under the
 * dataset name or "Untitled Profile" would otherwise leave it behind and block the next run). */
const before = new WeakMap<Page, Set<string>>();

async function ownProfiles(page: Page): Promise<{id: string; name: string}[]> {
  return page.evaluate(async () => (await grok.dapi.domains.table('mpo.profile').query().top(1000))
    .filter((row: any) => row.author_id === grok.shell.user.id).map((row: any) => ({id: row.id, name: row.name})));
}

async function deleteNewProfiles(page: Page): Promise<void> {
  const known = before.get(page);
  before.delete(page);
  if (!known)
    return;
  const created = (await ownProfiles(page)).filter((p) => !known.has(p.id));
  await deleteProfiles(page, created.map((p) => p.name));
}

const profileCount = (page: Page, name: string): Promise<number> => page.evaluate(async (n) =>
  (await grok.dapi.domains.table('mpo.profile').query().top(1000)).filter((row: any) => row.name === n).length, name);

export const noMpoProfile = Given('no MPO profile named {string} is on the server', async (page: Page, list: string) => {
  const names = namesOf(list);
  await deleteProfiles(page, names);
  if (!before.has(page)) {
    before.set(page, new Set((await ownProfiles(page)).map((p) => p.id)));
    atFeatureEnd(page, () => deleteNewProfiles(page));
  }
}, {tier: 'api', description: 'deletes the profiles (comma-separated) and their pMPO model files; when the feature ends, every profile the user saved during it goes too, with its model file'});

export const shippedMpoProfile = Given('the shipped MPO profile {string} is on the server', async (page: Page, name: string) => {
  if (await profileCount(page, name) === 0)
    await page.evaluate(() => grok.functions.call('Chem:seedMpoProfiles'));
  await expect.poll(() => profileCount(page, name), {message: `MPO profiles named "${name}" in the mpo domain table`}).toBe(1);
}, {tier: 'api', description: 'a profile of the package\'s mpo folder; on a stand that was never seeded, Chem:seedMpoProfiles (idempotent) adds the shipped ones, which stay'});

export const mpoProfileProperties = Then('the MPO profile {string} should have the properties {string}',
  (page: Page, name: string, list: string) =>
    expect.poll(() => page.evaluate(async (n) => {
      const rows = (await grok.dapi.domains.table('mpo.profile').query().top(1000)).filter((row: any) => row.name === n);
      return rows.length === 1 ? Object.keys(rows[0].properties ?? {}).join(', ') : `${rows.length} profiles named "${n}"`;
    }, name), {message: `the stored properties of the MPO profile "${name}"`}).toBe(namesOf(list).join(', ')),
{tier: 'api', description: 'the property names stored with the profile, in order (comma-separated)'});

export const mpoProfilesOnServer = Then('{int} MPO profile(s) named {string} should be on the server',
  (page: Page, count: number, name: string) =>
    expect.poll(() => profileCount(page, name), {message: `MPO profiles named "${name}" in the mpo domain table`}).toBe(count),
{tier: 'api', description: 'what the mpo domain table holds, not what the list draws'});

export const mpoProfileDescription = Then('the MPO profile {string} should have the description {string}',
  (page: Page, name: string, description: string) =>
    expect.poll(() => page.evaluate(async (n) => {
      const rows = (await grok.dapi.domains.table('mpo.profile').query().top(1000)).filter((row: any) => row.name === n);
      return rows.length === 1 ? rows[0].description ?? '' : `${rows.length} profiles named "${n}"`;
    }, name), {message: `the stored description of the MPO profile "${name}"`}).toBe(description),
{tier: 'api', description: 'the description stored with the profile'});

export const pmpoModelFile = Then('the pMPO model file of {string} should hold its name and the description {string}',
  (page: Page, name: string, description: string) =>
    expect.poll(() => page.evaluate(async (f) => {
      if (!await grok.dapi.files.exists(f))
        return 'no file';
      const model = JSON.parse(await grok.dapi.files.readAsText(f));
      return `${model.name} | ${model.description}`;
    }, pmpoFileOf(name)), {message: `the name and description in the pMPO model file ${pmpoFileOf(name)}`}).toBe(`${name} | ${description}`),
{tier: 'api', description: `the model file a data-driven save writes to ${PMPO_FOLDER}, named after the profile and holding its name and description`});

/** One pass at a person's pace, never retyped: the library's typing retypes until the field holds
 * the text, which a field that drops keys only on its first rename would survive. */
export const typeKeyByKey = When('user types {string} key by key into {element}', async (page: Page, text: string, target: ElementRef) => {
  const editor = await gestures.editorOf(page, target);
  await editor.click();
  await editor.press('ControlOrMeta+A');
  await page.keyboard.type(text, {delay: 120});
}, {tier: 'ui', description: 'clicks the field, selects its text and types once at 120 ms a key; the field is not read back here'});

/** A triple click: the profile title swallows the select-all key, so the text a user retypes is
 * selected with the pointer. */
export const selectTextOf = When('user selects the text of {element}', async (page: Page, target: ElementRef) => {
  const loc = await locate(page, target);
  await loc.click({clickCount: 3});
  await expect.poll(() => page.evaluate(() => window.getSelection()?.toString() ?? ''),
    {message: `the selected text of ${target.phrase}`}).toBe((await loc.innerText()).trim());
}, {tier: 'ui', description: 'a triple click on the element, checked by what the page reports as selected'});
