/* Diff Studio's own vocabulary. Everything else a feature here needs is platform vocabulary: the
   model's charts are real viewers on a table view (a Grid and a Line chart on a table named after
   the model), so the `viewers` tier drives them, and the Open model menu carries the platform's own
   item names (div-Library---Bioreactor), so its entries are plain menu items. */
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, savedScriptOf, takeErrors} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* The app's deep link per model of the library: /apps/DiffStudio/Library/<state>. The titles are
   what the Open model menu shows; the states are what the route takes. */
const LIBRARY: Record<string, string> = {
  'Chem reactions': 'chem-react',
  "Robertson's model": 'robertson',
  'Lotka-Volterra': 'lotka-volterra',
  'Fermentation': 'fermentation',
  'PK': 'pk',
  'PK-PD': 'pk-pd',
  'Acid production': 'ga-production',
  'Nimotuzumab': 'nimotuzumab',
  'Bioreactor': 'bioreactor',
  'Pollution': 'pollution',
};

/* Opening the app itself lands on its hub, which never builds the model UI, so a feature that is
   about a model asks for the model. The step is over when the app says it is ready: its ribbon
   widget is up and the platform's current view is the model (~7 s on dev, 2026-09-10), and the
   inputs the model declares are on the page. */
export const openLibraryModel = Given('user opens the {string} model of the Diff Studio library',
  async (page: Page, title: string) => {
    const state = LIBRARY[title];
    if (!state)
      throw new Error(`no "${title}" in the Diff Studio library; it has: ${Object.keys(LIBRARY).join(', ')}`);
    await page.goto(`/apps/DiffStudio/Library/${state}?params:`, {waitUntil: 'domcontentloaded', timeout: 180000});
    // what the document being left logs as its fetches are aborted, and what the shell logs while
    // booting again, is not the scenario's — the same floor the login step sets
    takeErrors(page);
    await expect(page.locator('.diff-studio-ribbon-widget').first(), 'the Diff Studio ribbon')
      .toBeVisible({timeout: 180000});
    await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')),
      {message: 'the current view', timeout: 60000}).toBe(title);
    await expect.poll(() => page.locator('[name^="input-host-"]').count(),
      {message: 'the inputs of the model', timeout: 60000}).toBeGreaterThan(2);
  }, {tier: 'ui', description: 'the app\'s own route for a library model; done when the ribbon, the view name and the model inputs are all there'});

/* Saving writes a new .ivp into System:AppData/DiffStudio/library under a name the app picks
   ("PK-PD(14).ivp" — the count is how many earlier runs left theirs behind), and announces it on
   the platform's event bus. The step notes what the folder held first and deletes exactly what the
   save added when the feature ends, so a run leaves the library as it found it. */
export const saveToLibrary = When('user saves the model to the Diff Studio library', async (page: Page) => {
  const folder = 'System:AppData/DiffStudio/library';
  const before: string[] = await page.evaluate(async (f) =>
    (await grok.dapi.files.list(f)).map((x: any) => String(x.name)), folder);
  await page.locator('.diff-studio-ribbon-save-to-model-catalog-icon').first().click();
  await expect.poll(async () => (await page.evaluate(async (f) =>
    (await grok.dapi.files.list(f)).map((x: any) => String(x.name)), folder)).length,
  {message: 'files in the Diff Studio library after the save', timeout: 60000}).toBeGreaterThan(before.length);
  atFeatureEnd(page, async () => {
    await page.evaluate(async ([f, known]) => {
      for (const file of await grok.dapi.files.list(f))
        if (!(known as string[]).includes(String(file.name)))
          await grok.dapi.files.delete(`${f}/${file.name}`).catch(() => undefined);
    }, [folder, before] as [string, string[]]);
  });
}, {tier: 'ui', description: 'the ribbon icon; the file it creates is deleted when the feature ends'});

/* The Model Hub is Compute2's catalog view, not a plain #app of the registry, so "user opens the
   … app" does not find it — the browse tree node runs Compute2:modelCatalog, and so does this. */
export const openModelHub = Given('user opens the Model Hub', async (page: Page) => {
  // A gallery appearing under the old pointer can open a tooltip over the next card.
  await page.mouse.move(0, 0);
  await page.evaluate(async () => {
    const view = await grok.functions.call('Compute2:modelCatalog', {});
    if (view?.root && !Array.from(grok.shell.views).some((v: any) => v.dart === view.dart))
      grok.shell.addView(view);
  });
  // the gallery is on the page, visible and EMPTY, for about six seconds after the call returns
  // (measured on dev, 2026-09-10), so waiting for the element is waiting for nothing
  await expect.poll(() => page.evaluate(() =>
    grok.shell.v?.root?.querySelectorAll('.grok-gallery-grid .d4-link-label').length ?? 0),
  {message: 'cards in the Model Hub gallery', timeout: 120000}).toBeGreaterThan(0);
}, {tier: 'api', description: 'the function the browse tree runs for Apps > Compute > Model Hub; done when the catalog has cards'});

/** A script's friendly name can also label a built-in model, so the card is found by the link the gallery
 * gives it, the saved script's grok name; the library's "user saves the script" remembers the script and
 * deletes it when the feature ends. */
function savedCard(page: Page) {
  const link = `/func/${savedScriptOf(page).nqName.replace(/:/g, '.')}`;
  return page.locator(`.grok-gallery-grid .d4-link-label[data-link=${JSON.stringify(link)}]`);
}

export const hubListsScript = Then('the Model Hub should list the saved script', async (page: Page) => {
  await expect(savedCard(page), `the card of the saved script (${savedScriptOf(page).name}) in the Model Hub`)
    .toBeVisible({timeout: 60000});
}, {tier: 'ui'});

export const openSavedScript = When('user opens the saved script from the Model Hub', async (page: Page) => {
  await savedCard(page).dblclick();
  await expect.poll(() => page.locator('[name^="input-host-"]').count(),
    {message: 'the inputs of the script the Model Hub opened', timeout: 120000}).toBeGreaterThan(0);
}, {tier: 'ui'});

/** Refresh must notice the deletion of this feature's script; same-named models stay intact. */
export const deleteSavedScript = When('the saved script is deleted on the server', async (page: Page) => {
  const {id, name} = savedScriptOf(page);
  await page.evaluate(async (savedId) => {
    const script = await grok.dapi.scripts.find(savedId).catch(() => null);
    if (script)
      await grok.dapi.scripts.delete(script);
  }, id);
  await expect.poll(() => page.evaluate(async (savedId) => (await grok.dapi.scripts.find(savedId).catch(() => null)) != null, id),
  {message: `"${name}" (${id}) among the scripts on the stand`, timeout: 60000}).toBe(false);
}, {tier: 'api', description: 'removed behind the back of the view, so the next Refresh has something to notice'});

export const hubDoesNotListScript = Then('the Model Hub should not list the saved script', async (page: Page) => {
  await expect(savedCard(page), `the card of the saved script (${savedScriptOf(page).name}) in the Model Hub`)
    .toHaveCount(0, {timeout: 60000});
}, {tier: 'ui'});

/** The same address again, from scratch — what pasting the link into a second tab does, minus the
 * tab. The platform's own step ends when the navigation does, and a model rebuilt from its address
 * takes longer than that: the claim after it used to read a page with no inputs on it yet. */
export const reopenModelAddress = When('user opens the model at the page address', async (page: Page) => {
  await page.goto(page.url(), {waitUntil: 'domcontentloaded', timeout: 180000});
  takeErrors(page);
  await expect.poll(() => page.locator('[name^="input-host-"]').count(),
    {message: 'the inputs of the model the address names', timeout: 180000}).toBeGreaterThan(2);
}, {tier: 'ui', description: 'done when the model the address names has its inputs on the page'});
