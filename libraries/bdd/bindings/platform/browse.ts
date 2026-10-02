/* The Browse panel beyond opening and clicking nodes (its Refresh is in steps.ts): the
   tree's own "children are there" state, favorites, and moving inside the app by address. Refresh and Find path fire
   `grok.events.onBrowseTreeRefreshed` once the tree is rebuilt and the path parsed
   (browse_panel.dart); a group fetching its children carries `data-state="loading"` on its host
   (tree_view.dart `loadChildren`). */
import {type Page} from '@playwright/test';
import {expect, pollMs} from '../../src/runtime/patience.js';
import {Given, Then, When} from '../../src/registry.js';
import {atFeatureEnd} from '../../src/runtime/harness.js';

declare const grok: any;
declare const DG: any;

const LOADING = '.grok-view-browse .d4-tree-view-group-host[data-state="loading"], .layout-browse .d4-tree-view-group-host[data-state="loading"]';

async function treeLoaded(page: Page): Promise<void> {
  await expect.poll(() => page.evaluate((sel) => document.querySelectorAll(sel).length, LOADING),
    {message: 'Browse tree groups still fetching their children', timeout: pollMs(30000)}).toBe(0);
}

export const browseTreeLoaded = Then('the browse tree should have finished loading', (page: Page) => treeLoaded(page),
  {description: 'no group of the tree carries data-state="loading" — the end anchor after anything that rebuilds it'});

/* Favorites belong to the account and outlive a feature: a feature that stars an entity takes it out
   of the favorites before it starts, and at its end leaves the favorites as it found them — a fixture
   it starred comes out, an entity the account had starred (the home share) goes back — the list read back. */
async function removeFavorite(page: Page, name: string): Promise<void> {
  await page.evaluate(async (n) => {
    for (const favorite of await grok.dapi.entities.getFavorites())
      if (favorite.friendlyName === n || favorite.name === n)
        await DG.Favorites.remove(favorite);
  }, name);
  await expect.poll(() => favoriteNames(page).then((names) => names.includes(name)),
    {message: `"${name}" among the account's favorites`, timeout: pollMs(30000)}).toBe(false);
}

function favoriteNames(page: Page): Promise<string[]> {
  return page.evaluate(async () => (await grok.dapi.entities.getFavorites())
    .flatMap((f: any) => [f.friendlyName, f.name]).filter((x: any) => !!x));
}

export const notFavorite = Given('{string} is not in favorites', async (page: Page, name: string) => {
  // the starred entities are kept in the page, which the feature keeps, to be starred again at its end
  const had = await page.evaluate(async (n) => {
    const starred = (await grok.dapi.entities.getFavorites()).filter((f: any) => f.friendlyName === n || f.name === n);
    ((window as any).__bddFavorites ??= {})[n] = starred;
    return starred.length > 0;
  }, name);
  atFeatureEnd(page, async () => {
    await removeFavorite(page, name);
    if (!had)
      return;
    await page.evaluate(async (n) => {
      for (const entity of (window as any).__bddFavorites?.[n] ?? [])
        await DG.Favorites.add(entity);
    }, name);
    await expect.poll(() => favoriteNames(page).then((names) => names.includes(name)),
      {message: `"${name}" back among the account's favorites`, timeout: pollMs(30000)}).toBe(true);
  });
  await removeFavorite(page, name);
}, {tier: 'api', description: 'the entity is taken out of the account\'s favorites now; at feature end the favorites are as the feature found them, the list read back'});

export const favoriteOnServer = Then('{string} should be in favorites on the server', async (page: Page, name: string) => {
  await expect.poll(() => favoriteNames(page).then((names) => names.includes(name)),
    {message: `"${name}" among the account's favorites on the server`, timeout: pollMs(30000)}).toBe(true);
}, {tier: 'api', description: 'what the server lists as the account\'s favorites, not what the tree draws'});

export const notFavoriteOnServer = Then('{string} should not be in favorites on the server', async (page: Page, name: string) => {
  await expect.poll(() => favoriteNames(page).then((names) => names.includes(name)),
    {message: `"${name}" among the account's favorites on the server`, timeout: pollMs(30000)}).toBe(false);
}, {tier: 'api'});

export const openAddress = When('user opens the address {string}', async (page: Page, path: string) => {
  await page.evaluate((p) => { grok.shell.route(p); }, path);
}, {tier: 'api', description: 'moves inside the running app to the address, as a link does (grok.shell.route) — no reload'});

/* Recent is what the account touched, read from the audit log, which the server writes after the
   gesture: a claim on the tree right after it reads the list as it was. The server's own list is the
   anchor; the tree is claimed after it. */
export const recentOnServer = Then('{string} should be among the recently used entities on the server', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(async (n) => (await grok.dapi.entities.getRecentEntities())
    .some((e: any) => e?.friendlyName === n || e?.name === n), name),
  {message: `"${name}" among the account's recently used entities`, timeout: pollMs(60000)}).toBe(true);
}, {tier: 'api', description: 'the server\'s Recent list (the audit log of the account), which My stuff > Recent shows when it is opened'});

export const counterMatchesFolder = Then('the gallery counter should show as many items as the {string} folder holds on the server', async (page: Page, path: string) => {
  const count = await page.evaluate(async (p) => (await grok.dapi.files.list(p, false)).length, path);
  if (count === 0)
    throw new Error(`the folder ${path} lists nothing on the server, so the counter claim could not fail`);
  await expect.poll(() => page.evaluate(() => {
    const counters = [...document.querySelectorAll('.grok-items-view-counts')]
      .filter((e) => (e as HTMLElement).offsetWidth > 0).map((e) => (e.textContent ?? '').trim());
    return counters.join(' | ');
  }), {message: `the gallery counter against the ${count} items of ${path}`}).toBe(String(count));
}, {tier: 'ui', description: 'the counter of the folder view equals the number of entries the files API lists in the folder (not recursive)'});

/* "Copy the URL and open it" (browse.md 5): the address the platform wrote for what is open, taken
   and later followed inside the app. */
const rememberedAddress = new WeakMap<Page, string>();

export const rememberAddress = When('user remembers the page address', async (page: Page) => {
  const address = await page.evaluate(() => location.pathname + location.search);
  if (address === '/' || address === '')
    throw new Error('the page address is the root: nothing open has an address to remember');
  if (!rememberedAddress.has(page))
    atFeatureEnd(page, async () => { rememberedAddress.delete(page); });
  rememberedAddress.set(page, address);
}, {tier: 'ui', description: 'the path and query the platform wrote into the address bar for what is open'});

export const openRememberedAddress = When('user opens the remembered address', async (page: Page) => {
  const address = rememberedAddress.get(page);
  if (!address)
    throw new Error('no address remembered: "user remembers the page address" first');
  await page.evaluate((p) => { grok.shell.route(p); }, address);
}, {tier: 'api', description: 'follows the remembered address inside the running app (grok.shell.route), as a pasted link does'});
