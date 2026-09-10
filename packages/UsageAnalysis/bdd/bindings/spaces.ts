/* Spaces as the Browse tree shows them. Everything general — the browse panel, the browse tree,
   the expanded state of a node — is platform vocabulary in the library; what stays here is the
   space view's own gallery and search, and the spaces a scenario leaves on the server. */
import {Locator, Page} from '@playwright/test';
import {element, Given, Then} from '@datagrok-libraries/bdd';
import {el, exactText, expect, gestures, locate} from '@datagrok-libraries/bdd/runtime';
import {atFeatureEnd} from '@datagrok-libraries/bdd/runtime';
import {sharingLogin} from '@datagrok-libraries/bdd/bindings/platform/steps';

declare const grok: any;

/* The gallery a space view shows, and its search — the same two elements in the Spaces list and
   inside a space. The cards themselves are links, so "X link in space gallery" names one without
   catching the identically-classed links of an open help pane. */
// the plain "gallery" is the platform's own element now (bindings/platform/elements.ts) — the same
// selector, so "X link in gallery" and "X link in space gallery" name the same cards
element('space gallery', {selector: '.grok-gallery-grid', aliases: ['space content']});
element('space search', {selector: '.grok-gallery-search-bar .ui-input-type-ahead'});
/* The access level in the sharing dialog is not a <select> any more (it was in the August spec):
   a Dart privilege selector with its current level as text and a popup behind the triangle. */
element('share access selector', {selector: '[name="div-share-selector"]'});

/* grok.dapi.spaces.filter('name = "…"') answers nothing on a stand that holds the space (probed
   2026-09-10 against dev, by grok name, friendly name and both), so a space is found in the list. */
async function deleteSpaces(page: Page, names: string[]): Promise<void> {
  const remaining = () => page.evaluate(async (wanted) => {
    const left: string[] = [];
    for (const space of await grok.dapi.spaces.list({pageSize: 1000}))
      if (wanted.includes(space.friendlyName) || wanted.includes(space.name)) {
        await grok.dapi.spaces.delete(space).catch(() => undefined);
        left.push(space.friendlyName ?? space.name);
      }
    return left;
  }, names);
  // the delete returns before the space is gone, and creating the same name meanwhile is refused as
  // a duplicate — so the step is over only once the server stops listing them
  await expect.poll(remaining, {message: `spaces still on the server under ${names.join(', ')}`, timeout: 60000}).toEqual([]);
}

export const noSpaceOnServer = Given('no space named {string} is on the server', async (page: Page, name: string) => {
  const names = name.split(',').map((n) => n.trim()).filter(Boolean);
  await deleteSpaces(page, names);
  atFeatureEnd(page, () => deleteSpaces(page, names));
}, {tier: 'api', description: 'deletes what an earlier run left under those names (comma-separated), and deletes them again when the feature ends'});

/* A space is listed once its save returns, and the save of a ROOT space is slow: 4.8 s alone and
   18 s with four features creating at once on a local stand (2026-09-10, POST /api/spaces in the
   trace's network log; a second root space in the same feature takes under half a second, so the
   server does one-time or serialized work on the first). The claim right after OK owns the same
   budget the dialog-close claim below does, or it fails while the dialog is still legitimately open. */
export const spacesOnServer = Then('{int} space(s) named {string} should be on the server', async (page: Page, count: number, name: string) => {
  await expect.poll(() => page.evaluate(async (n) =>
    (await grok.dapi.spaces.list({pageSize: 1000})).filter((s: any) => s.friendlyName === n || s.name === n).length,
  name), {message: `spaces the server holds under "${name}"`, timeout: 60000}).toBe(count);
}, {tier: 'api', description: 'what the server holds, not what the tree draws — the refusal of a duplicate is a space that was never created'});

/* Whom a space is shared with, read where the platform shows it. grok.dapi.permissions.get answers
   with the edit and view buckets only, and a share made through the dialog lands in neither — the
   Sharing pane calls it "has special permissions" — so the API cannot see it and the pane is the
   claim. The user appears there under the punctuation-stripped login ("a+b@x" as "abx"). Two
   expressions rather than one with "(not )": an optional literal is not a parameter, so a single
   step would always take the positive branch. */
async function sharingPane(page: Page): Promise<{pane: ReturnType<Page['locator']>; shown: RegExp}> {
  const login = sharingLogin();
  const header = page.locator('.grok-prop-panel [name="div-section--Sharing"]').first();
  await expect(header, 'the Sharing pane of the context panel').toBeVisible({timeout: 30000});
  if (await header.getAttribute('aria-expanded') !== 'true')
    await header.click();
  return {
    pane: page.locator('.grok-prop-panel .d4-pane-sharing').first(),
    shown: new RegExp(login.split('@')[0].replace(/[^a-z0-9]/gi, ''), 'i'),
  };
}

export const sharingPaneLists = Then('the sharing pane should list the sharing user', async (page: Page) => {
  const {pane, shown} = await sharingPane(page);
  await expect(pane).toContainText(shown);
}, {tier: 'ui', description: 'the Sharing section of the context panel, opened if it is closed'});

export const sharingPaneListsNot = Then('the sharing pane should not list the sharing user', async (page: Page) => {
  const {pane, shown} = await sharingPane(page);
  await expect(pane).not.toContainText(shown);
}, {tier: 'ui'});

/* Creating a space keeps the Create Space dialog open until the platform is done, and for a CHILD
   space that took 6-18 s on dev (2026-09-10) — the shared 15 s expect timeout sits inside that
   range, so the generic "should be hidden" passed or failed by luck. This claim owns its budget. */
export const createDialogCloses = Then('the Create Space dialog should close', async (page: Page) => {
  await expect(page.locator('.d4-dialog[name="dialog-Create-Space"]').filter({visible: true}),
    'the Create Space dialog').toHaveCount(0, {timeout: 60000});
}, {tier: 'ui', description: 'the platform closes it when the space is actually created'});

/* --- what the tree and a space view show -------------------------------------------------------
   Both refresh on a server round-trip the platform does not announce, and the tree's Spaces group
   is closed unless something opened it: the product reveals a space it has just created by opening
   the group, which on a loaded stand it sometimes does not do at all. So these four claims open the
   group themselves on every attempt — what a person does — and give the round-trip a minute. The
   negative pair needs the group open just as much: "absent" was true of a closed group whatever the
   server held (the Dart tree builds a group's children when it opens). */
const SPACES_GROUP = 'Spaces tree node inside browse tree';

async function spacesOpen(page: Page): Promise<void> {
  await gestures.setExpanded(page, el(SPACES_GROUP), true).catch(() => undefined);
}

function treeNode(page: Page, name: string): Promise<Locator> {
  return locate(page, el(`${name} tree node inside browse tree`));
}

async function treeShows(page: Page, name: string): Promise<number> {
  await spacesOpen(page);
  return (await treeNode(page, name)).filter({visible: true}).count();
}

export const treeShowsSpace = Then('the browse tree should show the {string} space', async (page: Page, name: string) => {
  await expect.poll(() => treeShows(page, name),
    {message: `"${name}" among the nodes of the browse tree`, timeout: 60000}).toBeGreaterThan(0);
}, {tier: 'ui', description: 'opens the Spaces group on each attempt, so the claim does not depend on the product revealing the space'});

export const treeHidesSpace = Then('the browse tree should not show the {string} space', async (page: Page, name: string) => {
  await expect.poll(() => treeShows(page, name),
    {message: `"${name}" among the nodes of the browse tree`, timeout: 60000}).toBe(0);
}, {tier: 'ui', description: 'with the group open, so the absence is the server’s answer rather than a closed group'});

function cards(page: Page, name: string): Locator {
  return page.locator('.grok-gallery-grid').locator('a, .d4-link-label')
    .filter({hasText: exactText(name)}).filter({visible: true});
}

export const spaceShowsCard = Then('the space should show the {string} card', async (page: Page, name: string) => {
  await expect.poll(() => cards(page, name).count(),
    {message: `a "${name}" card in the space`, timeout: 60000}).toBeGreaterThan(0);
}, {tier: 'ui', description: 'the view rebuilds its cards on a server round-trip it does not announce'});

export const spaceHidesCard = Then('the space should not show the {string} card', async (page: Page, name: string) => {
  await expect.poll(() => cards(page, name).count(),
    {message: `a "${name}" card in the space`, timeout: 60000}).toBe(0);
}, {tier: 'ui'});
