/* Spaces as the Browse tree shows them. Only two things are not platform vocabulary: the browse
   panel a bdd page keeps closed, and the spaces a scenario leaves on the server. */
import {expect, Page} from '@playwright/test';
import {element, Given, Then} from '@datagrok-libraries/bdd';
import {atFeatureEnd} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* The gallery a space view shows, and its search — the same two elements in the Spaces list and
   inside a space. The cards themselves are links, so "X link in space gallery" names one without
   catching the identically-classed links of an open help pane. */
element('space gallery', {selector: '.grok-gallery-grid', aliases: ['space content', 'gallery']});
element('space search', {selector: '.grok-gallery-search-bar .ui-input-type-ahead'});
/* The access level in the sharing dialog is not a <select> any more (it was in the August spec):
   a Dart privilege selector with its current level as text and a popup behind the triangle. */
element('share access selector', {selector: '[name="div-share-selector"]'});

/* grok.dapi.spaces.filter('name = "…"') answers nothing on a stand that holds the space (probed
   2026-09-10 against dev, by grok name, friendly name and both), so a space is found in the list. */
async function deleteSpaces(page: Page, names: string[]): Promise<void> {
  await page.evaluate(async (wanted) => {
    for (const space of await grok.dapi.spaces.list({pageSize: 1000}))
      if (wanted.includes(space.friendlyName) || wanted.includes(space.name))
        await grok.dapi.spaces.delete(space);
  }, names);
}

export const browsePanelOpen = Given('the browse panel is open', async (page: Page) => {
  await page.evaluate(() => {
    grok.shell.windows.simpleMode = false;
    grok.shell.windows.showBrowse = true;
  });
  await page.locator('[name="tree-Spaces"]').first().waitFor({timeout: 60000});
  atFeatureEnd(page, () => page.evaluate(() => {
    grok.shell.windows.simpleMode = true;
  }));
}, {description: 'the Browse tree, which a bdd page hides: "user is logged in" puts the shell in simple mode, and the feature puts it back at the end'});

export const noSpaceOnServer = Given('no space named {string} is on the server', async (page: Page, name: string) => {
  const names = name.split(',').map((n) => n.trim()).filter(Boolean);
  await deleteSpaces(page, names);
  atFeatureEnd(page, () => deleteSpaces(page, names));
}, {tier: 'api', description: 'deletes what an earlier run left under those names (comma-separated), and deletes them again when the feature ends'});

export const nodeExpanded = Given('the {string} tree node is expanded', async (page: Page, path: string) => {
  const key = path.replace(/\s*>\s*/g, '---').replace(/\s+/g, '-');
  const twistie = page.locator(`[name="tree-expander-${key}"]`).first();
  await twistie.waitFor({state: 'visible', timeout: 30000});
  if (!((await twistie.getAttribute('class')) ?? '').includes('d4-tree-view-tri-expanded'))
    await twistie.click();
  await expect(twistie).toHaveClass(/d4-tree-view-tri-expanded/);
}, {description: 'idempotent, unlike "user expands": a Dart tree node carries no aria-expanded, so the generic step toggles blindly and closes a group that is already open — this reads the twistie\'s own class'});

export const spacesOnServer = Then('{int} space(s) named {string} should be on the server', async (page: Page, count: number, name: string) => {
  await expect.poll(() => page.evaluate(async (n) =>
    (await grok.dapi.spaces.list({pageSize: 1000})).filter((s: any) => s.friendlyName === n || s.name === n).length,
  name), {message: `spaces the server holds under "${name}"`}).toBe(count);
}, {tier: 'api', description: 'what the server holds, not what the tree draws — the refusal of a duplicate is a space that was never created'});
