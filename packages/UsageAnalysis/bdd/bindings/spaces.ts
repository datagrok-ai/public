/* Spaces as the Browse tree shows them. Everything general — the browse panel, the browse tree,
   the expanded state of a node — is platform vocabulary in the library; what stays here is the
   space view's own gallery and search, and the spaces a scenario leaves on the server. */
import {expect, Page} from '@playwright/test';
import {element, Given, Then} from '@datagrok-libraries/bdd';
import {atFeatureEnd} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/* The gallery a space view shows, and its search — the same two elements in the Spaces list and
   inside a space. The cards themselves are links, so "X link in space gallery" names one without
   catching the identically-classed links of an open help pane. */
// the plain "gallery" is the platform's own element now (bindings/platform/elements.ts) — the same
// selector, so "X link in gallery" and "X link in space gallery" name the same cards
element('space gallery', {selector: '.grok-gallery-grid', aliases: ['space content']});
element('space search', {selector: '.grok-gallery-search-bar .ui-input-type-ahead'});
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
  // a listing the server refuses once (a stand under load) is one attempt, not the answer
  await expect.poll(() => page.evaluate(async (n) => {
    try {
      return (await grok.dapi.spaces.list({pageSize: 1000})).filter((s: any) => s.friendlyName === n || s.name === n).length;
    } catch (e) {
      return `the listing failed: ${String(e)}`;
    }
  }, name), {message: `spaces the server holds under "${name}"`, timeout: 60000}).toBe(count);
}, {tier: 'api', description: 'what the server holds, not what the tree draws — the refusal of a duplicate is a space that was never created'});
