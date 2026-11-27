/* Steps the "Files" how-to clips (features/howto/14-files) need and the library's vocabulary lacks: a file or a
   folder of My files the clip makes through the UI, removed when the feature starts and ends, and a file dropped from
   disk onto a folder of the file manager, which uploads it there (the library's drop step opens a file locally and
   waits for the "Incoming!" layer, which a folder view keeps away: `files_view.dart` prevents FILE_DRAG_ENTER), the
   text a file of My files holds, and the sharing user taken from the typeahead with its row on film. */
import {resolve} from 'path';
import type {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, expect, guideHop, locate, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

const homePath = (page: Page, name: string): Promise<string> =>
  page.evaluate((n) => `${grok.shell.user.project.name}:Home/${n}`, name);

async function removeFromHome(page: Page, path: string): Promise<void> {
  await page.evaluate(async (p) => {
    if (await grok.dapi.files.exists(p))
      await grok.dapi.files.delete(p);
  }, path);
  await expect.poll(() => page.evaluate((p) => grok.dapi.files.exists(p), path),
    {message: `the file ${path} on the server`, timeout: pollMs(30000)}).toBe(false);
}

export const noFileInHome = Given('no file {string} is in the user\'s files, now or when the feature ends', async (page: Page, name: string) => {
  const path = await homePath(page, name);
  atFeatureEnd(page, () => removeFromHome(page, path));
  await removeFromHome(page, path);
}, {tier: 'api', description: 'the file or folder of that name in My files of the signed-in account is deleted now and at feature end, read back gone'});

export const fileInHome = Then('the file {string} should be in the user\'s files on the server', async (page: Page, name: string) => {
  const path = await homePath(page, name);
  await expect.poll(() => page.evaluate((p) => grok.dapi.files.exists(p), path),
    {message: `the file ${path} on the server`, timeout: pollMs(30000)}).toBe(true);
}, {tier: 'api'});

export const noFileInHomeNow = Then('the file {string} should not be in the user\'s files on the server', async (page: Page, name: string) => {
  const path = await homePath(page, name);
  await expect.poll(() => page.evaluate((p) => grok.dapi.files.exists(p), path),
    {message: `the file ${path} on the server`, timeout: pollMs(30000)}).toBe(false);
}, {tier: 'api'});

export const fileInHomeHolds = Then('the file {string} in the user\'s files should hold the text {string}', async (page: Page, name: string, text: string) => {
  const path = await homePath(page, name);
  await expect.poll(() => page.evaluate(async (p) => await grok.dapi.files.exists(p) ? String(await grok.dapi.files.readAsText(p)) : '(no file)', path),
    {message: `the text of ${path} on the server`, timeout: pollMs(30000)}).toContain(text);
}, {tier: 'api', description: 'the file as the server keeps it holds the text'});

/* The library's `user picks the sharing user in …` clicks the typeahead's row with no stop of its own, so the video
   never shows the drop-down and the press lands on nothing (STATUS: dashboards agent). The same pick, with the row
   pictured and lit before the press. */
export const pickSharingUserShown = When('user picks the sharing user from the suggestions of {element}', async (page: Page, target: ElementRef) => {
  const login = process.env.DATAGROK_SHARING_LOGIN;
  if (!login)
    throw new Error('no second account: set DATAGROK_SHARING_LOGIN');
  const host = (await locate(page, target)).filter({visible: true}).first();
  const editor = await host.evaluate((e) => e.tagName === 'INPUT') ? host : host.locator('input').first();
  await editor.click();
  await editor.pressSequentially(login.split('@')[0]);
  const shown = new RegExp(login.split('@')[0].replace(/[^a-z0-9]/gi, ''), 'i');
  const row = page.locator('.d4-user-selector-drop-down tr').filter({has: page.locator('.d4-user-selector-user-name', {hasText: shown})}).first();
  await expect(row, `the "${login}" row of the user typeahead`).toBeVisible({timeout: pollMs(15000)});
  await guideHop(page, row);
  await row.click();
}, {tier: 'ui', description: 'types the login of DATAGROK_SHARING_LOGIN and takes its row from the typeahead, the row lit on film before the press'});

/* The browser's own drag of a file from disk (DevTools protocol, as in the library's drop step), released in the
   lower part of the element: a folder's gallery lists its items at the top, and a drop on a subfolder's item uploads
   into that subfolder. The upload runs in the task bar ("Files upload"): claim the file it makes. */
export const dropToUpload = When('user drops the {string} file of the project onto {element} to upload it', async (page: Page, file: string, target: ElementRef) => {
  const loc = (await locate(page, target)).filter({visible: true}).first();
  await expect(loc, target.phrase).toBeVisible({timeout: pollMs(15000)});
  const box = (await loc.boundingBox())!;
  const at = {x: box.x + box.width / 2, y: box.y + box.height * 0.75};
  // a stop on the folder as it is, then the pointer goes where the file is released: on film it arrives there as
  // the file shows up
  await guideHop(page, loc);
  await page.mouse.move(at.x, at.y, {steps: 6});
  const cdp = await page.context().newCDPSession(page);
  const data = {items: [], files: [resolve(process.env.BDD_ROOT ?? process.cwd(), file)], dragOperationsMask: 1};
  try {
    await cdp.send('Input.dispatchDragEvent', {type: 'dragEnter', ...at, data});
    await cdp.send('Input.dispatchDragEvent', {type: 'dragOver', ...at, data});
    await cdp.send('Input.dispatchDragEvent', {type: 'drop', ...at, data});
  }
  finally {
    await cdp.detach().catch(() => undefined);
  }
}, {tier: 'ui', description: 'the browser\'s drag of the file from disk, dropped in the lower part of the element (a folder\'s gallery or tree node uploads it there)'});
